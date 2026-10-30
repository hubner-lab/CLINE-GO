# manhattan_utils.R — shared helpers for Manhattan / Miami plots
# Usage: source("/pipeline/scripts/R/utils/manhattan_utils.R")
# Requires: dplyr, ggplot2, stringr, scales, scattermore loaded by the sourcing script

# Parse "METHOD:ADJUST:FILEPATH,..." string into a named list.
parse_assoc_files_str <- function(files_str) {
    assoc_list <- stringr::str_split(files_str, ',')[[1]]
    assoc_info <- list()
    for (item in assoc_list) {
        parts <- stringr::str_split(item, ':')[[1]]
        if (length(parts) == 3) {
            method   <- parts[1]
            adjust   <- parts[2]
            filepath <- parts[3]
            assoc_info[[method]] <- list(method = method, adjust = adjust, filepath = filepath)
        }
    }
    return(assoc_info)
}

# Compute cumulative genomic positions and chromosome midpoints.
# Returns list(data = df_with_pos_cum_and_log10p, chr_info = chr_info_df).
prepare_manhattan_data <- function(df, chr_col = "chr", pos_col = "pos", pval_col) {
    chr_order <- df[[chr_col]] %>% unique() %>%
        .[order(as.numeric(stringr::str_extract(., "\\d+")))]

    df <- df %>%
        dplyr::mutate(chr_f = factor(.data[[chr_col]], levels = chr_order))

    chr_lengths <- df %>%
        dplyr::group_by(chr_f) %>%
        dplyr::summarise(chr_len = max(.data[[pos_col]]), .groups = 'drop')

    # 2% gap between chromosomes
    chr_gap <- mean(chr_lengths$chr_len) * 0.02

    chr_lengths <- chr_lengths %>%
        dplyr::mutate(
            chr_idx = dplyr::row_number(),
            tot     = cumsum(as.numeric(chr_len)) - chr_len + (chr_idx - 1) * chr_gap,
            center  = tot + chr_len / 2
        ) %>%
        dplyr::select(-chr_idx)

    df <- df %>%
        dplyr::left_join(chr_lengths %>% dplyr::select(chr_f, tot), by = "chr_f") %>%
        dplyr::mutate(
            pos_cum = .data[[pos_col]] + tot,
            log10p  = -log10(.data[[pval_col]])
        )

    list(data = df, chr_info = chr_lengths)
}

# Shared ggplot2 theme for Manhattan and Miami plots.
theme_manhattan <- function() {
    ggplot2::theme_minimal() +
    ggplot2::theme(
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank(),
        panel.border       = ggplot2::element_blank(),
        axis.line.y        = ggplot2::element_line(color = "grey30", linewidth = 0.3),
        axis.ticks.y       = ggplot2::element_line(color = "grey30", linewidth = 0.3),
        axis.text.x        = ggplot2::element_text(angle = 60, hjust = 1, size = 8),
        axis.text.y        = ggplot2::element_text(size = 9),
        axis.title         = ggplot2::element_text(size = 10),
        plot.title         = ggplot2::element_text(size = 11, face = "bold", hjust = 0.5),
        plot.subtitle      = ggplot2::element_text(size = 9, hjust = 0.5, color = "grey40"),
        legend.position    = "none",
        plot.margin        = ggplot2::margin(10, 15, 10, 10)
    )
}

# Alternating dark/light blue chromosome color palette.
get_chr_colors <- function(n_chr) {
    rep(c("#2166AC", "#92C5DE"), length.out = n_chr)
}

# Okabe-Ito colorblind-safe palette for traits (up to 8); turbo for more.
get_trait_colors <- function(traits) {
    okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                   "#0072B2", "#D55E00", "#CC79A7", "#999999")
    if (length(traits) <= 8) {
        setNames(okabe_ito[1:length(traits)], traits)
    } else {
        setNames(viridis::turbo(length(traits)), traits)
    }
}

# Up to 10 distinct region highlight colors.
get_region_colors <- function(n_regions) {
    if (n_regions <= 10) {
        colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                    "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5")
        return(colors[1:n_regions])
    } else {
        return(scales::hue_pal()(n_regions))
    }
}

# Distinct point shapes for methods (up to 10).
get_method_shapes <- function(methods) {
    shapes <- c(16, 17, 15, 18, 3, 4, 8, 6, 9, 14)
    setNames(shapes[1:length(methods)], methods)
}

# scattermore layer with pre-computed per-chromosome colors.
# chr_colors: named vector (names = chromosome levels from chr_f column).
add_scatter_layer <- function(data, chr_colors, alpha = 0.5, pointsize = 10) {
    # scattermore's C routine (scatter_points_rgbwt) aborts on NA/NaN/Inf coordinates
    # ("NA/NaN/Inf in foreign function call"). Drop non-finite points: NA log10p occurs
    # for WZA windows with no valid SNPs for a trait; +Inf occurs when a p-value is
    # exactly 0 (-log10(0)). geom_point tolerates these silently, but scattermore does not.
    data <- data[is.finite(data$log10p) & is.finite(data$pos_cum), ]
    if (nrow(data) == 0) {
        return(ggplot2::geom_blank())
    }
    data$point_color <- chr_colors[as.character(data$chr_f)]
    scattermore::geom_scattermore(data = data,
                                  ggplot2::aes(x = pos_cum, y = log10p, color = point_color),
                                  alpha = alpha, pointsize = pointsize,
                                  pixels = c(3000, 1200), interpolate = FALSE)
}

# lambda_GC of one method x trait for its QQ panel (d97202), from the {method}_lambda_K{k}.tsv
# assoc_lambda.R writes beside the p-value table, never recomputed here. Stops when the table
# has no row, or several, for the pair: the plot rule declares the table as an input, so that
# means the table is not the one written for these p-values.
# @return list(lambda_gc, correction)
qq_lambda_row <- function(lambda_table, method, trait) {
    # A data.frame, not fread's data.table: inside a data.table's [ ], `method` and `trait`
    # would name its columns, so every row would match.
    lt  <- as.data.frame(lambda_table)
    hit <- lt[lt$method == method & lt$trait == trait, , drop = FALSE]
    if (nrow(hit) != 1L)
        stop("lambda table has ", nrow(hit), " rows for method ", method, " and trait ", trait,
             " (expected 1)")
    list(lambda_gc = as.numeric(hit$lambda_gc), correction = as.character(hit$correction))
}

# The QQ panel's lambda_GC, as a plotmath string. A method that rescales its p-values by its own
# inflation factor (correction genomic_control: LFMM; rdadapt_gif: RDA) has lambda 1 by
# construction, so its panel says that instead of printing 1.000, which would read as a
# calibration result. NA when there is no value to show.
qq_lambda_label <- function(lambda_gc, correction) {
    if (!identical(correction, "none")) return('lambda[GC] == 1 ~ "by construction"')
    if (is.na(lambda_gc)) return(NA_character_)
    sprintf('lambda[GC] == "%.3f"', lambda_gc)
}
