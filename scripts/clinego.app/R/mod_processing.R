

#' Processing QC tab server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @noRd
mod_processing_server <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        summary_data <- shiny::reactive({
            pd <- project_data()
            load_pipeline_summary(pd$name)
        })

        summary_val <- function(metric_key, default = "—") {
            shiny::reactive({
                s <- summary_data()
                if (nrow(s) == 0) return(default)
                row <- s[s$step == "processing" & s$metric == metric_key, ]
                if (nrow(row) == 0) default else as.character(row$value[1])
            })
        }

        # ── Value boxes (dynamic themes + conditional het/depth boxes) ─────────
        output$summary_boxes <- shiny::renderUI({
            pd  <- project_data()
            # Convert to data.frame to avoid data.table column-name shadowing:
            # inside data.table's [, bare names like `metric` resolve to the column
            # named `metric` in the table, not the function parameter.
            s <- as.data.frame(summary_data())

            sv <- function(metric, default = NA) {
                row <- s[s$step == "processing" & s$metric == metric, ]
                if (nrow(row) == 0) default else suppressWarnings(as.numeric(row$value[1]))
            }
            sv_str <- function(metric, default = "\u2014") {
                row <- s[s$step == "processing" & s$metric == metric, ]
                if (nrow(row) == 0) default else as.character(row$value[1])
            }

            # Samples — warning if >10% removed
            n_samp_raw <- sv("samples_total")
            n_samp_flt <- sv("samples_after_filtering")
            samp_theme <- if (!is.na(n_samp_raw) && !is.na(n_samp_flt) && n_samp_raw > 0) {
                if ((n_samp_raw - n_samp_flt) / n_samp_raw > 0.1) "warning" else "success"
            } else "primary"
            samp_box <- lab_value_box(
                "Samples",
                lab_arrow_value(sv_str("samples_total"), sv_str("samples_after_filtering")),
                samp_theme, "people-fill")

            # SNPs — warning if >50% removed
            n_snp_raw <- sv("snps_raw")
            n_snp_flt <- sv("snps_after_filtering")
            snp_theme <- if (!is.na(n_snp_raw) && !is.na(n_snp_flt) && n_snp_raw > 0) {
                if ((n_snp_raw - n_snp_flt) / n_snp_raw > 0.5) "warning" else "success"
            } else "info"
            snp_box <- lab_value_box(
                "SNPs", lab_arrow_value(sv_str("snps_raw"), sv_str("snps_after_filtering")),
                snp_theme, "database-fill")

            # SNPs LD-pruned — always info
            ld_box <- lab_value_box("SNPs (LD-pruned)", sv_str("snps_after_ld_pruning"),
                                    "info", "scissors")

            # Ti/Tv — color-coded by quality
            titv_val <- sv("titv_ratio")
            titv_theme <- if (is.na(titv_val)) "secondary"
                          else if (titv_val >= 2.0) "success"
                          else if (titv_val >= 1.5) "warning"
                          else "danger"
            titv_box <- lab_value_box("Ti/Tv Ratio", sv_str("titv_ratio"),
                                      titv_theme, "arrow-left-right")

            # Depth filter — conditional on config
            min_dp <- config_get(pd$config, "Filter", "min_depth", default = NULL)
            max_dp <- config_get(pd$config, "Filter", "max_depth", default = NULL)
            depth_box <- if (!is.null(min_dp) || !is.null(max_dp)) {
                dp_parts <- c(
                    if (!is.null(min_dp)) paste0("min=", min_dp),
                    if (!is.null(max_dp)) paste0("max=", max_dp)
                )
                lab_value_box("Depth Filter", paste(dp_parts, collapse = ", "),
                              "info", "layers")
            } else NULL

            extra_boxes <- Filter(Negate(is.null), list(depth_box))
            n_extra <- length(extra_boxes)
            width <- if (n_extra == 0) 1/4 else if (n_extra == 1) 1/5 else 1/6

            bslib::layout_column_wrap(
                width = width,
                fill  = FALSE,
                samp_box, snp_box, ld_box, titv_box,
                !!!extra_boxes
            )
        })

        # ── LEA conversion losses alert ────────────────────────────────────────
        # Converted from a full-width alert banner to a compact hover badge (same
        # wording, now in the tooltip) — tidy FYI status, not a page-wide banner.
        output$lea_losses_alert <- shiny::renderUI({
            pd <- project_data()
            p  <- removed_snps_path(pd$name)
            if (!nzchar(p) || !file.exists(p)) return(NULL)
            lines <- tryCatch(readLines(p, warn = FALSE), error = function(e) character(0))
            n <- length(lines)
            if (n == 0) return(NULL)
            htmltools::div(
                class = "d-flex justify-content-end mb-2",
                filter_note(
                    paste0(n, " removed"),
                    htmltools::p(
                        htmltools::tags$strong(n, if (n == 1) "SNP" else "SNPs",
                                           "removed during LEA format conversion"),
                        " (multi-allelic or invariant after VCF\u2192LFMM conversion). ",
                        "These are excluded from PCA, sNMF, and association analyses."
                    ),
                    class = "bg-secondary"
                )
            )
        })

        # Keep legacy outputs in case anything still references them (safe no-ops)
        output$samples_display <- shiny::renderUI(NULL)
        output$snps_display    <- shiny::renderUI(NULL)
        output$snps_ld         <- shiny::renderText("\u2014")
        output$titv            <- shiny::renderText("\u2014")

        # ── Image cards ────────────────────────────────────────────────────────
        shiny::observe({
            pd <- project_data()

            mod_image_card_server("attrition",
                path    = shiny::reactive(qc_plot_path(pd$name, "filtering_attrition.png")),
                title   = shiny::reactive("Filtering Attrition"),
                dl_name = shiny::reactive(paste0(pd$name, "_filtering_attrition")),
                suggestion = shiny::reactive("Run mode=processing to generate QC plots"),
                note    = shiny::reactive(filter_note(
                    "stages",
                    "SNP counts surviving each filter stage: MAF → missingness → LD pruning."
                ))
            )
            mod_image_card_server("sample_miss",
                path    = shiny::reactive(qc_plot_path(pd$name, "sample_missingness_distribution.png")),
                title   = shiny::reactive("Sample Missingness"),
                dl_name = shiny::reactive(paste0(pd$name, "_sample_missingness")),
                note    = shiny::reactive({
                    t <- config_get(pd$config, "Filter", "sample_miss", default = NULL)
                    if (is.null(t)) return(NULL)
                    filter_note(
                        paste0("miss > ", t),
                        paste0("Samples missing more than ", t, " of genotypes are dropped.")
                    )
                })
            )
            mod_image_card_server("het_miss",
                path    = shiny::reactive(qc_plot_path(pd$name, "het_vs_missingness.png")),
                title   = shiny::reactive("Heterozygosity vs Missingness"),
                dl_name = shiny::reactive(paste0(pd$name, "_het_vs_missingness")),
                note    = shiny::reactive({
                    sd_thresh <- config_get(pd$config, "Filter", "het_outlier_sd", default = NULL)
                    miss_thresh <- config_get(pd$config, "Filter", "sample_miss", default = NULL)
                    if (is.null(sd_thresh) && is.null(miss_thresh)) return(NULL)
                    # Since release 4c the plot shows what the het filter did when it ran;
                    # without it, a ±3 SD preview on the retained samples (plot_qc_processing.R).
                    txt <- if (!is.null(sd_thresh)) {
                        paste0("Red: samples the heterozygosity filter removed, beyond ±", sd_thresh,
                              " population SD of the samples that passed the missingness filter. ",
                              "The vertical line is the missingness cutoff.")
                    } else {
                        paste0("Red: retained samples beyond ±3 SD, shown only: no heterozygosity ",
                              "filter is set. The vertical line is the missingness cutoff.")
                    }
                    filter_note(
                        if (!is.null(sd_thresh)) paste0("±", sd_thresh, " SD") else "±3 SD preview",
                        txt
                    )
                })
            )
            mod_image_card_server("relatedness",
                path    = shiny::reactive(qc_plot_path(pd$name, "relatedness_distribution.png")),
                title   = shiny::reactive("Relatedness (IBS allele-sharing)"),
                dl_name = shiny::reactive(paste0(pd$name, "_relatedness_distribution")),
                note    = shiny::reactive({
                    thresh <- config_get(pd$config, "Filter", "relatedness", default = NULL)
                    if (is.null(thresh) || is.na(thresh)) return(NULL)
                    action <- config_get(pd$config, "Filter", "relatedness_action", default = "keep")
                    pairs   <- load_relatedness_pairs(pd$name)
                    removed <- load_relatedness_removed(pd$name)
                    n_pairs_above <- if (nrow(pairs) == 0) NA_integer_
                                      else sum(pairs$IBS > thresh, na.rm = TRUE)
                    n_would_remove <- if (nrow(removed) == 0) NA_integer_ else nrow(removed)
                    relatedness_note(thresh, action, n_pairs_above, n_would_remove)
                })
            )
            mod_image_card_server("relatedness_mds",
                path    = shiny::reactive(qc_plot_path(pd$name, "relatedness_mds.png")),
                title   = shiny::reactive("Relatedness MDS (IBS)"),
                dl_name = shiny::reactive(paste0(pd$name, "_relatedness_mds")),
                note    = shiny::reactive({
                    thresh <- config_get(pd$config, "Filter", "relatedness", default = NULL)
                    if (is.null(thresh) || is.na(thresh)) return(NULL)
                    removed <- load_relatedness_removed(pd$name)
                    n_would_remove <- if (nrow(removed) == 0) NA_integer_ else nrow(removed)
                    filter_note(
                        if (is.na(n_would_remove)) "0 discarded" else paste0(n_would_remove, " discarded"),
                        "Samples that cluster near an IBS of 1.0 (duplicates/clones) are marked red."
                    )
                })
            )
            mod_image_card_server("snp_miss",
                path    = shiny::reactive(qc_plot_path(pd$name, "snp_missingness_distribution.png")),
                title   = shiny::reactive("SNP Missingness"),
                dl_name = shiny::reactive(paste0(pd$name, "_snp_missingness")),
                note    = shiny::reactive({
                    t <- config_get(pd$config, "Filter", "snp_miss", default = NULL)
                    if (is.null(t)) return(NULL)
                    filter_note(
                        paste0("miss > ", t),
                        paste0("SNPs missing in more than ", t, " of samples are dropped.")
                    )
                })
            )
            mod_image_card_server("maf",
                path    = shiny::reactive(qc_plot_path(pd$name, "maf_distribution.png")),
                title   = shiny::reactive("MAF Distribution"),
                dl_name = shiny::reactive(paste0(pd$name, "_maf_distribution")),
                note    = shiny::reactive({
                    t <- config_get(pd$config, "Filter", "maf", default = NULL)
                    if (is.null(t)) return(NULL)
                    # What the MAF step itself removed, from filtering_summary.tsv's rows for
                    # plink's two SNP steps (--geno runs first, then --maf). The histogram's red
                    # bars count every raw SNP below the threshold, those --geno removed first
                    # included. A table from before release 4c has no "After SNP missingness
                    # filter" row, and the note then gives no count.
                    fs <- as.data.frame(load_filtering_summary(pd$name))
                    n_maf <- NA_integer_
                    if (nrow(fs) > 0) {
                        geno_row <- fs[fs$stage == "After SNP missingness filter", ]
                        maf_row  <- fs[fs$stage == "After MAF + SNP missingness filter", ]
                        if (nrow(geno_row) > 0 && nrow(maf_row) > 0)
                            n_maf <- as.integer(geno_row$n_snps[1]) - as.integer(maf_row$n_snps[1])
                    }
                    removed_txt <- if (is.na(n_maf)) "" else
                        paste0(n_maf, " SNP", if (n_maf != 1) "s" else "",
                               " removed by this filter. It runs after the SNP missingness ",
                               "filter, so the red bars can show more. ")
                    filter_note(
                        paste0("MAF < ", t),
                        paste0(removed_txt, "Variants with minor-allele frequency below ", t,
                              " are removed before analysis.")
                    )
                })
            )
            mod_image_card_server("snp_density",
                path    = shiny::reactive(qc_plot_path(pd$name, "snp_density_by_chr.png")),
                title   = shiny::reactive("SNP Density by Chromosome"),
                dl_name = shiny::reactive(paste0(pd$name, "_snp_density")),
                note    = shiny::reactive(help_note("snp_density"))
            )
        })

        # ── Filtering summary table ────────────────────────────────────────────
        output$filtering_table <- DT::renderDataTable({
            pd <- project_data()
            df <- as.data.frame(load_filtering_summary(pd$name))
            safe_datatable(
                df,
                options  = list(dom = "t", scrollX = TRUE),
                rownames = FALSE
            )
        })

        # ── Depth section (always shown — placeholder explains unavailable states) ──
        # Three states written by write_summary.R: not_provided / declared_empty / available.
        depth_status <- shiny::reactive({
            s <- summary_data()
            if (nrow(s) == 0) return("not_provided")
            row <- s[s$step == "processing" & s$metric == "depth_qc", ]
            if (nrow(row) == 0) return("not_provided")
            row$value[1]
        })

        output$depth_section <- shiny::renderUI({
            htmltools::tagList(
                mod_image_card_ui(ns("depth")),
                shiny::br()
            )
        })

        shiny::observe({
            pd <- project_data()
            mod_image_card_server("depth",
                path    = shiny::reactive(qc_plot_path(pd$name, "depth_distribution.png")),
                title   = shiny::reactive("Depth Distribution"),
                dl_name = shiny::reactive(paste0(pd$name, "_depth_distribution")),
                placeholder = shiny::reactive(switch(depth_status(),
                    not_provided   = "VCF has no FORMAT/DP field",
                    declared_empty = "FORMAT/DP is declared, but no genotype depth was read",
                    "Plot not available"
                )),
                # 7ca55e: until release 4c the pipeline read depth from the filtered VCF, which
                # plink writes without DP, so every project with DP landed in declared_empty and
                # this text blamed the data. Now depth comes from the input VCF.
                suggestion = shiny::reactive(switch(depth_status(),
                    not_provided   = "Depth QC needs per-genotype DP in the VCF.",
                    declared_empty = paste0("Depth comes from each genotype's DP in the input VCF. ",
                                            "Results made before 2026-10-30 read the filtered VCF, ",
                                            "which has none: re-run mode=processing."),
                    NULL
                )),
                note    = shiny::reactive({
                    min_d <- config_get(pd$config, "Filter", "min_depth", default = NULL)
                    max_d <- config_get(pd$config, "Filter", "max_depth", default = NULL)
                    if (is.null(min_d) && is.null(max_d)) return(NULL)
                    lo <- if (!is.null(min_d)) min_d else "–"
                    hi <- if (!is.null(max_d)) max_d else "–"
                    # depth_filter masks calls (bcftools filter -S .), it removes no site
                    filter_note(
                        paste0(lo, "–", hi, "×"),
                        paste0("Genotype calls with DP outside ", lo, "–", hi,
                              " are set to missing before the SNP filters. ",
                              "The plot shows depth before that masking.")
                    )
                })
            )
        })

        # ── Sample heterozygosity table ────────────────────────────────────────
        het_data <- shiny::reactive({
            pd <- project_data()
            as.data.frame(load_sample_heterozygosity(pd$name))
        })

        output$het_table <- DT::renderDataTable({
            safe_datatable(
                het_data(),
                options  = list(dom = "frtip", scrollX = TRUE, pageLength = 25),
                rownames = FALSE
            )
        })

        output$dl_het <- shiny::downloadHandler(
            filename = function() paste0("sample_heterozygosity_", project_data()$name, ".csv"),
            content  = function(file) utils::write.csv(het_data(), file, row.names = FALSE)
        )
    })
}

#' processing tab UI — dashboard layout.
#'
#' Promoted from scripts/layout_lab/ on 2026-08-29; the comments inside
#' carry the measurements behind each sizing decision.
#' @param id module namespace id
#' @noRd
mod_processing_ui <- function(id) {
    ns <- shiny::NS(id)

    lab_root("a", module = "processing",

        lab_kpi_row(ns),

        htmltools::div(
            class = "lab-overview",
            bslib::layout_columns(
                col_widths = c(5, 7),
                # layout_columns defaults to fill=TRUE, which would stretch the
                # columns and fight the fraction-based heights below.
                fill = FALSE, fillable = FALSE,
                gap = "0.75rem",

                # LEFT — hero figure over the short filtering table
                htmltools::div(
                    class = "lab-hero-col",
                    lab_hero(mod_image_card_ui(ns("attrition"))),
                    lab_filtering_summary_card(ns)
                ),

                # RIGHT — grouped small multiples
                htmltools::div(
                    class = "lab-multiples-col",

                    lab_section_header("Sample QC", icon = "people-fill"),
                    lab_thumb_grid(
                        cols = 4,
                        lab_thumb(mod_image_card_ui(ns("sample_miss"))),
                        lab_thumb(mod_image_card_ui(ns("het_miss"))),
                        lab_thumb(mod_image_card_ui(ns("relatedness"))),
                        lab_thumb(mod_image_card_ui(ns("relatedness_mds")))
                    ),

                    lab_section_header("SNP QC", icon = "database-fill"),
                    # Same 4-column rhythm as Sample QC so every tile is the same
                    # width and the two groups read as one grid. A 2-col grid made
                    # these tiles ~520px wide while the height cap held the image
                    # to ~280px, stranding it in empty space.
                    lab_thumb_grid(
                        cols = 4,
                        lab_thumb(mod_image_card_ui(ns("snp_miss"))),
                        lab_thumb(mod_image_card_ui(ns("maf"))),
                        # a chromosome strip chart is unreadable in a square tile
                        lab_thumb(mod_image_card_ui(ns("snp_density")), span = 2)
                    )
                )
            )
        ),

        # ── Below the fold: both collapsed -> outputs suspended ───────────────
        bslib::accordion(
            class = "lab-below-fold",
            open = FALSE, multiple = TRUE,
            bslib::accordion_panel(
                value = "depth",
                title = htmltools::tagList(
                    bsicons::bs_icon("layers"), " Depth Distribution ",
                    htmltools::span(class = "text-muted small",
                                    "(context — needs per-genotype DP)")
                ),
                shiny::uiOutput(ns("depth_section"))
            ),
            bslib::accordion_panel(
                value = "het",
                title = htmltools::tagList(bsicons::bs_icon("table"),
                                           " Sample Heterozygosity"),
                htmltools::div(
                    class = "d-flex justify-content-end mb-2",
                    lab_download_button(ns("dl_het"), "CSV")
                ),
                DT::DTOutput(ns("het_table"))
            )
        )
    )
}

