#' Interactive Region Explorer UI
#'
#' Replaces the static region dropdown in the Association tab.
#' Regions are auto-computed from the current sig SNPs on every filter/strategy change.
#' Clicking a region row shows genes (on-the-fly), and on-demand enrichment / regionplot
#' via processx subprocesses.
#'
#' @param id module namespace id
#' @noRd
mod_region_explorer_ui <- function(id) {
    ns <- shiny::NS(id)

    htmltools::tagList(

        # ── Region list card ─────────────────────────────────────────────────
        bslib::card(
            bslib::card_header(
                class = "d-flex justify-content-between align-items-center",
                shiny::uiOutput(ns("explorer_title"), inline = TRUE),
                shiny::uiOutput(ns("region_count_badge"), inline = TRUE)
            ),
            bslib::card_body(
                class = "p-2",
                DT::DTOutput(ns("region_table"), width = "100%")
            )
        ),

        # ── Region detail (conditional on selection) ─────────────────────────
        shiny::uiOutput(ns("region_detail_ui"))
    )
}

#' Interactive Region Explorer Server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @param module character: MOD_GEA or MOD_GWAS
#' @param interactive_sigsnps reactive data.table of current filtered sig SNPs
#' @return list with computed_regions (reactive) and selected_region_id (reactiveVal)
#' @noRd
mod_region_explorer_server <- function(id, project_data, module = MOD_GEA,
                                        interactive_sigsnps = shiny::reactive(NULL),
                                        region_distance     = shiny::reactive(2000000L),
                                        override_regions    = shiny::reactive(NULL),
                                        regime              = shiny::reactive(FALSE)) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # ── Computed regions ───────────────────────────────────────────────────
        computed_regions <- shiny::reactive({
            # Overlapping tab supplies pre-computed overlap regions directly
            ovr <- override_regions()
            if (!is.null(ovr)) return(ovr)
            snps <- interactive_sigsnps()
            if (is.null(snps) || nrow(snps) == 0) return(.empty_regions())
            compute_all_regions(snps, region_distance())
        })

        # ── Selected region state ──────────────────────────────────────────────
        selected_region_id <- shiny::reactiveVal(NULL)

        shiny::observeEvent(computed_regions(), {
            rid  <- selected_region_id()
            regs <- computed_regions()
            if (!is.null(rid) && (nrow(regs) == 0 || !rid %in% regs$region_id)) {
                selected_region_id(NULL)
            }
        }, ignoreInit = TRUE)

        # ── GFF genes (cached per project) ─────────────────────────────────────
        gff_genes <- shiny::reactive({
            pd <- project_data()
            load_gff_genes(pd$name, pd$config)
        })

        # ── Explorer title ──────────────────────────────────────────────────────
        output$explorer_title <- shiny::renderUI("Regions")

        # ── Region count badge ─────────────────────────────────────────────────
        output$region_count_badge <- shiny::renderUI({
            n   <- nrow(computed_regions())
            cls <- if (n == 0) "bg-secondary" else "bg-primary"
            htmltools::span(class = paste("badge rounded-pill", cls), n, "regions")
        })

        # ── Region list table ──────────────────────────────────────────────────
        output$region_table <- DT::renderDataTable({
            regs <- computed_regions()
            if (nrow(regs) == 0) {
                return(safe_datatable(
                    data.frame(Message = "No significant SNPs with current filters"),
                    options  = list(dom = "t"),
                    rownames = FALSE
                ))
            }

            rid <- selected_region_id()
            selected_row <- if (!is.null(rid)) which(regs$region_id == rid) else integer(0)

            disp <- data.frame(
                Region  = sapply(regs$region_id, format_region_id),
                Chr     = regs$chr,
                Length  = .format_bp(regs$length),
                `Sig SNPs` = regs$snp_count,
                Traits  = regs$traits,
                Methods = regs$methods,
                MinP    = signif(regs$min_pvalue, 3),
                stringsAsFactors = FALSE
            )

            safe_datatable(
                disp,
                selection  = list(mode = "single", selected = selected_row),
                options    = list(
                    dom        = "ftp",
                    pageLength = 10,
                    scrollX    = TRUE,
                    order      = list(list(3L, "desc"))
                ),
                rownames = FALSE
            )
        })

        shiny::observeEvent(input$region_table_rows_selected, {
            row  <- input$region_table_rows_selected
            regs <- computed_regions()
            if (length(row) == 1 && nrow(regs) >= row) {
                selected_region_id(regs$region_id[row])
            } else {
                selected_region_id(NULL)
            }
        }, ignoreNULL = FALSE)

        # ── Genes for selected region ──────────────────────────────────────────
        region_genes <- shiny::reactive({
            rid <- selected_region_id()
            shiny::req(rid)
            regs <- computed_regions()
            row  <- regs[region_id == rid]
            if (nrow(row) == 0) return(data.table::data.table())
            gff  <- gff_genes()
            if (is.null(gff) || nrow(gff) == 0) return(data.table::data.table())
            pd   <- project_data()
            plen <- config_get(pd$config, "GEA", "promoter_length", default = 10000L)
            find_genes_in_region(gff, row, as.integer(plen))
        })

        # ── SNPs for selected region ───────────────────────────────────────────
        region_snps <- shiny::reactive({
            rid  <- selected_region_id()
            shiny::req(rid)
            regs <- computed_regions()
            row  <- regs[region_id == rid]
            if (nrow(row) == 0) return(data.table::data.table())
            snps <- interactive_sigsnps()
            if (is.null(snps) || nrow(snps) == 0) return(data.table::data.table())
            snps[as.character(chr) == as.character(row$chr[1]) &
                 as.integer(pos) >= row$start[1] &
                 as.integer(pos) <= row$end[1]]
        })

        # Total SNP count in region (all SNPs in filtered VCF, not just sig).
        # Haplotype analysis operates on every variant in the region — this is
        # the relevant number for the min-SNPs check, not the sig-SNP count.
        #
        # T2.4: Pure read from vcf_snp_count_cache (reactiveValues); the blocking awk
        # command runs in a background process (see observe after is_running definition).
        # Returns NA_integer_ while pending — haplotype UI shows "unavailable" and
        # re-renders automatically once the async count arrives.
        region_total_snps <- shiny::reactive({
            row <- current_region_row()
            if (is.null(row)) return(NA_integer_)
            rid    <- row$region_id[1]
            cached <- vcf_snp_count_cache[[rid]]
            if (is.null(cached) || identical(cached, "pending")) return(NA_integer_)
            as.integer(cached)
        })

        # ── Enrichment state ───────────────────────────────────────────────────
        enrichment_cache <- shiny::reactiveValues()
        regionplot_cache <- shiny::reactiveValues()
        hap_scan_cache   <- shiny::reactiveValues()  # keyed by region_id___tag
        hap_viz_cache    <- shiny::reactiveValues()  # keyed by region_id___tag

        # Session-level cache for VCF region SNP counts (keyed by region_id).
        # "pending" = background awk in flight; integer = done; NA_integer_ = VCF missing.
        # Populated by the async launch observe and read by region_total_snps().
        vcf_snp_count_cache <- shiny::reactiveValues()

        # Handles for background processes: keyed by cache key
        # Each entry: list(process, log_file, type, ...) or NULL when done
        subprocess_handles <- shiny::reactiveValues()

        # Persists user-chosen epsilon across renderUI re-renders (NULL = use saved/config default)
        user_epsilon <- shiny::reactiveVal(NULL)

        shiny::observeEvent(input$hap_epsilon_selected, {
            user_epsilon(input$hap_epsilon_selected)
        }, ignoreInit = TRUE)

        # Per-region exploratory params (persisted to disk in region_params.json)
        # Source of truth for "modified" badges and input defaults
        region_params <- shiny::reactiveVal(list(global = list(), regions = list()))
        shiny::observe({
            pd <- project_data()
            shiny::req(pd)
            region_params(read_region_params(pd$name))
        })

        # Persists user-chosen metadata type across renderUI re-renders (NULL = use config default)
        user_meta_type <- shiny::reactiveVal(NULL)

        shiny::observeEvent(input$hap_meta_type, {
            user_meta_type(input$hap_meta_type)
        }, ignoreInit = TRUE)

        # ── Haplotype tag for this module ──────────────────────────────────────
        # expected_hap_tag: derived from user's chosen metadata type (or config default).
        # Stable across renderUI re-renders because it uses user_meta_type() not disk state.
        expected_hap_tag <- shiny::reactive({
            src <- .module_to_hap_source(module)
            if (is.null(src)) return(NULL)
            pd   <- project_data()
            meta <- user_meta_type() %||%
                    config_get(pd$config, "haplotype", "scan", "metadata_type", default = "site")
            paste0(meta, "_", src)
        })

        # ── Derived keys (region + selected trait) ────────────────────────────
        current_region_row <- shiny::reactive({
            rid  <- selected_region_id()
            shiny::req(rid)
            regs <- computed_regions()
            row  <- regs[region_id == rid]
            if (nrow(row) == 0) return(NULL)
            row
        })

        # Derive trait automatically from region (first trait); no user selector needed
        # since GO enrichment is per-region (same gene set regardless of trait).
        current_trait <- shiny::reactive({
            row <- current_region_row()
            if (is.null(row) || nrow(row) == 0) return(NULL)
            traits <- trimws(strsplit(as.character(row$traits[1]), ",")[[1]])
            traits <- traits[nzchar(traits)]
            if (length(traits) == 0) return(NULL)
            traits[1]
        })

        enrich_key <- shiny::reactive({
            rid <- selected_region_id()
            if (is.null(rid)) return(NULL)
            rid
        })

        # combo_hash: deterministic 8-char hash of sorted trait/method pairs + regime.
        # Different filter state or regime → different hash → different cache key and disk path.
        combo_hash <- shiny::reactive({
            snps <- region_snps()
            if (is.null(snps) || nrow(snps) == 0) return(NULL)
            arg <- .build_trait_method_arg(snps)
            if (is.null(arg)) return(NULL)
            regime_str <- if (isTRUE(regime())) "wza" else "snp"
            substr(digest::digest(paste0(arg, "__", regime_str), algo = "md5"), 1, 8)
        })

        rplot_key <- shiny::reactive({
            rid  <- selected_region_id()
            hash <- combo_hash()
            if (is.null(rid) || is.null(hash)) return(NULL)
            paste0(rid, "___", hash, "___rplot")
        })

        hap_key <- shiny::reactive({
            rid <- selected_region_id()
            tag <- expected_hap_tag()
            if (is.null(rid) || is.null(tag)) return(NULL)
            paste0(rid, "___", tag)
        })

        # ── Disk-cache hydration: load persisted regionplot on region/filter change ─
        shiny::observe({
            key <- rplot_key()
            if (is.null(key) || !is.null(regionplot_cache[[key]]) || is_running(key)) return()
            rid  <- selected_region_id()
            hash <- combo_hash()
            pd   <- project_data()
            if (is.null(rid) || is.null(hash) || is.null(pd)) return()
            region_safe <- gsub("[:\\-]", "_", rid)
            disk_dir  <- ondemand_regionplot_dir(pd$name, module, hash)
            disk_path <- file.path(disk_dir,
                                   paste0("regionplot_custom_", region_safe, ".png"))
            if (file_ok(disk_path)) {
                regionplot_cache[[key]] <- list(path = disk_path)
            }
        })

        # ── Disk-cache hydration: load persisted enrichment on trait/region change ─
        shiny::observe({
            key <- enrich_key()
            if (is.null(key) || !is.null(enrichment_cache[[key]]) || is_running(key)) return()
            rid <- selected_region_id()
            tr  <- current_trait()
            pd  <- project_data()
            if (is.null(rid) || is.null(tr) || is.null(pd)) return()
            tsv <- enrichment_table_path(pd$name, module, tr, rid)
            if (!file_ok(tsv)) return()
            tbl <- tryCatch(
                data.table::fread(tsv, sep = "\t", header = TRUE),
                error = function(e) NULL
            )
            if (is.null(tbl) || nrow(tbl) == 0) return()
            dotplot_p  <- enrichment_plot_path(pd$name, module, tr, rid, "dotplot")
            emapplot_p <- enrichment_plot_path(pd$name, module, tr, rid, "emapplot")
            enrichment_cache[[key]] <- list(
                table         = tbl,
                dotplot_path  = if (file_ok(dotplot_p))  dotplot_p  else NULL,
                emapplot_path = if (file_ok(emapplot_p)) emapplot_p else NULL
            )
        })

        # ── Disk-cache hydration: load pre-computed haplotype scan results ────
        # Uses exact region_id matching (like enrichment/regionplot) — no fuzzy overlap.
        # This prevents stale results from a different region_distance being loaded.
        shiny::observe({
            key <- hap_key()
            if (is.null(key) || !is.null(hap_scan_cache[[key]]) || is_running(paste0(key, "__scan"))) return()
            rid <- selected_region_id()
            tag <- expected_hap_tag()
            pd  <- project_data()
            if (is.null(rid) || is.null(tag) || is.null(pd)) return()
            # Exact match only: check if HapObject or clustree exists for this exact region_id
            inter_dir <- hap_intermediate_dir(pd$name, tag)
            plots_dir <- hap_scan_plots_dir(pd$name, tag)
            hapobj    <- file.path(inter_dir, paste0("Region_", rid, "_HapObject.qs"))
            mg_png    <- file.path(plots_dir, paste0("Region_", rid, "_clustree_MG.png"))
            if (file_ok(hapobj) || file_ok(mg_png)) {
                fake_handle <- list(
                    region_id = rid,
                    plots_dir = plots_dir,
                    inter_dir = inter_dir
                )
                hap_scan_cache[[key]] <- .collect_hap_scan_result(fake_handle)
            }
        })

        # ── Disk-cache hydration: load pre-computed haplotype viz results ────
        # Exact region_id matching only — no fuzzy overlap.
        shiny::observe({
            key <- hap_key()
            if (is.null(key) || !is.null(hap_viz_cache[[key]]) || is_running(paste0(key, "__viz"))) return()
            rid <- selected_region_id()
            tag <- expected_hap_tag()
            pd  <- project_data()
            if (is.null(rid) || is.null(tag) || is.null(pd)) return()
            # Exact match: look for crosshap_viz PNGs with this exact region_id
            viz_dir <- hap_viz_plots_dir(pd$name, tag)
            viz_files <- Sys.glob(file.path(viz_dir, paste0("Region_", rid, "_crosshap_viz_*.png")))
            if (length(viz_files) > 0) {
                traits <- sub("\\.png$", "", sub(paste0("^Region_", rid, "_crosshap_viz_"), "", basename(viz_files)))
                hap_viz_cache[[key]] <- list(traits = traits, matched_rid = rid)
            }
        })

        # ── Helper: is a subprocess currently running for this key? ────────────
        is_running <- function(key) {
            if (is.null(key)) return(FALSE)
            h <- subprocess_handles[[key]]
            !is.null(h) && h$process$is_alive()
        }

        # ── VCF SNP count: launch async awk subprocess on region selection (T2.4) ──
        # Decoupled from region_total_snps() (which is a pure read reactive) so the
        # session never blocks on the ~2M-line VCF scan. Memoized per region_id for
        # the session lifetime — re-selecting the same region is instant.
        shiny::observe({
            row <- current_region_row()
            if (is.null(row)) return()
            rid        <- row$region_id[1]
            handle_key <- paste0("vcf_count_", rid)

            # Already cached or in flight — nothing to do
            if (!is.null(vcf_snp_count_cache[[rid]]) || is_running(handle_key)) return()

            pd  <- project_data()
            vcf <- hap_filtered_vcf_path(pd$name, pd$config)
            if (is.null(vcf) || !nzchar(vcf) || !file_ok(vcf)) {
                vcf_snp_count_cache[[rid]] <- NA_integer_
                return()
            }

            # Mark as pending so parallel reactive reads don't re-launch
            vcf_snp_count_cache[[rid]] <- "pending"

            # Use a temp file for output — mirrors the file-redirect pattern used by
            # all other subprocess handlers (proved pattern; avoids pipe-buffer uncertainty).
            out_file <- tempfile("vcf_count_", fileext = ".txt")
            cmd <- sprintf(
                "awk -v c=%s -v s=%d -v e=%d 'BEGIN{n=0} !/^#/ && $1==c && $2>=s && $2<=e {n++} END{print n}' %s > %s",
                shQuote(as.character(row$chr[1])),
                as.integer(row$start[1]), as.integer(row$end[1]),
                shQuote(vcf),
                shQuote(out_file)
            )

            subprocess_handles[[handle_key]] <- list(
                process     = processx::process$new("bash", c("-c", cmd)),
                type        = "vcf_count",
                region_id   = rid,
                output_file = out_file
            )
        })

        # ── Polling observer: checks background processes every second ──────────
        shiny::observe({
            handles_list <- shiny::reactiveValuesToList(subprocess_handles)
            active <- Filter(Negate(is.null), handles_list)
            if (length(active) == 0) return()

            shiny::invalidateLater(1000)

            for (key in names(active)) {
                h <- active[[key]]
                if (h$process$is_alive()) next  # still running

                status <- tryCatch(h$process$get_exit_status(), error = function(e) 1L)

                if (h$type == "enrichment") {
                    result <- tryCatch(
                        .collect_enrichment_result(h),
                        error = function(e) list(error = conditionMessage(e))
                    )
                    if (!is.null(result$error) && status != 0L) {
                        log_lines <- tryCatch(readLines(h$log_file, warn = FALSE), error = function(e) character())
                        result$log_tail <- paste(utils::tail(log_lines, 10), collapse = "\n")
                    }
                    enrichment_cache[[key]] <- result

                } else if (h$type == "regionplot") {
                    result <- tryCatch(
                        .collect_regionplot_result(h),
                        error = function(e) list(error = conditionMessage(e))
                    )
                    if (!is.null(result$error) && status != 0L) {
                        log_lines <- tryCatch(readLines(h$log_file, warn = FALSE), error = function(e) character())
                        result$log_tail <- paste(utils::tail(log_lines, 10), collapse = "\n")
                    }
                    regionplot_cache[[key]] <- result

                } else if (h$type == "hap_scan") {
                    result <- tryCatch(
                        .collect_hap_scan_result(h),
                        error = function(e) list(error = conditionMessage(e))
                    )
                    if (!is.null(result$error) && status != 0L) {
                        log_lines <- tryCatch(readLines(h$log_file, warn = FALSE), error = function(e) character())
                        result$log_tail <- paste(utils::tail(log_lines, 10), collapse = "\n")
                    }
                    # Create scan_done.flag so find_haplotype_tags() can discover this scan on restart
                    if (is.null(result$error)) {
                        flag_dir  <- file.path(dirname(dirname(h$inter_dir)), "flags")
                        flag_path <- file.path(flag_dir,
                                               paste0("haplotype_", h$tag, "_scan_done.flag"))
                        dir.create(flag_dir, recursive = TRUE, showWarnings = FALSE)
                        if (!file.exists(flag_path))
                            file.create(flag_path, showWarnings = FALSE)
                        # Persist scan params to region_params.json
                        if (!is.null(h$params)) {
                            rp <- region_params()
                            rp <- set_region_param(rp, module, h$region_id, "hap_scan", h$params)
                            save_region_params(project_data()$name, rp)
                            region_params(rp)
                        }
                    }
                    # key for hap_scan is stored with __scan suffix
                    base_key <- sub("__scan$", "", key)
                    hap_scan_cache[[base_key]] <- result

                } else if (h$type == "hap_viz") {
                    result <- tryCatch(
                        .collect_hap_viz_result(h),
                        error = function(e) list(error = conditionMessage(e))
                    )
                    if (!is.null(result$error) && status != 0L) {
                        log_lines <- tryCatch(readLines(h$log_file, warn = FALSE), error = function(e) character())
                        result$log_tail <- paste(utils::tail(log_lines, 10), collapse = "\n")
                    }
                    # Persist viz params to region_params.json
                    if (is.null(result$error) && !is.null(h$params)) {
                        rp <- region_params()
                        rp <- set_region_param(rp, module, h$region_id, "hap_viz", h$params)
                        save_region_params(project_data()$name, rp)
                        region_params(rp)
                    }
                    base_key <- sub("__viz$", "", key)
                    hap_viz_cache[[base_key]] <- result

                } else if (h$type == "vcf_count") {
                    # Read awk output from temp file (file-redirect pattern,
                    # mirrors other subprocess handlers — no pipe-buffer uncertainty)
                    lines <- tryCatch(
                        readLines(h$output_file, warn = FALSE),
                        error = function(e) character(0)
                    )
                    unlink(h$output_file)  # clean up temp file
                    count <- suppressWarnings(as.integer(trimws(lines[1])))
                    vcf_snp_count_cache[[h$region_id]] <-
                        if (length(lines) > 0 && !is.na(count)) count else NA_integer_
                }

                subprocess_handles[[key]] <- NULL
            }
        })

        # ── Detail panel (conditional) ─────────────────────────────────────────
        output$region_detail_ui <- shiny::renderUI({
            rid <- selected_region_id()
            if (is.null(rid)) return(NULL)
            htmltools::tagList(
                shiny::uiOutput(ns("detail_info_bar")),
                bslib::accordion(
                    id       = ns("detail_sections"),
                    open     = c("genes", "enrichment", "regionplot", "haplotype"),
                    multiple = TRUE,

                    bslib::accordion_panel(
                        "Genes in Region",
                        value = "genes",
                        icon  = bsicons::bs_icon("list-ul"),
                        DT::DTOutput(ns("genes_table"))
                    ),

                    bslib::accordion_panel(
                        "SNPs in Region",
                        value = "snps",
                        icon  = bsicons::bs_icon("dot"),
                        DT::DTOutput(ns("snps_table"))
                    ),

                    bslib::accordion_panel(
                        "GO Enrichment",
                        value = "enrichment",
                        icon  = bsicons::bs_icon("bar-chart"),
                        shiny::uiOutput(ns("enrichment_ui"))
                    ),

                    bslib::accordion_panel(
                        "Regional Manhattan",
                        value = "regionplot",
                        icon  = bsicons::bs_icon("graph-up"),
                        shiny::uiOutput(ns("regionplot_ui"))
                    ),

                    if (!is.null(.module_to_hap_source(module)))
                        bslib::accordion_panel(
                            "Haplotype Analysis",
                            value = "haplotype",
                            icon  = bsicons::bs_icon("diagram-3"),
                            shiny::uiOutput(ns("haplotype_ui"))
                        )
                )
            )
        })

        # ── Detail info bar ────────────────────────────────────────────────────
        output$detail_info_bar <- shiny::renderUI({
            row <- current_region_row()
            if (is.null(row)) return(NULL)
            genes <- region_genes()
            n_g <- if (nrow(genes) > 0) nrow(genes) else NULL
            region_info_bar(row$region_id[1], n_snps = row$snp_count[1], n_genes = n_g)
        })

        # ── Genes table ────────────────────────────────────────────────────────
        output$genes_table <- DT::renderDataTable({
            shiny::req(selected_region_id())
            genes <- region_genes()
            if (nrow(genes) == 0) {
                return(safe_datatable(
                    data.frame(Message = "No genes found (check GFF config)"),
                    options = list(dom = "t"), rownames = FALSE))
            }
            exclude_cols <- c("attributes", "source", "feature", "score", "strand", "phase")
            show_cols <- setdiff(names(genes), exclude_cols)
            safe_datatable(
                as.data.frame(genes[, show_cols, with = FALSE]),
                extensions = "Buttons",
                options    = list(dom = "Bfrtip", buttons = list("csv"),
                                  scrollX = TRUE, pageLength = 10),
                rownames = FALSE
            )
        })

        # ── SNPs table ─────────────────────────────────────────────────────────
        output$snps_table <- DT::renderDataTable({
            shiny::req(selected_region_id())
            snps <- region_snps()
            if (nrow(snps) == 0) {
                return(safe_datatable(data.frame(Message = "No SNPs in region"),
                    options = list(dom = "t"), rownames = FALSE))
            }
            keep <- intersect(c("SNPID", "chr", "pos", "pvalue", "method", "trait"), names(snps))
            safe_datatable(
                as.data.frame(snps[, keep, with = FALSE]),
                extensions = "Buttons",
                options    = list(dom = "Bfrtip", buttons = list("csv"),
                                  scrollX = TRUE, pageLength = 15),
                rownames = FALSE
            )
        })

        # ── Enrichment UI ──────────────────────────────────────────────────────
        output$enrichment_ui <- shiny::renderUI({
            rid <- selected_region_id()
            if (is.null(rid)) return(NULL)
            key    <- enrich_key()
            cached <- if (!is.null(key)) enrichment_cache[[key]] else NULL
            running <- is_running(key)

            if (running) {
                return(htmltools::div(
                    class = "pipeline-rule-item d-flex align-items-center gap-2 py-1",
                    bsicons::bs_icon("arrow-repeat",
                        class = "text-warning spin-icon flex-shrink-0", size = "0.85em"),
                    htmltools::span(class = "small", "Running GO enrichment...")
                ))
            }

            if (!is.null(cached)) {
                if (!is.null(cached$error)) {
                    return(htmltools::div(
                        class = "text-muted small py-2",
                        bsicons::bs_icon("exclamation-triangle"), " ", cached$error,
                        htmltools::br(),
                        shiny::actionButton(ns("run_enrichment"),
                            "Re-run Enrichment",
                            class = "btn-sm btn-outline-secondary mt-2",
                            icon  = shiny::icon("rotate-right"))
                    ))
                }
                return(htmltools::tagList(
                    bslib::layout_column_wrap(
                        width = 1 / 2,
                        if (!is.null(cached$dotplot_path) && file_ok(cached$dotplot_path))
                            mod_image_card_ui(ns("enrich_dotplot"))
                        else
                            bslib::card(bslib::card_body(
                                plot_placeholder("Dotplot not available"))),
                        if (!is.null(cached$emapplot_path) && file_ok(cached$emapplot_path))
                            mod_image_card_ui(ns("enrich_emapplot"))
                        else
                            bslib::card(bslib::card_body(
                                plot_placeholder("Emapplot not available (\u22652 GO terms needed)")))
                    ),
                    DT::DTOutput(ns("enrich_table"))
                ))
            }

            # Show warnings that block enrichment (config error or no annotations)
            genes_df <- region_genes()
            pd_enr   <- project_data()
            go_field <- config_get(pd_enr$config, "GFF", "go_field", default = NULL)
            go_status <- if (nrow(genes_df) == 0) {
                NULL
            } else if (is.null(go_field) || !go_field %in% names(genes_df)) {
                htmltools::div(
                    class = "alert alert-warning small d-flex gap-2 align-items-start mb-2 py-2",
                    bsicons::bs_icon("exclamation-triangle-fill", class = "flex-shrink-0 mt-1"),
                    htmltools::div(
                        "No GO annotation column found. Check ",
                        htmltools::tags$code("GFF.go_field"), " in config."
                    )
                )
            } else {
                n_total <- nrow(genes_df)
                n_go    <- sum(!is.na(genes_df[[go_field]]) & nzchar(genes_df[[go_field]]))
                if (n_go == 0) {
                    htmltools::div(
                        class = "alert alert-warning small d-flex gap-2 align-items-start mb-2 py-2",
                        bsicons::bs_icon("exclamation-triangle-fill", class = "flex-shrink-0 mt-1"),
                        htmltools::div(
                            "No GO annotations found in ", n_total, " region gene",
                            if (n_total != 1) "s", ". ",
                            "Check ", htmltools::tags$code("GFF.go_field"), " config."
                        )
                    )
                } else {
                    NULL  # annotations present — no need to say anything before running
                }
            }

            htmltools::tagList(
                go_status,
                shiny::actionButton(ns("run_enrichment"),
                    "Run Enrichment",
                    class = "btn-sm btn-outline-primary",
                    icon  = shiny::icon("play"))
            )
        })

        # Image card servers (re-initialized when new enrichment result cached)
        shiny::observe({
            key    <- enrich_key()
            cached <- if (!is.null(key)) enrichment_cache[[key]] else NULL
            if (is.null(cached) || !is.null(cached$error)) return()
            mod_image_card_server("enrich_dotplot",
                path    = shiny::reactive(cached$dotplot_path),
                title   = shiny::reactive("GO Dotplot"),
                dl_name = shiny::reactive("enrichment_dotplot"),
                note    = shiny::reactive(help_note("enrich_dotplot"))
            )
            mod_image_card_server("enrich_emapplot",
                path    = shiny::reactive(cached$emapplot_path),
                title   = shiny::reactive("GO Emapplot"),
                dl_name = shiny::reactive("enrichment_emapplot"),
                note    = shiny::reactive(help_note("enrich_emapplot"))
            )
        })

        output$enrich_table <- DT::renderDataTable({
            key    <- enrich_key()
            cached <- if (!is.null(key)) enrichment_cache[[key]] else NULL
            if (is.null(cached) || !is.null(cached$error) || is.null(cached$table))
                return(safe_datatable(data.frame(), options = list(dom = "t"), rownames = FALSE))
            dt_format_pvals(safe_datatable(
                as.data.frame(cached$table),
                extensions = "Buttons",
                options    = list(dom = "Bfrtip", buttons = list("csv"),
                                  scrollX = TRUE, pageLength = 10),
                rownames = FALSE
            ))
        })

        shiny::observeEvent(input$run_enrichment, {
            key <- enrich_key()
            if (is.null(key) || is_running(key)) return()
            row   <- current_region_row()
            if (is.null(row)) return()
            tr    <- current_trait()
            if (is.null(tr)) return()

            handle <- launch_enrichment_subprocess(
                region_genes(), row$region_id[1], tr, project_data(), module
            )
            if (!is.null(handle$error)) {
                enrichment_cache[[key]] <- handle   # store error immediately
            } else {
                subprocess_handles[[key]] <- handle # start polling
            }
        })

        # ── Regionplot UI ──────────────────────────────────────────────────────
        output$regionplot_ui <- shiny::renderUI({
            rid    <- selected_region_id()
            if (is.null(rid)) return(NULL)
            key    <- rplot_key()
            cached <- if (!is.null(key)) regionplot_cache[[key]] else NULL
            running <- is_running(key)

            if (running) {
                return(htmltools::div(
                    class = "pipeline-rule-item d-flex align-items-center gap-2 py-1",
                    bsicons::bs_icon("arrow-repeat",
                        class = "text-warning spin-icon flex-shrink-0", size = "0.85em"),
                    htmltools::span(class = "small", "Generating regional Manhattan...")
                ))
            }

            if (!is.null(cached)) {
                if (!is.null(cached$error)) {
                    return(htmltools::div(
                        class = "text-muted small py-2",
                        bsicons::bs_icon("exclamation-triangle"), " ", cached$error,
                        htmltools::br(),
                        shiny::actionButton(ns("run_regionplot"),
                            "Retry Region Plot",
                            class = "btn-sm btn-outline-secondary mt-2",
                            icon  = shiny::icon("rotate-right"))
                    ))
                }
                path <- cached$path
                if (!is.null(path) && file_ok(path)) {
                    return(mod_image_card_ui(ns("regionplot_img")))
                }
            }

            shiny::actionButton(ns("run_regionplot"),
                "Generate Region Plot",
                class = "btn-sm btn-outline-primary",
                icon  = shiny::icon("chart-line"))
        })

        shiny::observe({
            key    <- rplot_key()
            cached <- if (!is.null(key)) regionplot_cache[[key]] else NULL
            if (is.null(cached) || !is.null(cached$error)) return()
            path <- cached$path
            if (is.null(path) || !file_ok(path)) return()
            mod_image_card_server("regionplot_img",
                path    = shiny::reactive(path),
                title   = shiny::reactive("Regional Manhattan"),
                dl_name = shiny::reactive("regionplot"),
                note    = shiny::reactive(help_note("regionplot_img"))
            )
        })

        shiny::observeEvent(input$run_regionplot, {
            key <- rplot_key()
            if (is.null(key) || is_running(key)) return()
            row  <- current_region_row()
            if (is.null(row)) return()
            snps <- region_snps()
            if (is.null(snps) || nrow(snps) == 0) return()

            handle <- launch_regionplot_subprocess(row, snps, project_data(), module,
                                                     regime = if (isTRUE(regime())) "wza" else "snp")
            if (!is.null(handle$error)) {
                regionplot_cache[[key]] <- handle   # store error immediately
            } else {
                subprocess_handles[[key]] <- handle # start polling
            }
        })

        # ── Haplotype UI ───────────────────────────────────────────────────────
        output$haplotype_ui <- shiny::renderUI({
            rid <- selected_region_id()
            if (is.null(rid)) return(NULL)
            key      <- hap_key()
            tag      <- expected_hap_tag()
            scan_res <- if (!is.null(key)) hap_scan_cache[[key]] else NULL
            viz_res  <- if (!is.null(key)) hap_viz_cache[[key]]  else NULL
            scan_run <- is_running(paste0(key %||% "", "__scan"))
            viz_run  <- is_running(paste0(key %||% "", "__viz"))

            # Config defaults for parameter inputs
            pd <- project_data()
            cfg_eps   <- config_get(pd$config, "haplotype", "epsilon_selected", default = NULL)
            cfg_range <- paste(config_get(pd$config, "haplotype", "scan",
                                          "epsilon_range", default = "0.3,0.5,0.7,0.9"),
                               collapse = ",")
            if (is.list(cfg_range)) cfg_range <- paste(unlist(cfg_range), collapse = ",")
            cfg_mgmin  <- config_get(pd$config, "haplotype", "scan", "min_group_size",  default = 50L)
            cfg_minhap <- config_get(pd$config, "haplotype", "scan", "min_haplotype_size", default = 15L)
            cfg_minsnp <- config_get(pd$config, "haplotype", "scan", "min_snps", default = 3L)
            cfg_meta   <- config_get(pd$config, "haplotype", "scan", "metadata_type", default = "site")

            # Saved params from previous computations for this region (source of truth for defaults + badges)
            saved_scan <- get_region_param(region_params(), module, rid, "hap_scan")
            saved_viz  <- get_region_param(region_params(), module, rid, "hap_viz")
            # Use saved params as defaults; fall back to config defaults for new regions
            use_range  <- saved_scan$epsilon_range %||% cfg_range
            use_mgmin  <- saved_scan$mgmin         %||% cfg_mgmin
            use_minhap <- saved_scan$minhap        %||% cfg_minhap
            use_minsnp <- saved_scan$min_snps      %||% cfg_minsnp
            use_meta   <- saved_scan$meta_type     %||% cfg_meta

            # ── Sub-box 1: Scan ──────────────────────────────────────────────
            # Shared input form — used in initial, error, and re-run states
            # Modified badges compare against saved params (what produced current output)
            .scan_modified <- function(field, current) {
                ref <- if (!is.null(saved_scan)) saved_scan[[field]] else NULL
                if (is.null(ref)) return(NULL)  # no saved run yet — no badge
                changed <- tryCatch(as.character(current) != as.character(ref), error = function(e) FALSE)
                if (isTRUE(changed))
                    htmltools::span(class = "badge bg-warning text-dark small ms-1", "modified")
            }
            # H4: SNP sufficiency warning — counts ALL SNPs in the region
            # from the filtered VCF (haplotype uses all variants, not just sig).
            n_region_snps <- region_total_snps()
            min_snp_thr   <- as.integer(use_minsnp)
            snp_warning <- if (is.na(n_region_snps)) {
                htmltools::div(
                    class = "alert alert-info small d-flex gap-2 align-items-start mb-2 py-2",
                    bsicons::bs_icon("info-circle-fill", class = "flex-shrink-0 mt-1"),
                    htmltools::div(
                        "Total SNP count in region unavailable ",
                        "(filtered VCF not found — run mode=processing)."
                    )
                )
            } else if (n_region_snps < min_snp_thr) {
                htmltools::div(
                    class = "alert alert-warning small d-flex gap-2 align-items-start mb-2 py-2",
                    bsicons::bs_icon("exclamation-triangle-fill", class = "flex-shrink-0 mt-1"),
                    htmltools::div(
                        "Region has ", htmltools::tags$strong(n_region_snps), " total SNPs",
                        " (minimum: ", min_snp_thr, "). ",
                        "Haplotype scan may fail. Lower ",
                        htmltools::tags$code("Min SNPs in region"), " or choose a different region."
                    )
                )
            } else if (n_region_snps < 2L * min_snp_thr) {
                htmltools::div(
                    class = "alert alert-info small d-flex gap-2 align-items-start mb-2 py-2",
                    bsicons::bs_icon("info-circle-fill", class = "flex-shrink-0 mt-1"),
                    htmltools::div(
                        "Region has ", htmltools::tags$strong(n_region_snps), " total SNPs",
                        " (minimum: ", min_snp_thr, "). Results may be limited."
                    )
                )
            } else NULL

            scan_inputs <- htmltools::tagList(
                snp_warning,
                bslib::layout_column_wrap(
                    width = 1/2,
                    htmltools::div(
                        shiny::textInput(ns("hap_epsilon_range"), "Epsilon values",
                            value = use_range, placeholder = "e.g. 0.3,0.5,0.7,0.9"),
                        .scan_modified("epsilon_range", use_range)
                    ),
                    htmltools::div(
                        shiny::numericInput(ns("hap_mgmin"), "Min group size (MGmin)",
                            value = as.integer(use_mgmin), min = 1L, step = 1L),
                        .scan_modified("mgmin", as.integer(use_mgmin))
                    )
                ),
                bslib::layout_column_wrap(
                    width = 1/2,
                    htmltools::div(
                        shiny::numericInput(ns("hap_minhap"), "Min haplotype size",
                            value = as.integer(use_minhap), min = 1L, step = 1L),
                        .scan_modified("minhap", as.integer(use_minhap))
                    ),
                    htmltools::div(
                        shiny::numericInput(ns("hap_min_snps"), "Min SNPs in region",
                            value = as.integer(use_minsnp), min = 1L, step = 1L),
                        .scan_modified("min_snps", as.integer(use_minsnp))
                    )
                ),
                htmltools::div(
                    shiny::selectInput(ns("hap_meta_type"), "Metadata grouping",
                        choices  = c("site", paste0("cluster_K", pd$k_best)),
                        selected = user_meta_type() %||% use_meta, width = "220px"),
                    .scan_modified("meta_type", user_meta_type() %||% use_meta)
                ),
                shiny::actionButton(ns("reset_scan_defaults"), "Reset to defaults",
                    class = "btn-sm btn-link text-muted px-0",
                    icon  = shiny::icon("arrow-rotate-left")),
                htmltools::tags$details(
                    htmltools::tags$summary(
                        class = "text-muted small mb-1",
                        bsicons::bs_icon("info-circle"), " How to choose epsilon"
                    ),
                    htmltools::p(class = "text-muted small ms-3",
                        "Run the scan to explore haplotype clustering across epsilon values. ",
                        "Review the clustree plot \u2014 look for stable clusters where lines ",
                        "don\u2019t cross. Fewer crossings = more stable partitioning. ",
                        "Then select your chosen epsilon in the Visualization box below."
                    )
                )
            )

            scan_content <- if (scan_run) {
                htmltools::div(
                    class = "pipeline-rule-item d-flex align-items-center gap-2 py-1",
                    bsicons::bs_icon("arrow-repeat",
                        class = "text-warning spin-icon flex-shrink-0", size = "0.85em"),
                    htmltools::span(class = "small", "Running haplotype scan...")
                )
            } else if (!is.null(scan_res) && isTRUE(scan_res$skip_clustree)) {
                # Scan succeeded but only 1 epsilon worked — clustree not possible
                eps_hint <- scan_res$epsilon_hint
                htmltools::tagList(
                    htmltools::div(
                        class = "alert alert-info small mb-3",
                        bsicons::bs_icon("info-circle"), " ",
                        "Scan complete \u2014 only 1 epsilon produced marker groups ",
                        "(clustree requires \u22652). ",
                        if (!is.null(eps_hint))
                            htmltools::tagList(
                                "Proceed to Step\u00a02 with epsilon\u00a0=\u00a0",
                                htmltools::strong(eps_hint), ".")
                        else
                            "Proceed to Step 2 Visualization."
                    ),
                    scan_inputs,
                    shiny::actionButton(ns("run_hap_scan"), "Re-run Scan",
                        class = "btn-sm btn-outline-secondary mt-1",
                        icon  = shiny::icon("rotate-right"))
                )
            } else if (!is.null(scan_res) && is.null(scan_res$error)) {
                htmltools::tagList(
                    shiny::uiOutput(ns("hap_clustree_ui")),
                    htmltools::hr(),
                    htmltools::p(class = "text-muted small fw-semibold mb-1",
                        bsicons::bs_icon("arrow-repeat"), " Re-run scan with different parameters"),
                    scan_inputs,
                    shiny::actionButton(ns("run_hap_scan"), "Re-run Scan",
                        class = "btn-sm btn-outline-secondary mt-1",
                        icon  = shiny::icon("rotate-right"))
                )
            } else if (!is.null(scan_res) && !is.null(scan_res$error)) {
                htmltools::tagList(
                    htmltools::div(
                        class = "alert alert-warning small mb-3",
                        bsicons::bs_icon("exclamation-triangle"), " ", scan_res$error,
                        if (!is.null(scan_res$log_tail) && nzchar(scan_res$log_tail))
                            htmltools::tagList(
                                htmltools::br(),
                                htmltools::tags$details(
                                    htmltools::tags$summary(
                                        class = "text-muted",
                                        style = "cursor:pointer",
                                        "Show log"
                                    ),
                                    htmltools::tags$pre(
                                        class = "text-muted mb-0 mt-1",
                                        style = "white-space:pre-wrap;font-size:0.72em;max-height:120px;overflow-y:auto",
                                        scan_res$log_tail
                                    )
                                )
                            )
                    ),
                    scan_inputs,
                    shiny::actionButton(ns("run_hap_scan"), "Retry Scan",
                        class = "btn-sm btn-outline-secondary",
                        icon  = shiny::icon("rotate-right"))
                )
            } else {
                htmltools::tagList(
                    scan_inputs,
                    shiny::actionButton(ns("run_hap_scan"), "Run Scan",
                        class = "btn-sm btn-outline-primary",
                        icon  = shiny::icon("magnifying-glass"))
                )
            }

            # ── Sub-box 2: Visualization ────────────────────────────────────
            # Always show epsilon form + Run button; show plots below when available
            viz_status <- if (viz_run) {
                htmltools::div(
                    class = "pipeline-rule-item d-flex align-items-center gap-2 py-1 mb-2",
                    bsicons::bs_icon("arrow-repeat",
                        class = "text-warning spin-icon flex-shrink-0", size = "0.85em"),
                    htmltools::span(class = "small", "Running haplotype visualization...")
                )
            } else if (!is.null(viz_res) && !is.null(viz_res$error)) {
                htmltools::div(
                    class = "alert alert-warning small mb-2",
                    bsicons::bs_icon("exclamation-triangle"), " ", viz_res$error,
                    if (!is.null(viz_res$log_tail) && nzchar(viz_res$log_tail))
                        htmltools::tagList(
                            htmltools::br(),
                            htmltools::tags$details(
                                htmltools::tags$summary(
                                    class = "text-muted",
                                    style = "cursor:pointer",
                                    "Show log"
                                ),
                                htmltools::tags$pre(
                                    class = "text-muted mb-0 mt-1",
                                    style = "white-space:pre-wrap;font-size:0.72em;max-height:120px;overflow-y:auto",
                                    viz_res$log_tail
                                )
                            )
                        )
                )
            } else NULL

            # Epsilon dropdown: choices from successful scan epsilons; fall back to saved/config default
            scan_epsilons <- if (!is.null(scan_res) && !is.null(scan_res$epsilons) && length(scan_res$epsilons) > 0)
                                 sort(unique(scan_res$epsilons))
                             else numeric(0)
            eps_stored   <- user_epsilon()
            eps_default  <- if (!is.null(cfg_eps) && !is.na(suppressWarnings(as.numeric(cfg_eps))))
                                as.numeric(cfg_eps) else NA_real_
            # Prefer saved viz epsilon, then user-stored, then config default
            eps_saved    <- if (!is.null(saved_viz$epsilon_selected))
                                as.numeric(saved_viz$epsilon_selected) else NA_real_
            eps_choices  <- if (length(scan_epsilons) > 0) scan_epsilons else
                            if (!is.na(eps_saved)) eps_saved else
                            if (!is.na(eps_default)) eps_default else numeric(0)
            eps_value    <- eps_stored %||% (if (!is.na(eps_saved)) eps_saved else NULL) %||% eps_default
            # Select: prefer stored value if among choices, else saved, else first
            eps_selected <- if (!is.null(eps_value) && eps_value %in% eps_choices) eps_value else
                            if (!is.na(eps_saved) && eps_saved %in% eps_choices) eps_saved else
                            if (length(eps_choices) > 0) eps_choices[1] else NULL
            # "modified": current selection differs from saved viz params
            eps_modified <- !is.null(saved_viz$epsilon_selected) && !is.null(eps_selected) &&
                            !is.na(eps_selected) &&
                            as.numeric(eps_selected) != as.numeric(saved_viz$epsilon_selected)

            viz_content <- htmltools::tagList(
                htmltools::div(
                    class = "d-flex align-items-center gap-2 flex-wrap",
                    if (length(eps_choices) > 0)
                        shiny::selectInput(ns("hap_epsilon_selected"), "Epsilon (from scan)",
                            choices = eps_choices, selected = eps_selected, width = "160px")
                    else
                        shiny::numericInput(ns("hap_epsilon_selected"), "Epsilon",
                            value = eps_value %||% NA_real_, min = 0.01, step = 0.1, width = "160px"),
                    if (isTRUE(eps_modified))
                        htmltools::span(class = "badge bg-warning text-dark small mt-3", "modified")
                ),
                if (length(scan_epsilons) > 0)
                    htmltools::p(class = "text-muted small",
                        bsicons::bs_icon("arrow-up"), " Select epsilon from successful scan runs above.")
                else
                    htmltools::p(class = "text-muted small",
                        bsicons::bs_icon("arrow-up"), " Run the scan first to get epsilon choices."),
                shiny::actionButton(ns("run_hap_viz"),
                    if (viz_run) "Running..." else "Run Visualization",
                    class = paste("btn-sm", if (viz_run) "btn-outline-secondary disabled" else "btn-outline-primary"),
                    icon  = shiny::icon("chart-bar")),
                viz_status,
                # Viz plots always present — servers use reactive paths (placeholder when empty)
                if (!is.null(viz_res) && is.null(viz_res$error)) {
                    htmltools::tagList(
                        htmltools::hr(),
                        mod_image_card_ui(ns("hap_crosshap_viz")),
                        bslib::layout_column_wrap(
                            width = 1/2,
                            mod_image_card_ui(ns("hap_boxplot")),
                            mod_image_card_ui(ns("hap_piemap"))
                        ),
                        htmltools::hr(),
                        bslib::accordion(
                            open     = FALSE,
                            multiple = TRUE,
                            bslib::accordion_panel(
                                "Haplotype Assignments",
                                value = "hap_assignments",
                                icon  = bsicons::bs_icon("table"),
                                DT::DTOutput(ns("hap_tbl_assignments"))
                            ),
                            bslib::accordion_panel(
                                "Haplotype Frequencies",
                                value = "hap_frequencies",
                                icon  = bsicons::bs_icon("bar-chart"),
                                DT::DTOutput(ns("hap_tbl_frequencies"))
                            )
                        )
                    )
                }
            )

            # Collect all available traits from scan + viz results
            all_traits <- unique(c(
                if (!is.null(scan_res) && !is.null(scan_res$traits)) scan_res$traits else character(0),
                if (!is.null(viz_res)  && !is.null(viz_res$traits))  viz_res$traits  else character(0)
            ))
            trait_selector <- if (length(all_traits) > 1) {
                htmltools::div(
                    class = "mb-3",
                    shiny::selectInput(ns("hap_trait"), "Trait",
                        choices = all_traits, width = "220px")
                )
            } else if (length(all_traits) == 1) {
                # Single trait — no selector needed, but set the value
                shiny::tagList(shiny::tags$script(
                    shiny::HTML(sprintf(
                        "$(document).ready(function(){ Shiny.setInputValue('%s', '%s'); });",
                        ns("hap_trait"), all_traits[1]
                    ))
                ))
            } else NULL

            htmltools::tagList(
                trait_selector,
                bslib::layout_column_wrap(
                    width = 1/2,
                    bslib::card(
                        bslib::card_header(
                            bsicons::bs_icon("search"), " Step 1: Scan"
                        ),
                        bslib::card_body(scan_content)
                    ),
                    bslib::card(
                        bslib::card_header(
                            bsicons::bs_icon("palette"), " Step 2: Visualization"
                        ),
                        bslib::card_body(viz_content)
                    )
                )
            )
        })

        # Clustree images UI container (rendered inside scan success state)
        output$hap_clustree_ui <- shiny::renderUI({
            bslib::layout_column_wrap(
                width = 1/2,
                mod_image_card_ui(ns("hap_clustree_mg")),
                mod_image_card_ui(ns("hap_clustree_hap"))
            )
        })

        # Image servers for clustree plots (re-initialized when scan cache or trait updates)
        shiny::observe({
            key      <- hap_key()
            scan_res <- if (!is.null(key)) hap_scan_cache[[key]] else NULL
            if (is.null(scan_res) || !is.null(scan_res$error) || isTRUE(scan_res$skip_clustree)) return()

            # Pick per-trait or default paths
            trait <- input$hap_trait
            use_trait <- !is.null(trait) && length(scan_res$traits) > 0 &&
                         trait %in% scan_res$traits
            if (use_trait) {
                rid   <- scan_res$region_id
                pdir  <- scan_res$plots_dir
                mg_p  <- file.path(pdir, paste0("Region_", rid, "_clustree_MG_",  trait, ".png"))
                hap_p <- file.path(pdir, paste0("Region_", rid, "_clustree_hap_", trait, ".png"))
                lbl   <- paste0(" (", trait, ")")
            } else {
                mg_p  <- scan_res$mg_path
                hap_p <- scan_res$hap_path
                lbl   <- ""
            }

            if (!is.null(mg_p) && file_ok(mg_p))
                mod_image_card_server("hap_clustree_mg",
                    path    = shiny::reactive(mg_p),
                    title   = shiny::reactive(paste0("Clustree \u2014 Marker Groups", lbl)),
                    dl_name = shiny::reactive(paste0("clustree_MG", if (use_trait) paste0("_", trait) else "")),
                    note    = shiny::reactive(help_note("hap_clustree_mg")))
            if (!is.null(hap_p) && file_ok(hap_p))
                mod_image_card_server("hap_clustree_hap",
                    path    = shiny::reactive(hap_p),
                    title   = shiny::reactive(paste0("Clustree \u2014 Haplotypes", lbl)),
                    dl_name = shiny::reactive(paste0("clustree_hap", if (use_trait) paste0("_", trait) else "")),
                    note    = shiny::reactive(help_note("hap_clustree_hap")))
        })

        # Reactive viz image paths — update when trait/region/viz results change
        hap_viz_trait_r <- shiny::reactive({
            key     <- hap_key()
            viz_res <- if (!is.null(key)) hap_viz_cache[[key]] else NULL
            if (is.null(viz_res) || !is.null(viz_res$error)) return(NULL)
            input$hap_trait %||% viz_res$traits[1]
        })

        hap_viz_path_r <- function(path_fn) shiny::reactive({
            trait   <- hap_viz_trait_r()
            rid     <- selected_region_id()  # always exact — no fuzzy matching
            tag     <- expected_hap_tag()
            pd      <- project_data()
            if (is.null(trait) || is.null(rid) || is.null(tag) || is.null(pd)) return(NULL)
            path_fn(pd$name, tag, rid, trait)
        })

        mod_image_card_server("hap_crosshap_viz",
            path    = hap_viz_path_r(crosshap_viz_path),
            title   = shiny::reactive({
                t <- hap_viz_trait_r()
                if (!is.null(t)) paste0("Haplotype Visualization (", t, ")") else "Haplotype Visualization"
            }),
            dl_name = shiny::reactive(paste0("crosshap_", selected_region_id() %||% "region")),
            note    = shiny::reactive(help_note("hap_crosshap_viz")))

        mod_image_card_server("hap_boxplot",
            path    = hap_viz_path_r(hap_boxplot_path),
            title   = shiny::reactive({
                t <- hap_viz_trait_r()
                if (!is.null(t)) paste0("Boxplot (", t, ")") else "Haplotype Phenotype Boxplot"
            }),
            dl_name = shiny::reactive(paste0("hap_boxplot_", selected_region_id() %||% "region")),
            note    = shiny::reactive(help_note("hap_boxplot")))

        mod_image_card_server("hap_piemap",
            path    = hap_viz_path_r(hap_piemap_path),
            title   = shiny::reactive({
                t <- hap_viz_trait_r()
                if (!is.null(t)) paste0("Piemap (", t, ")") else "Haplotype Piemap"
            }),
            dl_name = shiny::reactive(paste0("hap_piemap_", selected_region_id() %||% "region")),
            note    = shiny::reactive(help_note("hap_piemap")))

        # ── Haplotype data tables (Assignments + Frequencies) ──────────────────
        hap_table_rid_r <- shiny::reactive({
            selected_region_id()  # always exact — no fuzzy matching
        })

        output$hap_tbl_assignments <- DT::renderDataTable({
            rid <- shiny::req(hap_table_rid_r())
            tag <- shiny::req(expected_hap_tag())
            pd  <- project_data()
            p   <- hap_assignments_path(pd$name, tag)
            dt  <- if (file_ok(p)) data.table::fread(p, colClasses = c(region_id = "character", sample = "character")) else data.table::data.table()
            if (nrow(dt) > 0 && "region_id" %in% names(dt)) dt <- dt[dt$region_id == rid, ]
            safe_datatable(as.data.frame(dt),
                extensions = "Buttons",
                options    = list(dom = "Bfrtip", buttons = list("csv"),
                                  scrollX = TRUE, pageLength = 20),
                rownames   = FALSE)
        })

        output$hap_tbl_frequencies <- DT::renderDataTable({
            rid <- shiny::req(hap_table_rid_r())
            tag <- shiny::req(expected_hap_tag())
            pd  <- project_data()
            p   <- hap_frequencies_path(pd$name, tag)
            dt  <- if (file_ok(p)) data.table::fread(p, colClasses = c(region_id = "character")) else data.table::data.table()
            if (nrow(dt) > 0 && "region_id" %in% names(dt)) dt <- dt[dt$region_id == rid, ]
            safe_datatable(as.data.frame(dt),
                extensions = "Buttons",
                options    = list(dom = "Bfrtip", buttons = list("csv"),
                                  scrollX = TRUE, pageLength = 20),
                rownames   = FALSE)
        })

        # ── Run Scan button ────────────────────────────────────────────────────
        shiny::observeEvent(input$run_hap_scan, {
            row <- current_region_row()
            if (is.null(row)) return()
            src <- .module_to_hap_source(module)
            if (is.null(src)) return()
            meta <- input$hap_meta_type %||% user_meta_type() %||% "site"
            tag  <- paste0(meta, "_", src)
            # Use a key consistent with expected_hap_tag (which also constructs from config)
            key  <- hap_key()
            if (is.null(key) || is_running(paste0(key, "__scan"))) return()

            params <- list(
                epsilon_range = input$hap_epsilon_range %||% "0.3,0.5,0.7,0.9",
                mgmin         = input$hap_mgmin  %||% 50L,
                minhap        = input$hap_minhap %||% 15L,
                min_snps      = input$hap_min_snps %||% 3L,
                meta_type     = meta
            )
            handle <- launch_hap_scan_subprocess(row, project_data(), tag, params)
            if (!is.null(handle$error)) {
                hap_scan_cache[[key]] <- handle
            } else {
                subprocess_handles[[paste0(key, "__scan")]] <- handle
            }
        })

        # ── Run Viz button ─────────────────────────────────────────────────────
        shiny::observeEvent(input$run_hap_viz, {
            row <- current_region_row()
            if (is.null(row)) return()
            src <- .module_to_hap_source(module)
            if (is.null(src)) return()
            meta <- input$hap_meta_type %||% user_meta_type() %||% "site"
            tag  <- paste0(meta, "_", src)
            key  <- hap_key()
            if (is.null(key) || is_running(paste0(key, "__viz"))) return()

            eps <- input$hap_epsilon_selected
            params <- list(
                epsilon_selected = eps,
                meta_type        = meta
            )
            handle <- launch_hap_viz_subprocess(row, project_data(), tag, params)
            if (!is.null(handle$error)) {
                hap_viz_cache[[key]] <- handle
            } else {
                subprocess_handles[[paste0(key, "__viz")]] <- handle
            }
        })

        # ── Reset scan params to config defaults ───────────────────────────────
        shiny::observeEvent(input$reset_scan_defaults, {
            pd <- project_data()
            shiny::updateTextInput(session, "hap_epsilon_range",
                value = paste(config_get(pd$config, "haplotype", "scan",
                                          "epsilon_range", default = "0.3,0.5,0.7,0.9"),
                              collapse = ","))
            shiny::updateNumericInput(session, "hap_mgmin",
                value = as.integer(config_get(pd$config, "haplotype", "scan",
                                               "min_group_size", default = 50L)))
            shiny::updateNumericInput(session, "hap_minhap",
                value = as.integer(config_get(pd$config, "haplotype", "scan",
                                               "min_haplotype_size", default = 15L)))
            shiny::updateNumericInput(session, "hap_min_snps",
                value = as.integer(config_get(pd$config, "haplotype", "scan",
                                               "min_snps", default = 3L)))
            shiny::updateSelectInput(session, "hap_meta_type",
                selected = config_get(pd$config, "haplotype", "scan",
                                       "metadata_type", default = "site"))
        })

        # ── Return value for parent module ─────────────────────────────────────
        list(
            computed_regions   = computed_regions,
            selected_region_id = selected_region_id
        )
    })
}

# ── Helper ────────────────────────────────────────────────────────────────────

#' Format a bp distance as Mb/kb/bp string
#' @noRd
.format_bp <- function(bp) {
    bp <- as.numeric(bp)
    ifelse(bp >= 1e6, paste0(round(bp / 1e6, 1), " Mb"),
    ifelse(bp >= 1e3, paste0(round(bp / 1e3,  1), " kb"),
                      paste0(bp, " bp")))
}
