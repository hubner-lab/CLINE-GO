#' Per-cell significance rule helpers.
#'
#' A "rule" is list(type=, value=). The MASTER rule is the threshold bar's
#' scalar (type, value). An OVERRIDE map is "trait::method" -> rule for the
#' subset of cells that deviate from the fallback. A cell absent from the map
#' falls back to its method's REGISTRY rule when that method's
#' significance_family differs from "univariate_pvalue" (today: RDA only —
#' rdadapt emits one joint test per SNP, so the univariate master threshold
#' has no standing there, gea.py:125-135), else it follows master.
#'
#' Row and column popups do NOT introduce a separate precedence tier — they
#' are bulk writers that expand into N individual "trait::method" entries
#' (one per method in a trait's row, or one per trait in a method's column).
#' Precedence is therefore exactly: cell override > registry-or-master.
#' @noRd

#' Deterministic, order-independent cache key for an override map.
#'
#' MUST be sorted: an unsorted list key produces cache misses that look like a
#' perf regression (bindCache in mod_gea.R, digest() in compute_method_sigsnps_cached).
#' Folds in registry_defaults so editing workflow/methods/gea.py's
#' adjust_default/threshold_default (a pipeline-code change, not a Shiny one)
#' still invalidates caches built under the old registry rule.
#' @noRd
threshold_overrides_key <- function(overrides, registry_defaults = list()) {
    ov_part <- if (length(overrides) == 0) "" else {
        ov <- overrides[order(names(overrides))]
        paste(vapply(names(ov), function(k) sprintf(
            "%s=%s:%s", k, ov[[k]]$type,
            format(as.numeric(ov[[k]]$value), digits = 15, scientific = TRUE)),
            character(1)), collapse = "|")
    }
    reg_part <- if (length(registry_defaults) == 0) "" else {
        rd <- registry_defaults[order(names(registry_defaults))]
        paste(vapply(names(rd), function(m) sprintf(
            "%s=%s:%s:%s", m, rd[[m]]$adjust %||% "", rd[[m]]$threshold %||% "",
            rd[[m]]$family %||% ""), character(1)), collapse = "|")
    }
    paste(ov_part, reg_part, sep = "||")
}

#' Resolve the effective rule for one (method, trait) cell.
#'
#' @param method   Method name (e.g. "EMMAX")
#' @param trait    Trait name, or NULL for a method-only lookup (no cell-level
#'   override can apply — used only by the save_snp_set provenance summary,
#'   which reports one rule per method, not per cell)
#' @param overrides Named list "trait::method" -> list(type=,value=)
#' @param master_type  Master threshold type
#' @param master_value Master threshold value
#' @param registry_defaults Named list method -> list(adjust=,threshold=,family=)
#'   from gea_method_significance_defaults(). A method whose family differs
#'   from "univariate_pvalue" is pinned to its registry rule unless a cell
#'   override exists for it.
#' @return list(type=, value=, source=) — source one of "override", "registry", "master"
#' @noRd
effective_rule_for <- function(method, trait = NULL, overrides = list(),
                               master_type, master_value,
                               registry_defaults = list()) {
    if (!is.null(trait)) {
        r <- overrides[[paste0(trait, "::", method)]]
        if (!is.null(r)) {
            return(list(type = r$type %||% master_type,
                       value = as.numeric(r$value %||% master_value),
                       source = "override"))
        }
    }
    rd <- registry_defaults[[method]]
    if (!is.null(rd) && !identical(rd$family %||% "univariate_pvalue", "univariate_pvalue")) {
        return(list(type = rd$adjust, value = as.numeric(rd$threshold), source = "registry"))
    }
    list(type = master_type, value = as.numeric(master_value), source = "master")
}

#' Normalise an override map read back from region_params.json.
#'
#' read_json(simplifyVector = FALSE) yields nested lists; drop malformed
#' entries, coerce value to numeric, and migrate legacy bare-method keys
#' (pre-rework: METHOD -> rule, applied to every trait) into per-cell
#' "trait::method" keys for every CURRENT trait — preserves the override
#' instead of silently losing it on first load after the rework. A key whose
#' trait or method no longer exists in the project is dropped rather than
#' piling up as unreachable dead weight in region_params.json.
#'
#' @param raw     Parsed JSON list (or NULL)
#' @param traits  Character vector of current trait names (migration + pruning)
#' @param methods Character vector of current method names (migration + pruning)
#' @return Named list "trait::method" -> list(type=, value=). list() for
#'   absent/empty/all-invalid input.
#' @noRd
normalize_threshold_overrides <- function(raw, traits = character(0), methods = character(0)) {
    if (is.null(raw) || length(raw) == 0 || is.null(names(raw))) return(list())
    out <- list()
    for (k in names(raw)) {
        e  <- raw[[k]]
        ty <- tryCatch(as.character(e$type), error = function(...) NA_character_)
        if (length(ty) != 1) ty <- NA_character_
        # e$value on an entry missing "value" returns NULL (no error) — coerce
        # NULL/absent to a scalar NA_real_ before is.na(), or as.numeric(NULL)
        # (a zero-length vector) breaks the length-1 `if (... || is.na(va))`
        # condition below instead of taking the "skip this entry" branch.
        va_raw <- tryCatch(e$value, error = function(...) NULL)
        va <- if (is.null(va_raw)) NA_real_ else suppressWarnings(as.numeric(va_raw))
        if (length(va) != 1) va <- NA_real_
        if (is.na(ty) || !nzchar(ty) || is.na(va)) next
        if (!threshold_value_valid_for_type(ty, va)) next   # fct_combine.R

        if (grepl("::", k, fixed = TRUE)) {
            parts <- strsplit(k, "::", fixed = TRUE)[[1]]
            if (length(parts) != 2) next
            tr <- parts[1]; m <- parts[2]
            if (length(traits) > 0 && !(tr %in% traits)) next
            if (length(methods) > 0 && !(m %in% methods)) next
            out[[k]] <- list(type = ty, value = va)
        } else {
            # Legacy bare-method key — expand to every current trait so the
            # rule still applies everywhere it used to.
            if (length(methods) > 0 && !(k %in% methods)) next
            for (tr in traits) out[[paste0(tr, "::", k)]] <- list(type = ty, value = va)
        }
    }
    out
}
