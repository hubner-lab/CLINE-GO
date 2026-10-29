# gapit.R's copy of each model's native GAPIT files into GAPIT_native_output/<model>/,
# quick tier (finding fd5085). The copy selects files by a name pattern; GAPIT itself
# needs a genotype panel, so the pattern is taken from the script's own list.files() call
# and run here on the names GAPIT 3.5 writes (GAPIT.<kind>.<model>.<trait>.<ext>, as in
# SIMDATA's GWAS/GAPIT_native_output/BLINK/). The unanchored "GAPIT.*<model>" sent CMLM's,
# ECMLM's and MLMM's files into MLM's folder, and ECMLM's into CMLM's.

# The `pattern` argument of gapit.R's `native_files <- list.files(...)`, unevaluated.
native_pattern_expr <- function() {
    src <- readLines(file.path(.clinego_root, "scripts", "gapit.R"))
    i <- grep("native_files <- list.files(", src, fixed = TRUE)
    expect_length(i, 1L)
    for (j in i:min(length(src), i + 5L)) {
        call <- tryCatch(parse(text = paste(src[i:j], collapse = "\n"))[[1]],
                         error = function(e) NULL)
        if (!is.null(call)) return(call[[3]]$pattern)
    }
    stop("gapit.R: the native_files call does not parse")
}

# Every GAPIT model of the method registry (workflow/methods/gea.py).
GAPIT_MODELS <- c("BLINK", "FarmCPU", "MLM", "MLMM", "GLM", "CMLM", "ECMLM", "SUPER")

gapit_file_names <- function(models, traits = c("height", "bio_1")) {
    kinds <- c("Association.GWAS_Results.%s.%s.csv", "Association.GWAS_StdErr.%s.%s.csv",
               "Association.Manhattan_Chro.%s.%s.pdf", "Association.Manhattan_Geno.%s.%s.pdf",
               "Association.QQ.%s.%s.pdf", "Phenotype.Distribution_Significantmarkers.%s.%s.pdf")
    grid <- expand.grid(kind = kinds, model = models, trait = traits, stringsAsFactors = FALSE)
    paste0("GAPIT.", sprintf(grid$kind, grid$model, grid$trait))
}

test_that("each model's native folder receives that model's GAPIT files and no other's", {
    pattern_expr <- native_pattern_expr()
    files <- c(gapit_file_names(GAPIT_MODELS),
               "GAPIT.Genotype.PCA_eigenvalues.pdf", "GAPIT.Association.Filter_GWAS_results.csv")
    for (model in GAPIT_MODELS) {
        got <- grep(eval(pattern_expr, list(model = model)), files, value = TRUE)
        expect_setequal(got, gapit_file_names(model))
    }
})
