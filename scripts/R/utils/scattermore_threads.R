# scattermore_threads.R — run scattermore's point kernel on the job's own threads.
# Usage: source("/pipeline/scripts/R/utils/scattermore_threads.R") in every script that draws
# geom_scattermore(). Requires: scattermore installed (it need not be attached).
#
# scattermore 1.2's scattermore(), which geom_scattermore() draws through, blurs the points with
# apply_kernel_rgbwt(threads = 0) whenever pointsize is not 0, and 0 means one thread per host
# core (std::thread::hardware_concurrency()). No environment variable reaches it: not
# OMP_NUM_THREADS, not data.table's cap, not docker --cpus. On the 192-core lab every Manhattan
# save ran 194 threads inside a 1-thread job, and 48 containers drove the load to ~1,900
# (de7512). The kernel's default is set here to the job's thread count, which Snakemake exports
# as OMP_NUM_THREADS (1 when unset). Each thread blurs its own band of rows, so the raster is the
# same for any count.
local({
    threads <- suppressWarnings(as.integer(Sys.getenv("OMP_NUM_THREADS", "1")))
    if (is.na(threads) || threads < 1L) threads <- 1L
    kernel <- get("apply_kernel_rgbwt", envir = asNamespace("scattermore"))
    if (!"threads" %in% names(formals(kernel))) {
        stop("scattermore_threads.R: scattermore ", as.character(utils::packageVersion("scattermore")),
             "'s apply_kernel_rgbwt() has no `threads` argument; find how this version sizes its ",
             "thread pool before changing this file")
    }
    formals(kernel)$threads <- threads
    utils::assignInNamespace("apply_kernel_rgbwt", kernel, ns = "scattermore")
})
