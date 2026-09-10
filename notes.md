# CLINE-GO — session notes (USE-only mode; I run the app/pipeline, do not modify code)

## Current state
- [x] WBDC_final GEA mode finished correctly (51/51 steps, outputs valid, all K5).
- [!] Shiny app showed no GEA figures → cause: stale app-level `bindCache` held `k_best=6` while outputs are K5; every GEA path resolved to non-existent `_K6_` files. Browser reload can't clear it (cache is process-level). **Fix = restart container.**
- [x] FIXED (user authorized code change): Per-Method Details crash. See [BUG_per_method_manhattan_region_id.md](docs/BUG_per_method_manhattan_region_id.md). Four edits, all parse-OK:
  - `fct_manhattan.R` — `region_id` now optional in `build_manhattan_plotly` (conditional hover line + `customdata = NULL` when absent). Root-cause fix, protects every caller.
  - `mod_manhattan_overlay.R` — `renderPlotly` wrapped in `tryCatch` → fail-soft placeholder instead of killing the app.
  - `app_server.R` — `project_data` bindCache now `cache="session"` (browser reload picks up external config/output changes; fixes stale k_best).
  - `utils_helpers.R` — `resolve_k_best()` emits a `warning()` when falling back to max(K) instead of failing silently.
- [ ] VERIFY in browser: GEA/GWAS → Per-Method Details renders (no "Region:" line, no crash) on WBDC_final (K5); Combined/region detail/QQ unchanged. User tests manually.
- [x] WBDC_final GWAS: 2 WZA Manhattan rules failed (assoc_wza_manhattan Tiller_Number/BLINK + _combined) — WZA window p=0 → -log10=Inf → scattermore crash. Diagnosed: see [BUG_wza_manhattan_scattermore_inf.md](docs/BUG_wza_manhattan_scattermore_inf.md).
- [x] FIXED (user authorized code change). 3 files: `compute_wza.R` (pnorm lower.tail=FALSE + `pmax(p,.Machine$double.xmin)` floor → no exact-0), `manhattan_utils.R::add_scatter_layer` (drop non-finite log10p before scattermore), `plot_manhattan_combined.R` (df_background filter is.finite). Toy-validated: both plot scripts now survive p=0+NA windows (exit 0).
- [x] WBDC_final gwas re-run DONE: `-R assoc_wza` 35/35 steps, exit 0. Both failed plots now exist; regen WZA tables have no 0/NA/neg. The underflowed window now reads true p=8.48e-18 (-log10=17.07) — `1-pnorm` was just cancellation, `lower.tail=FALSE` recovers it; floor never triggered. NOTE: snakemake needs `-e USER=cline-go` (getpass fails on uid 1020).
