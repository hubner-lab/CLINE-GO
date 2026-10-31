#' The theme's colours, the one place they are written (I1): app_theme() passes them to
#' bslib, the SCSS reads them as Sass variables ($primary, ...) and R code that draws
#' (plotly lines, ggplot fills) reads them from here.
#' @noRd
APP_COLORS <- c(
    primary   = "#1B7A6E",   # deep teal-green
    secondary = "#4A5568",
    success   = "#38A169",
    info      = "#3182CE",
    warning   = "#D69E2E",
    danger    = "#E53E3E",
    navbar    = "#1A2332",
    navbar_fg = "#E2E8F0"
)

#' Build the Precision Genomics bslib theme
#'
#' Dark mode only (A4 I2): app_ui.R sets data-bs-theme="dark", so Bootstrap draws the
#' page from its dark-mode body colours; `bg`/`fg` below are the light pair Bootstrap
#' needs to derive its colour maps, never shown.
#' @noRd
app_theme <- function() {
    bslib::bs_theme(
        version = 5,
        preset  = "shiny",

        bg        = "#FAFBFC",
        fg        = "#2D3748",
        primary   = APP_COLORS[["primary"]],
        secondary = APP_COLORS[["secondary"]],
        success   = APP_COLORS[["success"]],
        info      = APP_COLORS[["info"]],
        warning   = APP_COLORS[["warning"]],
        danger    = APP_COLORS[["danger"]],

        # Typography
        base_font    = bslib::font_google("Source Sans 3"),
        heading_font = bslib::font_google("Inter", wght = c(500, 600, 700)),
        code_font    = bslib::font_google("JetBrains Mono"),
        font_scale   = 0.95,

        # Component variables
        "navbar-bg"          = APP_COLORS[["navbar"]],
        "navbar-dark-color"  = APP_COLORS[["navbar_fg"]],
        "card-border-radius" = "0.5rem",
        "card-border-color"  = APP_COLORS[["navbar_fg"]],
        "border-radius"      = "0.375rem",
        "accordion-border-color" = APP_COLORS[["navbar_fg"]]
    ) |>
        bslib::bs_add_rules(sass::sass_file(app_sys("app/www/custom.scss"))) |>
        # Dashboard layout LAST: it refines what custom.scss establishes, and
        # this order is the one every layout measurement was taken against.
        bslib::bs_add_rules(sass::sass_file(app_sys("app/www/dashboard.scss")))
}
