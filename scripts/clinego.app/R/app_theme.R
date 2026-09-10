#' Build the Precision Genomics bslib theme
#' @noRd
app_theme <- function() {
    bslib::bs_theme(
        version = 5,
        preset  = "shiny",

        # Colors — deep teal-green scientific palette
        bg        = "#FAFBFC",
        fg        = "#2D3748",
        primary   = "#1B7A6E",
        secondary = "#4A5568",
        success   = "#38A169",
        info      = "#3182CE",
        warning   = "#D69E2E",
        danger    = "#E53E3E",

        # Typography
        base_font    = bslib::font_google("Source Sans 3"),
        heading_font = bslib::font_google("Inter", wght = c(500, 600, 700)),
        code_font    = bslib::font_google("JetBrains Mono"),
        font_scale   = 0.95,

        # Component variables
        "navbar-bg"          = "#1A2332",
        "navbar-dark-color"  = "#E2E8F0",
        "card-border-radius" = "0.5rem",
        "card-border-color"  = "#E2E8F0",
        "border-radius"      = "0.375rem",
        "accordion-border-color" = "#E2E8F0"
    ) |>
        bslib::bs_add_rules(sass::sass_file(app_sys("app/www/custom.scss"))) |>
        # Dashboard layout LAST: it refines what custom.scss establishes, and
        # this order is the one every layout measurement was taken against.
        bslib::bs_add_rules(sass::sass_file(app_sys("app/www/dashboard.scss")))
}
