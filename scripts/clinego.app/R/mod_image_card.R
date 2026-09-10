#' Generic image card module UI
#'
#' Displays a static PNG/SVG image in a card with full-screen expansion
#' and SVG/PNG download buttons.
#'
#' @param id module namespace id
#' @param height card body height (default "auto")
#' @noRd
mod_image_card_ui <- function(id, height = "auto") {
    ns <- shiny::NS(id)
    bslib::card(
        full_screen = TRUE,
        bslib::card_header(
            class = "d-flex justify-content-between align-items-center",
            htmltools::div(
                class = "d-flex align-items-center gap-2",
                shiny::textOutput(ns("title"), inline = TRUE),
                shiny::uiOutput(ns("note"), inline = TRUE)
            ),
            bslib::popover(
                trigger = bsicons::bs_icon("download", title = "Download plot"),
                title = "Download",
                shiny::downloadButton(ns("dl_svg"), "SVG",
                                      class = "btn-sm btn-outline-secondary"),
                shiny::downloadButton(ns("dl_png"), "PNG",
                                      class = "btn-sm btn-outline-secondary")
            )
        ),
        bslib::card_body(
            fillable = FALSE,
            class = "p-2 text-center",
            shiny::uiOutput(ns("image_or_placeholder"))
        )
    )
}

#' Generic image card module server
#'
#' @param id module namespace id
#' @param path reactive returning file path (PNG or SVG); NULL = placeholder
#' @param title reactive string for card header title
#' @param dl_name reactive string used as download filename base (no extension)
#' @param placeholder reactive string shown when path is unavailable
#' @param suggestion reactive string shown as smaller hint text below placeholder message
#' @param note reactive returning a small htmltools tag (e.g. a colored badge with a hover
#'   tooltip) shown next to the title, or NULL to show nothing. Used for compact data-derived
#'   status (e.g. relatedness pair counts) instead of a static caption baked into the plot.
#' @noRd
mod_image_card_server <- function(id, path, title = shiny::reactive("Plot"),
                                   dl_name = shiny::reactive("plot"),
                                   placeholder = shiny::reactive("Plot not available"),
                                   suggestion  = shiny::reactive(NULL),
                                   note        = shiny::reactive(NULL)) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        output$title <- shiny::renderText(title())
        output$note  <- shiny::renderUI(note())

        output$image_or_placeholder <- shiny::renderUI({
            p <- path()
            if (file_ok(p)) {
                shiny::imageOutput(ns("img"), height = "auto",
                                   width = "100%", inline = FALSE)
            } else {
                plot_placeholder(placeholder(), suggestion())
            }
        })

        output$img <- shiny::renderImage({
            p <- shiny::req(path())
            shiny::validate(shiny::need(file_ok(p), "File not found"))
            ext <- tolower(tools::file_ext(p))
            mime <- switch(ext, svg = "image/svg+xml", "image/png")
            list(src = p, contentType = mime, width = "100%", alt = title())
        }, deleteFile = FALSE)

        # Download handlers — try SVG first, fall back to PNG
        output$dl_svg <- shiny::downloadHandler(
            filename = function() paste0(dl_name(), ".svg"),
            content  = function(file) {
                p <- path()
                svg_path <- sub("\\.png$", ".svg", p, ignore.case = TRUE)
                src <- if (file_ok(svg_path)) svg_path else p
                if (file_ok(src)) file.copy(src, file) else writeLines("", file)
            }
        )

        output$dl_png <- shiny::downloadHandler(
            filename = function() paste0(dl_name(), ".png"),
            content  = function(file) {
                p <- path()
                png_path <- sub("\\.svg$", ".png", p, ignore.case = TRUE)
                src <- if (file_ok(png_path)) png_path else p
                if (file_ok(src)) file.copy(src, file) else writeLines("", file)
            }
        )
    })
}
