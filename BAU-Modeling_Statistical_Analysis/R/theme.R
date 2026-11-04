bau_colours <- c(
  bag = "#285F8F",
  bau = "#C1644A",
  nmr = "#2F7F78",
  proteomic = "#80669D",
  text = "#202426",
  grid = "#D8DDDF"
)

paper_theme <- function(base_size = 7, base_family = "Helvetica") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(colour = bau_colours[["text"]]),
      axis.title = ggplot2::element_text(colour = bau_colours[["text"]]),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      legend.title = ggplot2::element_blank()
    )
}

save_panel <- function(plot, path, width = 90, height = 65) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(path, ".svg"), plot,
    device = svglite::svglite,
    width = width, height = height, units = "mm", bg = "white"
  )
  ggplot2::ggsave(paste0(path, ".pdf"), plot,
    device = grDevices::cairo_pdf,
    width = width, height = height, units = "mm", bg = "white"
  )
}
