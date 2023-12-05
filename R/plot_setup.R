theme_custom <-
  function(base_size = 15, base_family = "Open Sans", base_line_size = base_size/22, 
          base_rect_size = base_size/22) {
    ggplot2::theme_minimal(
      base_size = base_size, base_family = base_family, 
      base_line_size = base_line_size, base_rect_size = base_rect_size
    ) %+replace% 
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey96", linewidth = .8),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 12)),
      axis.title.y = ggplot2::element_text(angle = 90, margin = ggplot2::margin(r = 12)), 
      legend.title = ggplot2::element_text(face = "bold", hjust = 0),
      strip.text = ggplot2::element_text(face = "bold", size = rel(1), 
                                         margin = ggplot2::margin(t = 12, b = 12)),
      #strip.background = ggplot2::element_rect(fill = "grey96", color = "grey96"),
      plot.title = ggplot2::element_text(face = "bold", size = rel(1.3), hjust = 0,
                                         margin = ggplot2::margin(b = 12)),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      panel.spacing = grid::unit(2.5, "lines"),
      plot.margin = ggplot2::margin(rep(.5, 4)),
      complete = TRUE
    )
  }

ggplot2::theme_set(theme_custom())

theme_map <-
  function(base_size = 15, base_family = "Open Sans", base_line_size = base_size/22, 
           base_rect_size = base_size/22) {
    ggplot2::theme_void(
      base_size = base_size, base_family = base_family, 
      base_line_size = base_line_size, base_rect_size = base_rect_size
    ) %+replace% 
      ggplot2::theme(
        legend.position = "top", 
        legend.justification = "left", 
        axis.text = ggplot2::element_text(color = "black", size = rel(.8), 
                                          margin = ggplot2::margin(rep(5, 4))),
        axis.ticks.length = grid::unit(.4, "lines"), 
        axis.ticks = ggplot2::element_line(color = "grey75"),
        plot.tag = ggplot2::element_text(face = "bold"),
        plot.margin = ggplot2::margin(rep(5, 4)),
        complete = TRUE
      )
  }


## color palette for Berlin and Brandenburg regions
colors_regions <- c('Berlin' = "#e7b800", 'Brandenburg' = "#2e6c61")

## color palette for effects plots showing trends per season
colors_seasons <- c('spring' = "#29a089", 'S_autumn' = "#b26b2f", 'winter' = "#886ed1")
