#' penguin_plot UI
#' 
#' @param id Unique id for module instance.
#' 
#' @keywords internal
penguin_plotUI <- function(id){
	ns <- NS(id)

	tagList(
		h2("Plots of penguins data"),
    plotOutput(ns("bill_size_scatterplot")),
    plotOutput(ns("species_barplot"))
	)
}

#' penguin_plot Server
#' 
#' @param id Unique id for module instance.
#' 
#' @keywords internal
#' @import ggplot2
penguin_plot_server <- function(id, penguins){
	moduleServer(
		id,
		function(
			input, 
			output, 
			session
			){
				
				ns <- session$ns
				send_message <- make_send_message(session)

        species_palette <- c("Adelie" = "darkorange",
                             "Chinstrap" = "darkorchid",
                             "Gentoo" = "cyan4")
        output$bill_size_scatterplot <- renderPlot({
          ggplot(penguins) +
            geom_point(aes(x=bill_length_mm, y=bill_depth_mm, color=species),
                      size=3, alpha=0.8) +
            scale_color_manual(values=species_palette) +
            xlab("Bill length (mm)") +
            ylab("Bill depth (mm)") +
            coord_fixed(ratio=1)
        })

        output$species_barplot <- renderPlot({
          ggplot(penguins) +
            geom_bar(aes(x=species, fill=species)) +
            scale_fill_manual(values=species_palette, guide=FALSE)
        })

		}
	)
}

# UI
# penguin_plotUI('id')

# server
# penguin_plot_server('id')
