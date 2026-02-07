#' penguin_table UI
#' 
#' @param id Unique id for module instance.
#' 
#' @keywords internal
penguin_tableUI <- function(id){
	ns <- NS(id)

	tagList(
		h2("What does the penguins data look like?"),
    reactable::reactableOutput(ns("penguins_reactable"))
	)
}

#' penguin_table Server
#' 
#' @param id Unique id for module instance.
#' 
#' @keywords internal
penguin_table_server <- function(id, penguins){
	moduleServer(
		id,
		function(
			input, 
			output, 
			session
			){
				
				ns <- session$ns
				send_message <- make_send_message(session)

				output$penguins_reactable <- reactable::renderReactable(reactable::reactable(penguins()))
		}
	)
}

# UI
# penguin_tableUI('id')

# server
# penguin_table_server('id')
