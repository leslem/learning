library(shiny)

ui <- fluidPage(
	'Hello World',
	sliderInput(inputId='num', label='Choose a number',
				value=25, min=1, max=100),
	plotOutput('hist')
)

server <- function(input, output) {
	output$hist <- renderPlot({
		hist(rnorm(input$num), main='Histogram')
	})
	
}

shinyApp(ui = ui, server = server)
