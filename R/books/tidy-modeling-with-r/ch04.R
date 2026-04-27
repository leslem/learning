library(dplyr)
library(ggplot2)

data(ames, package = "modeldata")

glimpse(ames)
ggplot(ames, aes(x = Sale_Price)) +
    geom_histogram(bins = 50, col = "white") +
    scale_x_log10()

ames <- ames |>
    mutate(Sale_Price = log10(Sale_Price))
