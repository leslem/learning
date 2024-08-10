source("~/devel/set_proxy.R")

renv::settings$snapshot.type("all")
renv::init(bare=TRUE)
my_repos <- c(
    'CRAN'="https://cran.rstudio.com/"
)
options(repos=my_repos)
options()$repos


# Install packages from CRAN
renv::install(c(
    "dplyr",
    "ggplot2",
    "glue",
    "knitr",
    "markdown",
    "rmarkdown"
), prompt=FALSE)

renv::snapshot(type="all", prompt=FALSE)

renv::restore(prompt=FALSE)
