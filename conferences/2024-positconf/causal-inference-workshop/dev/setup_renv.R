source("~/devel/dotfiles/set_proxy.R")
source("~/devel/dotfiles/unset_proxy.R")

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

install.packages("pak")
pak::pak("r-causal/causalworkshop")
causalworkshop::install_workshop("~/devel/learning/conferences/2024-positconf/causal-inference-workshop")
