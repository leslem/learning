## Following the [Quick Start guide](https://leprechaun.opifex.org/#/guide/quick-start)

# already created the project directory and set up an renv
renv::install(c("usethis", "shiny", "bslib", "devtools", "attachment", "leprechaun"), prompt=F)
renv::snapshot(type="all", prompt=FALSE)

usethis::create_package(".", rstudio=FALSE)
# Complained about my project inside a project, but I went ahead anyway

leprechaun::scaffold(ui="navbarPage", bs_version = 5)

devtools::document()
devtools::load_all()
run()

# Add some stuff to the app
renv::install("palmerpenguins", prompt=F)
# Add the data to the server
usethis::use_package("palmerpenguins")
renv::snapshot(type="all", prompt=FALSE)

# I ran into trouble here because I hadn't renamed my package, so it got loaded as "leprechaun"
# I updated the DESCRIPTION file and then reloaded my R session
# I later ran into a problem because I hadn't updated the package name in the zzz.R file
# Create a module
leprechaun::add_module(name="penguin_table")
renv::install("reactable", prompt=F)
renv::snapshot(type="all", prompt=FALSE)

