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
usethis::use_package("palmerpenguins")
renv::snapshot(type="all", prompt=FALSE)


# Create a module
leprechaun::add_module()