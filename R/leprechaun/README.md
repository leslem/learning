# Learning the `leprechaun` shiny framework

## Comparison to `golem`
- Does not make `leprechaun` a dependency of the shiny app
- Approach is to make most things optional, not included by default
- `leprechaun` is harder to type and I'm already mildly irritated by it
- This sounds more like the `usethis` approach than `golem` is, but also might put more of the maintenance burden on the app developer (vs. `golem` maintenance which is done by ThinkR)
- Does not put shiny or bslib in the namespace, so you have to do `shiny::h1()` and `bslib::card()`for everything unless you change the package options yourself
- Not many issues = not many users

## General impressions
- I like the default deps (shiny, bslib, htmltools, pkgload)
- Potential gotcha: the template files created have tabs and not spaces

## Discussion from book club
- `leprechaun` is good for lightweight shiny-app-as-a-package without extra dependencies, but is not ready for production apps
- `rhino` is good if you're forced to use R but want to pretend that it's python
    - `tapyr` is their equivalent for Shiny for Python
- `golem` is good if you're okay with an opnionated framework, and it **is** production ready
- we should check out `periscope` as well
