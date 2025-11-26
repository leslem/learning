# Shiny UI testing

## Appsilon comparison of Cypress and `shinytest2`
https://www.appsilon.com/post/shinytest2-vs-cypress-e2e-testing
- shinytest2
    - snapshot tests are useful when the ui doesn't change often
    - doesn't emulate typing, but instead uses "the mechanics of shiny" e.g. sets inputs directly, doesn't emulate someone clicking or typing into inputs
    - this can bypass custom JS or other parts of client side code
    - you can only use CSS selectors
    - you can only interact with clicks (no scroll or text input)
    - much of the snapshot based test can change with a code update (tests your implementation, not your functionality)
- cypress
    - a node.js package; not specifically made for R/shiny
    - more stable than Selenium
    - tests must be written in javascript

## Jumping Rivers series on e2e testing with `shinytest2`
https://www.jumpingrivers.com/blog/end-to-end-testing-shinytest2-part-1/
https://www.jumpingrivers.com/blog/end-to-end-testing-shinytest2-part-2/

- automated interactive testing frameworks
- usually include webdriver + assertion library
- Existing options by language:
    - JS: cypress, puppeteer, playwright
    - many: selenium
    - R: shinytest2 (and shinytest, but that's deprecated now)
- Main caveats for e2e tests:
    - Very slow
    - Harder to write
    - More fragile
    - Flaky b/c of external dependencies
- Best shinytest2 tutorial I've seen; not focused on using the test recorder
- Shows how to make an R6 object for testing your specific app

## Appsilon: how to write tests with shiny::testServer similar to shinytest2
https://www.appsilon.com/post/how-to-write-tests-with-shiny-testserver
- shiny::testServer() lets you test reactivity without running the entire app
- use it with testthat

```
shiny::testServer(server, {
    session$setInputs(x=1)
    expect_equal(myreactive(), 1)
    expect_equal(output$txt, "x is 1")
})
```
- you can use shiny::testServer for testing modules too

## Jakub Sobolewski: Choosing the Best Library for Acceptance Testing Shiny Apps: {shinytest2} vs Cypress
https://jakubsobolewski.com/blog/which-testing-library-for-shiny/
- shinytest2
    - easier to set up and start, everything written in R
    - `shiny::exportTestValues()` is great
    - dependent on implementation details
    - snapshot tests are encouraged too much
    - hard to see what's going on to write tests
    - very slow (needs to init Chrome session and R session)
    - the inputs you set with `app$set_inputs()` are the server-facing inputs, not user-facing inputs
- cypress
    - complete separation of tests from app code; can't access Shiny server
    - better selectors
    - need to write tests in JS
    - Set up can be difficult
    
## playwright

- Iain mentioned this as something that shiny for python devs used at posit::conf
- the documentation looks nice
- the getting started info is more approachable than cypress
- not trying to sell you a bunch of services like cypress
- I followed this example project to try it out on a toy shiny app:
    - https://github.com/kamilsi/shiny-playwright
    - I got it working
    - Just had to install first (following slide instructions)
    - Had to add `renv::load()` to the config webServer argument to launch the app via Rscript
    - Had to increase some timeouts even further
- I already have nvm, node, npm installed in my home dir at ~/.nvm and ~/.npm from the data harmonization agent setup
    - I installed nvm via macports and then used nvm to install node v20
    
    sudo port install nvm
    sudo nvm install v20
    echo "source /opt/local/share/nvm/init-nvm.sh" >> ~/.zshrc
    
    - I also have nodejs22 and npm10 installed via macports, for no reason I can discern    
    - More info on node, nvm, npm: https://stackoverflow.com/questions/32660993/difference-between-npm-and-nvm
    - This guide that seems very obviously LLM-generated is actually distressingly good: https://www.w3tutorials.net/blog/macports-nodejs/
    - My options for node version management are:
        1. install nvm via Macports and then use nvm to install whichever versions I want
            - this is what I already did for the harmonization agent
        2. install nvm outside of Macports
        3. individually install each nodejs version via macports
    - I'll uninstall the npm10 and nodejs22 ports and then proceed with option 1
       - I'll install node v22 for setting up playwright testing
    

```
# to install and set up playwright for the toy shiny app
npm init -y
sudo npm install --save-dev @playwright/test
sudo npx playwright install
```

`npm init -y`

> npm init -y is a non‑interactive shortcut for initializing a new Node.js project. Normally npm init launches a prompt asking for fields like package name, version, entry point, author, license, and scripts. The -y (or --yes) flag automatically accepts all default values and immediately writes a package.json file using npm’s built‑in defaults, speeding up setup when you don’t need custom metadata yet or are scripting project creation.
> The generated package.json becomes the manifest for the project: it records dependencies, scripts (like test or build commands), and metadata used by tooling, CI, and publishing. After running this command you can manually edit package.json to adjust fields (e.g. change "main", add "type": "module", insert "scripts" like "start" or "test"). Using -y is helpful in automation (scaffolding tools, demos, throwaway prototypes) because it avoids user prompts, but for production packages it is wise to review and refine the defaults to ensure clarity and correctness (accurate license, semantic versioning, repository links).

`npm install --save-dev @playwright/test`

> npm install --save-dev @playwright/test installs the Playwright test runner package and records it under the devDependencies section of your package.json. Using --save-dev signals that this library is only needed during development and CI (for writing and running tests) and is not required in production runtime builds. The command downloads the package (and its transitive dependencies) into node_modules and updates package-lock.json to lock exact versions, ensuring reproducible test environments.
> The specific package @playwright/test provides the test runner, assertion API (expect), fixtures, parallelization, reporters, and auto‑waiting logic tailored for end‑to‑end and component style browser tests. After installation you typically (a) add a script like "test": "playwright test" to package.json, and (b) run npx playwright install (or npx playwright install chromium firefox webkit) to fetch the actual browser binaries Playwright controls. Keeping it as a dev dependency keeps production deploys lighter while still enabling full-featured local and CI browser testing.

- playwright's default is to use TypeScript, "a superset of JavaScript that uses static typing"
- https://playwright.dev/docs/best-practices
- A customized test class using a "fixture": https://playwright.dev/docs/test-fixtures

## Playwright guidelines from Houseful

https://www.houseful.blog/posts/2023/playwright-standards/
