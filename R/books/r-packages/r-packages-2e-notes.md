# R Packages, 2nd edition

https://r-pkgs.org/

## 13 Testing basics

- Testing workflow:
    - One time in the package: run `usethis::use_testthat(3)` to use edition 3 of `testthat`
    - Create a test with `usethis::use_test("suffix-for-test-file")`
    - Write `test_that("name", {code})` calls in the test file
    - Run tests
        - Develop test code interactively with `devtools::load_all(); test_that("test-name", {code to test})` interactively
        - Run one test file with `testthat::test_file("tests/tesththat/test-suffix.R")` or `devtools::test_active_file()`
        - Run all tests with `devtools::test()`
        - Finally, run `devtools::check()`
- Each test may have multiple expectations, but one test should cover a single unit of functionality
- Snapshot tests
    - New to `testthat` 3e
    - Uses `waldo` package
    - Only function when run non-interactively
    - Advantages: moves the static test data values into a readable md format that can be updated more easily than if in quoted strings etc.
    - `testthat::snapshot_review('diff')` lets you review the changes and accept or skip in a local shiny app
    - if you want to use an error output in a snapshot you have to use `error=TRUE`, because errors are not allowed in snapshots by default
    - default is `cran=FALSE`, which means the test won't run on CRAN's servers; but this is due to the assumption that other non-snapshot tests are doing the checks for correctness
    - `variant` arg lets you specify different snapshots for different conditions (e.g. OS, dependency version, etc.)    

## Designing by your test suite

### What to test
- Test the external interface, not the internal implementation
- Test each thing once and only once
- Avoid testing simple, obviously correct code
- Always write a test when you discover a bug
- You can use `devtools::test_coverage_active_file()` and `devtools::test_coverage()` to generate coverage data

### High-level principles for testing
- Self-sufficiency
    - Your `tests/testthat/test-*.R` files should mostly contain calls to `test_that()` and minimal other top-level code
    - This may mean repeating code in multiple `test_that()` calls
- Self-containment
    - Watch out for things that affect the filesystem, the search path, or global options like `options()`, `par()`, and env vars
    - If you **have to** make such changes in the course of a test, you need to do clean up after
    - Try using `withr::local_*()` or `withr::with_*()` functions, or just use `withr::defer()` (the nicer equivalent of `on.exit()`)
    - If you use `withr` in tests, put it in `Suggests`
    - In `testthat` 3e, `testthat::local_reproducible_output()` is used implicitly in every `test_that()` call, which "suppresses colored output, turns off fancy quotes, sets the console width, and sets LC_COLLATE = "C""
- Plan on test failure
    - The `devtools` setup lets you run one test at a time to fix things; this is one of the best ways to debug a test interactively
    - This is one of the big reasons for self-sufficient/self-contained test code
- Repetition is OK
    - It's more important to have obvious/understandable test code than DRY
    - See this blog post: https://mtlynch.io/good-developers-bad-tests/
    - If you really need to make something that's big and cumbersome to repeat, use a test fixture
- Remove tension between interactive and automated testing
    - Don't set up anything for interactive convenience that will compromise the health of automated tests
    - You WILL NOT NEED any `library()` calls in the test files!!!
    - Don't use `source()` in tests! (probably)
- Many tidyverse packages might violate these principles because `testthat` is evolving and it's been a long time

### Files relevant to testing
- Test helpers (e.g. factories) can be defined in `R/some_file.R` and used in tests
- Don't edit anything in `tests/testthat.R`; it's run in R CMD check, but not used in most other testing or interactively
- `tests/testthat/helper*.R>` gets executed by `devtools::load_all()`
    - It's your choice whether you put test helper functions here or in the R dir
    - This is a "good location for setup code that is needed for its side effects"
    - In an API wrapper, this is also a good place for authenticating with test creds
- Setup files in `tests/testthat/setup*.R`
    - Treated like a helper, EXCEPT not executed by `devtools::load_all()`
    - Often contain corresponding tear down code
- Any other files in the `tests/testthat` dir won't be automatically executed
- Storing test data
    - Suggestion here is to haver `tests/testthat/fixtures` and store the code to create the data and the data itself in this dir
    - In your tests, access these files with `readRDS(testthat::test_path("fixtures", "data-object.rds"))`
- Only ever write files to the session temp dir during tests! (and still make sure to clean up after)
    - `withr::local_tempfile()` is the best way to do this
        - The temp file's lifetime is tied to the "local" environment
    - If you need a specific file name, use `withr::local_tempdir()` first as part of the path containing the specific file


## Advanced testing techniques

### Test fixtures
- For when it's not practical to make your test entirely self-sufficient and -contained
1. Create `useful_thing`s with a helper function
    - If it's fiddly or takes several lines of code to create a thing
    - This is basically a factory like I used in my Django tests
    - Convention: `new_useful_thing <- function(){}`
2. Create (and destroy) a "local" `useful_thing`
    - When creating the useful thing has side effects (on file system, remote resource, R options, env vars, etc.)
    - Need to make sure the helper cleans up afterwards
    - Write a custom function in the style of `withr`
    - Convention: `local_useful_thing <- function(..., env=parent.frame()){}`
3. Store a concrete `useful_thing` persistently
    - For when it's costly to create (time, memory, or it doesn't need to be new for each test)
    - Create the thing once, store it as a static test fixture, and load it in the tests that need it

### Building your own testing tools
- Helper defined inside a test
- Custom expectations
    - a nice example from `usethis`:

```{r}
expect_proj_file <- function(...) {
  expect_true(file_exists(proj_path(...)))
}
```

- `testthat` creates the env var `TESTTHAT` and sets it to true when testing
    - You can use `testthat::is_testing()` as a shortcut for checking this env var
- The package-under-test is available in `TESTTHAT_PKG` and the shortcut `testthat::testing_package()`

### When testing gets hard
1. Skipping a test
    - Use `testthat::skip()` to write a custom skipper function
    - e.g. `skip_if_no_api()`
    - Several built in skipper functions
        - `skip_if()`
        - `skip_if_not_installed(pkg)`
        - `skip_if_offline()`
        - `skip_on_cran()`
        - `skip_on_os("windows")`
    - Be careful about skipping too much    
2. Mocking an external service
    - "mocking is when we replace something that’s complicated or unreliable or out of our control with something simpler, that’s fully within our control"
    - Usually an external service or a function that reports about session state
    - "Our main advice about mocking is to avoid it if you can"
    - Mostly if the API requires authentication that's hard to do in tests, or if it often has downtime that will break your tests
3. Dealing with secrets
    - Main advice is to try to design the package so that most of it can be tested without the creds
    - Then test the authentication parts in an env that provides secure env vars (e.g. GitHub Actions)

### Special considerations for CRAN packages
- You can skip if needed
- Tests need to run quickly (prefer < 1 min)
- Beware of reproducibility on CRAN systems
- "There is basically no latitude for a test that's 'just flaky'"
    - Test external service connections elsewhere and skip on cran
- Snapshot tests are skipped on CRAN!
- On CRAN you can only write in the R session's temp dir
- Do not start external software (PDF viewers, browsers, etc.)
- Turn off clipboard functionality in tests on CRAN

## References mentioned a lot
- Martin Fowler
- [Software Engineering at Google](https://abseil.io/resources/swe-book)
