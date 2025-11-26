# Software testing

## Research on mocking (or not)

- Two Scoops of Django:
    - "Use Mock to keep unit tests from touching the world"
    - external APIs, emails, webhooks, etc.
    - A mock will "fake a response from an external API"
    - "Don't rely on fixtures" i.e. static test data in json format
    - Instead use a "factory" to produce fake data for testing
- Test Driven Development in Python
    - "Using mocks to test external dependencies"
    - example is emails
    - "Monkey patching" is mocking manually
    - e.g. replace `send_mail()` with a fake version at runtime; `fake_send_mail()` will return the appropriate parts of the email, e.g. subject line, from email, and to email
    - "Tests using mocks end up tightly coupled to the implementation"
    - "It's better to test behavior, not implementation details"
    - Reasons to use mocks:
        1. Isolate from external side effects (e.g. creating file, sending emails, making API calls)
        2. To reduce code duplication across tests
    - Types of tests
        - Functional (acceptance test, or end-to-end test): track a user story; from the POV of a user; from the outside in
        - Unit: from the POV of a programmer; from the inside out
        - Functional tests drive what code you develop; unit tests drive how you implement the code
        - Integration test: checks that the code you control is integrated correctly with some external system that you don't control
- [Martin Fowler "Mocks aren't Stubs"](https://martinfowler.com/articles/mocksArentStubs.html)
    - This replaces the style of testing where you set up, test, and then tear down ("state verification")
    - Another style of mocking mocks objects that are internal to the codebase ("behavior verification")
    - A mock is just one type of a "test double" -- a stand in to use during testing
        - dummy: passed around but not used, just to fill in the param as needed
        - fake: has working implementation, but use a shortcut to make it
        - stub: provide canned answers to a call made during test
        - spy: stub that also records information based on how it was called
        - mock: object pre-programmed with expectations that form a specification of the calls it's expected to receive
    - mocks do behavior verification and the others all do state verification
    - classical TDD: use real objects if possible and a double if the real thing is awkward
    - mockist TDD: always use a mock for any object with interesting behavior
    - "Although the various mock frameworks were designed with mockist testing in mind, many classicists find them useful for creating doubles"

### `testthat` mocking article

https://testthat.r-lib.org/articles/mocking.html

