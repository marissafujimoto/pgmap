# CONTRIBUTING

This document describes how to contribute to this R package.

## GitHub Workflow

- File an issue with anything you see that needs to be considered, added or fixed.
- Assign yourself to an issue that you'd like to work on. If you have any questions about how to go about the issue you've been assigned to, communicate your questions on that issue!
- Clone this repository locally and create a new branch to work on this issue from.
- When you feel you have enough to discuss, push the branch and open a pull request. When you do this, you will see a series of checks begin to run:

## Automatic checks

- `R-CMD-checks` - these checks make sure that the R package files are built properly. It will also rebuild the vignette underneath 4 different operating systems:
  - windows
  - macOS
  - ubuntu release
  - ubuntu dev
- `docker-build.yml` - handles the docker image associated with this pipeline. Any changes to the docker image will be tested by this GitHub Action by rebuilding it.
- `style-code.yml` - this will style your code and commit any changes to the pull request


## Testing

Whenever appropriate, new functions should be added to:
  - the `vignette/getting-started.Rmd` if they are a part of the essential workflow.
  - a unit test using testthat. You can create a new test using `usethis::use_test("name")`. then open up that file in the `tests/testthat` folder and write a test there. See this chapter for how to write tests https://r-pkgs.org/testing-basics.html


## Docker image

For reproducibility purposes, this repository has a Docker image. And software requirements should be added to this Docker image as development occurs.
