# CONTRIBUTING

This document describes how to contribute to this python package.

## GitHub Workflow

- File an issue with anything you see that needs to be considered, added or fixed.
- Assign yourself to an issue that you'd like to work on. If you have any questions about how to go about the issue you've been assigned to, communicate your questions on that issue!
- Clone this repository locally and create a new branch to work on this issue from.
- When you feel you have enough to discuss, push the branch and open a pull request. When you do this, you will see a series of checks begin to run:

## Automatic checks

- `test-package.yml` - this will run tests and report the results.

## Setting Up Your Development Environment

- Follow the instructions to install from source in the [README](README.md).
- Use the editor of your choice to make changes.
- Build your changes and test with `pip install . && python -m tests`.
- If adding new dependencies, specify them in `requirements.txt` and `pip install -r requirements.txt` to update your environment. Considering using an environment manager like [conda](https://docs.conda.io/projects/conda/) to manage your environments.

## Testing

For new code you should write tests which cover both happy case and edge cases with python's built in `unittest` library. Examples of this are included in the `tests` folder.
