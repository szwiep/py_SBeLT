# How to contribute

You want to contribute? That's awesome. Thank you! Things can only get better from here. 

## Testing

The model is tested with Python's [unittest](https://docs.python.org/3/library/unittest.html) framework. Please write new tests (where applicable) for any code you create. The repository is linked to Circleci's CI tool so any commited code will be built and run automatically and test results will be visible both on the branch and when pull requests to master are made. Locally you can test your changes by calling the unittest module on the test files.

## Submitting Changes

Once you're ready to submit code, please create a pull request that clearly outlines the changes you've made. Additionally, when committing code please use 
clear commit messages describing the feature being committed. This helps keep track of changes to the project and makes reviews much easier!

## Coding Conventions

- Functions should all have docstrings in accordance with [PEP-0257](https://www.python.org/dev/peps/pep-0257/)
- If you think the code isn't simple enough to understand alone, write comments!
- Use tabs for spacing
- Use snake_case for function and variable names and camelCase for class names
- Avoid putting logic in `run.py`
- Keep variable/function/class names in line with the vocabulary described in `README.md`
