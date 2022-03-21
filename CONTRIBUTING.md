# How to contribute

You want to contribute? That's awesome. Thank you! Things can only get better from here. 

## Testing

The model is tested with Python's [unittest](https://docs.python.org/3/library/unittest.html) framework. Please write new tests (where applicable) for any code you create. The repository is linked to Circleci's CI tool so any commited code will be built and run automatically and test results will be visible both on the branch and when pull requests are made.

To locally test the model and/or any changes you make to the source code, set your working directory to the root of the project `py_SBeLT/` and run the following

```bash
python3 -m unittest src/tests/test_*
```

## Submitting Changes

Once you're ready to submit code, please create a pull request that clearly outlines the changes you've made. Additionally, when committing code please use 
clear commit messages describing the feature being committed. This helps keep track of changes to the project and makes reviews much easier!

## Coding Conventions

- Docstrings are written in accordance with [Google Python Style Guide](http://google.github.io/styleguide/pyguide.html)
- If you think the code isn't simple enough to understand alone, write comments!
- Use tabs for spacing (I know, I know...)
- Any logic related to the model should exist in `logic.py`
- Keep variable/function/class names in line with the vocabulary described in `docs/NOMECULTURE.md`
