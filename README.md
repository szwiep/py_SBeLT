# Py_SBeLT

![gif of 500 model runs](figures/Cropped_ModelGif.gif)

Rivers transport sediment particles. Individual particles can exhibit transport behavior that differs significantly when compared to other particles. py_SBeLT provides a simple Python framework to numerically examine how individual particle motions in rivers combine to produce rates of transport that can be measured at one of a number of downstream points. The model can be used for basic research, and the model's relatively straightforward set-up makes it an effective and efficient teaching tool to help students build intuition about river transport of sediment particles.

## Installation

### Quick Installation

```bash
pip install sbelt
```

### Installation from Source

Clone the `py_SBeLT` GitHub repository

```bash
git clone https://github.com/szwiep/py_SBeLT.git
```

Then set your working directory to `py_SBeLT/` and build the project

```bash
 cd py_SBeLT/
 python setup.py build_ext --inplace
 pip install -e .
```

## Documentation

Documentation, including Jupyter Notebooks outlining basic usage and data storage methods, API documentation, and model nomeculture/vocabulary, can be found in the repository's `docs/` directory. Additional information regarding the theoritical motivation of the model can be found in the `paper/paper.md` and `THEORY.md` files.

## Examples

## Attribution

TBD
