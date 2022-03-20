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

## Getting Started

Users can work through the Jupyter Notebooks provided to gain a better understanding of py_SBeLT's basic usage, potential, and data storage methods. Either launch the binder instance (), clone the repository, or download the notebooks directly to get started.

If notebook's aren't your thing, simply run:

```bash
sbelt-run
```

or

```bash
from sbelt import sbelt_runner
sbelt_runner.run()
```

To get started. For help, reach out with questions to the repository owner `szwiep` and reference the documenation in `docs/` and `paper/`! 


## Documentation

Documentation, including Jupyter Notebooks, API documentation, default parameters, and model nomeculture, can be found in the repository's `docs/` directory. Additional information regarding the theoritical motivation of the model can be found in the `paper/paper.md` and `THEORY.md` files.





## Attribution

TBD
