[![status](https://joss.theoj.org/papers/d7b9cc16b87e8875ec7115a22e1413fe/status.svg)](https://joss.theoj.org/papers/d7b9cc16b87e8875ec7115a22e1413fe)
[![szwiep](https://circleci.com/gh/szwiep/py_SBeLT.svg?style=svg)](https://app.circleci.com/pipelines/github/szwiep/py_SBeLT?filter=all)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/szwiep/py_SBeLT/master?labpath=docs%2Fnotebooks%2F)
# Py_SBeLT

|![gif of 500 model runs](figures/Cropped_ModelGif.gif)
|:--:|
| *Figure 1. A .gif of 500 iterations from a `py_SBeLT`run using default parameters. |

Rivers transport sediment particles. Individual particles can exhibit transport behavior that differs significantly when compared to other particles. py_SBeLT--Stochastic Bed Load Transport--provides a simple Python framework to numerically examine how individual particle motions in rivers combine to produce rates of transport that can be measured at one of a number of downstream points. The model can be used for basic research, and the model's relatively straightforward set-up makes it an effective and efficient teaching tool to help students build intuition about river transport of sediment particles.

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

### Testing your installation

If you've installed from source, you can test the installation by setting your working directory to `py_SBeLT/` 
and running the following

```bash
python -m unittest src/tests/test_*
```

## Getting Started

Users can work through the Jupyter Notebooks provided to gain a better understanding of py_SBeLT's basic usage, potential, and data storage methods. Either launch the binder instance ([![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/szwiep/py_SBeLT/master?labpath=docs%2Fnotebooks%2F)), clone the repository, or download the notebooks directly to get started.

If notebook's aren't your thing, simply run:

```bash
sbelt-run
```

or

```bash
from sbelt import sbelt_runner
sbelt_runner.run()
```

For help, reach out with questions to the repository owner `szwiep` and reference the documenation in `docs/` and `paper/`! 


## Documentation

Documentation, including Jupyter Notebooks, API documentation, [default parameters](https://github.com/szwiep/py_SBeLT/blob/master/docs/DEFAULT_PARAMS.md), and [model nomenclature](https://github.com/szwiep/py_SBeLT/blob/master/docs/NOMENCLATURE.md), can be found in the repository's `docs/` directory. Additional information regarding the theoritical motivation of the model can be found in the [`paper/paper.md`](https://github.com/szwiep/py_SBeLT/blob/master/paper/paper.md) and [`THEORY.md`](https://github.com/szwiep/py_SBeLT/blob/master/THEORY.md) files.

Two Jupyter Notebooks describe [basic usage](https://github.com/szwiep/py_SBeLT/blob/master/docs/notebooks/basic_usage_sbelt.ipynb) and [the structure of the output hdf5 format file](data_storage_sbelt.ipynb).

The API documentation is in HTML format. These files can either be downloaded and viewed directly in your browser or can be viewed using the [GitHub HTML preview project](https://htmlpreview.github.io/). 
  - [sbelt_runner](https://htmlpreview.github.io/?https://github.com/szwiep/py_SBeLT/blob/update_docs/docs/API/sbelt_runner.html)
  - [plotting](https://htmlpreview.github.io/?https://github.com/szwiep/py_SBeLT/blob/update_docs/docs/API/plotting.html)
  - [logic](https://htmlpreview.github.io/?https://github.com/szwiep/py_SBeLT/blob/update_docs/docs/API/logic.html)
  - [test_logic](https://htmlpreview.github.io/?https://github.com/szwiep/py_SBeLT/blob/update_docs/docs/API/test_logic.html)
  - [utils](https://htmlpreview.github.io/?https://github.com/szwiep/py_SBeLT/blob/update_docs/docs/API/utils.html)


## Attribution and Citation

If you use `Simframe`, please remember to cite (to be updated later)[]().

```
@article{pySBelt,
  doi = {},
  url = {},
  year = {2022},
  publisher = {The Open Journal},
  volume = {},
  number = {},
  pages = {},
  author = {Sarah Zwiep, Shawn Chartrand},
  title = {pySBeLT: A Python software for stochastic sediment transport under rarefied conditions},
  journal = {Journal of Open Source Software}
}

```

## Ackowledgements

`py_SBeLT` has received funding from NSERC Undergraduate Student Research Awards Program and Simon Fraser University.

`py_SBeLT` was developed at the [Simon Fraser University](https://www.sfu.ca/) within the [School of Environmental Science](https://www.sfu.ca/evsc.html).
