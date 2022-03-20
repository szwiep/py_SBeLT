""" Set up gently borrowed from:
https://github.com/pypa/sampleproject/blob/main/setup.py
See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / 'README.md').read_text(encoding='utf-8')

setup(
    name='sbelt',  # Required

    # Versions should comply with PEP 440:
    # https://www.python.org/dev/peps/pep-0440/
    version='1.0.1',  # Required

    # A one-line description of what your project does.
    description='A Markov-type numerical model of sediment particle transport in rivers',  
    
    long_description=long_description,  
    
    long_description_content_type='text/markdown',  
    
    url='https://github.com/szwiep/py_SBeLT',  

    author='S. Zwiep', 
     
    author_email='szwiep@sfu.ca',  

    keywords='earth science, hydrology, sediment transport',  

    # When your source code is in a subdirectory under the project root, e.g.
    # `src/`, it is necessary to specify the `package_dir` argument.
    package_dir={'': 'src'},  # Optional

    # You can just specify package directories manually here if your project is
    # simple. Or you can use find_packages().
    #
    # Alternatively, if you just want to distribute a single Python file, use
    # the `py_modules` argument instead as follows, which will expect a file
    # called `my_module.py` to exist:
    #
    #   py_modules=["my_module"],
    #
    packages=['sbelt.plots', 'sbelt', 'tests'],  # Required

    # Specify which Python versions you support. In contrast to the
    # 'Programming Language' classifiers above, 'pip install' will check this
    # and refuse to install the project if the version does not match. See
    # https://packaging.python.org/guides/distributing-packages-using-setuptools/#python-requires
    python_requires='>=3.6, <4',

    # This field lists other packages that your project depends on to run.
    # Any package you put here will be installed by pip when your project is
    # installed, so they must be valid existing projects.
    #
    # For an analysis of "install_requires" vs pip's requirements files see:
    # https://packaging.python.org/discussions/install-requires-vs-requirements/
    install_requires=['numpy', 
                        'scipy',
                        'h5py',
                        'tqdm',
                        'matplotlib'],  

    # For example, the following would provide a command called `sample` which
    # executes the function `main` from this package when invoked:
    entry_points={  
        'console_scripts': [
            'sbelt-run=sbelt.sbelt_runner:run',
        ],
    },
)
