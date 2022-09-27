# dendroplot
Python analysis scripts for applying dendrograms to LMC ALMA data, and visualizing the results.  Still under active development, but currently available are the following:

`from dendroplot import run_dendro`

`from dendroplot.analysis import calc_phys_props, find_clusters, get_response_width`

`from dendroplot.lte import lte, add_ltemass`

`from dendroplot.plotting import pltprops, colorcode, comp_props`

To get started, please see the sample notebook in `example` and the various scripts in `published`.

## Required packages

- [numpy](https://numpy.org)
- [matplotlib](https://matplotlib.org)
- [astropy](https://astropy.org)
- [astrodendro](https://github.com/tonywong94/astrodendro) (fork by Tony Wong)
- [scimes](https://github.com/tonywong94/SCIMES) (fork by Tony Wong)
- [kapteyn](https://www.astro.rug.nl/software/kapteyn/)

## Installation

Install the package directly from PyPI (https://pypi.org/project/dendroplot/) using

    pip install dendroplot

or, if you prefer to keep the latest source code handy, by cloning this Github repository and running

    pip install .

in the directory containing `setup.py`.