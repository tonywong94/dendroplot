# Scripts used in [Wong et al. (2019)](https://doi.org/10.3847/1538-4357/ab46ba)

To reproduce the main results of this paper, you should download from the [data repository](https://doi.org/10.13012/B2IDB-7090706_V1) this file:

- lmc_cyc4_data.zip

and unpack it within the `images` subdirectory.  The data repository is missing the SAGE 8µm images, which can be downloaded [from this link](https://mmwave.astro.illinois.edu/almalmc/lmc_cyc4_sage.zip) and also unpacked in the `images` subdirectory. The main scripts to run have file names starting with a number and will reproduce most of the LTE results in `lte` and the dendrogram results in `dendro`.  Additional scripts can largely be run using the output in `dendro`.  More details about the tables and figures in the paper are given below.

- Table 2: the data are generated by `2-doscimes_3p5.py` which writes a few TeX files.

- Figures 1-2: Not available.

- Figures 3-4: these are generated by `0-mom0plot_3p5.py`.

- Figure 5: this is generated by `2-doscimes_3p5.py`.  The filename is `A439_12_dendrogram.pdf`.  Note that this script may have a memory leak and you may need to break the clouds into separate lists and re-run the script to complete it.

- Figures 6-7: these are generated by `3-doplots_3p5.py`.

- Figures 8-11: these are generated by `4-do_comp_3p5.py`.

- Figure 12: these are generated by `rdv_intercept.py`.  The N113 data are not included here though.

- Figure 13: this is generated by `energy.py`.  Note that the caption in the paper is wrong: the x-axis is the CO-based mass of each structure.
