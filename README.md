
# mtor-pancan

Analysis of TCGA/IMPACT datasets, and mTOR mutation status across cancer types

## Visualizations

The most up-to-date output of this workflow is [available online](http://www.cs.toronto.edu/~charrigan/plots/a0320a8/)

## To reproduce

All plots are generated with [open access data](https://portal.gdc.cancer.gov/). 
These can be individually reproduced by knitting `.Rmd` files, or generated en-mass using the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow provided. 
[MutSigCV](https://software.broadinstitute.org/cancer/cga/mutsig) is used to generate p-values for significantly mutated genes.
