
# mtor-rcc

Analysis of TCGA datasets, and mTOR mutations across cancer types

# Dependencies

* Snakemake
* R
* [MutSigCV](https://software.broadinstitute.org/cancer/cga/mutsig)

The dev environment can be reproduced using `conda env create -f mtor-rcc.yml`. 

Additional R packages required can be installed via BiocManager:

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("maftools")
BiocManager::install("TCGAbiolinks")
BiocManager::install("PoisonAlien/TCGAmutations")
```

# Run 

Run the analysis with the following:

`conda activate mtor-rcc`
`snakemake -s Snakefile -j1`

The flag `-j1` tells snakemake to parallelize over one job. In `Snakefile`, [taskset](https://man7.org/linux/man-pages/man1/taskset.1.html) is used to restrict the cores that MutSigCV is allowed to run on (27,28,29,30) - you may wish to customize this to your machine/compute environment.

# Results

* MutSigCV results are by default output to the directory `mutsigcv_results/<sample set>/<source>` 
* Plots are output to `plots/<sample set>/<source>`

