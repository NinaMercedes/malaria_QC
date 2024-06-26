# malaria_QC
Calculating coverage and making coverage plots for P falciparum data:

'create_coverage_plots_samtoolsdepth.py' uses samtools to loop over all cram files and calculate coverage for each file. Samtools and bash are then used to make (1kb) windows calculations of coverage to plot on a coverage plot which is made using 'chromo_coverage.py'. 'create_coverage_plots_samtoolsdepth_bam.py' does the same thing but uses bam files instead (e.g. for nanopore)!

For use reccomend installing the following conda environment.

```
conda create -n cov_plots_malaria
conda activate cov_plots_malaria
conda install python=3.7
conda install matplotlib
conda install samtools
conda install R
conda install conda-forge::r-dplyr
conda install conda-forge::r-data.table
conda install conda-forge::r-readr
```

We are often interested in other QC metrics such as % of genome with coverage> 5. This can be calculated using the 'parse_coverage.R' script or 'get_region_cov.R' for specific regions.

Credit to Sophie Moss, Holly Acford Palmer, Matthew Higgins and Emilia Manko for code development.
