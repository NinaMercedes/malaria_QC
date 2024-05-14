# malaria_QC
Calculating coverage and making coverage plots for P falciparum data:

'create_coverage_plots_samtoolsdepth.py' uses samtools to loop over all cram files and calculate coverage for each file. Samtools and bash are then used to make (1kb) windows calculations of coverage to plot on a coverage plot which is made using 'chromo_coverage.py'.

We are often interested in other QC metrics such as % of genome with coverage> 5. This can be calculated using the 'parse_coverage.R' script or 'get_region_cov.R' for specific regions.

