# Using samtools depth instead to be able to keep mitochondria in, think these were being filtered out by sambamba when making chromo.cov csv fles

import subprocess
import pandas as pd
import numpy as np

# Execute the first bash command to go through bam files and use sambamba to calculate depth of coverage in windows across the genome
# windows for plasmodium were 1000, but for anopheles I could use 10x this 10,000, because the genome is 10x the size
# samtools depth ur1001_Combined.mkdup.bam > ur1001_samtools_depth.csv
subprocess.run('for f in *bam ; do samtools depth "$f" > "$f.samtools_depth.csv" ; done', shell=True)

# take samtools_depth.csv and convert it to be in windows instead of each position

subprocess.run('''
for f in *bam ; do
    samtools depth "$f" | \
    awk 'BEGIN {OFS="\\t"; sum=0; count=0; prevChrom=""; start=0; end=0}
         {
           if ($1 != prevChrom || NR == 1) {
             if (NR > 1) {
               printf "%s\\t%d\\t%d\\t%.2f\\n", prevChrom, start, end, sum/count;
             }
             start = $2; 
             end = $2 + 999; # changed this to 999 for 1000 base pair windows
             sum = $3; 
             count = 1;
             prevChrom = $1;
           } else {
             if ($2 <= end) {
               sum += $3; 
               count++;
             } else {
               printf "%s\\t%d\\t%d\\t%.2f\\n", prevChrom, start, end, sum/count;
               start = end + 1; 
               end = start + 999; # changed this to 999 for 1000 base pair windows 
               sum = $3; 
               count = 1;
             }
           }
         }
         END {printf "%s\\t%d\\t%d\\t%.2f\\n", prevChrom, start, end, sum/count}' >> "$f.samtools_windowed_depth.csv"
done
''' , shell = True)


# Filter out contigs from the windowed depth files
# Read the file containing sample names
with open('samples.txt', 'r') as file:
    sample_names = file.read()
print(sample_names)

# Loop through each sample name to filter out contigs

# Read the file containing sample names
with open('samples.txt', 'r') as file:
    sample_names = file.read().splitlines()
print('sample names opened')

# Loop through each sample name to filter out contigs
for sample_name in sample_names:
  # Construct the file paths
  input_file = sample_name + '.bam.samtools_windowed_depth.csv'
  output_file = sample_name + '.bam.samtools_windowed_depth_filtered.csv'
  # Read the input file
  df = pd.read_csv(input_file, sep='\t', error_bad_lines=False, warn_bad_lines=True, header=None)
  # Add the header
  df.columns = ['# chrom', 'chromStart', 'chromEnd', 'meanCoverage']
  # Filter the data
  sub = df[df['# chrom'].isin(['Pf3D7_01_v3', 'Pf3D7_02_v3', 'Pf3D7_03_v3', 'Pf3D7_04_v3', 'Pf3D7_05_v3', 'Pf3D7_06_v3', 'Pf3D7_07_v3', 'Pf3D7_08_v3',       'Pf3D7_09_v3', 'Pf3D7_10_v3', 'Pf3D7_11_v3', 'Pf3D7_12_v3','Pf3D7_13_v3', 'Pf3D7_14_v3', 'Pf3D7_API_v3', 'Pf_M76611'])]
  # Export the filtered data to a new file
  sub.to_csv(output_file, sep='\t', index=None)

# Loop through each sample name to get mean coverage across entire chromosome and say if it meets 5x on each!!
QC=[]
mean_cov = []
for sample_name in sample_names:
  # Construct the file paths
  input_file = sample_name + '.bam.samtools_windowed_depth_filtered.csv'
  output_file1 = 'mean_windowed_depth_filtered_all.csv'
  output_file2 = 'Brazil_pass_QC.csv'
  # Read the input file
  df = pd.read_csv(input_file, sep='\t', error_bad_lines=False, warn_bad_lines=True)
  # Add the header
  #df.columns = ['# chrom', 'chromStart', 'chromEnd', 'meanCoverage']
  # Filter the data
  cov = df.groupby(by=['# chrom'])['meanCoverage'].mean()
  mean_cov.append(pd.DataFrame({'# CHROM': cov.index, 'coverage': cov, 'sample': sample_name}, index=None))
  if (all(cov[0:14]>5)):
    QC.append(pd.DataFrame({'sample': [sample_name], 'QC': ['TRUE']}))
  else:
    QC.append(pd.DataFrame({'sample': [sample_name], 'QC': ['FALSE']}))


mean_cov_df = pd.concat(mean_cov, ignore_index=True)
QC_df = pd.concat(QC, ignore_index=True)
QC_df = pd.DataFrame(QC_df)
# Export the filtered data to a new file
mean_cov_df.to_csv(output_file1,  index=None)  
QC_df.to_csv(output_file2, index=None)


# Use the other python script to create figures of coverage over the genome
subprocess.run('cat samples.txt | parallel -j 30 "python chromo_coverage.py --sample {} --csvfile {}.bam.samtools_windowed_depth_filtered.csv" > chromocoverage.txt', shell=True)
