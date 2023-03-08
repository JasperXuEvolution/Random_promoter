#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=XXX
#SBATCH --mail-user=XXX
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=15g 
#SBATCH --time=72:00:00
#SBATCH --account=XXX
#SBATCH --partition=standard


# your job commands below
LP='XXX'

#load modules
module load Bioinformatics
module load fastqc/0.11.8
module load adapterremoval/2.3.1
module load bartender/1.1
module load singularity

# Step 1 QC
mkdir -p "$LP/QC"
fastqc -o $LP/QC $LP/Sequencing_data/19186FL-07-01_S9_L004_R1_001.fastq.gz $LP/Sequencing_data/19186FL-07-01_S9_L004_R2_001.fastq.gz


# Step 2 Adapter removal
mkdir -p "$LP/Merged"
AdapterRemoval --file1 $LP/Sequencing_data/19186FL-07-01_S9_L004_R1_001.fastq.gz  --file2 $LP/Sequencing_data/19186FL-07-01_S9_L004_R2_001.fastq.gz \
--adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCTCGTATGCCGTCTTCTGCTTG --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTGTGTAGATCTCGGTGGTCGCCGTATCATT \
--basename $LP/Merged/YPD_New_N6 --trimns --trimqualities --collapse 
# Combine truncated and intact collapsed reads into one 
cat $LP/Merged/YPD_New_N6.collapsed $LP/Merged/YPD_New_N6.collapsed.truncated > $LP/Merged/YPD_New_N6_total

# Step 3 Length filtering

# I request collapsed sequence to be more than 80 bp. The target length is 85 bp.
awk 'BEGIN {RS="\n@";FS="\n"} {if (length($2) >=80) {print "@"$0} }' $LP/Merged/YPD_New_N6_total > $LP/Filtered_data/YPD_New_N6_filtered_S1.fastq


# Step 4 Barcode extraction

mkdir -p "$LP/Filtered_data"
python YPD_new_N6_Barcode_extract.py --a $LP/Filtered_data/YPD_New_N6_filtered_S1.fastq \
--o $LP/Filtered_data/YPD_new_N6_Q20 --q 20

# Step 5 Clustering barcode and UMI, and combine results
mkdir -p "$LP/Bartender"
bartender_single_com -z -1 -d 2 -l 5 -f $LP/Filtered_data/YPD_new_N6_Q20_barcode \
-o $LP/Bartender/YPD_new_N6_Q20_barcode_d2

bartender_single_com -z -1 -d 2 -l 5 -f $LP/Filtered_data/YPD_new_N6_Q20_UMI \
-o $LP/Bartender/YPD_new_N6_Q20_UMI
mkdir -p "$LP/Cluster_analysis"
python YPD_new_N6_Barcode_analysis_bartender_data.py --a1 $LP/Bartender/YPD_new_N6_Q20_barcode_d2_barcode.csv --a2 $LP/Bartender/YPD_new_N6_Q20_barcode_d2_cluster.csv\
 --a3 $LP/Bartender/YPD_new_N6_Q20_UMI_barcode.csv --a4 $LP/Bartender/YPD_new_N6_Q20_UMI_cluster.csv --a5 $LP/Filtered_data/YPD_new_N6_Q20_information_list \
 --a6 $LP/Filtered_data/YPD_new_N6_Q20_barcode --o $LP/Cluster_analysis/YPD_new_N6_Q20_total