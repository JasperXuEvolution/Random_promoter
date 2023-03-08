#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=XX
#SBATCH --mail-user=XXX
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=15g 
#SBATCH --time=72:00:00
#SBATCH --account=lsa1
#SBATCH --partition=standard

# Home directory
LP='XXXX'

#load modules
module load Bioinformatics
module load fastqc/0.11.8
module load adapterremoval/2.3.1
module load bartender/1.1
module load singularity
module load cd-hit/4.7

# Step 1 QC

mkdir -p "$LP/QC"
fastqc -o $LP/QC $LP/Sequencing_data/19186FL-07-02_S10_L005_R1_001.fastq.gz $LP/Sequencing_data/19186FL-07-02_S10_L005_R2_001.fastq.gz


# Step 2 Adapter removal
mkdir -p "$LP/Merged"
AdapterRemoval --file1 $LP/Sequencing_data/19186FL-07-02_S10_L005_R1_001.fastq.gz  --file2 $LP/Sequencing_data/19186FL-07-02_S10_L005_R2_001.fastq.gz \
--adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAACCGCGGATCTCGTATGCCGTCTTCTGCTTG --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGCGCTAGGTGTAGATCTCGGTGGTCGCCGTATCATT \
--basename $LP/Merged/N6_Linkage --trimns --trimqualities --collapse 
cat $LP/Merged/N6_Linkage.collapsed $LP/Merged/N6_Linkage.collapsed.truncated > $LP/Merged/N6_Linkage_total
#Step_3 extract promoter and barcode
mkdir -p "$LP/Filtered_data"
python Total_lib_extraction_March_V2.py --a $LP/Merged/N6_Linkage_total --o $LP/Filtered_data/Linkage_lib_20b_120p_q10 --q 10

#Step_4_1 cluster barcode using bartender
python Barcode_deduplex.py --a $LP/Filtered_data/Linkage_lib_20b_q10_barcode \
--b $LP/Filtered_data/Linkage_barcode_info --c $LP/Filtered_data/Linkage_barcode_bartender_input

mkdir -p "$LP/Bartender"
bartender_single_com -z -1 -d 2 -l 5 -f $LP/Filtered_data/Linkage_barcode_bartender_input \
-o $LP/Bartender/Linkage_barcode_clustered

#Step_4_2 cluster promoter using cd-hit-est
mkdir -p "$LP/Cluster_cdhit"
python Promoter_deduplex.py --a $LP/Filtered_data/Linkage_lib_20b_q10_promoter \
--b $LP/Filtered_data/Linkage_promoter_info --c $LP/Filtered_data/Linkage_promoter_cdhit_input

cd-hit-est -i $LP/Filtered_data/Linkage_promoter_cdhit_input \
-o $LP/Cluster_cdhit/Linkage_promoter_clustered -d 30 -n 10 -g 1 -p 1 -r 0 -c 0.95 -uS 0.05 -uL 0.05 -mismatch -1 -M 20000 -T 0


#Step_5 combine barcode and promoter 

mkdir -p "$LP/Linkage_summary"
python New_barcode_promoter_combine.py --a1 $LP/Filtered_data/Linkage_lib_20b_q10_barcode \
--a2 $LP/Bartender/Linkage_barcode_clustered_barcode.csv \
--a3 $LP/Bartender/Linkage_barcode_clustered_cluster.csv \
--a4 $LP/Filtered_data/Linkage_promoter_info \
--a5 0.8 \
--a6 $LP/Filtered_data/Linkage_barcode_info \
--a7 $LP/Cluster_cdhit/Linkage_promoter_clustered.clstr \
--o1 $LP/Linkage_summary/linkage_barcode_promoter_combined \
--o2 $LP/Linkage_summary/Linkage_barcode \
--o3 $LP/Linkage_summary/Linkage_promoter
