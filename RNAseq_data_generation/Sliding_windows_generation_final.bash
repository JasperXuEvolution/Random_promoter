#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=XXX
#SBATCH --mail-user=XXX
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10g 
#SBATCH --time=24:00:00
#SBATCH --account=lsa1
#SBATCH --partition=standard


# your job commands below
#for Q6
LP='XXXX'
module load Bioinformatics

module load samtools

module load bedtools2/2.29.2



# Step_7 generate different tilling windows
cat $LP/bed_file/gene_merged.bed $LP/annotation_file/total_3UTR_50_raw.bed $LP/annotation_file/total_5UTR_50_raw.bed \
$LP/annotation_file/RNA_50bp_raw.bed $LP/bed_file/Centromere_merged.bed $LP/bed_file/Telomere_merged.bed \
$LP/bed_file/ARS_merged.bed $LP/bed_file/LTR_merged.bed $LP/bed_file/LTR_retrotransposon_merged.bed \
| bedtools sort | grep -v 'Mito' > $LP/bed_file/New_filter/Filter_N1.bed

cat $LP/bed_file/New_filter/Filter_N1.bed $LP/bed_file/SUT.bed | bedtools sort | grep -v 'Mito'> $LP/bed_file/New_filter/Filter_N2.bed


cat $LP/bed_file/gene_merged.bed $LP/annotation_file/total_3UTR_100_raw.bed $LP/annotation_file/total_5UTR_100_raw.bed \
$LP/annotation_file/RNA_50bp_raw.bed $LP/bed_file/Centromere_merged.bed $LP/bed_file/Telomere_merged.bed \
$LP/bed_file/ARS_merged.bed $LP/bed_file/LTR_merged.bed $LP/bed_file/LTR_retrotransposon_merged.bed \
| bedtools sort | grep -v 'Mito' > $LP/bed_file/New_filter/Filter_N3.bed

cat $LP/bed_file/New_filter/Filter_N3.bed $LP/bed_file/SUT.bed | bedtools sort | grep -v 'Mito'> $LP/bed_file/New_filter/Filter_N4.bed


cat $LP/bed_file/gene_merged.bed $LP/annotation_file/total_3UTR_200_raw.bed $LP/annotation_file/total_5UTR_200_raw.bed \
$LP/annotation_file/RNA_50bp_raw.bed $LP/bed_file/Centromere_merged.bed $LP/bed_file/Telomere_merged.bed \
$LP/bed_file/ARS_merged.bed $LP/bed_file/LTR_merged.bed $LP/bed_file/LTR_retrotransposon_merged.bed \
| bedtools sort | grep -v 'Mito' > $LP/bed_file/New_filter/Filter_N5.bed

cat $LP/bed_file/New_filter/Filter_N5.bed $LP/bed_file/SUT.bed | bedtools sort | grep -v 'Mito'> $LP/bed_file/New_filter/Filter_N6.bed


cat $LP/bed_file/gene_merged.bed $LP/annotation_file/total_3UTR_400_raw.bed $LP/annotation_file/total_5UTR_400_raw.bed \
$LP/annotation_file/RNA_50bp_raw.bed $LP/bed_file/Centromere_merged.bed $LP/bed_file/Telomere_merged.bed \
$LP/bed_file/ARS_merged.bed $LP/bed_file/LTR_merged.bed $LP/bed_file/LTR_retrotransposon_merged.bed \
| bedtools sort | grep -v 'Mito' > $LP/bed_file/New_filter/Filter_N7.bed

cat $LP/bed_file/New_filter/Filter_N7.bed $LP/bed_file/SUT.bed | bedtools sort | grep -v 'Mito'> $LP/bed_file/New_filter/Filter_N8.bed


cat $LP/bed_file/gene_merged.bed $LP/annotation_file/total_3UTR_800_raw.bed $LP/annotation_file/total_5UTR_800_raw.bed \
$LP/annotation_file/RNA_50bp_raw.bed $LP/bed_file/Centromere_merged.bed $LP/bed_file/Telomere_merged.bed \
$LP/bed_file/ARS_merged.bed $LP/bed_file/LTR_merged.bed $LP/bed_file/LTR_retrotransposon_merged.bed \
| bedtools sort | grep -v 'Mito' > $LP/bed_file/New_filter/Filter_N9.bed

cat $LP/bed_file/New_filter/Filter_N9.bed $LP/bed_file/SUT.bed | bedtools sort | grep -v 'Mito'> $LP/bed_file/New_filter/Filter_N10.bed

cat $LP/bed_file/gene_merged.bed $LP/annotation_file/total_3UTR_1200_raw.bed $LP/annotation_file/total_5UTR_1200_raw.bed \
$LP/annotation_file/RNA_50bp_raw.bed $LP/bed_file/Centromere_merged.bed $LP/bed_file/Telomere_merged.bed \
$LP/bed_file/ARS_merged.bed $LP/bed_file/LTR_merged.bed $LP/bed_file/LTR_retrotransposon_merged.bed \
| bedtools sort | grep -v 'Mito' > $LP/bed_file/New_filter/Filter_N11.bed

cat $LP/bed_file/New_filter/Filter_N11.bed $LP/bed_file/SUT.bed | bedtools sort | grep -v 'Mito'> $LP/bed_file/New_filter/Filter_N12.bed

cat $LP/bed_file/gene_merged.bed $LP/annotation_file/total_3UTR_1600_raw.bed $LP/annotation_file/total_5UTR_1600_raw.bed \
$LP/annotation_file/RNA_50bp_raw.bed $LP/bed_file/Centromere_merged.bed $LP/bed_file/Telomere_merged.bed \
$LP/bed_file/ARS_merged.bed $LP/bed_file/LTR_merged.bed $LP/bed_file/LTR_retrotransposon_merged.bed \
| bedtools sort | grep -v 'Mito' > $LP/bed_file/New_filter/Filter_N13.bed

cat $LP/bed_file/New_filter/Filter_N13.bed $LP/bed_file/SUT.bed | bedtools sort | grep -v 'Mito'> $LP/bed_file/New_filter/Filter_N14.bed

#intergenic region without mitochondria
bedtools complement -i $LP/bed_file/New_filter/Filter_N1.bed -g $LP/Yeast_genome/yeast_genome_sorted | grep -v 'Mito' > $LP/bed_file/New_filter/intergenic_N1.bed
bedtools complement -i $LP/bed_file/New_filter/Filter_N2.bed -g $LP/Yeast_genome/yeast_genome_sorted | grep -v 'Mito' > $LP/bed_file/New_filter/intergenic_N2.bed
bedtools complement -i $LP/bed_file/New_filter/Filter_N3.bed -g $LP/Yeast_genome/yeast_genome_sorted | grep -v 'Mito' > $LP/bed_file/New_filter/intergenic_N3.bed
bedtools complement -i $LP/bed_file/New_filter/Filter_N4.bed -g $LP/Yeast_genome/yeast_genome_sorted | grep -v 'Mito' > $LP/bed_file/New_filter/intergenic_N4.bed
bedtools complement -i $LP/bed_file/New_filter/Filter_N5.bed -g $LP/Yeast_genome/yeast_genome_sorted | grep -v 'Mito' > $LP/bed_file/New_filter/intergenic_N5.bed
bedtools complement -i $LP/bed_file/New_filter/Filter_N6.bed -g $LP/Yeast_genome/yeast_genome_sorted | grep -v 'Mito' > $LP/bed_file/New_filter/intergenic_N6.bed
bedtools complement -i $LP/bed_file/New_filter/Filter_N7.bed -g $LP/Yeast_genome/yeast_genome_sorted | grep -v 'Mito' > $LP/bed_file/New_filter/intergenic_N7.bed
bedtools complement -i $LP/bed_file/New_filter/Filter_N8.bed -g $LP/Yeast_genome/yeast_genome_sorted | grep -v 'Mito' > $LP/bed_file/New_filter/intergenic_N8.bed
bedtools complement -i $LP/bed_file/New_filter/Filter_N9.bed -g $LP/Yeast_genome/yeast_genome_sorted | grep -v 'Mito' > $LP/bed_file/New_filter/intergenic_N9.bed
bedtools complement -i $LP/bed_file/New_filter/Filter_N10.bed -g $LP/Yeast_genome/yeast_genome_sorted | grep -v 'Mito' > $LP/bed_file/New_filter/intergenic_N10.bed
bedtools complement -i $LP/bed_file/New_filter/Filter_N11.bed -g $LP/Yeast_genome/yeast_genome_sorted | grep -v 'Mito' > $LP/bed_file/New_filter/intergenic_N11.bed
bedtools complement -i $LP/bed_file/New_filter/Filter_N12.bed -g $LP/Yeast_genome/yeast_genome_sorted | grep -v 'Mito' > $LP/bed_file/New_filter/intergenic_N12.bed
bedtools complement -i $LP/bed_file/New_filter/Filter_N13.bed -g $LP/Yeast_genome/yeast_genome_sorted | grep -v 'Mito' > $LP/bed_file/New_filter/intergenic_N13.bed
bedtools complement -i $LP/bed_file/New_filter/Filter_N14.bed -g $LP/Yeast_genome/yeast_genome_sorted | grep -v 'Mito' > $LP/bed_file/New_filter/intergenic_N14.bed
# mkdir $LP/UTR_test/SW
# mkdir $LP/UTR_test/expression
########################################################################################################################################################
# intergenic 
# sliding window with step size of 10
bedtools makewindows -b $LP/bed_file/New_filter/intergenic_N1.bed -w 20 -s 10 > $LP/UTR_test/SW/SW_intergenic_N1_w20_s10.bed
bedtools makewindows -b $LP/bed_file/New_filter/intergenic_N2.bed -w 20 -s 10 > $LP/UTR_test/SW/SW_intergenic_N2_w20_s10.bed
bedtools makewindows -b $LP/bed_file/New_filter/intergenic_N3.bed -w 20 -s 10 > $LP/UTR_test/SW/SW_intergenic_N3_w20_s10.bed
bedtools makewindows -b $LP/bed_file/New_filter/intergenic_N4.bed -w 20 -s 10 > $LP/UTR_test/SW/SW_intergenic_N4_w20_s10.bed
bedtools makewindows -b $LP/bed_file/New_filter/intergenic_N5.bed -w 20 -s 10 > $LP/UTR_test/SW/SW_intergenic_N5_w20_s10.bed
bedtools makewindows -b $LP/bed_file/New_filter/intergenic_N6.bed -w 20 -s 10 > $LP/UTR_test/SW/SW_intergenic_N6_w20_s10.bed
bedtools makewindows -b $LP/bed_file/New_filter/intergenic_N7.bed -w 20 -s 10 > $LP/UTR_test/SW/SW_intergenic_N7_w20_s10.bed
bedtools makewindows -b $LP/bed_file/New_filter/intergenic_N8.bed -w 20 -s 10 > $LP/UTR_test/SW/SW_intergenic_N8_w20_s10.bed
bedtools makewindows -b $LP/bed_file/New_filter/intergenic_N9.bed -w 20 -s 10 > $LP/UTR_test/SW/SW_intergenic_N9_w20_s10.bed
bedtools makewindows -b $LP/bed_file/New_filter/intergenic_N10.bed -w 20 -s 10 > $LP/UTR_test/SW/SW_intergenic_N10_w20_s10.bed
bedtools makewindows -b $LP/bed_file/New_filter/intergenic_N11.bed -w 20 -s 10 > $LP/UTR_test/SW/SW_intergenic_N11_w20_s10.bed
bedtools makewindows -b $LP/bed_file/New_filter/intergenic_N12.bed -w 20 -s 10 > $LP/UTR_test/SW/SW_intergenic_N12_w20_s10.bed
bedtools makewindows -b $LP/bed_file/New_filter/intergenic_N13.bed -w 20 -s 10 > $LP/UTR_test/SW/SW_intergenic_N13_w20_s10.bed
bedtools makewindows -b $LP/bed_file/New_filter/intergenic_N14.bed -w 20 -s 10 > $LP/UTR_test/SW/SW_intergenic_N14_w20_s10.bed

# gene
bedtools makewindows -b $LP/bed_file/gene_merged.bed -w 20 -s 10 > $LP/UTR_test/SW/SW_genic_w20_s10.bed




# sliding window with step size of 20
bedtools makewindows -b $LP/bed_file/New_filter/intergenic_N1.bed -w 20 -s 20 > $LP/UTR_test/SW/SW_intergenic_N1_w20_s20.bed
bedtools makewindows -b $LP/bed_file/New_filter/intergenic_N2.bed -w 20 -s 20 > $LP/UTR_test/SW/SW_intergenic_N2_w20_s20.bed
bedtools makewindows -b $LP/bed_file/New_filter/intergenic_N3.bed -w 20 -s 20 > $LP/UTR_test/SW/SW_intergenic_N3_w20_s20.bed
bedtools makewindows -b $LP/bed_file/New_filter/intergenic_N4.bed -w 20 -s 20 > $LP/UTR_test/SW/SW_intergenic_N4_w20_s20.bed
bedtools makewindows -b $LP/bed_file/New_filter/intergenic_N5.bed -w 20 -s 20 > $LP/UTR_test/SW/SW_intergenic_N5_w20_s20.bed
bedtools makewindows -b $LP/bed_file/New_filter/intergenic_N6.bed -w 20 -s 20 > $LP/UTR_test/SW/SW_intergenic_N6_w20_s20.bed
bedtools makewindows -b $LP/bed_file/New_filter/intergenic_N7.bed -w 20 -s 20 > $LP/UTR_test/SW/SW_intergenic_N7_w20_s20.bed
bedtools makewindows -b $LP/bed_file/New_filter/intergenic_N8.bed -w 20 -s 20 > $LP/UTR_test/SW/SW_intergenic_N8_w20_s20.bed
bedtools makewindows -b $LP/bed_file/New_filter/intergenic_N9.bed -w 20 -s 20 > $LP/UTR_test/SW/SW_intergenic_N9_w20_s20.bed
bedtools makewindows -b $LP/bed_file/New_filter/intergenic_N10.bed -w 20 -s 20 > $LP/UTR_test/SW/SW_intergenic_N10_w20_s20.bed
bedtools makewindows -b $LP/bed_file/New_filter/intergenic_N11.bed -w 20 -s 20 > $LP/UTR_test/SW/SW_intergenic_N11_w20_s20.bed
bedtools makewindows -b $LP/bed_file/New_filter/intergenic_N12.bed -w 20 -s 20 > $LP/UTR_test/SW/SW_intergenic_N12_w20_s20.bed
bedtools makewindows -b $LP/bed_file/New_filter/intergenic_N13.bed -w 20 -s 20 > $LP/UTR_test/SW/SW_intergenic_N13_w20_s20.bed
bedtools makewindows -b $LP/bed_file/New_filter/intergenic_N14.bed -w 20 -s 20 > $LP/UTR_test/SW/SW_intergenic_N14_w20_s20.bed

# gene
bedtools makewindows -b $LP/bed_file/gene_merged.bed -w 20 -s 20 > $LP/UTR_test/SW/SW_genic_w20_s20.bed
