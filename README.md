# Random_promoter
Scrips for analyzing the random promoter project
## 1. Data
* Sequencing data for expression and linkage is deposited to NCBI (PRJNA876017 [https://www.ncbi.nlm.nih.gov/bioproject/876017])
* RNA-Seq data for is from PRJNA728585 [https://www.ncbi.nlm.nih.gov/bioproject/PRJNA728585] and PRJNA392312 [https://www.ncbi.nlm.nih.gov/bioproject/PRJNA392312]
* Active promoter information is in the Data folder
## 2. Scrips to extract linkage information
* Refer to the folder Linkage_information_generation
* Linkage_information_generation.bash contains all the procedures.
* For using the linkage information, please refer to YPD_RD_expression_total_method.ipynb

## 3. Scrips to extract expression information
* Refer to the folder expression_information_generation
* Expression_information_generation.bash contains all the procdures for extracting the expression information.
* YPD_expression_analysis.r contains the initial filtering and combining steps.
* YPD_RD_expression_total_method.ipynb contains the steps for further filtering and combines the promtoer information.

## 4. Expression data generation for random promoters
* Refer to the folder expression analysis
* RD_data_processing.ipynb is for processing random promoter expression data

## 5. RNAseq data pipelines
* Refer to the folder RNAseq_data_generation

## 5. Plot figures
* Refer to the folder Figures

