# Random_promoter
Scrips for analyzing the random promoter project
## 1. Data
* Sequencing data for expression and linkage is deposited at: https://www.ncbi.nlm.nih.gov/bioproject/876017
* Active promoters with expression and transcriptional factor binding site information is in file YPD_active_promoter_info_df.csv and SCD_active_promoter_info_df.csv
## 2. Scrips to extract linkage information
* Refer to the folder Linkage_information_generation
* Linkage_information_generation.bash contains all the procedures.
* For using the linkage information, please refer to YPD_RD_expression_total_method.ipynb

## 3. Scrips to extract expression information
* Refer to the folder expression_information_generation
* Expression_information_generation.bash contains all the procdures for extracting the expression information.
* YPD_expression_analysis.r contains the initial filtering and combining steps.
* YPD_RD_expression_total_method.ipynb contains the steps for further filtering and combines the promtoer information.


