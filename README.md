# Random_promoter
Scrips for analyzing the data from "Chance promoter activities illuminate the origins of eukaryotic intergenic transcriptions"
## 1. Data
* Sequencing data for expression and linkage is deposited to NCBI (PRJNA876017 [https://www.ncbi.nlm.nih.gov/bioproject/876017])
* RNA-Seq data for is from PRJNA728585 [https://www.ncbi.nlm.nih.gov/bioproject/PRJNA728585] and PRJNA392312 [https://www.ncbi.nlm.nih.gov/bioproject/PRJNA392312]
* Active promoter information is in the folder Data 
* All the intermediate data used for this study can be found at figshare [https://figshare.com/articles/dataset/Intermediate_data_for_Chance_promoter_activities_illuminate_the_origins_of_eukaryotic_intergenic_transcriptions_/22231603]
* Source data are in the folder data

## 2. Scrips to extract linkage information
* Refer to the folder Linkage_information_generation
* Linkage_information_generation.bash contains the pipeline
* For using the linkage information, please refer to YPD_RD_expression_total_method.ipynb

## 3. Scrips to extract expression information
* Refer to the folder expression_information_generation
* Expression_information_generation.bash contains the pipeline

## 4. Expression data generation for random promoters
* Refer to the folder expression analysis
* RD_data_processing_SCD.ipynb and RD_data_processing_YPD.ipynb contains the initial filtering and combining steps
* RD_data_processing_Step2.ipynb contains the steps for further filtering and compared expression of random promoter to controls

## 5. RNAseq data pipelines
* Refer to the folder RNAseq_data_generation
* YPD_data_pipeline.bash/SCD_data_pipeline.bash are the pipeline files

## 6. Motif analysis
* Refer to the folder Motif analysis
## 7. Plot figures
* Refer to the folder Figures
* Source_data_plotting_main_to_figS10.ipynb and Source_data_plotting_main_to_figS11_20.ipynb can plot all the figures in the paper using the source_data.zip
* Figs_final.ipynb can generate all the plots from the intermediate data

