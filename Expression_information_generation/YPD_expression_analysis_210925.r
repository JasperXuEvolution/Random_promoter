##clear data
rm(list=ls())
#Load library and source file####
library(GGally)
library(dplyr)
library(ggplot2)
library(VennDiagram)
library(grid)
library(gridExtra)
#set working directory and read data####
parent.dir="C:/Users/lumia/Dropbox (University of Michigan)/1. Haiqing Xu/1. Project/2. Random Promoter Projects/New_Analysis/New_YPD_analysis/Input"
setwd(parent.dir)
S1=read.csv(file="YPD_new_N6_Q20_total_S1_b.csv", sep=",",header = TRUE,stringsAsFactors = FALSE)
S2=read.csv(file="YPD_new_N6_Q20_total_S2_b.csv", sep=",",header = TRUE,stringsAsFactors = FALSE)
S3=read.csv(file="YPD_new_N6_Q20_total_S3_b.csv", sep=",",header = TRUE,stringsAsFactors = FALSE)

n_control<-c(1769576,2281214,2698104,4741707,2328337,1871693,3375381,2102647)
p_control<-c(2140910,3304444,3556204,3944581,2442367,3375101,1981999)

#Adding the zero RNA data as well####
#input data
S1_D=read.csv(file="YPD_new_N6_Q20_total_S1_d.csv", sep=",",header = TRUE,stringsAsFactors = FALSE)
S2_D=read.csv(file="YPD_new_N6_Q20_total_S2_d.csv", sep=",",header = TRUE,stringsAsFactors = FALSE)
S3_D=read.csv(file="YPD_new_N6_Q20_total_S3_d.csv", sep=",",header = TRUE,stringsAsFactors = FALSE)
S1_D['RNA_absolute_count']=rep(0,length(S1_D$barcode_cluster))
S1_D['RNA_relative_count']=rep(0,length(S1_D$barcode_cluster))
S2_D['RNA_absolute_count']=rep(0,length(S2_D$barcode_cluster))
S2_D['RNA_relative_count']=rep(0,length(S2_D$barcode_cluster))
S3_D['RNA_absolute_count']=rep(0,length(S3_D$barcode_cluster))
S3_D['RNA_relative_count']=rep(0,length(S3_D$barcode_cluster))
S1=S1[,1:6]
S2=S2[,1:6]
S3=S3[,1:6]
S1=rbind(S1,S1_D)
S2=rbind(S2,S2_D)
S3=rbind(S3,S3_D)

#calculate total count in each library


DNA_S1=read.csv(file="YPD_new_N6_Q20_total_DNA_S1.csv", sep=";",header = TRUE)
DNA_S2=read.csv(file="YPD_new_N6_Q20_total_DNA_S2.csv", sep=";",header = TRUE)
DNA_S3=read.csv(file="YPD_new_N6_Q20_total_DNA_S3.csv", sep=";",header = TRUE)
RNA_S1=read.csv(file="YPD_new_N6_Q20_total_RNA_S1.csv", sep=";",header = TRUE)
RNA_S2=read.csv(file="YPD_new_N6_Q20_total_RNA_S2.csv", sep=";",header = TRUE)
RNA_S3=read.csv(file="YPD_new_N6_Q20_total_RNA_S3.csv", sep=";",header = TRUE)

TC_DNA_S1=sum(DNA_S1$absolute_count)
TC_DNA_S2=sum(DNA_S2$absolute_count)
TC_DNA_S3=sum(DNA_S3$absolute_count)
TC_RNA_S1=sum(RNA_S1$absolute_count)
TC_RNA_S2=sum(RNA_S2$absolute_count)
TC_RNA_S3=sum(RNA_S3$absolute_count)

TC_DNA_S1_r=sum(DNA_S1$relative_count)
TC_DNA_S2_r=sum(DNA_S2$relative_count)
TC_DNA_S3_r=sum(DNA_S3$relative_count)
TC_RNA_S1_r=sum(RNA_S1$relative_count)
TC_RNA_S2_r=sum(RNA_S2$relative_count)
TC_RNA_S3_r=sum(RNA_S3$relative_count)




S1["N_DNA_absolute_count"]=S1["DNA_absolute_count"]/TC_DNA_S1
S1["N_RNA_absolute_count"]=S1["RNA_absolute_count"]/TC_RNA_S1
S1["N_DNA_relative_count"]=S1["DNA_relative_count"]/TC_DNA_S1_r
S1["N_RNA_relative_count"]=S1["RNA_relative_count"]/TC_RNA_S1_r

S2["N_DNA_absolute_count"]=S2["DNA_absolute_count"]/TC_DNA_S2
S2["N_RNA_absolute_count"]=S2["RNA_absolute_count"]/TC_RNA_S2
S2["N_DNA_relative_count"]=S2["DNA_relative_count"]/TC_DNA_S2_r
S2["N_RNA_relative_count"]=S2["RNA_relative_count"]/TC_RNA_S2_r

S3["N_DNA_absolute_count"]=S3["DNA_absolute_count"]/TC_DNA_S3
S3["N_RNA_absolute_count"]=S3["RNA_absolute_count"]/TC_RNA_S3
S3["N_DNA_relative_count"]=S3["DNA_relative_count"]/TC_DNA_S3_r
S3["N_RNA_relative_count"]=S3["RNA_relative_count"]/TC_RNA_S3_r

name_list<-union(union(S1$barcode_cluster,S2$barcode_cluster),S3$barcode_cluster)
Final_data=data.frame()


haha=merge(S1, S2, by='representative_sequence', all = TRUE)
haha1=merge(haha, S3, by='representative_sequence', all = TRUE)

haha1=haha1[,-c(11,20)]
#replace NA with 0
haha1[is.na(haha1)] <- 0
#haha1["DNA_absolute_count"]
names(haha1)[3:26]=c("DNA_absolute_count_S1","DNA_relative_count_S1","RNA_absolute_count_S1","RNA_relative_count_S1",
                     "N_DNA_absolute_count_S1","N_RNA_absolute_count_S1","N_DNA_relative_count_S1","N_RNA_relative_count_S1", 
                     "DNA_absolute_count_S2","DNA_relative_count_S2","RNA_absolute_count_S2","RNA_relative_count_S2",
                     "N_DNA_absolute_count_S2","N_RNA_absolute_count_S2","N_DNA_relative_count_S2","N_RNA_relative_count_S2",
                     "DNA_absolute_count_S3","DNA_relative_count_S3","RNA_absolute_count_S3","RNA_relative_count_S3",
                     "N_DNA_absolute_count_S3","N_RNA_absolute_count_S3","N_DNA_relative_count_S3","N_RNA_relative_count_S3")
haha1["DNA_absolute_count_Total"]=haha1[,3]+haha1[,11]+haha1[,19]
haha1["DNA_relative_count_Total"]=haha1[,4]+haha1[,12]+haha1[,20]
haha1["RNA_absolute_count_Total"]=haha1[,5]+haha1[,13]+haha1[,21]
haha1["RNA_relative_count_Total"]=haha1[,6]+haha1[,14]+haha1[,22]
haha1["N_DNA_absolute_count_Total"]=haha1[,7]+haha1[,15]+haha1[,23]
haha1["N_RNA_absolute_count_Total"]=haha1[,8]+haha1[,16]+haha1[,24]
haha1["N_DNA_relative_count_Total"]=haha1[,9]+haha1[,17]+haha1[,25]
haha1["N_RNA_relative_count_Total"]=haha1[,10]+haha1[,18]+haha1[,26]

Final_nc=haha1[haha1$barcode_cluster %in% n_control,]
Final_pc=haha1[haha1$barcode_cluster %in% p_control,]
haha1<-haha1[!haha1$barcode_cluster %in% union(n_control,p_control),]
#Output analyzed data
write.csv(haha1,'C:/Users/lumia/Dropbox (University of Michigan)/1. Haiqing Xu/1. Project/2. Random Promoter Projects/New_Analysis/New_YPD_analysis/Output/YPD_final_sum_up_data_V2.csv')
write.csv(Final_nc,'C:/Users/lumia/Dropbox (University of Michigan)/1. Haiqing Xu/1. Project/2. Random Promoter Projects/New_Analysis/New_YPD_analysis/Output/YPD_final_negative_control_data_V2.csv')
write.csv(Final_pc,'C:/Users/lumia/Dropbox (University of Michigan)/1. Haiqing Xu/1. Project/2. Random Promoter Projects/New_Analysis/New_YPD_analysis/Output/YPD_final_positive_control_data_V2.csv')


