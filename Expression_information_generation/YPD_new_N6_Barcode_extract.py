#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Here I import all the package that I potentially will use to analyse my data
import sys
from Bio import SeqIO
from fuzzysearch import find_near_matches
import itertools
from Bio import Align 
import argparse
import statistics as ST
aligner = Align.PairwiseAligner()
aligner.open_gap_score = -0.5
aligner.extend_gap_score = -0.5


# In[2]:


parser = argparse.ArgumentParser(description='A function to extract barcode, UMI and other information!')
parser.add_argument("--a", required=True, help="This is the input file of reads")
parser.add_argument("--o", required=True, help="This is the name of output file")
parser.add_argument("--q", required=True, help="This is quality cutoff")

args = parser.parse_args()
input_file=args.a
output_name=args.o
Q_cutoff=int(args.q)

# In[3]:


def Find_neighbour(ref,input_seq,d,method):
    if method ==2:# 1 is using alignment, 2 is using hamming distance
        for x in ref:
            if hamming_distance(x,input_seq)<=d:
                return(x)
    else:
        for x in ref:
            if aligner.score(x,input_seq)>=d:
                return(x)


# In[4]:


ref_list=["CTAGCGCT","TCGATATC","CGTCTGCG","CGCTATGT","ACGCACCT","GTATGTTC"]


# In[5]:


M_index_List={"CGCTATGT":"RNA_S1","TCGATATC":"RNA_S2","CGTCTGCG":"RNA_S3","CTAGCGCT":"DNA_S1","ACGCACCT":"DNA_S2","GTATGTTC":"DNA_S3"}


# In[6]:


constant_2="TGGCTGAACCTAGTTTTGCCC"
constant_3="CTCTCGGAACATAGCAGTTT"


# In[ ]:


output1=output_name+"_information_list"
output2=output_name+"_barcode"
output3=output_name+"_UMI"


# In[28]:

#This is requiring the error for each position to be less than 25%
file_a = open(output1,'w')
file_b = open(output2,'w')
file_c = open(output3,'w')
tc=1
number_pass_length_filter=1
total_average_quality=0
with open(input_file, "r") as handle:
    for record in SeqIO.parse(handle, "fastq"):
        item_2=record.seq
        if (len(item_2)>=82) and (len(item_2)<=88):
            rb,UMI_1,UMI_2,Sample_condition=[],[],[],[]
            start1,end1=[],[]
            #first I will check if they have the multiplexing barcode for the first 16 barcode (the expected length is 6+8=14)
            if item_2[6:14] in M_index_List.keys():# First, I check the multiplexing barcode
                UMI_1=item_2[0:6]# This is the forward barcode
                Sample_condition=M_index_List[item_2[6:14]]
                if constant_2 in item_2:
                    start1=item_2.find(constant_2)+len(constant_2)
                    if constant_3 in item_2:
                        end1=item_2.find(constant_3)
                        UMI_2=item_2[end1+24-len(item_2):]
                        if len(UMI_2)<6:
                            UMI_2=UMI_2+"N"*(6-len(UMI_2))
                        if len(UMI_2)>10:
                            UMI_2=item_2[-6:]
                        UMI_2=UMI_2[0:6]
                        rb=str(item_2[start1:end1])
                        quality=record[start1:end1].letter_annotations["phred_quality"]
                        total_average_quality+=ST.mean(quality)
                        if len(rb)==20:
                            number_pass_length_filter+=1
                            if min(quality)>Q_cutoff:
                                stringtowrite_b="{},{}\n".format(rb,str(tc),)
                                file_b.write(stringtowrite_b)# this for writing the fasta file that containg the barcode
                                stringtowrite_c="{},{}\n".format(str(UMI_1+UMI_2),str(tc))
                                file_c.write(stringtowrite_c)# this for writing the fasta file that containg the UMI
                                line="{},{},{},{},{}\n".format(str(tc),UMI_1,UMI_2,rb,Sample_condition)
                                file_a.write(line) # this is for generating the csv file containing all the information
                                tc+=1
                    else:
                        ts3=item_2[-32:-8] #this is target region for finding match  
                        Test_c=find_near_matches(constant_3, ts3, max_l_dist=3)
                        if Test_c:
                            end1=-32+Test_c[0].start
                            UMI_2=item_2[end1+24:]
                            if len(UMI_2)<6:
                                UMI_2=UMI_2+"N"*(6-len(UMI_2))
                            if len(UMI_2)>10:
                                UMI_2=item_2[-6:]
                            UMI_2=UMI_2[0:6]
                            rb=str(item_2[start1:end1])
                            quality=record[start1:end1].letter_annotations["phred_quality"]
                            total_average_quality+=ST.mean(quality)
                            if len(rb)==20:
                                number_pass_length_filter+=1
                                if min(quality)>Q_cutoff:
                                    stringtowrite_b="{},{}\n".format(rb,str(tc),)
                                    file_b.write(stringtowrite_b)# this for writing the fasta file that containg the barcode
                                    stringtowrite_c="{},{}\n".format(str(UMI_1+UMI_2),str(tc))
                                    file_c.write(stringtowrite_c)# this for writing the fasta file that containg the UMI
                                    line="{},{},{},{},{}\n".format(str(tc),UMI_1,UMI_2,rb,Sample_condition)
                                    file_a.write(line)
                                    tc+=1
                else:
                    ts2=item_2[12:37]  #this is target region for finding match
                    Test_b=find_near_matches(constant_2, ts2, max_l_dist=3)
                    if Test_b:
                        start1=12+Test_b[0].end
                        if constant_3 in item_2:
                            end1=item_2.find(constant_3)
                            UMI_2=item_2[end1+24-len(item_2):]
                            if len(UMI_2)<6:
                                UMI_2=UMI_2+"N"*(6-len(UMI_2))
                            if len(UMI_2)>10:
                                UMI_2=item_2[-6:]
                            UMI_2=UMI_2[0:6]
                            rb=str(item_2[start1:end1])
                            quality=record[start1:end1].letter_annotations["phred_quality"]
                            total_average_quality+=ST.mean(quality)
                            if len(rb)==20:
                                number_pass_length_filter+=1
                                if min(quality)>Q_cutoff:                                
                                    stringtowrite_b="{},{}\n".format(rb,str(tc),)
                                    file_b.write(stringtowrite_b)# this for writing the fasta file that containg the barcode
                                    stringtowrite_c="{},{}\n".format(str(UMI_1+UMI_2),str(tc))
                                    file_c.write(stringtowrite_c)# this for writing the fasta file that containg the UMI	
                                    line="{},{},{},{},{}\n".format(str(tc),UMI_1,UMI_2,rb,Sample_condition)
                                    file_a.write(line)

                                    tc+=1
                        else:
                            ts3=item_2[-32:-8] #this is target region for finding match  
                            Test_c=find_near_matches(constant_3, ts3, max_l_dist=3)
                            if Test_c:
                                end1=-32+Test_c[0].start
                                UMI_2=item_2[end1+24:]
                                if len(UMI_2)<6:
                                    UMI_2=UMI_2+"N"*(6-len(UMI_2))
                                if len(UMI_2)>10:
                                    UMI_2=item_2[-6:]
                                UMI_2=UMI_2[0:6]
                                rb=str(item_2[start1:end1])
                                quality=record[start1:end1].letter_annotations["phred_quality"]
                                total_average_quality+=ST.mean(quality)
                                if len(rb)==20:
                                    number_pass_length_filter+=1
                                    if min(quality)>Q_cutoff:                                    
                                        stringtowrite_b="{},{}\n".format(rb,str(tc),)
                                        file_b.write(stringtowrite_b)# this for writing the fasta file that containg the barcode
                                        stringtowrite_c="{},{}\n".format(str(UMI_1+UMI_2),str(tc))
                                        file_c.write(stringtowrite_c)# this for writing the fasta file that containg the UMI
                                        line="{},{},{},{},{}\n".format(str(tc),UMI_1,UMI_2,rb,Sample_condition)
                                        file_a.write(line)                                
                                        tc+=1

            elif Find_neighbour(ref_list,item_2[6:14],6,1):# When there is no exact match 
                ts1=item_2[4:16]
                Test_a=find_near_matches(Find_neighbour(ref_list,item_2[6:14],6,1), ts1, max_l_dist=2)
                UMI_1=item_2[0:Test_a[0].start+4] 
                if len(UMI_1)<6:
                    UMI_1="N"*(6-len(UMI_1))+UMI_1
                if len(UMI_1)>10:
                    UMI_1=item_2[0:6]
                UMI_2=UMI_1[-6:]
                Sample_condition=M_index_List[Find_neighbour(ref_list,item_2[6:14],6,1)]
                if constant_2 in item_2:
                    start1=item_2.find(constant_2)+len(constant_2)
                    if constant_3 in item_2:
                        end1=item_2.find(constant_3)
                        UMI_2=item_2[end1+24-len(item_2):]  
                        if len(UMI_2)<6:
                            UMI_2=UMI_2+"N"*(6-len(UMI_2))   
                        if len(UMI_2)>10:
                            UMI_2=item_2[-6:]
                        UMI_2=UMI_2[0:6]
                        rb=str(item_2[start1:end1])
                        quality=record[start1:end1].letter_annotations["phred_quality"]
                        total_average_quality+=ST.mean(quality)
                        if len(rb)==20:
                            number_pass_length_filter+=1
                            if min(quality)>Q_cutoff:                            
                                stringtowrite_b="{},{}\n".format(rb,str(tc),)
                                file_b.write(stringtowrite_b)# this for writing the fasta file that containg the barcode
                                stringtowrite_c="{},{}\n".format(str(UMI_1+UMI_2),str(tc))
                                file_c.write(stringtowrite_c)# this for writing the fasta file that containg the UMI			
                                line="{},{},{},{},{}\n".format(str(tc),UMI_1,UMI_2,rb,Sample_condition)
                                file_a.write(line)                  
                                tc+=1
                    else:
                        ts3=item_2[-32:-8] #this is target region for finding match  
                        Test_c=find_near_matches(constant_3, ts3, max_l_dist=3)
                        if Test_c:
                            end1=-32+Test_c[0].start
                            UMI_2=item_2[end1+24:]
                            if len(UMI_2)<6:
                                UMI_2=UMI_2+"N"*(6-len(UMI_2))
                            if len(UMI_2)>10:
                                UMI_2=item_2[-6:]
                            UMI_2=UMI_2[0:6]
                            rb=str(item_2[start1:end1])
                            quality=record[start1:end1].letter_annotations["phred_quality"]
                            total_average_quality+=ST.mean(quality)
                            if len(rb)==20:
                                number_pass_length_filter+=1
                                if min(quality)>Q_cutoff:                                
                                    stringtowrite_b="{},{}\n".format(rb,str(tc),)
                                    file_b.write(stringtowrite_b)# this for writing the fasta file that containg the barcode
                                    stringtowrite_c="{},{}\n".format(str(UMI_1+UMI_2),str(tc))
                                    file_c.write(stringtowrite_c)# this for writing the fasta file that containg the UMI		
                                    line="{},{},{},{},{}\n".format(str(tc),UMI_1,UMI_2,rb,Sample_condition)
                                    file_a.write(line)                            
                                    tc+=1
                else:
                    ts2=item_2[12:37]  #this is target region for finding match
                    Test_b=find_near_matches(constant_2, ts2, max_l_dist=3)
                    if Test_b:
                        start1=12+Test_b[0].end
                        if constant_3 in item_2:
                            end1=item_2.find(constant_3)
                            UMI_2=item_2[end1+24-len(item_2):]
                            if len(UMI_2)<6:
                                UMI_2=UMI_2+"N"*(6-len(UMI_2))
                            if len(UMI_2)>10:
                                UMI_2=item_2[-6:]
                            UMI_2=UMI_2[0:6]
                            rb=str(item_2[start1:end1])
                            quality=record[start1:end1].letter_annotations["phred_quality"]
                            total_average_quality+=ST.mean(quality)
                            if len(rb)==20:
                                number_pass_length_filter+=1
                                if min(quality)>Q_cutoff:                                
                                    stringtowrite_b="{},{}\n".format(rb,str(tc),)
                                    file_b.write(stringtowrite_b)# this for writing the fasta file that containg the barcode
                                    stringtowrite_c="{},{}\n".format(str(UMI_1+UMI_2),str(tc))
                                    file_c.write(stringtowrite_c)# this for writing the fasta file that containg the UMI		
                                    line="{},{},{},{},{}\n".format(str(tc),UMI_1,UMI_2,rb,Sample_condition)
                                    file_a.write(line)                            
                                    tc+=1
                        else:
                            ts3=item_2[-32:-8] #this is target region for finding match  
                            Test_c=find_near_matches(constant_3, ts3, max_l_dist=3)
                            if Test_c:
                                end1=-32+Test_c[0].start
                                UMI_2=item_2[end1+24:]
                                if len(UMI_2)<6:
                                    UMI_2=UMI_2+"N"*(6-len(UMI_2))
                                if len(UMI_2)>10:
                                    UMI_2=item_2[-6:]
                                UMI_2=UMI_2[0:6]
                                rb=str(item_2[start1:end1])
                                quality=record[start1:end1].letter_annotations["phred_quality"]
                                total_average_quality+=ST.mean(quality)
                                if len(rb)==20:
                                    number_pass_length_filter+=1
                                    if min(quality)>Q_cutoff:                                    
                                        stringtowrite_b="{},{}\n".format(rb,str(tc),)
                                        file_b.write(stringtowrite_b)# this for writing the fasta file that containg the barcode
                                        stringtowrite_c="{},{}\n".format(str(UMI_1+UMI_2),str(tc))
                                        file_c.write(stringtowrite_c)# this for writing the fasta file that containg the UMI	
                                        line="{},{},{},{},{}\n".format(str(tc),UMI_1,UMI_2,rb,Sample_condition)
                                        file_a.write(line)

                                        tc+=1

file_a.close()
file_b.close()
file_c.close()

# In[22]:
# In[22]:
