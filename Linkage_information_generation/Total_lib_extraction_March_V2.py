#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
from Bio import SeqIO
from fuzzysearch import find_near_matches
import itertools
from Bio import Align 
import argparse
import statistics as ST


# In[ ]:


import time
t = time.localtime()
current_time = time.strftime("%H:%M:%S", t)
print(current_time)


# In[98]:


parser = argparse.ArgumentParser(description='A function to extract promoter and barcode!')
parser.add_argument("--a", required=True, help="This is the input file of reads")
parser.add_argument("--o", required=True, help="This is the name of output file")
parser.add_argument("--q", required=True, help="This is quality cutoff")

args = parser.parse_args()
input_file=args.a
output_name=args.o
Q_cutoff=int(args.q)


# In[161]:


constant_1="TAGCGCGACGTTTTCTTCCT"
constant_2="TGGCTGAACCTAGTTTTGCCC"
constant_3="CTCTCGGAACATAGCAGTTT"


# In[267]:


#input_file="test2.txt"


# record.letter_annotations["phred_quality"]

# In[272]:


#end2


# In[268]:


#Q_cutoff=10
Alligned_Seq ={}
for record in SeqIO.parse(input_file, "fastq"):
    x=record.seq
    if constant_1 in x:
        start1=x.find(constant_1)+len(constant_1)
        if constant_2 in x:
            end1=x.find(constant_2)
            start2=end1+len(constant_2)
            if constant_3 in x:
                end2=x.find(constant_3)
                rp=x[start1:end1]
                rb=x[start2:end2]
                quality=record[start2:end2].letter_annotations["phred_quality"]
                if (len(rp) ==120) and (len(rb)==20): #trimmed according to the length of barcode and promoter
                    if min(quality)>Q_cutoff:#the quality cutoff
                        if rp+rb in Alligned_Seq.keys():
                            Alligned_Seq[rp+rb]+=1
                        else:
                            Alligned_Seq[rp+rb]=1
            else:
                ts3=x[-34:-6] #this is target region for finding match  
                Test_c=find_near_matches(constant_3, ts3, max_l_dist=3)
                if Test_c:
                    end2=-34+Test_c[0].start
                    rp=x[start1:end1]
                    rb=x[start2:end2]
                    quality=record[start2:end2].letter_annotations["phred_quality"]
                    if (len(rp) ==120) and (len(rb)==20): #trimmed according to the length of barcode and promoter
                        if min(quality)>Q_cutoff:#the quality cutoff
                            if rp+rb in Alligned_Seq.keys():
                                Alligned_Seq[rp+rb]+=1
                            else:
                                Alligned_Seq[rp+rb]=1

        else:
            ts2=x[-75:-46]  #this is target region for finding match
            Test_b=find_near_matches(constant_2, ts2, max_l_dist=3)
            if Test_b:
                end1=-75+Test_b[0].start
                start2=-75+Test_b[0].end
                if constant_3 in x:
                    end2=x.find(constant_3)
                    rp=x[start1:end1]
                    rb=x[start2:end2]
                    quality=record[start2:end2].letter_annotations["phred_quality"]
                    if (len(rp) ==120) and (len(rb)==20): #trimmed according to the length of barcode and promoter
                        if min(quality)>Q_cutoff:#the quality cutoff
                            if rp+rb in Alligned_Seq.keys():
                                Alligned_Seq[rp+rb]+=1
                            else:
                                Alligned_Seq[rp+rb]=1

                else:
                    ts3=x[-34:-6] #this is target region for finding match  
                    Test_c=find_near_matches(constant_3, ts3, max_l_dist=3)
                    if Test_c:
                        end2=-34+Test_c[0].start
                        rp=x[start1:end1]
                        rb=x[start2:end2]
                        quality=record[start2:end2].letter_annotations["phred_quality"]
                        if (len(rp) ==120) and (len(rb)==20): #trimmed according to the length of barcode and promoter
                            if min(quality)>Q_cutoff:#the quality cutoff
                                if rp+rb in Alligned_Seq.keys():
                                    Alligned_Seq[rp+rb]+=1
                                else:
                                    Alligned_Seq[rp+rb]=1
    else:
        ts1=x[5:33]
        Test_a=find_near_matches(constant_1, ts1, max_l_dist=3)
        if Test_a:
            start1=Test_a[0].end+5
            if constant_2 in x:
                end1=x.find(constant_2)
                start2=end1+len(constant_2)
                if constant_3 in x:
                    end2=x.find(constant_3)
                    rp=x[start1:end1]
                    rb=x[start2:end2]
                    quality=record[start2:end2].letter_annotations["phred_quality"]
                    if (len(rp) ==120) and (len(rb)==20): #trimmed according to the length of barcode and promoter
                        if min(quality)>Q_cutoff:#the quality cutoff
                            if rp+rb in Alligned_Seq.keys():
                                Alligned_Seq[rp+rb]+=1
                            else:
                                Alligned_Seq[rp+rb]=1

                else:
                    ts3=x[-34:-6] #this is target region for finding match  
                    Test_c=find_near_matches(constant_3, ts3, max_l_dist=3)
                    if Test_c:
                        end2=-34+Test_c[0].start
                        rp=x[start1:end1]
                        rb=x[start2:end2]
                        quality=record[start2:end2].letter_annotations["phred_quality"]
                        if (len(rp) ==120) and (len(rb)==20): #trimmed according to the length of barcode and promoter
                            if min(quality)>Q_cutoff:#the quality cutoff
                                if rp+rb in Alligned_Seq.keys():
                                    Alligned_Seq[rp+rb]+=1
                                else:
                                    Alligned_Seq[rp+rb]=1

            else:
                ts2=x[-75:-46]  #this is target region for finding match
                Test_b=find_near_matches(constant_2, ts2, max_l_dist=3)
                if Test_b:
                    end1=-75+Test_b[0].start
                    start2=-75+Test_b[0].end
                    if constant_3 in x:
                        end2=x.find(constant_3)
                        rp=x[start1:end1]
                        rb=x[start2:end2]
                        quality=record[start2:end2].letter_annotations["phred_quality"]
                        if (len(rp) ==120) and (len(rb)==20): #trimmed according to the length of barcode and promoter
                            if min(quality)>Q_cutoff:#the quality cutoff
                                if rp+rb in Alligned_Seq.keys():
                                    Alligned_Seq[rp+rb]+=1
                                else:
                                    Alligned_Seq[rp+rb]=1
                    else:
                        ts3=x[-34:-6] #this is target region for finding match  
                        Test_c=find_near_matches(constant_3, ts3, max_l_dist=3)
                        if Test_c:
                            end2=-34+Test_c[0].start
                            rp=x[start1:end1]
                            rb=x[start2:end2]
                            quality=record[start2:end2].letter_annotations["phred_quality"]
                            if (len(rp) ==120) and (len(rb)==20): #trimmed according to the length of barcode and promoter
                                if min(quality)>Q_cutoff:#the quality cutoff
                                    if rp+rb in Alligned_Seq.keys():
                                        Alligned_Seq[rp+rb]+=1
                                    else:
                                        Alligned_Seq[rp+rb]=1


# In[266]:


#len(Alligned_Seq.keys())


# In[258]:


#output_name="whynot"


# In[259]:


ouputname_p=output_name+"_promoter"
ouputname_b=output_name+"_barcode"


# In[260]:


file_p = open(ouputname_p,'w')
file_b = open(ouputname_b,'w')
tc=1
for key, value in Alligned_Seq.items():
    stringtowrite_b="{},{},{}\n".format(str(key)[-20:],str(tc),str(value))
    file_b.write(stringtowrite_b)
    stringtowrite_p="{},{},{}\n".format(str(key)[:-20],str(tc),str(value))
    file_p.write(stringtowrite_p)
    tc+=1
file_p.close()
file_b.close()


# In[ ]:
t = time.localtime()
current_time = time.strftime("%H:%M:%S", t)
print(current_time)

