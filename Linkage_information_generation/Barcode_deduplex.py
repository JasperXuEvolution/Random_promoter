#!/usr/bin/env python
# coding: utf-8

# In[18]:


import sys
import itertools
import argparse
import statistics as ST
import time
import os
import pandas as pd


# In[ ]:


parser = argparse.ArgumentParser(description='A function to extract barcode information!')
parser.add_argument("--a", required=True, help="This is the input file of reads")
parser.add_argument("--b", required=True, help="This is the name of output file for recording")
parser.add_argument("--c", required=True, help="This is the name of output file for cdhit")

args = parser.parse_args()
input_file=args.a
output_file_a=args.b
output_file_b=args.c


# In[21]:


df=pd.read_csv(input_file, sep=',',header=None)

# In[23]:

file_b=open(output_file_b,'w') #This is the output for cdhit
barcode_dic={}
#This is a dictionary, key is the barcode sequence and value is a list of two item. First is an integer of total count of this barcode and second is a list of ID 
for index, row in df.iterrows():
    temp_barcode=row[0]
    temp_ID=row[1]
    temp_count=row[2]
    for length in range(0,temp_count): #I duplex sequnece so it will be used a input for bartender
        string_to_write="{},{}\n".format(temp_barcode,str(temp_ID)+"_"+str(length))
        file_b.write(string_to_write)
    if temp_barcode in barcode_dic.keys():
        barcode_dic[temp_barcode][0]+=temp_count
        barcode_dic[temp_barcode][1].append(temp_ID)
    else:
        barcode_dic[temp_barcode]=[temp_count,[temp_ID]]


# ### Output for recording and cdhit

# In[40]:


# In[41]:


file_a=open(output_file_a,'w') #This is the output for recording
string_to_write="Barcode_ID,Sequence,Total_count,ID_list\n"
file_a.write(string_to_write)



# In[42]:


tc=1
for keys,values in barcode_dic.items():
    string_to_write="{},{},{},{}\n".format('B'+str(tc),keys,str(values[0]),':'.join([str(x) for x in values[1]]))
    file_a.write(string_to_write)
    tc+=1



# In[43]:


file_a.close()
file_b.close()

