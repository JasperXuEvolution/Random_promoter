#!/usr/bin/env python
# coding: utf-8

# In[17]:



import numpy as np
import pandas as pd
import random
from Bio import SeqIO
import copy
import sys
import argparse
import time
import ast 
import copy
import os
import itertools

# In[18]:


t = time.localtime()
current_time = time.strftime("%H:%M:%S", t)
print(current_time)


# In[ ]:


parser = argparse.ArgumentParser(description='A function to combine the promoter and barcode information')
parser.add_argument("--a1", required=True, help="This is the input file of barcode and sequence id")
parser.add_argument("--a2", required=True, help="This is the output file of bartender_barcode for 20 bp barcode")
parser.add_argument("--a3", required=True, help="This is the output file of bartender_cluster for 20 bp barcode")
parser.add_argument("--a4", required=True, help="This is the file containing the deduplexed promoter")
parser.add_argument("--a5", required=True, help="This is the cutoff for major promoter cluster")
parser.add_argument("--a6", required=True, help="This is the file containing the deduplexed barcode")
parser.add_argument("--a7", required=True, help="This is the file of cdhit result")
parser.add_argument("--o1", required=True, help="This is the name of output file")
parser.add_argument("--o2", required=True, help="This is the name of output file for barcode clustering")
parser.add_argument("--o3", required=True, help="This is the name of output file for promoter clustering")
args = parser.parse_args()


# In[4]:


seq_id_count=args.a1
barcode_bartender_barcode=args.a2
barcode_bartender_cluster=args.a3
promoter_info=args.a4
barcode_info=args.a6
promoter_cluster_cdhit=args.a7


output1=args.o1+"_shared"
output2=args.o1+"_Summary"
output3=args.o2+"_deduplex_final"
output4=args.o2+"_cluster_final"

output5=args.o3+"_deduplex_final"
output6=args.o3+"_cluster_final"

cutoff=float(args.a5)



# In[19]:


# seq_id_count="test_barcode"
# barcode_bartender_barcode="test_barcode_bartender_barcode.csv"
# barcode_bartender_cluster="test_barcode_bartender_cluster.csv"
# promoter_ref="why"
# output1="shared"
# output2="Final"
# cutoff=0.8


# In[20]:


seqID_to_count={} #The most original ID for a promoter+barcode combination, or I call it seq_ID, as key and count as the value
with open(seq_id_count,'r') as handler:
    temp=handler.readline().strip()
    while temp:
        temp_list=temp.split(',')
        seq_id=temp_list[1]
        seq_count=int(temp_list[2])
        seq=temp_list[0]
        seqID_to_count[seq_id]=seq_count
#         seq_to_seqID[seq]=seq_id
        temp=handler.readline().strip()

# In[21]:
df_1=pd.read_csv(barcode_info, sep=',',header=0)
seqID_to_deduplex_barcode_ID={} # the key is the seqID, the value is the corresponding barcode ID 
represent_barcode_seq_to_seqID={} # the key is the deduplex barcode sequence, the value is a list of seqID 
for index, row in df_1.iterrows():
    temp_list=[] # a list of seq ID
    for x in row['ID_list'].split(':'):
        temp_list.append(x)
        seqID_to_deduplex_barcode_ID[x]=row['Barcode_ID']
    represent_barcode_seq_to_seqID[row['Sequence']]=temp_list


b_cluster_to_seqID={}#barcode cluster ID as key and a list of seqID as value
b_cluster_to_repSeq={}#barcode cluster ID as key and the representative sequence as value
with open (barcode_bartender_cluster,'r') as handler:
    temp=handler.readline().strip()
    temp=handler.readline().strip()
    while temp:
        b_cluster=temp.split(',')[0]#This is the barcode cluster ID
        repSeq=temp.split(',')[1] #This is the central sequence
        b_cluster_to_repSeq[b_cluster]=repSeq
        temp=handler.readline().strip()
        
with open(barcode_bartender_barcode,'r') as handler:
    temp=handler.readline().strip()
    temp=handler.readline().strip()
    while temp:
        if temp.split(',')[2] in b_cluster_to_seqID.keys():
            b_cluster_to_seqID[temp.split(',')[2]]+=represent_barcode_seq_to_seqID.get(temp.split(',')[0])
        else:
            b_cluster_to_seqID[temp.split(',')[2]]=represent_barcode_seq_to_seqID.get(temp.split(',')[0])
        temp=handler.readline().strip()  

b_cluster_to_barcodeID={} #barcode cluster ID as key and deduplex barcode ID as value
b_cluster_to_total_count={}
for key,value in b_cluster_to_seqID.items():
    temp_count=0
    b_cluster_to_barcodeID[key]=[]
    for x in value:
        temp_bID=seqID_to_deduplex_barcode_ID.get(x)
        temp_count+=seqID_to_count.get(x)
        if temp_bID not in b_cluster_to_barcodeID[key]:
            b_cluster_to_barcodeID[key].append(temp_bID)
    b_cluster_to_total_count[key]=temp_count

barcodeID_to_b_cluster={}
for key,value in b_cluster_to_barcodeID.items():
    for x in value:
        barcodeID_to_b_cluster[x]=key



df_1["Barcode_cluster_ID"]=df_1['Barcode_ID'].apply(lambda x: barcodeID_to_b_cluster.get(x))
df_1["Cluster_total_count"]=df_1['Barcode_cluster_ID'].apply(lambda x: b_cluster_to_total_count.get(x))
df_1["Cluster_central_sequence"]=df_1['Barcode_cluster_ID'].apply(lambda x: b_cluster_to_repSeq.get(x))
df_1=df_1.sort_values(by=['Cluster_total_count', 'Barcode_cluster_ID'],ascending=False)
df_1=df_1[['Barcode_ID','Sequence','Total_count','Barcode_cluster_ID','Cluster_total_count','Cluster_central_sequence','ID_list']]
df_1.to_csv(output3,index=False)


data = {'Barcode_cluster_ID':list(b_cluster_to_total_count.keys()), 
        'Cluster_total_count':list(b_cluster_to_total_count.values())} 

df_2=pd.DataFrame(data) 
df_2['Barcode_ID']=df_2['Barcode_cluster_ID'].apply(lambda x: ':'.join([str(x) for x in b_cluster_to_barcodeID.get(x)]))
df_2["Cluster_central_sequence"]=df_2['Barcode_cluster_ID'].apply(lambda x: b_cluster_to_repSeq.get(x))
df_2=df_2.sort_values(by=['Cluster_total_count', 'Barcode_cluster_ID'], ascending=False)
df_2=df_2[['Barcode_cluster_ID','Cluster_total_count','Cluster_central_sequence','Barcode_ID']]
df_2.to_csv(output4,index=False)



del represent_barcode_seq_to_seqID, df_1, df_2, b_cluster_to_barcodeID, b_cluster_to_total_count, barcodeID_to_b_cluster


# In[22]:


print("Barcode dictionary is finished.")
t = time.localtime()
current_time = time.strftime("%H:%M:%S", t)
print(current_time)


# In[23]:


df_3=pd.read_csv(promoter_info, sep=',',header=0)# read the promoter information

Promoter_ID_to_seq={} # the promoter deduplex ID as key and the [sequence,count] as value 
Promoter_ID_to_seqid={} # the promoter deduplex ID as key, a list of orginal ID as value 
for index, row in df_3.iterrows():
    temp_ID=row['Promoter_ID']
    temp_sequence=row['Sequence']
    temp_count=row['Total_count']
    Promoter_ID_to_seq[temp_ID]=[temp_sequence,int(temp_count)]
    Promoter_ID_to_seqid[temp_ID]=row['ID_list'].split(':')



# In[64]:


p_cluster_file=open(promoter_cluster_cdhit,'r')
p_cluster_to_seqID={} #promoter cluster ID as key and a list of seqID as value
p_cluster_to_repSeq={} #promoter cluster ID as key and the representative sequence as value
p_cluster_to_size={} #promoter cluster ID as key and the cluster size as value

Promoter_cluster_summary={}#this is summary dict for promoter_cluster, key is the promoter cluster ID

# parse through the cluster file and store the cluster name + sequences in the dictionary
p_cluster_groups = (x[1] for x in itertools.groupby(p_cluster_file, key=lambda line: line[0] == '>'))
for cluster in p_cluster_groups:
    seqs=[] #promoter_deduplex_ID
    main_seq=[]
    name = str(int(cluster.__next__().strip().split(' ')[1])+1)
    for seq in p_cluster_groups.__next__():
        temp=seq.split('>')[1].split('...')[0]
        seqs.append(temp)
    total_count=0 #This for recording the total count of seq in this promoter cluster
    ID_list=[] # a list of all the original ID
    largest_cout=0# This is for recording the largest count
    representative_seq=None
    temp_count=0
    for seq_id in seqs:
        if Promoter_ID_to_seq.get(seq_id):
            temp_count=Promoter_ID_to_seq.get(seq_id)[1]
            ID_list+=Promoter_ID_to_seqid.get(seq_id)
            if temp_count>largest_cout:
                largest_cout=temp_count
                representative_seq=Promoter_ID_to_seq.get(seq_id)[0]
            total_count+=temp_count
    if total_count!=0:
        ratio=round(largest_cout/total_count,4)
        Promoter_cluster_summary[name]=[representative_seq,total_count,ratio,ID_list,seqs]  
        p_cluster_to_seqID[name]=ID_list
        p_cluster_to_repSeq[name]=representative_seq
        p_cluster_to_size[name]=total_count

seqID_to_p_cluster={} # seqID as key and promoter cluster ID as value

for key,values in Promoter_cluster_summary.items():
    for x in values[3]:
        seqID_to_p_cluster[x]=key

promoterID_to_p_cluster={} #promoter ID as key, promoter cluster ID as value
for key,value in Promoter_cluster_summary.items():
    for x in value[4]:
        promoterID_to_p_cluster[x]=key
        
        
df_3["Promoter_cluster_ID"]=df_3['Promoter_ID'].apply(lambda x: promoterID_to_p_cluster.get(x))
df_3["Cluster_total_count"]=df_3['Promoter_cluster_ID'].apply(lambda x: p_cluster_to_size.get(x))
df_3["Cluster_central_sequence"]=df_3['Promoter_cluster_ID'].apply(lambda x: p_cluster_to_repSeq.get(x))
df_3=df_3.sort_values(by=['Cluster_total_count', 'Promoter_cluster_ID'],ascending=False)
df_3=df_3[['Promoter_ID','Sequence','Total_count','Promoter_cluster_ID','Cluster_total_count','Cluster_central_sequence','ID_list']]
df_3.to_csv(output5,index=False)

temp_list=[x[4] for x in list(Promoter_cluster_summary.values())]
ll_new=[]
for x in temp_list:
    ll_new.append(':'.join([str(y) for y in x]))
data = {'Promoter_cluster_ID':list(Promoter_cluster_summary.keys()), 
        'Cluster_total_count':[x[1] for x in list(Promoter_cluster_summary.values())],
        "Cluster_central_sequence":[x[0] for x in list(Promoter_cluster_summary.values())],
        "Major ratio":[x[2] for x in list(Promoter_cluster_summary.values())],
        'Promoter_ID':ll_new,
        } 
df_4=pd.DataFrame(data) 
df_4=df_4.sort_values(by=['Cluster_total_count', 'Promoter_cluster_ID'], ascending=False)
df_4.to_csv(output6,index=False)

del Promoter_cluster_summary, df_3, df_4, Promoter_ID_to_seq,Promoter_ID_to_seqid



print("promoter dictionary is finished.")
t = time.localtime()
current_time = time.strftime("%H:%M:%S", t)
print(current_time)


# In[26]:


#This is a function of calculating hamming distance
def hamming_distance(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


# In[27]:


file_a=open(output1,'w')
file_a.write('Barcode_cluster,Promoter_cluster:count\n')
b_cluster_to_seqID_copy=copy.deepcopy(b_cluster_to_seqID)
for key,value in b_cluster_to_seqID_copy.items():
    temp_dic={} #This is the temporary dictionary with promoter cluster ID as key, [count],[seqID]
    for x in value: #loop over all the seqID stored in the value
        temp_p_cluster=seqID_to_p_cluster[x]
        if temp_p_cluster in temp_dic.keys(): #store the promoter cluster name as key
            temp_dic[temp_p_cluster][0]+=seqID_to_count[x]
            temp_dic[temp_p_cluster][1].append(x) #store the seqID as the second items in value
        else:
            s1=p_cluster_to_repSeq[temp_p_cluster]
            r_record=0 # This is the counting machine
            for p_cluster in temp_dic.keys():
                s2=p_cluster_to_repSeq[p_cluster]
                if hamming_distance(s1[-117:], s2[-117:])<20: #if the hamming distance between a existing prmoter cluster and a new promoter clustser is less than 20.
                    c1=p_cluster_to_size[p_cluster]
                    c2=p_cluster_to_size[temp_p_cluster]
                    if c1>c2:
                        #the new cluster is merged to the order cluster 
                        temp_dic[p_cluster][0]+=seqID_to_count[x]
                        temp_dic[p_cluster][1].append(x)
                        # i update the promoter cluster dictionary
                        p_cluster_to_seqID[p_cluster]+=p_cluster_to_seqID[temp_p_cluster] #merge this p_cluster into existing one 
                        for l in p_cluster_to_seqID[temp_p_cluster]:
                            seqID_to_p_cluster[l]=p_cluster
                        del p_cluster_to_seqID[temp_p_cluster] #delete the p_cluster
                    else:
                        #print ("work")
                        #print (p_cluster)
                        #print (x)
                        #the old cluster is merged to the new cluster 
                        temp_dic[temp_p_cluster]=[seqID_to_count[x],[x]]
                        temp_dic[temp_p_cluster][0]+=temp_dic[p_cluster][0]
                        temp_dic[temp_p_cluster][1]+=temp_dic[p_cluster][1]
                        del temp_dic[p_cluster]
                        p_cluster_to_seqID[temp_p_cluster]+=p_cluster_to_seqID[p_cluster]
                        for l in p_cluster_to_seqID[p_cluster]:
                            seqID_to_p_cluster[l]=temp_p_cluster
                        del p_cluster_to_seqID[p_cluster]
                    r_record=1
                    break
            if r_record==0: #if not within 20 distance, then there is case of multiple promoter with same barcode 
                temp_dic[temp_p_cluster]=[seqID_to_count[x],[x]]
    if len(temp_dic.keys())>1:
        count_list=[a[0] for a in temp_dic.values()]
        p_cluster_list=[x for x in temp_dic.keys()]
        position=count_list.index(max(count_list))
        main_p_cluster= p_cluster_list[position]
        if max(count_list)<cutoff*sum(count_list):# if there is no major promoter with more than 80%, I will record this shared barcode event and delete this barcode in the original dictionary
            file_a.write(key+',')
            for key1,value1 in temp_dic.items():
                stringtowrite="{}:{};".format(key1,value1[0])
                file_a.write(stringtowrite)
            del b_cluster_to_seqID[key]
            file_a.write('\n')
        else: #if there is a major promoter for this barcode. I will delete the seqID corresponding to minor promoter in barcode ditionary
            a=set(b_cluster_to_seqID[key])
            b=set(temp_dic[main_p_cluster][1])
            overlap=[x for x in (a&b)]
            b_cluster_to_seqID[key]=overlap


# In[28]:


file_a=open(output2,'w')
promoter_set=set() #This is for recording promoter set
nc=0 #count total number of barcode promoter combination
stringtowrite="{}:{}:{}:{}:{}:{}\n".format('P_Cluster_ID','B_Clustser_ID','P_Cluster_seq','B_Cluster_seq','Total_count','seq.ID')
file_a.write(stringtowrite)
for key,value in b_cluster_to_seqID.items():
    nc+=1
    temp_dict={} 
    total_count=0#this is a dictionary to store promoter cluster as key, and number of count as value
    for x in value: #loop over all the seqID stored in the value
        temp_count=seqID_to_count[x]
        total_count+=temp_count
    pn=seqID_to_p_cluster[value[0]]#promoter cluster name
    promoter_set.update({pn})
    bn=key#barcode cluster name
    ps=p_cluster_to_repSeq[pn] #promoter sequence
    bs=b_cluster_to_repSeq[bn] #barcode sequence
    stringtowrite="{},{},{},{},{},{}\n".format(pn,bn,ps,bs,total_count,':'.join([str(x) for x in value]))
    file_a.write(stringtowrite)
file_a.close()
print ("There are totally %i promoter clusters after merging data."%(len(promoter_set)))
print ("There are totally %i promoter and barcode combination after merging data."%(nc))


# In[29]:


t = time.localtime()
current_time = time.strftime("%H:%M:%S", t)
print(current_time)


# In[ ]:




