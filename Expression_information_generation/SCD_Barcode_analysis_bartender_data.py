#!/usr/bin/env python
# coding: utf-8

# In[13]:


#Here I import all the package that I potentially will use to analyse my data
import sys
import statistics as ST
import argparse
import numpy


# In[14]:





parser = argparse.ArgumentParser(description='A function to analyse extracted barcode by its condition and UMI')
parser.add_argument("--a1", required=True, help="This is the input file of barcode cluster")
parser.add_argument("--a2", required=True, help="This is the input file of representative barcode")
parser.add_argument("--a3", required=True, help="This is the input file of unique UMI and its corresponding cluster")
parser.add_argument("--a4", required=True, help="This is the input file of representative UMI")
parser.add_argument("--a5", required=True, help="This is the input file of information summary")
parser.add_argument("--a6", required=True, help="This is the input barcode cluster ID")
parser.add_argument("--o", required=True, help="This is the name of output file")
args = parser.parse_args()
barcode_unique=args.a1
barcode_central=args.a2
UMI_unique=args.a3
UMI_central=args.a4
info_summary=args.a5
Seq_ID=args.a6
output_name=args.o
output1=info_summary+"_UMI"
output2=output_name+'barcode_cluster_final'

UMI_cluster_file={}
# parse through the cluster file and store the cluster name + sequences in the dictionary
with open(UMI_unique,'r') as handler:
    temp=handler.readline().strip()
    temp=handler.readline().strip()
    while temp:
        UMI=temp.split(',')[0] #record the sequence of unique UMI 
        Cluster_ID=temp.split(',')[2] #record cluster ID
        UMI_cluster_file[UMI]=Cluster_ID
        temp=handler.readline().strip()

#this list store the the representative barcode information
UMI_central_sequence=[]
with open(UMI_central,'r') as handler:
    temp=handler.readline().strip()
    temp=handler.readline().strip()
    while temp:
        UMI_central_sequence.append(temp.split(',')[1])
        temp=handler.readline().strip()

        
limit_N=0
file_x=open(output1,'w') #I generate a new file with UMI identified 
with open(info_summary,'r') as handler:
    temp=handler.readline().strip()
    while temp:
        cl=UMI_cluster_file.get(temp.split(',')[1]+temp.split(',')[2]) 
        central_UMI=UMI_central_sequence[int(cl)-1]#This is the sequence of central barcode
        stringtowrite="{},{},{}\n".format(temp,cl,central_UMI)
        file_x.write(stringtowrite)
        limit_N+=1
        temp=handler.readline().strip()  
file_x.close()
del UMI_cluster_file
del UMI_central_sequence
print ("Finished UMI")

barcode_to_ID_dic={}
with open(Seq_ID,'r') as handler:
    temp=handler.readline().strip()
    while temp:
        if temp.split(',')[0] in barcode_to_ID_dic.keys(): #using barcode sequence as key, the a list of read ID as value
            barcode_to_ID_dic[temp.split(',')[0]].append(int(temp.split(',')[1]))
        else:
            barcode_to_ID_dic[temp.split(',')[0]]=[int(temp.split(',')[1])]
        temp=handler.readline().strip()
        
        
        
b_cluster_dic={} #creat a dict, key is the cluster ID and the value are a lsit of read ID
with open(barcode_unique,'r') as handler:
    temp=handler.readline().strip()
    temp=handler.readline().strip()
    while temp:
        if temp.split(',')[2] in b_cluster_dic.keys():
            b_cluster_dic[temp.split(',')[2]]+=barcode_to_ID_dic.get(temp.split(',')[0])
        else:
            b_cluster_dic[temp.split(',')[2]]=barcode_to_ID_dic.get(temp.split(',')[0])
        temp=handler.readline().strip()
del barcode_to_ID_dic


#this step add the central sequence to the dictionary 
file_a=open(output2,'w')
with open (barcode_central,'r') as handler:
    temp=handler.readline().strip() 
    stringtowrite="{},{}\n".format(temp,'read_ID')
    file_a.write(stringtowrite)
    temp=handler.readline().strip() 
    while temp:
        seq_id_list=b_cluster_dic.get(temp.split(',')[0])
        stringtowrite="{};{}\n".format(temp,seq_id_list)
        file_a.write(stringtowrite)
        temp=handler.readline().strip() 
file_a.close()
del b_cluster_dic
print ("barcode_cluster_final is generated")




types='U7,int32'
my_data = numpy.zeros(shape=(limit_N,1),dtype=types)
n=0
with open(output1,'r') as handler:
    line=handler.readline().strip()
    while line:
        temp=line.split(',')
        umi_number=temp[5]
        ar=(temp[4],umi_number)
        my_data[n,]=ar
        n+=1
        line=handler.readline().strip()
print ("reading data_frame finished")

Condition_name=['DNA_S1','DNA_S2','DNA_S3','RNA_S1','RNA_S2','RNA_S3']
o1=output_name+'_DNA_S1.csv'
o2=output_name+'_DNA_S2.csv'
o3=output_name+'_DNA_S3.csv'
o4=output_name+'_RNA_S1.csv'
o5=output_name+'_RNA_S2.csv'
o6=output_name+'_RNA_S3.csv'

o7=output_name+'_S1_b.csv'
o8=output_name+'_S1_d.csv'
o9=output_name+'_S1_r.csv'

o10=output_name+'_S2_b.csv'
o11=output_name+'_S2_d.csv'
o12=output_name+'_S2_r.csv'

o13=output_name+'_S3_b.csv'
o14=output_name+'_S3_d.csv'
o15=output_name+'_S3_r.csv'

o16=output_name+'_long'+'_DNA_S1.csv'
o17=output_name+'_long'+'_DNA_S2.csv'
o18=output_name+'_long'+'_DNA_S3.csv'
o19=output_name+'_long'+'_RNA_S1.csv'
o20=output_name+'_long'+'_RNA_S2.csv'
o21=output_name+'_long'+'_RNA_S3.csv'




file_a=open(o1,'w')
file_b=open(o2,'w')
file_c=open(o3,'w')
file_d=open(o4,'w')
file_e=open(o5,'w')
file_f=open(o6,'w')
Summary_list=[file_a,file_b,file_c,file_d,file_e,file_f]
file_al=open(o16,'w')
file_bl=open(o17,'w')
file_cl=open(o18,'w')
file_dl=open(o19,'w')
file_el=open(o20,'w')
file_fl=open(o21,'w')
Summary_list_l=[file_al,file_bl,file_cl,file_dl,file_el,file_fl]



for l in range(0,6):
        stringtowrite="{};{};{};{};{};{}\n".format('barcode_cluster','absolute_count','ID_for_sequence_in_cluster','relative_count','UMI_cluster_names','representative_sequence')
        # semicolon as delimiter because I have comma in the list
        Summary_list_l[l].write(stringtowrite) 
        stringtowrite="{};{};{}\n".format('barcode_cluster','absolute_count','relative_count')
        Summary_list[l].write(stringtowrite)
file_x_b=open(o7,'w') #This is the file storing barcode in sample one that has both DNA reads and RNA reads for sample1
file_x_d=open(o8,'w') #This is the file storing barcode in sample one that has both DNA reads only for sample1
file_x_r=open(o9,'w') #This is the file storing barcode in sample one that has both RNA reads only for sample1
file_y_b=open(o10,'w') #This is the file storing barcode in sample one that has both DNA reads and RNA reads for sample2
file_y_d=open(o11,'w') #This is the file storing barcode in sample one that has both DNA reads only for sample2
file_y_r=open(o12,'w') #This is the file storing barcode in sample one that has both RNA reads only for sample2
file_z_b=open(o13,'w') #This is the file storing barcode in sample one that has both DNA reads and RNA reads for sample3
file_z_d=open(o14,'w') #This is the file storing barcode in sample one that has both DNA reads only for sample3
file_z_r=open(o15,'w') #This is the file storing barcode in sample one that has both RNA reads only for sample3
both_list=[file_x_b,file_y_b,file_z_b]    
d_list=[file_x_d,file_y_d,file_z_d]
r_list=[file_x_r,file_y_r,file_z_r]
for l in range (0,3):
    stringtowrite="{},{},{},{},{},{},{},{}\n".format('barcode_cluster','representative_sequence','DNA_absolute_count','DNA_relative_count','RNA_absolute_count','RNA_relative_count','absolute_ratio','relative_ratio')#barcode cluster name,absolute count,all the id in this barcode cluster,
    both_list[l].write(stringtowrite)
    stringtowrite="{},{},{},{}\n".format('barcode_cluster','representative_sequence','DNA_absolute_count','DNA_relative_count')
    d_list[l].write(stringtowrite)
    stringtowrite="{},{},{},{}\n".format('barcode_cluster','representative_sequence','RNA_absolute_count','RNA_relative_count')
    r_list[l].write(stringtowrite)
with open(output2,'r') as handler:
    temp=handler.readline().strip()
    temp=handler.readline().strip()
    while temp:
        name=temp.split(',')[0] #This is the barcode cluster name
        ID=temp.split(';')[1].split('[')[1].split(']')[0].split(',')
        central_seq=temp.split(',')[1]
        dic1={} #This is a dictionary with sample condition as key, and correponding id as value
        for i in ID:
            if my_data[int(i)-1,][0][0] in dic1.keys():
                dic1[my_data[int(i)-1,][0][0]][1]+=1
                dic1[my_data[int(i)-1,][0][0]][0].append(i)
            else:
                dic1[my_data[int(i)-1,][0][0]]=[[i],1]
        for key,value in dic1.items():
            dic2={} # This is a dictionary with UMI cluster name as key.
            for sub_id in value[0]:
                if my_data[int(sub_id)-1,][0][1] in dic2.keys():
                    dic2[my_data[int(sub_id)-1,][0][1]]+=1
                else:
                    dic2[my_data[int(sub_id)-1,][0][1]]=1
            value.append(list(dic2.keys()))
            value.append(len(value[2]))
        for j in range(0,6):
            if dic1.get(Condition_name[j]):
                stringtowrite="{};{};{};{};{};{}\n".format(name,dic1[Condition_name[j]][1],dic1[Condition_name[j]][0],dic1[Condition_name[j]][3],dic1[Condition_name[j]][2],central_seq)#barcode cluster name,absolute count,all the id in this barcode cluster,
                # semicolon as delimiter because I have comma in the list
                Summary_list_l[j].write(stringtowrite) 
                stringtowrite="{};{};{}\n".format(name,dic1[Condition_name[j]][1],dic1[Condition_name[j]][3])#barcode cluster name,absolute count,all the id in this barcode cluster,
                Summary_list[j].write(stringtowrite) 
        for k in range(0,3):
            if dic1.get(Condition_name[k]) and dic1.get(Condition_name[k+3]): #both reads for this barcode presence
                ratio_1=dic1[Condition_name[k+3]][1]/dic1[Condition_name[k]][1]
                ratio_2=dic1[Condition_name[k+3]][3]/dic1[Condition_name[k]][3]
                stringtowrite="{},{},{},{},{},{},{},{}\n".format(name,central_seq,dic1[Condition_name[k]][1],dic1[Condition_name[k]][3],\
                                                          dic1[Condition_name[k+3]][1],dic1[Condition_name[k+3]][3],ratio_1,ratio_2)
                #barcode cluster name,cluster central sequence, absolute count for DNA, relative count for DNA, absoulte count for RNA, relative count for RNA, ratio based on absolute count, ratio based on relative count
                both_list[k].write(stringtowrite)
            if dic1.get(Condition_name[k]) and (not dic1.get(Condition_name[k+3])):
                stringtowrite="{},{},{},{}\n".format(name,central_seq,dic1[Condition_name[k]][1],dic1[Condition_name[k]][3])
                #barcode cluster name,cluster central sequence, absolute count for DNA, relative count for DNA, 
                d_list[k].write(stringtowrite)
            if not dic1.get(Condition_name[k]) and dic1.get(Condition_name[k+3]):
                stringtowrite="{},{},{},{}\n".format(name,central_seq,dic1[Condition_name[k+3]][1],dic1[Condition_name[k+3]][3])
                #barcode cluster name,cluster central sequence, absolute count for DNA, relative count for DNA, 
                r_list[k].write(stringtowrite)
        temp=handler.readline().strip()
file_a.close()
file_b.close()
file_c.close()
file_d.close()
file_e.close()
file_f.close()
file_x_b.close() 
file_x_d.close() 
file_x_r.close() 
file_y_b.close() 
file_y_d.close() 
file_y_r.close() 
file_z_b.close() 
file_z_d.close() 
file_z_r.close() 
file_al.close()
file_bl.close()
file_cl.close()
file_dl.close()
file_el.close()
file_fl.close()
