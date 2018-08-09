# -*- coding: utf-8 -*-
"""
last edit: 06/22/17
Created on Fri Oct 21 15:51:46 2016

transform read counts to a table
rows: mutants
columns: different conditions (plasmid, transfection, selection conditions)
mutants not observed: assign count=0

input: 
count_aa (output file from Mapper4.py)
sequence_list_aa (output file from Mapper3.py)
sequence_list_number (output file from Mapper3.py)

output:
count_aa_table

@author: lei
"""
#oct 2016
#samples = ['plasmid','transfection','TPV0_D3','TPV100_D3','TPV316_D3','TPV1000_D3','TPV0_D9','TPV100_D9','TPV316_D9','TPV1000_D9']
#infile2 = open('/mnt/hgfs/analysis_sunlab/analysis_protease/count_aa','r')  
#outfile = open('/mnt/hgfs/analysis_sunlab/analysis_protease/count_aa_table','w')
#dec 2016 TPV
#samples = ['transfection','TPV0_D7','TPV100_D7','TPV316_D7','TPV1000_D7']
#infile2 = open('/mnt/hgfs/analysis_sunlab/analysis_protease/dec2016_TPV/count_aa','r')  
#outfile = open('/mnt/hgfs/analysis_sunlab/analysis_protease/dec2016_TPV/count_aa_table','w')
#dec 2016 DRV
#samples = ['DRV0_D3','DRV3_D3','DRV10_D3','DRV100_D3','DRV0_D6','DRV3_D6','DRV10_D6','DRV100_D6','DRV0_D9','DRV3_D9','DRV10_D9','DRV100_D9']
#infile2 = open('/mnt/hgfs/analysis_sunlab/analysis_protease/dec2016_DRV/count_aa','r')  
#outfile = open('/mnt/hgfs/analysis_sunlab/analysis_protease/dec2016_DRV/count_aa_table','w')
#June 2017
samples = ['tr_R1','tr_R2','tr_R3','nd_R1','nd_R2','nd_R3','TPV_R1','TPV_R2','TPV_R3','DRV_R1','DRV_R2','DRV_R3','plasmid']
infile2 = open('/mnt/hgfs/analysis_sunlab/analysis_HIV_protease/june2017/count_aa','r')  
outfile = open('/mnt/hgfs/analysis_sunlab/analysis_HIV_protease/june2017/count_aa_table','w')

#list of sequences
infile =  open('/mnt/hgfs/analysis_sunlab/analysis_HIV_protease/sequence_list_aa','r')
sequence_list_aa=[] 
for line in infile:
    genotype = line.strip()
    sequence_list_aa.append(genotype)    
infile3 =  open('/mnt/hgfs/analysis_sunlab/analysis_HIV_protease/sequence_list_number','r')
sequence_list_number=[] 
for line in infile3:
    genotype_number = line.strip()
    sequence_list_number.append(genotype_number)   

#iterate the list of combinatorial mutants
#if a genotype is found in count_aa, assign read counts; if not found, count=0

count_table={}
for genotype in sequence_list_aa:
    count_table[genotype] = {}
    for sample in samples:
        count_table[genotype][sample] = 0
        
#read file: count_aa             
for line in infile2:
    [sample, muts_aa, count] = line.strip().split('\t')
    count_table[muts_aa][sample] = count

#write to file
#header
header = ['sequence_aa','sequence_num']
for sample in samples:
    header.extend([sample]) 
outfile.write("\t".join(header)+"\n")

for i in range(0,len(sequence_list_aa)):
    genotype = sequence_list_aa[i]
    genotype_number = sequence_list_number[i]
    out = [genotype, genotype_number]
    for sample in samples:
        out.append(str(count_table[genotype][sample]))
    outfile.write("\t".join(out)+"\n")

infile.close()
infile2.close()
outfile.close()