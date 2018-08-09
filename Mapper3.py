# -*- coding: utf-8 -*-
"""
last edit: 12/14/16
notes: no need to re-run this script for new data

Created on Thu Oct 20 16:31:12 2016

HIV protease: downstream analysis
generate a list of all possible aa mutants in the combinatorial library

@author: lei
"""

# generate all combinations 
outfile =  open('/mnt/hgfs/analysis_sunlab/analysis_protease/sequence_list_aa','w') 
outfile2 =  open('/mnt/hgfs/analysis_sunlab/analysis_protease/sequence_list_number','w') 

import itertools
mutation_library = ['L10F','V32I','M46I','I47V','I50V',['I54M','I54L'],'T74P','L76V',['V82T','V82F'],'I84V','L90M']
#sequence is a 11-site string of 0/1/2: e.g. 11111211211 
#residue 10 on the left. 0:no mutation, 1/2: aa mutation
#iterate
sequence_list_tuple = list(itertools.product(range(2),range(2),range(2),range(2),range(2),range(3),range(2),range(2),range(3),range(2),range(2)))
for seq in sequence_list_tuple:
    muts_aa=[]
    muts=[]
    for i in range(0,len(seq)):
        muts.append(str(seq[i]))        
        if seq[i]!=0:
            if i==5 or i==8: #residue 54 or 82
                muts_aa.append(mutation_library[i][seq[i]-1]) 
            else:
                muts_aa.append(mutation_library[i])
                
    if len(muts_aa)==0:
        muts_seq_aa = 'WT'
    else:
        muts_seq_aa = '_'.join(muts_aa)                        
    muts_seq=''.join(muts)    
    
    #write to file: sequence_list_aa, sequence_list
    outfile.write(muts_seq_aa + '\n')
    outfile2.write(muts_seq + '\n')    
    
    
outfile.close()
outfile2.close()









