"""
last edit: 06/22/17

summary: this script processes MiSeq data, HIV protease (Oct 2016 run)
this script performs the following analysis
1) filter mutants: only mutations in the library are allowed 
further filtering in downstream analysis is required
2) count mutants for each condition/barcode

input: 
mutation_called (output file from Mapper1.py)

output: 
count
count_summary

run time: ~300s (~6M mapped reads)

@author: lei
"""


#!/usr/bin/env python
import os
import sys
import glob
import time

#specify input/output files
#Oct 2016 run
#infile = '/home/lei/Data/LD10042016-40246275/mutation_called' 
#outfile = open('/mnt/hgfs/analysis_sunlab/analysis_protease/count','w') 
#outfile2 = open('/mnt/hgfs/analysis_sunlab/analysis_protease/count_summary','w')#summary
#Dec 2016 run, sample 1(TPV)
#infile = '/home/lei/Data/RS12122016/RS-Sample1-41647632/mutation_called' 
#outfile = open('/mnt/hgfs/analysis_sunlab/analysis_protease/dec2016_TPV/count','w') 
#outfile2 = open('/mnt/hgfs/analysis_sunlab/analysis_protease/dec2016_TPV/count_summary','w')#summary
#Dec 2016 run, sample 2(DRV)
#infile = '/home/lei/Data/RS-Sample2-41681208/mutation_called' 
#outfile = open('/mnt/hgfs/analysis_sunlab/analysis_protease/dec2016_DRV/count','w') 
#outfile2 = open('/mnt/hgfs/analysis_sunlab/analysis_protease/dec2016_DRV/count_summary','w')#summary
#June 2017 run
infile = '/home/lei/Data/RS06072017-40849810/LD060117-49370361/mutation_called'
outfile = open('/mnt/hgfs/analysis_sunlab/analysis_HIV_protease/june2017/count','w') 
outfile2 = open('/mnt/hgfs/analysis_sunlab/analysis_HIV_protease/june2017/count_summary','w')#summary

start_time = time.time()
#%%
#READ BARCODE FILE (demultiplex)
#output: barcodes
reffolder = '/mnt/hgfs/analysis_sunlab/analysis_HIV_protease/reference/'
#Oct 2016 run
#reffilename = reffolder +'barcode_miseq'
#Dec 2016 run: sample 1 (TPV)
#reffilename = reffolder +'barcode_miseq_dec2016_TPV'
#Dec 2016 run: sample 2 (DRV)
#reffilename = reffolder +'barcode_miseq_dec2016_DRV'
#June 2017 run
reffilename = reffolder +'barcode_miseq_june2017'    

reffile   = open(reffilename,'r') # Barcode in sequencing library preparation: demultiplex 
barcodes = {}
for line in reffile.xreadlines():
    if '>' in line:
        ID = line.rstrip().replace('>','')    #3bp barcode
    else:
        barcodes[ID] = line.rstrip()     #specified library  
reffile.close()      
#barcode_list = barcodes.keys()

#%%
#initialize: for each condition/barcode
count_all={}
count_wt={}
count_mutant={}
count_others={}
mutdict={}
for barcode in barcodes.keys():
    count_all[barcode] = 0
    count_wt[barcode] =0 #WT
    count_mutant[barcode] = 0 #mutants in the library
    count_others[barcode] = 0 #mutants not in the library
    #initialize mutant list
    mutdict[barcode] = {}
    mutdict[barcode]['WT'] = 0

#%%    
#  mutpos  = [0,66,110,111,120,132,134,192,198,216,217,222,240]  
#L10F: C0T
#V32I: G66A
#M46I: G110C
#I47V: A111G
#I50V: A120G
#I54M/L: A132C, C134G
#T74P: A192C
#L76V: T198G
#V82T/F: G216A-T217C, G216T
#I84V: A222G
#L90M: T240A  
mutation_list = ['C0T','G66A','G110C','A111G','A120G','A132C','C134G','A192C','T198G','G216A','T217C','G216T','A222G','T240A']

#read file: mutation       
inhandle = open(infile,'r')

count = 0
#FOR LOOP STARTS HERE
for line in inhandle:
    count +=1 
#    if count>100000:
    if count%10000 == 0:      
        print count
        print "---time elapsed: %s seconds ---" % (time.time() - start_time)    
#        break #debug
        
    [bc, mut_called] = line.strip().split('\t')
#    if bc not in barcodes.keys(): #this is not necessary. debug Mapper1.py
#        continue 
    count_all[bc] += 1
    
#    print bc
#    print mut_called
#  bc_index = barcode_list.index(bc) 
  #check WT
    if mut_called == 'WT': 
        mutdict[bc]['WT'] += 1
        count_wt[bc] += 1
    else:  #if not WT, check if there are mutations not supposed to be in the library 
        muts = mut_called.rstrip().rsplit('_')[1:] #skip element 0 because the string starts with '_'
        flag = 0
        for mut in muts:
            if mut not in mutation_list: 
                flag = 1
        if flag == 0: #if the mutant is in the library, add to count list
            count_mutant[bc] += 1   
            #muts.sort(key = lambda s: int(s[1:-1])) #sort mutation by position
            muts = '_'.join(muts)
            #print muts
            if muts not in mutdict[bc].keys(): #if the mutant is not already in the list
                mutdict[bc][muts] = 0
            mutdict[bc][muts] += 1       
        else:
            count_others[bc] += 1           

#FOR LOOP ENDS HERE  
inhandle.close()

#%%
for bc in mutdict.keys():
    for mut in mutdict[bc].keys():
        outfile.write(barcodes[bc]+'\t'+ mut+'\t'+str(mutdict[bc][mut])+'\n')
    #write summary
    outfile2.write(barcodes[bc] +'\t'+'all'+'\t' +str(count_all[bc])+'\n')
    outfile2.write(barcodes[bc] +'\t'+'WT' +'\t' +str(count_wt[bc])+'\n')
    outfile2.write(barcodes[bc] +'\t'+'mutlib' +'\t' +str(count_mutant[bc])+'\n')
    outfile2.write(barcodes[bc] +'\t'+'others' +'\t' +str(count_others[bc])+'\n')

outfile.close()
outfile2.close()



