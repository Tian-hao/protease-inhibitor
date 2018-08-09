# -*- coding: utf-8 -*-
"""
last edit: 06/22/17
Created on Fri Oct 21 15:20:22 2016

HIV protease: downstream analysis
transform nt mutation to aa mutation

input: count (output file from Mapper2.py)
output: count_aa

@author: lei
"""

# convert to aa mutation
map_aa = {'C0T':'L10F','G66A':'V32I','G110C':'M46I','A111G':'I47V','A120G':'I50V',
               'A132C':'I54M','C134G':'I54L',  #residue 54
               'A192C':'T74P','T198G':'L76V', 
               'G216AT217C':'V82T','G216T':'V82F', #residue 82
               'A222G':'I84V','T240A':'L90M'}
#read file: count
#Oct 2016 run
#infile =  open('/mnt/hgfs/analysis_sunlab/analysis_protease/count','r')            
#outfile = open('/mnt/hgfs/analysis_sunlab/analysis_protease/count_aa','w')
#Dec 2016 run: TPV 
#infile =  open('/mnt/hgfs/analysis_sunlab/analysis_protease/dec2016_TPV/count','r')            
#outfile = open('/mnt/hgfs/analysis_sunlab/analysis_protease/dec2016_TPV/count_aa','w')
#Dec 2016 run: DRV
#infile =  open('/mnt/hgfs/analysis_sunlab/analysis_protease/dec2016_DRV/count','r')            
#outfile = open('/mnt/hgfs/analysis_sunlab/analysis_protease/dec2016_DRV/count_aa','w')
#June 2017 run
infile =  open('/mnt/hgfs/analysis_sunlab/analysis_HIV_protease/june2017/count','r')            
outfile = open('/mnt/hgfs/analysis_sunlab/analysis_HIV_protease/june2017/count_aa','w')

count_remove = 0
count_all = 0
for line in infile:
    count_all += 1
    [sample, mut_called, count] = line.strip().split('\t')
    if mut_called == 'WT':
        muts_seq_aa = 'WT'
    else:
        muts = mut_called.rstrip().rsplit('_')             
        #handle special case: A132C and C134 should not appear together
        if ('A132C' in muts) and ('C134G' in muts):
            count_remove +=1
            continue
        #handle special case: G216A and T217C should appear together
        if ('G216A' in muts) and ('T217C' in muts):
            index_tmp = muts.index('G216A')
            muts[index_tmp] = 'G216AT217C' #combine into one entry
            muts.remove('T217C')
        elif ('G216A' in muts) != ('T217C' in muts):
#            if sample == 'plasmid':
#                print mut_called, count #debug
            count_remove += 1
            continue       
        
        muts_aa = []
        for mut in muts:
            muts_aa.append(map_aa[mut])
        muts_seq_aa = '_'.join(muts_aa)    
    outfile.write(sample+'\t'+ muts_seq_aa +'\t'+ count +'\n')    #write to file: count_aa
    
print count_all
print count_remove

infile.close()
outfile.close()