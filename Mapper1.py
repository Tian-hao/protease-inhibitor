"""
last edit: 06/22/17

summary: this script processes MiSeq data, HIV protease 
1) Oct 2016 run
2) Dec 2016 run 1&2
3) June 2017 run

this script performs the following analysis
1) quality control: ID and barcodes of paired-end reads match
2) mapping: identify forward/reverse read, check read length
3) call mutations: compare to reference sequence, a mutation is called if it is found on both PE reads and quality score >=30

input: 
MiSeq data, FASTQ files (R1 and R2)
reference sequence (in the script)
barcode file

output:
1) summary: # of reads
2) mutation: for successfully mapped reads, record barcode and mutations 

run time: ~2200s (~20M reads)

@author: lei
"""

#!/usr/bin/env python
import os
import sys
import string
import operator
import time
from itertools import imap
from Bio import SeqIO

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(imap(operator.ne, str1, str2))

def rc(seq):
    complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    rcseq = seq.translate(complements)[::-1]
    return rcseq

def mapread(seq,primerF,primerR):
    seq_priF = seq[0:12]
    if hamming(primerF,seq_priF) < 2: return ['F'] #at most 1 mutation from the reference. 
    seq_priR = seq[0:12]
    if hamming(primerR,seq_priR) < 2: return ['R']
    return 'no'

start_time = time.time()
def main():
    #input files: FASTQ
    #  infile1 = open('/u/scratch/y/yushendu/OY05022016/split/'+sys.argv[1])
    #  infile2 = open('/u/scratch/y/yushendu/OY05022016/split/'+sys.argv[1].replace('_R1_','_R2_'))
    #output files
    #  outfile = open('/u/scratch/y/yushendu/OY05022016/pro/'+sys.argv[1]+'.m','w')
    #  numfile = open('/u/scratch/y/yushendu/OY05022016/pro/'+sys.argv[1]+'.count','w')

    #Oct 2016 run
#    infile1 = open('/home/lei/Data/LD10042016-40246275/LD10042016_S1_L001_R1_001.fastq')
#    infile2 = open('/home/lei/Data/LD10042016-40246275/LD10042016_S1_L001_R2_001.fastq')
#    outfile = open('/home/lei/Data/LD10042016-40246275/mutation_called','w')
#    outfile2 = open('/home/lei/Data/LD10042016-40246275/mutation_summary','w')
    #Dec 2016 run: sample 1 (TPV)
#    infile1 = open('/home/lei/Data/RS12122016/RS-Sample1-41647632/RS-Sample1_S1_L001_R1_001.fastq')
#    infile2 = open('/home/lei/Data/RS12122016/RS-Sample1-41647632/RS-Sample1_S1_L001_R2_001.fastq')
#    outfile = open('/home/lei/Data/RS12122016/RS-Sample1-41647632/mutation_called','w')
#    outfile2 = open('/home/lei/Data/RS12122016/RS-Sample1-41647632/mutation_summary','w')
    #Dec 2016 run: sample 2 (DRV)
#    infile1 = open('/home/lei/Data/RS-Sample2-41681208/RS-Sample2_S1_L001_R1_001.fastq')
#    infile2 = open('/home/lei/Data/RS-Sample2-41681208/RS-Sample2_S1_L001_R2_001.fastq')
#    outfile = open('/home/lei/Data/RS-Sample2-41681208/mutation_called','w')
#    outfile2 = open('/home/lei/Data/RS-Sample2-41681208/mutation_summary','w')
    #June 2017 run
    infile1 = open('/home/lei/Data/RS06072017-40849810/LD060117-49370361/LD060117_S1_L001_R1_001.fastq')
    infile2 = open('/home/lei/Data/RS06072017-40849810/LD060117-49370361/LD060117_S1_L001_R2_001.fastq')
    outfile = open('/home/lei/Data/RS06072017-40849810/LD060117-49370361/mutation_called','w')
    outfile2 = open('/home/lei/Data/RS06072017-40849810/LD060117-49370361/mutation_summary','w')    
    
    
    #primer sequences
    #MiSeq, after BpmI sites
    #  primerF = 'TCTTTGGCAGCGACCC'
    #  primerR = 'TGCAGCCAATCTGAGT'
      #reference sequence: HIV protease aa #10-90, 243 bp
    refseq = 'CTCGTCACAATAAAGATAGGGGGGCAATTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGCGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTG'
      #the primer is cut for MiSeq samples, take the first 12 bp for mapping forward/reverse reads
    primerF = 'CTCGTCACAATA'
    primerR = rc('AGAAATCTGTTG')

    #READ BARCODE FILE (demultiplex)
    #output: barcodes  
    reffolder = '/mnt/hgfs/analysis_sunlab/analysis_HIV_protease/reference/'
    #Oct 2016 run
#    reffilename = reffolder +'barcode_miseq' #October 2016 run
    #Dec 2016 run: sample 1 (TPV)
#    reffilename = reffolder +'barcode_miseq_dec2016_TPV'
    #Dec 2016 run: sample 2 (DRV)
#    reffilename = reffolder +'barcode_miseq_dec2016_DRV'
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

    #track read counts
    count_run = 0 # all reads
    count_hiv = 0 # reads with correct barcodes for HIV protease
    count_others = 0 # reads with correct barcodes for other samples   
    count_mapped = 0 #reads correctly mapped to HIV protease, with the right length     

    #initilize: for each condition/barcode
    count_wt = 0 # wild type reads
        
#    mutdict  = {}  #initilize the list of observed mutants
    amp_len=243 #length of amplicon
    
    #read FASTQ files
    handle1  = SeqIO.parse(infile1,'fastq')
    handle2  = SeqIO.parse(infile2,'fastq')
    
    ##FOR LOOP STARTS HERE
    for record1 in handle1:
        #debug
        if count_run%10000 == 0:
#            break
            print count_run
            print count_hiv
#            print count_others
            print count_mapped
            print "---time elapsed: %s seconds ---" % (time.time() - start_time)
        count_run += 1
        record2 = handle2.next()

        seq1  = str(record1.seq)
        seq2  = str(record2.seq)
        qua1  = record1.letter_annotations["phred_quality"]
        qua2  = record2.letter_annotations["phred_quality"] 
    
        #demultiplex: barcode
        seq1_bc = seq1[0:3]    
        seq2_bc = seq2[0:3]
#        print seq1_bc #debug
#        print seq2_bc
        if record1.id != record2.id or seq1_bc != seq2_bc: 
#        if record1.id != record2.id or hamming(seq1_bc,seq2_bc)>1: 
        #note: set to perfect match of barcodes to see how many reads are left
        #note: many reads in R2 starts with N
            continue    
        else:
            if seq1_bc in barcodes.keys():
                count_hiv += 1                 
#                print count_hiv #debug
            else:
#                if seq1_bc == 'TAC' or seq1_bc =='CCG': #    Oct 2016 run, barcode_flu:'TAC', 'CCG'
                count_others += 1
                continue
                    
        #map reads: return forwared/reverse. no offset
        #amplicon: remove barcode + 'T'
        seq1_amp = seq1[4::]
        seq2_amp = seq2[4::]
        qua1_amp = qua1[4::]
        qua2_amp = qua2[4::]
    
        info1 = mapread(seq1_amp,primerF,primerR)
        info2 = mapread(seq2_amp,primerF,primerR)
        if info1 == 'no' or info2 == 'no' or info1 == info2 or len(seq1_amp)<amp_len or len(seq2_amp)<amp_len: 
            continue               
        count_mapped += 1
        
        if info1[0] == 'F':
          readF = seq1_amp[0:amp_len]
          qualF = qua1_amp[0:amp_len]
          #reverse
          readR = rc(seq2_amp[0:amp_len])
          qualR = list(reversed(qua2_amp))[0:amp_len]
        if info1[0] == 'R':
          readR = rc(seq1_amp[0:amp_len])
          qualR = list(reversed(qua2_amp))[0:amp_len]
          readF = seq2_amp[0:amp_len]
          qualF = qua2_amp[0:amp_len]
      
        #debug
#        print readF
#        print readR
#        print qualF
#        print qualR
      
        muts = '' #the string that lists all mutations in a read, such as: C0T_G20A
        for i in range(0,len(readF)):
            if readF[i] != refseq[i] and readF[i]==readR[i] and int(qualF[i]) >= 30 and int(qualR[i])>=30:
                mut = refseq[i]+str(i)+readF[i]
                muts += ('_'+mut)
            
        #write to file (mutation): 1) barcode; 2) mutations or WT
        #parse barcode/condition in downstream analysis          
        if muts != '': 
            outfile.write(seq1_bc + '\t' + muts +'\n')
        else: 
            outfile.write(seq1_bc + '\t' + 'WT' +'\n')
            count_wt += 1
    ##FOR LOOP ENDS HERE
    
    #write to file (summary)     
    outfile2.write('all reads:\t'+str(count_run)+'\n')  
    outfile2.write('HIV protease reads:\t'+str(count_hiv)+'\n')
    outfile2.write('other reads:\t'+str(count_others)+'\n')  
    outfile2.write('correctly mapped reads:\t'+str(count_mapped)+'\n')
    #WT for all conditions/barcodes 
    outfile2.write('wild type reads:\t'+str(count_wt)+'\n')
    
    outfile.close()
    outfile2.close()
    infile1.close()
    infile2.close()

if __name__ == '__main__':
  main()
