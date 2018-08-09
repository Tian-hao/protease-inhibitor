#!/usr/bin/env python
from Bio import SeqIO
from collections import Counter
from scipy.stats import entropy

def main():
  seqdict = {}
  infile = open('analysis/HIV-1_pol_noD.fasta')
  for record in SeqIO.parse(infile,'fasta'):
    seq = str(record.seq)
    polseq = seq[138:248]
    #polseq = seq[201:406]
    seqdict[str(record.id)] = polseq
  refname = 'B.FR.1983.HXB2-LAI-IIIB-BRU.K03455'
  refseq = seqdict[refname]
  newseqdict = {}
  pos = 0
  for i in range(len(refseq)):
    if refseq[i] == '-': continue
    pos += 1
    newseqdict[pos] = {}
    for name in seqdict:
      newseqdict[pos][name] = seqdict[name][i]
  entdict = {}
  for pos in newseqdict:
    entdict[pos] = entropy(Counter(newseqdict[pos].values()).values())
  outfile = open('analysis/diversity_noD.txt','w')
  for pos in range(len(refseq.replace('-',''))):
    outfile.write(str(pos+1)+'\t'+refseq.replace('-','')[pos]+'\t'+str(entdict[pos+1])+'\n')
  infile.close()
  outfile.close()


if __name__ == '__main__':
  main()
