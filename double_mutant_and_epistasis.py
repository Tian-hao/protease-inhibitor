#!/usr/bin/env python
from numpy import std
from math import log

def main():
  infile = open('analysis/fitness_filter.txt')
  header = infile.readline().rstrip().rsplit('\t')
  liblist = header[1:]
  sindict = {}
  doudict = {}
  for line in infile:
    line = line.rstrip().rsplit('\t')
    mut = line[0]
    mut_num = len(mut.rsplit('_'))
    if mut == 'WT': continue
    if mut_num > 2: continue
    if mut_num == 1: 
      sindict[mut] = {}
      for i,fit in enumerate(line[1:]): 
        sindict[mut][liblist[i]] = float(fit)
    if mut_num == 2:
      doudict[mut] = {}
      for i,fit in enumerate(line[1:]): 
        doudict[mut][liblist[i]] = float(fit)
  infile.close()

  for mut in sindict:
    fitlist = [log(x,10) for x in sindict[mut].values()]
    sindict[mut] = [sum(fitlist)/3,std(fitlist)]
  outfile = open('analysis/fitness_predict_double.txt','w')
  outfile.write('mutations\tfitness\tstd\tpredicted\tpredstd\n')
  for muts in doudict:
    fitdict = doudict[muts]
    mean = []
    for fit in fitdict.values():
      if fit != 0: mean.append(log(fit,10))
    sd     = std(mean)
    mean   = sum(mean)/len(mean)
    doudict[muts] = mean
    outfile.write(muts+'\t'+str(mean)+'\t'+str(sd))
    muts   = muts.rsplit('_')
    predm  = sum([sindict[mut][0] for mut in muts])
    predsd = sum([sindict[mut][1] for mut in muts])
    outfile.write('\t'+str(predm)+'\t'+str(predsd)+'\n')
  outfile.close()

  outfile = open('analysis/epistasis_double.txt','w')
  outfile.write('mutations\tepistasis\n')
  for muts in doudict:
    fit = doudict[muts]
    muts = muts.rsplit('_')
    epi = fit
    for mut in muts:
      epi -= sindict[mut][0]
    outfile.write('_'.join(muts)+'\t'+str(epi)+'\n')
  outfile.close()

if __name__ == '__main__':
  main()
