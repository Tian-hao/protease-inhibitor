#!/usr/bin/env python
from math import log
import numpy as np

def main():
  valfile = open('validation/hiseq/relative_fitness_designed.txt')
  valdict = {}
  alldict = {}
  seqdict = {}
  valstd  = {}
  allstd  = {}
  seqstd  = {}
  header  = valfile.readline()
  for line in valfile:
    line = line.rstrip().rsplit('\t')
    mut  = line[1]
    lib  = line[0].rsplit('_')[0]
    fit  = log(float(line[2]),10)
    if mut not in alldict:
      alldict[mut] = []
      valdict[mut] = []
    if lib == 'mixA':
      alldict[mut].append(fit)
    else:
      valdict[mut].append(fit)
  for mut in alldict:
    if len(alldict[mut]) > 0:
      print alldict[mut]
      allstd[mut] = np.std(alldict[mut])
      alldict[mut] = sum(alldict[mut])/len(alldict[mut])
    else:
      alldict[mut] = 'NA'
    if len(valdict[mut]) > 0:
      valstd[mut] = np.std(valdict[mut])
      valdict[mut] = sum(valdict[mut])/len(valdict[mut])
    else:
      valdict[mut] = 'NA'
  valfile.close()
  seqfile = open('analysis/fitness_single.txt')
  header = seqfile.readline()
  for line in seqfile:
    line = line.rstrip().rsplit('\t')
    mut  = line[0]
    fit  = float(line[1])
    std  = float(line[2])
    seqdict[mut] = fit
    seqstd[mut]  = std
  seqfile.close()

  outfile = open('analysis/validation_fitness_correlation.txt','w')
  outfile.write('mut\tfitness\tfitness_std\tvalidation\tvalidation_std\tmixed_validation\tmixed_validation_std\n')
  for mut in seqdict:
    outfile.write(mut+'\t'+str(seqdict[mut])+'\t'+str(seqstd[mut])+'\t')
    if mut in valdict:
      outfile.write(str(valdict[mut])+'\t'+str(valstd[mut])+'\t')
    else: 
      outfile.write('NA\tNA\t')
    if mut in alldict:
      outfile.write(str(alldict[mut])+'\t'+str(allstd[mut])+'\n')
    else:
      outfile.write('NA\tNA\n')
  outfile.close()

if __name__ == '__main__':
  main()
