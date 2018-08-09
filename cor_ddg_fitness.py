#!/usr/bin/env python
import numpy as np
from math import log

def main():
  infile = open('analysis/fitness_filter.txt')
  header = infile.readline()
  fitdict = {}
  stddict = {}
  for line in infile:
    line = line.rstrip().rsplit('\t')
    mut  = line[0]
    fitlist = []
    for fit in line[1:]:
      fit = float(fit)
      if fit < 0.0001: fit = 0.0001
      fitlist.append(log(fit,10))
    fit = sum(fitlist)/len(fitlist)
    stddict[mut] = np.std(fitlist)
    fitdict[mut] = fit
  infile.close()

  infile = open('analysis/old_fitness_TPV_ddG.csv')
  header = infile.readline()
  ddgdict = {}
  for line in infile:
    line = line.rstrip().rsplit(',')
    mut  = line[0]
    ddg  = float(line[-1])
    ddgdict[mut] = ddg
  infile.close()

  outfile = open('analysis/correlation_ddg_fitness.txt','w')
  outfile.write('mut\tfitness\tfitness_std\tddg\tmut_number\n')
  for mut in fitdict:
    if mut not in ddgdict: continue
    outfile.write(mut+'\t'+str(fitdict[mut])+'\t'+str(stddict[mut])+'\t'+str(ddgdict[mut])+'\t'+str(len(mut.rsplit('_')))+'\n')
  outfile.close()

if __name__ == '__main__':
  main()
