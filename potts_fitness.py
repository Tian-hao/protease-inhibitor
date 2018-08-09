#!/usr/bin/env python
import numpy as np
from math import log

def readfit():
  infile = open('analysis/fitness_filter.txt')
  header = infile.readline().rstrip().rsplit('\t')
  fitdict = {}
  stddict = {}
  for line in infile:
    line = line.rstrip().rsplit('\t')
    mut  = line[0]
    fitlist = []
    for fit in line[1:]:
      if float(fit) > 0:
        fitlist.append(log(float(fit),10))
    if len(fitlist) == 0:
      fitdict[mut] = -4
      stddict[mut] = 0
    else:
      fitdict[mut] = sum(fitlist)/len(fitlist)
      stddict[mut] = np.std(fitlist)
  infile.close()
  return fitdict,stddict

def main():
  fitdict,stddict = readfit()
  infile = open('code/mutation_protease_energies.csv','rU')
  header = infile.readline().rstrip().rsplit(',')
  line1 = infile.readline().rstrip().rsplit(',')
  outfile = open('analysis/potts_fitness.txt','w')
  outfile.write('mut\tepotts\tfitness\tstd\n')
  epotts_wt = float(line1[-1])
  for line in infile:
    records = line.rstrip().rsplit(',')
    mut = records[1]
    if len(mut.rsplit('_')) <= 12 and mut in fitdict:
      epotts = float(records[-1])-epotts_wt
      outfile.write(mut+'\t'+str(epotts)+'\t'+str(fitdict[mut])+'\t'+str(stddict[mut])+'\n')
  infile.close()
  outfile.close()

if __name__ == '__main__':
  main()
