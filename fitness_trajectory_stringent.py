#!/usr/bin/env python
#if one replicate is 0, then the mutant fitness is 0.0001
from itertools import combinations as comb
from math import log

def readfit():
  fitdict = {}
  sindict = {}
  infile = open('analysis/fitness_filter.txt')
  header = infile.readline()
  for line in infile:
    line = line.rstrip().rsplit('\t')
    muts = line[0]
    if muts == 'WT': continue
    fitlist = [float(x) for x in line[1:]]
    mean = []
    for fit in fitlist:
      if fit != 0:
        mean.append(log(fit,10))
    if len(mean) < 3:
      fitdict[muts] = -4
    else:
      fitdict[muts] = sum(mean)/len(mean)
    muts = muts.rsplit('_')
    if len(muts) == 1:
      sindict[muts[0]] = sum([log(float(x),10) for x in line[1:]])/3
  infile.close()
  return fitdict,sindict

def main():
  fitdict,sindict = readfit()
  mutlist = sindict.keys()
  
  outfile = open('analysis/stringent_trajectory_fitness.txt','w')
  tradict = calc_trajectory(fitdict,mutlist)
  for i in range(2,5):
    for geno in tradict[i]:
      for child_geno in tradict[i][geno]:
        outfile.write(str(i-1)+'\t'+str(fitdict[child_geno])+'\t'+str(i)+'\t'+str(fitdict[geno])+'\n')
  for mut in sindict:
    outfile.write('0\t0\t1\t'+str(sindict[mut])+'\n')
  outfile.close()

def calc_trajectory(fitdict,mutlist):
  tradict = {}
  #trajectory dict[mutation number][genotype] = [child genotypes]
  tradict[1] = {}
  for mut in mutlist:
    tradict[1][mut] = ['WT']
  for i in range(2,5):
    tradict[i] = {}
    for geno in fitdict:
      muts = geno.rsplit('_')
      if len(muts) != i: continue
      for mut in muts:
        child_geno = '_'.join([x for x in muts if x != mut])
        if child_geno in tradict[i-1]:
	  if geno not in tradict[i]: tradict[i][geno] = []
	  tradict[i][geno].append(child_geno)
  return tradict



if __name__ == '__main__':
  main()
