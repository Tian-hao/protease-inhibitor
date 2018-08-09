#!/usr/bin/env python
import sys
from math import log
from math import factorial as fac
from itertools import combinations as comb
#it is actually count_viable
def count_lethal(fitdict,sindict,thre):
  olethal = {}
  plethal = {}
  depth   = make_depth(sindict)
  for i in range(1,10):
    olethal[i] = 0
    plethal[i] = 0
  for muts in fitdict:
    mutlist = muts.rsplit('_')
    fit = fitdict[muts]
    if fit > thre: 
      olethal[len(mutlist)] += 1
    pfit = 0
    for mut in mutlist:
      pfit += sindict[mut]
    if pfit > thre:
      plethal[len(mutlist)] += 1
  for i in range(1,10):
    olethal[i] /= float(depth[i])
    plethal[i] /= float(depth[i])
  return olethal, plethal

def make_depth(sindict):
  depdict = {}
  for i in range(1,10):
    depth = 0
    for muts in comb(sindict.keys(),i):
      if len(set([mut[1:-1] for mut in muts])) == i:
        depth += 1
    depdict[i] = depth
  return depdict

def readfit():
  infile = open('analysis/fitness_filter.txt')
  header = infile.readline().rstrip().rsplit('\t')
  fitdict = {}
  sindict = {}
  for line in infile:
    line = line.rstrip().rsplit('\t')
    muts = line[0]
    if muts == 'WT': continue
    fitlist = [float(x) for x in line[1:]]
    mean = []
    for fit in fitlist:
      if fit != 0:
        mean.append(log(fit,10))
    if len(mean) == 0:
      continue
      #fitdict[muts] = -4
    else:
      fitdict[muts] = sum(mean)/len(mean)
    muts = muts.rsplit('_')
    if len(muts) == 1:
      sindict[muts[0]] = sum([log(float(x),10) for x in line[1:]])/3
  infile.close()
  return fitdict,sindict

def main():
  fitdict,sindict = readfit()
  thre = -2
  olethal2,plethal2 = count_lethal(fitdict,sindict,thre)
  thre = -4
  olethal4,plethal4 = count_lethal(fitdict,sindict,thre)
  outfile = open('analysis/lethal_fraction.txt','w')
  outfile.write('Num\tviable2\tpredviable2\tviable4\tpredviable4\n')
  for i in range(1,10):
    outfile.write(str(i)+'\t'+str(olethal2[i])+'\t'+str(plethal2[i])+'\t'+str(olethal4[i])+'\t'+str(plethal4[i])+'\n')
  outfile.close()

if __name__ == '__main__':
  main()
