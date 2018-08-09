#!/usr/bin/env python
from fitness_trajectory_stringent import readfit
from fitness_trajectory_stringent import calc_trajectory
from math import factorial as fac
from itertools import combinations as comb

def main():
  fitdict,sindict = readfit()
  mutlist = sindict.keys()
  tradict = calc_trajectory(fitdict,mutlist)
  accdict = {}
  #accessible dict[mutation number][frequency of accessible path] = genotype count
  pathdict = {}
  #pathdict[mutation number][geno] = [path1,path2,...]  path=[subgeno1,subgeno2...]
  pathdict[1] = {}
  for i in range(2,5):
    pathdict[i] = {}
    for geno in tradict[i]:
      pathdict[i][geno] = []
      for subgeno in tradict[i][geno]:
        if subgeno not in pathdict[i-1]:
	  pathdict[i][geno].append([subgeno])
	else:
          pathdict[i][geno].extend([[subgeno]+subpath for subpath in pathdict[i-1][subgeno]])
  for i in range(2,5):
    accdict[i] = {}
    for geno in pathdict[i]:
      totalpath = fac(i)
      accpath = 0
      for path in pathdict[i][geno]:
        _acc = 1
	for j in range(len(path)-1):
	  subgeno = path[j]
	  if fitdict[subgeno] > fitdict[path[j+1]]:
	    _acc = 0
	if fitdict[geno] > fitdict[path[0]]:
	  _acc = 0
	if _acc == 1:
	  accpath += 1
      acc_freq = float(accpath)/totalpath
      if acc_freq not in accdict[i]:
        accdict[i][acc_freq] = 0
      accdict[i][acc_freq] += 1
  outfile = open('analysis/stringent_accessible_trajectory.txt','w')
  outfile.write('hamming_distance\taccessible_frequency\tgenotype_count\n')
  for i in range(2,5):
    for freq in accdict[i]:
      outfile.write(str(i)+'\t'+str(freq)+'\t'+str(accdict[i][freq])+'\n')
  outfile.close()

if __name__ == '__main__':
  main()
