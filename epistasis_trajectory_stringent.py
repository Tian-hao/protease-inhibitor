#!/usr/bin/env python
from fitness_trajectory_stringent import readfit
from fitness_trajectory_stringent import calc_trajectory
from numpy import std

def main():
  fitdict,sindict = readfit()
  mutlist = sindict.keys()
  tradict = calc_trajectory(fitdict,mutlist)
  mutlist = ['L10F','V32I','M46I','I47V','I50V','I54M','I54L','T74P','L76V','V82T','V82F','I84V','L90M']
  newfitdict = {}
  outfile = open('analysis/stringent_trajectory_epistasis.txt','w')
  stdfile = open('analysis/stringent_trajectory_epistasis_std.txt','w')
  for col,mut in enumerate(mutlist):
    epidict = {}
    epidict[1] = 0
    for i in range(2,5):
      for geno in tradict[i]:
        if mut in geno:
	  muts = geno.rsplit('_')
	  if len(muts) > 1:
	    muts.remove(mut)
	    newgeno = '_'.join(muts)
	  else:
	    continue
	  if newgeno in tradict[i-1]:
	    newfitdict[newgeno] = fitdict[geno]-fitdict[newgeno]-fitdict[mut]
    for i in range(2,5):
      mean = []
      for newgeno in newfitdict:
        if len(newgeno.rsplit('_')) == i-1:
          mean.append(newfitdict[newgeno])
      epidict[i] = sum(mean)/len(mean)
      stdfile.write(mut+'\t'+str(i)+'\t'+str(epidict[i])+'\t'+str(std(mean))+'\t'+str(col+1)+'\n')
      outfile.write(mut+'\t'+str(i-1)+'\t'+str(epidict[i-1])+'\t'+str(i)+'\t'+str(epidict[i])+'\t'+str(col+1)+'\n')
  outfile.close()
  stdfile.close()

if __name__ == '__main__':
  main()
