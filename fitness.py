#!/usr/bin/env python
from scipy.stats.mstats import gmean
from math import log
from itertools import product
from numpy import std

def main():
  infile = open('code/count_aa_table')
  countdict = {}
  header = infile.readline().rstrip().rsplit('\t')
  liblist = header[2:]
  for line in infile:
    line = line.rstrip().rsplit('\t')
    mut  = line[0]
    countdict[mut] = {}
    for i, count in enumerate(line[2:]):
      countdict[mut][liblist[i]] = float(count)
  infile.close()

  freqdict = {}
  depthdict = {}
  for lib in liblist:
    depthdict[lib] = 0
  for mut in countdict:
    for lib in liblist:
      depthdict[lib] += countdict[mut][lib]
  for mut in countdict:
    freqdict[mut] = {}
    for lib in liblist:
      freqdict[mut][lib] = countdict[mut][lib] / depthdict[lib]

  outlier = []
  trlibs  = ['tr_R1','tr_R2','tr_R3']
  for mut in freqdict:
    for lib in trlibs:
      if freqdict[mut][lib] < 5e-5:
        outlier.append(mut)
        break
    if freqdict[mut][lib] < 5e-5: continue
    for lib1,lib2 in product(trlibs,trlibs):
      if freqdict[mut][lib1] / freqdict[mut][lib2] > 10:
	outlier.append(mut)
	break
  
  fitdict = {}
  fitlibs = ['nd_R1','nd_R2','nd_R3']
  for mut in freqdict:
    if mut in outlier: 
      continue
    fitdict[mut] = {}
    for lib in fitlibs:
      reflib = 'tr_'+lib.rsplit('_')[1]
      fitdict[mut][lib] = freqdict[mut][lib]/freqdict[mut][reflib]

  outfile = open('analysis/fitness_3mut.txt','w')
  outfile.write('mutations\t'+'\t'.join(fitlibs)+'\tcolor\n')
  for mut in fitdict:
    if mut == 'WT': continue
    mut = mut.rsplit('_')
    if len(mut) > 3: continue
    outfile.write('_'.join(mut))
    for lib in fitlibs:
      outfile.write('\t'+str(fitdict['_'.join(mut)][lib]/fitdict['WT'][lib]))
    outfile.write('\t'+str(len(mut))+'\n')
  outfile.close()

  fitlist = []
  outfile = open('analysis/fitness_single.txt','w')
  outfile.write('mutation\tmean\tstd\t'+'\t'.join(fitlibs)+'\n')
  for mut in fitdict:
    if mut == 'WT' or '_' in mut: continue
    fitlist.append(mut)
  fitlist = sorted(fitlist,key=lambda mut:int(mut[1:-1]))
  for mut in fitlist:
    mean = []
    for lib in fitlibs:
      mean.append(log(fitdict[mut][lib]/fitdict['WT'][lib],10))
    sd = std(mean)
    mm = sum(mean)/len(mean)
    outfile.write(mut+'\t'+str(mm)+'\t'+str(sd)+'\t'+'\t'.join([str(x) for x in mean])+'\n')
  outfile.close()
  
  outfile = open('analysis/fitness_filter.txt','w')
  outfile.write('mutations\t'+'\t'.join(fitlibs)+'\n')
  for mut in fitdict:
    outfile.write(mut)
    for lib in fitlibs:
      outfile.write('\t'+str(fitdict[mut][lib]/fitdict['WT'][lib]))
    outfile.write('\n')
  outfile.close()

  outfile = open('analysis/fitness_one2four.txt','w')
  outfile.write('mutation\t'+'\t'.join(fitlibs)+'\tmut_number\n')
  for mut in fitdict:
    if mut == 'WT': continue
    mut = mut.rsplit('_')
    outfile.write('_'.join(mut))
    for lib in fitlibs:
      outfile.write('\t'+str(fitdict['_'.join(mut)][lib]/fitdict['WT'][lib]))
    outfile.write('\t'+str(len(mut))+'\n')
  outfile.close() 

  outfile = open('analysis/fitness_mut_number.txt','w')
  outfile.write('mutation\tAverageFitness\tMutNumber\n')
  for mut in fitdict:
    if mut == 'WT': continue
    outfile.write(mut)
    fitlist = []
    for lib in fitlibs:
      if fitdict[mut][lib] >= 0.0001:
        fitlist.append(fitdict[mut][lib]/fitdict['WT'][lib])
      else: 
        fitlist.append(0.0001)
    fit = log(gmean(fitlist),10)
    mut = mut.rsplit('_')
    outfile.write('\t'+str(fit)+'\t'+str(len(mut))+'\n')
  outfile.close()


if __name__ == '__main__':
  main()
