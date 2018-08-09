#!/usr/bin/env python
from random import sample

def main():
  infile = open('analysis/potts_fitness.txt')
  header = infile.readline()
  sindict = {}
  doudict = {}
  for line in infile:
    line = line.rstrip().rsplit('\t')
    mut  = line[0]
    pots = float(line[1])
    mutlist = mut.rsplit('_')
    if len(mutlist) > 2: continue
    if len(mutlist) == 2: doudict[mut] = pots
    if len(mutlist) == 1: sindict[mut] = pots
  infile.close()

  predict = {}
  for mut in doudict:
    mutlist = mut.rsplit('_')
    predict[mut] = sindict[mutlist[0]]+sindict[mutlist[1]]

  outfile = open('analysis/potts_double_prediction.txt','w')
  outfile.write('muts\tobserved_energy\tpredicted_energy\n')
  for mut in doudict:
    outfile.write(mut+'\t'+str(doudict[mut])+'\t'+str(predict[mut])+'\n')
  outfile.close()

  infile = open('analysis/Potts_energy_for_all_single_and_double.txt')
  header = infile.readline()
  sindict = {}
  doudict = {}
  for line in infile:
    line = line.rstrip().rsplit('\t')
    mut  = line[0]
    pots = float(line[1])
    mutlist = mut.rsplit('_')
    if len(mutlist) == 1: sindict[mut] = pots
    if len(mutlist) == 2: doudict[mut] = pots
  infile.close()

  outfile = open('analysis/potts_double_prediction_all.txt','w')
  outfile.write('muts\tobserved_energy\tpredicted_energy\n')
  for muts in sample(doudict,1000):
    mutlist = muts.rsplit('_')
    predict = sindict[mutlist[0]]+sindict[mutlist[1]]
    outfile.write(muts+'\t'+str(doudict[muts])+'\t'+str(predict)+'\n')
  outfile.close()

if __name__ == '__main__':
  main()
