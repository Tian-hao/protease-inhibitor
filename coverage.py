#!/usr/bin/env python
import operator as op
from itertools import product

def choose(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, xrange(n, n-r, -1), 1)
    denom = reduce(op.mul, xrange(1, r+1), 1)
    return numer//denom

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

  outfile = open('analysis/coverage.txt','w')
  outfile.write('library\t'+'\t'.join([str(x) for x in range(1,12)])+'\tTotal\n')
  outfile.write('expected\t')
  total = 0
  for i in range(1,12):
    if i == 1: 
      n = 9 + 2*2; 
    if i > 1 and i <= 9:
      n = choose(9,i) + choose(9,i-1) * 4 + choose(9,i-2) * 4
    if i == 10:
      n = 4 + choose(9,8) * 4
    if i == 11: 
      n = 4
    total += n
    outfile.write(str(n)+'\t')
  outfile.write(str(total)+'\n')
  
  trlibs  = ['tr_R1','tr_R2','tr_R3']
  coverdict = {}
  for lib in trlibs:
    coverdict[lib] = {}
    for i in range(1,12):
      coverdict[lib][i] = 0
    coverdict[lib]['total'] = 0
  coverdict['overlap'] = {}
  for i in range(1,12):
    coverdict['overlap'][i] = 0
  coverdict['overlap']['total'] = 0
  for mut in freqdict:
    if mut == 'WT': continue
    mutnum = len(mut.rsplit('_'))
    for lib in trlibs:
      if freqdict[mut][lib] >= 5e-5:
        coverdict[lib][mutnum] += 1
  for lib in trlibs:
    for i in range(1,12):
      coverdict[lib]['total'] += coverdict[lib][i]

  outlier = []
  for mut in countdict:
    for lib in trlibs:
      if freqdict[mut][lib] < 5e-5:
        outlier.append(mut)
        break
    if freqdict[mut][lib] < 5e-5: continue
    for lib1, lib2 in product(trlibs,trlibs):
      if freqdict[mut][lib1] / freqdict[mut][lib2] > 10:
        outlier.append(mut)
        break
  for mut in freqdict:
    if mut == 'WT': continue
    mutnum = len(mut.rsplit('_'))
    if mut in outlier: continue
    coverdict['overlap'][mutnum] += 1
  for i in range(1,12):
    coverdict['overlap']['total'] += coverdict['overlap'][i]

  for lib in trlibs:
    outfile.write(lib)
    for i in range(1,12):
      outfile.write('\t'+str(coverdict[lib][i]))
    outfile.write('\t'+str(coverdict[lib]['total'])+'\n')
  outfile.write('overlap')
  for i in range(1,12):
    outfile.write('\t'+str(coverdict['overlap'][i]))
  outfile.write('\t'+str(coverdict['overlap']['total'])+'\n')
  outfile.close()

if __name__ == '__main__':
  main()
