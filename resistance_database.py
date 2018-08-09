#!/usr/bin/env python
from scipy.stats.mstats import gmean
def main():
  infile = open('analysis/PI_DataSet.txt')
  header = infile.readline().rstrip().rsplit()
  tpvdict = {}
  drvdict = {}
  mutlist = ['10F','32I','46I','47V','50V','54L','54M','74P','76V','82T','82F','84V','90M']
  for line in infile:
    line = line.rstrip().rsplit()
    muts = []
    for i,mut in enumerate(line[9:]):
      if str(i+1)+mut in mutlist:
        muts.append(str(i+1)+mut)
    if len(muts) == 1:
      mut = muts[0]
      if mut not in tpvdict:
        tpvdict[mut] = []
	drvdict[mut] = []
      if line[7] != 'NA':
        tpvdict[mut].append(float(line[7]))
      if line[8] != 'NA':
        drvdict[mut].append(float(line[8]))
  print tpvdict
  print drvdict
  infile.close()
  for mut in tpvdict:
    tpvdict[mut] = gmean(tpvdict[mut])
    drvdict[mut] = gmean(drvdict[mut])
  outfile = open('analysis/resistance_single.txt','w')
  outfile.write('mutation\tTPV\tDRV\n')
  for mut in tpvdict:
    outfile.write(mut+'\t'+str(tpvdict[mut])+'\t'+str(drvdict[mut])+'\n')
  outfile.close()
        

if __name__ == '__main__':
  main()
