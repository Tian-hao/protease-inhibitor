#!/usr/bin/env python
#sequencing depth

def main():
  infile = open('code/count_aa_table')
  header = infile.readline().rstrip().rsplit('\t')
  liblist = header[2:]
  depdict = {}
  wtdict  = {}
  for lib in liblist:
    depdict[lib] = 0
    wtdict[lib] = 0
  for line in infile:
    line = line.rstrip().rsplit('\t')
    for i,count in enumerate(line[2:]):
      if line[0] == 'WT':
        wtdict[liblist[i]] += int(count)
      depdict[liblist[i]] += int(count)
  infile.close()
  outfile = open('analysis/depth.txt','w')
  outfile.write('library\twild-type\tdepth\n')
  for lib in liblist:
    outfile.write(lib+'\t'+str(wtdict[lib])+'\t'+str(depdict[lib])+'\n')
  outfile.close()
    

if __name__ == '__main__':
  main()
