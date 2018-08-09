#!/usr/bin/env python

def main():
  poslist = ['10','32','46','47','50','54','54','74','76','82','82','84','90']
  infile = open('analysis/dist.txt')
  outfile = open('analysis/dist_category.txt','w')
  outfile.write('pos1\tpos2\tdist\tcategory\n')
  for line in infile:
    line = line.rstrip().rsplit('\t')
    pos1 = line[0]
    pos2 = line[1]
    dist = line[2]
    if pos1 == pos2: continue
    if pos1 in poslist and pos2 in poslist:
      cate = '1'
    elif pos1 in poslist or pos2 in poslist:
      cate = '2'
    else:
      cate = '3'
    outfile.write(pos1+'\t'+pos2+'\t'+dist+'\t'+cate+'\n')
  outfile.close()
  infile.close()


if __name__ == '__main__':
  main()
