from pymol import cmd

cmd.load('ref/4ll3.pdb')
# open dist.txt for writing
f=open('analysis/dist.txt','w')

for res1 in range(1,100):
  for res2 in range(1,100):
# calculate the distance and store it in dst
    dst=cmd.distance('tmp','A/'+str(res1)+'/ca','A/'+str(res2)+'/ca')

# write the formatted value of the distance (dst)
# to the output file
    f.write(str(res1)+'\t'+str(res2)+'\t'+"%8.3f\n"%dst)

# close the output file.
f.close()
