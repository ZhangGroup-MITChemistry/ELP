import sys
import os
import math

oldfile='prod_100.gro'
old=open(oldfile,'r')
newfile=sys.argv[1]+'_100.gro'
new=open(newfile,'w')
line=old.readline()
while line:
    line2=line.replace("ION     CL","CLA    CLA")
    line3 = line2.replace("ION     NA", "SOD    SOD")
    new.write(line3)
    line=old.readline()
new.close()
old.close()