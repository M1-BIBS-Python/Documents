#!/usr/bin/env python

import random, math, numpy, string, sys, glob, os
from pylab import *
import structureTools, ForceField3

indir = sys.argv[1]

filelist = glob.glob("%s/*DP.pdb"%(indir))

for pdb in filelist :
    num = pdb.split("/")[-1].split("_")[-2]
    curpdb = "rec_nat_lig_%s.pdb"%(num)
    os.system("cat Rec_natif.pdb %s > %s"%(pdb, curpdb))
    os.system("python ../CompEner.py -pdb %s -chain1 B -chain2 D")
    os.system("cat Eners.out ener.out > Eners.out")
