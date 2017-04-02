#!/usr/bin/env python

import random, math, numpy, string, sys
from pylab import *
import structureTools, ForceField3



########################################################################################################
#
#                   FUNCTIONS
#
########################################################################################################


def usage():
    print """

     obligatory:
     ===========

     -pdb     -> pdb file


  """



#####################################################################################
#
#                       MAIN
#
#####################################################################################

def compInter(dPDB, threshold, mode, chain1, chain2) :

    nbresInter = 0
    structureTools.initBfactor(dPDB)
    print dPDB["chains"]
    
    # verifie que la chaine1 existe
    if not chain1 in dPDB["chains"] :
        print "la chaine %s n'existe pas"%(chain1)
        sys.exit()
    # verifie que la chaine2 existe
    if not chain2 in dPDB["chains"] :
        print "la chaine %s n'existe pas"%(chain2)
        sys.exit()
        
    print "dealing chains %s and %s"%(chain1, chain2)
    
    for resi in dPDB[chain1]["reslist"] :
        for resj in dPDB[chain2]["reslist"] :
            dist = structureTools.computeDist_dico(dPDB[chain1][resi], dPDB[chain2][resj], mode = mode)
            if dist <= threshold :# means, the two residues belong to the interface
                if dPDB[chain1][resi]["bfactor"] == 0 :
                    nbresInter +=1
                    dPDB[chain1][resi]["bfactor"] = 1
                if dPDB[chain2][resj]["bfactor"] == 0 :
                    nbresInter +=1
                    dPDB[chain2][resj]["bfactor"] = 1
                print "bfactor ", dPDB[chain1][resi]["bfactor"]
                print "%s %s VS %s %s = %s"%(chain1, resi, chain2, resj, dist)


    return nbresInter                       



def compEner(dPDB, chain1, chain2) :

    nbresInter = 0
    #structureTools.initBfactor(dPDB)
    #print dPDB["chains"]
    
    # verifie que la chaine1 existe
    if not chain1 in dPDB["chains"] :
        print "la chaine %s n'existe pas"%(chain1)
        sys.exit()
    # verifie que la chaine2 existe
    if not chain2 in dPDB["chains"] :
        print "la chaine %s n'existe pas"%(chain2)
        sys.exit()
        
    print "dealing chains %s and %s"%(chain1, chain2)

    # DEBUG
    cmt = 0
    cmti = 0
    cmtj = 0
    Eij = 0
    
    # computes ener
    for resi in dPDB[chain1]["reslist"] :
        #print "resi ", resi 
        for atomi in dPDB[chain1][resi]["atomlist"] :
            cmti+=1
            cmtj = 0
            coordi = [dPDB[chain1][resi][atomi]["x"], dPDB[chain1][resi][atomi]["y"], dPDB[chain1][resi][atomi]["z"]]

            for resj in dPDB[chain2]["reslist"] :
                #print "resj" , resj
                for atomj in dPDB[chain2][resj]["atomlist"] :
                    cmt+=1
                    cmtj+=1
                    coordj = [dPDB[chain2][resj][atomj]["x"], dPDB[chain2][resj][atomj]["y"], dPDB[chain2][resj][atomj]["z"]]

                    # computes dist, Aij & Bij
                    
                    Rij = structureTools.distancePoints(coordi, coordj)
                    #print "%s %s %s VS %s %s %s"%(chain1, resi, atomi, chain2, resj, atomj)
                    #print dPDB[chain1][resi][atomi]
                    #print dPDB[chain2][resj][atomj]
                    Aij = computeAij(dPDB[chain1][resi][atomi], dPDB[chain2][resj][atomj])
                    Bij = computeBij(dPDB[chain1][resi][atomi], dPDB[chain2][resj][atomj])
                    #print "%s %s %s VS %s %s %s= %s %s %s"%(chain1, resi, atomi, chain2, resj, atomj, Rij, Aij, Bij)
                    coulij =  (332.0522*dPDB[chain1][resi][atomi]["charge"]*dPDB[chain2][resj][atomj]["charge"])/(20*Rij)
                    #print "%s %s %s VS %s %s %s= Rij %s Aij %s Bij %s coulij %s"%(chain1, resi, atomi, chain2, resj, atomj, Rij, Aij, Bij, coulij)

                    # computes ENER
                    Etmp = float(Aij)/pow(Rij,8) - float(Bij)/pow(Rij, 6) + coulij
                    #print "Etmp ", Etmp
                    Eij = Eij + Etmp
                    
    #print cmti, cmtj, cmt

    return Eij                       

def computeAij(d_atomi, d_atomj) :

    Aij = math.sqrt(d_atomi["epsilon"]*d_atomj["epsilon"])*(pow((d_atomi["vdw"]+d_atomj["vdw"]),8))

    return Aij

def computeBij(d_atomi, d_atomj) :
    Bij = 2*math.sqrt(d_atomi["epsilon"]*d_atomj["epsilon"])*(pow((d_atomi["vdw"]+d_atomj["vdw"]),6))

    return Bij


def computeChargesVdW(dPDB, dcharge, dvdw, depsilon, chain):

    first = True
    
    for resi in dPDB[chain]["reslist"] :
        
        # means this is the Nter residue
        if first and (dPDB[chain][resi].has_key("H1") or dPDB[chain][resi].has_key("H2") or dPDB[chain][resi].has_key("H3")) :

                #print "ENTER IN TEST Nterminal "
                #print dPDB[chain][resi]
                
                for atomi in dPDB[chain][resi]["atomlist"] :
                    if atomi == "H1" :
                        dPDB[chain][resi]["H1"]["charge"] = 0.1984
                        dPDB[chain][resi]["H1"]["vdw"] = 0.6
                        dPDB[chain][resi]["H1"]["epsilon"] = 0.0157


                    elif atomi == "H2" :
                        dPDB[chain][resi]["H2"]["charge"] = 0.1984
                        dPDB[chain][resi]["H2"]["vdw"] = 0.6
                        dPDB[chain][resi]["H2"]["epsilon"] = 0.0157

                    elif atomi == "H3" :
                        dPDB[chain][resi]["H3"]["charge"] = 0.1984   
                        dPDB[chain][resi]["H3"]["vdw"] = 0.6
                        dPDB[chain][resi]["H3"]["epsilon"] = 0.0157

                    elif atomi == "N" :
                        dPDB[chain][resi]["N"]["charge"] = 0.1592
                        dPDB[chain][resi]["N"]["vdw"] = 1.875
                        dPDB[chain][resi]["N"]["epsilon"] =  0.17

                    elif atomi == "CA" :
                        dPDB[chain][resi]["CA"]["charge"] = 0.0221
                        dPDB[chain][resi]["CA"]["vdw"] = 1.9080
                        dPDB[chain][resi]["CA"]["epsilon"] = 0.1094

                    elif atomi == "HA" :
                        dPDB[chain][resi]["HA"]["charge"] = 0.116
                        dPDB[chain][resi]["HA"]["vdw"] = 1.1
                        dPDB[chain][resi]["HA"]["epsilon"] = 0.0157

                    else:
                        dPDB[chain][resi][atomi]["charge"] = dcharge[dPDB[chain][resi]["resname"]][atomi]
                        dPDB[chain][resi][atomi]["vdw"] =  dvdw[dPDB[chain][resi]["resname"]][atomi]
                        dPDB[chain][resi][atomi]["epsilon"] =  depsilon[dPDB[chain][resi]["resname"]][atomi]
                        
                    
                first = False
            
        elif first == False and (dPDB[chain][resi].has_key("OXT")):

            for atomi in dPDB[chain][resi]["atomlist"] :

                if atomi == "CA" :
                    dPDB[chain][resi]["CA"]["charge"] = -0.2493
                    dPDB[chain][resi]["CA"]["vdw"] = 1.9080
                    dPDB[chain][resi]["CA"]["epsilon"] = 0.1094

                elif atomi == "C" :
                    dPDB[chain][resi]["C"]["charge"] = 0.7231
                    dPDB[chain][resi]["C"]["vdw"] = 1.9080
                    dPDB[chain][resi]["C"]["epsilon"] = 0.0860
                    
                elif atomi == "O" :
                    dPDB[chain][resi]["O"]["charge"] = -0.7855
                    dPDB[chain][resi]["O"]["vdw"] = 1.6612
                    dPDB[chain][resi]["O"]["epsilon"] = 0.2100

                elif atomi == "OXT" :
                    dPDB[chain][resi]["OXT"]["charge"] = -0.7855
                    dPDB[chain][resi]["OXT"]["vdw"] = 1.6612
                    dPDB[chain][resi]["OXT"]["epsilon"] = 0.2100
                    
                else:
                    dPDB[chain][resi][atomi]["charge"] = dcharge[dPDB[chain][resi]["resname"]][atomi]
                    dPDB[chain][resi][atomi]["vdw"] =  dvdw[dPDB[chain][resi]["resname"]][atomi]
                    dPDB[chain][resi][atomi]["epsilon"] =  depsilon[dPDB[chain][resi]["resname"]][atomi]

        
        else :
            for atomi in dPDB[chain][resi]["atomlist"] :
                #print "res %s atom %s resname %s"%(resi, atomi, dPDB[chain][resi]["resname"])

                if dPDB[chain][resi]["resname"] == "HIS" and dPDB[chain][resi].has_key("HD1") :
                    # means this is a HID
                    dPDB[chain][resi]["resname"] = "HID"
                    
                elif dPDB[chain][resi]["resname"] == "HIS" and dPDB[chain][resi].has_key("HE2") :
                    # means this is a HIE
                    dPDB[chain][resi]["resname"] = "HIE"
                    
                dPDB[chain][resi][atomi]["charge"] = dcharge[dPDB[chain][resi]["resname"]][atomi]
                dPDB[chain][resi][atomi]["vdw"] =  dvdw[dPDB[chain][resi]["resname"]][atomi]
                dPDB[chain][resi][atomi]["epsilon"] =  depsilon[dPDB[chain][resi]["resname"]][atomi]

    #print "dico N ", dPDB["B"][dPDB["B"]["reslist"][0]]["N"]
    #print "dico CA res 2", dPDB["B"][dPDB["B"]["reslist"][1]]["CA"]
    #sys.exit()
    #return dPDB



# Get Arguments
#===============

try:
    infile = sys.argv[sys.argv.index("-pdb")+1]
    print "pdb to treat:", infile
except:    
    usage()
    print "ERROR: please, enter the name of the pdb input"
    sys.exit()

try:
    threshold = float(sys.argv[sys.argv.index("-th")+1])
except:    
    threshold = 5.0

try:
    chain1 = sys.argv[sys.argv.index("-chain1")+1]
except:    
    chain1 = "A"


try:
    chain2 = sys.argv[sys.argv.index("-chain2")+1]
except:    
    chain2 = "B"

    
try:
    mode = sys.argv[sys.argv.index("-mode")+1]
    print "the mode of calculation for the distance is:", mode
except:    
    mode = "atom"

# Computes distances between the two proteins
#============================================

# get param

dcharge = ForceField3.chargePDB()
dvdw, depsilon = ForceField3.epsilon_vdw_PDB()

# parses the pdb file
dPDB = structureTools.parsePDBMultiChains(infile)
#print "the distance threshold is:", threshold

print "nb of chains: %s"%(len(dPDB["chains"]))
#print dPDB[dPDB["chains"][0]]["reslist"]
print "nb res chain %s: %s"%(chain1, len(dPDB[chain1]["reslist"]))
print "nb res chain %s: %s"%(chain2, len(dPDB[chain2]["reslist"]))


computeChargesVdW(dPDB, dcharge, dvdw, depsilon, chain1)
computeChargesVdW(dPDB, dcharge, dvdw, depsilon, chain2)

#print dPDB[dPDB["chains"][0]][dPDB[dPDB["chains"][0]]["reslist"][0]]

#print "dico N after loop", dPDB["B"][dPDB["B"]["reslist"][0]]["N"]
#print "dealing chain %s and chain %s"%(chain1, chain2)


#sys.exit()
Ener = compEner(dPDB, chain1, chain2)

out = open("ener.out","w")
out.write("%s %s\n"%(infile, Ener))
out.close()
print "Ener ", Ener

sys.exit()

#####
### computes charges and vdw of each atom




# computes interface

nbresInter = compInter(dPDB, threshold, mode, chain1, chain2)

for chain in dPDB["chains"]:
    for res in dPDB[chain]["reslist"] :
        print "res %s bfactor %s"%(res, dPDB[chain][res]["bfactor"])

print "nb of residues in the interface:",nbresInter

structureTools.writePDB(dPDB, filout = "interface.pdb", bfactor = True) 
