#!/usr/bin/env python

import random, math, numpy, string, sys
from pylab import *





def computeDist_dico(d_res1, d_res2, mode = "atom") :
    """res1, res2 are dico corresponding to residue 1 and residue 2 respectively """

    if mode == "atom" :
        minval = 1000000
        for atom1 in d_res1["atomlist"] :
            coord1 = [d_res1[atom1]["x"], d_res1[atom1]["y"], d_res1[atom1]["z"]]
            for atom2 in d_res2["atomlist"] :
                coord2 = [d_res2[atom2]["x"], d_res2[atom2]["y"], d_res2[atom2]["z"]]
                dist = distancePoints((coord1[0], coord1[1], coord1[2]),(coord2[0],coord2[1], coord2[2]))
                if minval > dist :
                    minval = dist

    elif mode == "center" : # computes the distance between the CM of the 2 given residues
        dPDBtmp = {}
        dPDBtmp["reslist"] = ["res1", "res2"]
        dPDBtmp["res1"] = d_res1
        dPDBtmp["res2"] = d_res2

        centerMassOfResidue(dPDBtmp)
        minval = distancePoints((dPDBtmp["res1"]["XCM"],dPDBtmp["res1"]["YCM"],dPDBtmp["res1"]["ZCM"]),(dPDBtmp["res2"]["XCM"],dPDBtmp["res2"]["YCM"],dPDBtmp["res2"]["ZCM"]))
        
    return minval


def extractContactResidues(matdist, seuil) :
    """from a distance matrix (matdist), returns pairs of residues in contacts (seuil) in a list of lists """

    contacts = []
    for i in range(len(matdist[0])) :
        for j in range (i+1, len(matdist[0])) :
            if matdist[i][j] <= seuil :
                      contacts.append([i, j])

    return contacts

def extractContactResidues2(matdist, seuil) :
    """from a distance matrix (matdist), returns pairs of residues in contacts (seuil) in a list of lists """

    contacts = []
    for i in range(len(matdist[0])) :
        for j in range(len(matdist[0])) :
            if matdist[i][j] <= seuil :
                      contacts.append([i, j])

    return contacts





def distancePoints((x1,y1,z1),(x2,y2,z2)):
    """Computes the distance between the two sets of coordinates
       input: 2 tuples with the corresponding coordinates 
       output: distance"""

    x = (x1-x2)
    y = (y1-y2)
    z = (z1-z2)
    return math.sqrt(x*x+y*y+z*z)



def centerMassOfResidue(dPDB, all = True, reslist = False):
    """Calculates the center of mass of each residue contained in dPDB (all = True & reslist = False) or a 
       subset of residues given in the residue list (["12_A", "13_A", "27_A"])"""

    if all == True :
        reslist = dPDB["reslist"]
    
        
    for res in reslist :        
        x = y = z = 0.0
        
        # looping over the current residue atoms
        for atom in dPDB[res]["atomlist"] :
            x +=dPDB[res][atom]["x"]
            y +=dPDB[res][atom]["y"]
            z +=dPDB[res][atom]["z"]
            
        Xcm = float(x)/len(dPDB[res]["atomlist"]) 
        Ycm = float(y)/len(dPDB[res]["atomlist"])
        Zcm = float(z)/len(dPDB[res]["atomlist"])
        dPDB[res]["XCM"] = Xcm
        dPDB[res]["YCM"] = Ycm
        dPDB[res]["ZCM"] = Zcm
        


def parsePDBMultiChains(infile) :

    # lecture du fichier PDB 
    f = open(infile, "r")
    lines = f.readlines()
    f.close()


    # var init
    chaine = True
    firstline = True
    prevres = None
    dPDB = {}
    dPDB["chains"] = []
    
    # parcoure le PDB   
    for line in lines :
        if line[0:4] == "ATOM" :
            chain = line[21]
            if not chain in dPDB["chains"] :
                dPDB["chains"].append(chain)
                dPDB[chain] = {}
                dPDB[chain]["reslist"] = []
            curres = "%s"%(line[22:26]).strip()
            if not curres in dPDB[chain]["reslist"] :
                dPDB[chain]["reslist"].append(curres)
                dPDB[chain][curres] = {}
                dPDB[chain][curres]["resname"] = string.strip(line[17:20])
                dPDB[chain][curres]["atomlist"] = []
            atomtype = string.strip(line[12:16])
            dPDB[chain][curres]["atomlist"].append(atomtype)
            dPDB[chain][curres][atomtype] = {}
            #print "cures ", curres
            #print dPDB[chain][curres]
 
            dPDB[chain][curres][atomtype]["x"] = float(line[30:38])
            dPDB[chain][curres][atomtype]["y"] = float(line[38:46])
            dPDB[chain][curres][atomtype]["z"] = float(line[46:54])
            dPDB[chain][curres][atomtype]["id"] = line[6:11].strip()

    return dPDB



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

    blablabla

  """



#####################################################################################
#
#                       MAIN
#
#####################################################################################

def compInter(dPDB, threshold, mode) :

    nbresInter = 0
    structureTools.initBfactor(dPDB)
    print dPDB["chains"]
    
    for chainidi in range(len(dPDB["chains"])-1) :
        chaini = dPDB["chains"][chainidi]
        chainidj = chainidi + 1
        chainj = dPDB["chains"][chainidj]
        print "dealing chains %s and %s"%(chaini, chainj)
        for chainidj in range(0, len(dPDB["chains"])):
            for resi in dPDB[chaini]["reslist"] :
                for resj in dPDB[chainj]["reslist"] :
                    dist = structureTools.computeDist_dico(dPDB[chaini][resi], dPDB[chainj][resj], mode = mode)
                    if dist <= threshold :# means, the two residues belong to the interface
                        if dPDB[chaini][resi]["bfactor"] == 0 :
                            nbresInter +=1
                        if dPDB[chainj][resj]["bfactor"] == 0 :
                            nbresInter +=1
                        dPDB[chaini][resi]["bfactor"] = 1
                        dPDB[chainj][resj]["bfactor"] = 1
                        #print "bfactor ", dPDB[chaini][resi]["bfactor"]
                        #print "%s %s VS %s %s = %s"%(chaini, resi, chainj, resj, dist)


    return nbresInter                       

#!/usr/bin/env python




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


     optional:
     =========

     -seuil   -> threshold to define a contact (in Angstrom)
     
     -mode    -> if mode = 'atom', compute the distance between all the atoms of the two 
                 residues and return the smallest distance.
                 if  mode = 'center', compute the distance between the two centers of mass
                 of the two residues and return it. (default = 'atom').

    -o        -> name of the output (default: distance_matrix.eps)

    -atomtype -> list that specifies the type of atoms that we want to consider in the distance calculation
                 example: ["CA", "C", "O", "N"], in this case, the distance will be compute only on CA, C, O and N
                 atoms. (default: ["all"]). Default value must be applied if the mode = "center".

  """


def computeContactMatrix(d_coords, mode) :
    """ input : list of lists which contains the coords of every atoms
    output : distance matrix
    """
    nbres = len(d_coords["reslist"])
    distmat = numpy.zeros((nbres,nbres))
    
    for i in range(nbres) :
        j = i + 1
        d_coordi = d_coords[d_coords["reslist"][i]]
        while j < nbres :
            d_coordj = d_coords[d_coords["reslist"][j]]
            dij = computeDist_dico(d_coordi, d_coordj, mode)
            distmat[i][j], distmat[j][i] = dij, dij
            j +=1

    return distmat 

def computeContactMatrix2(d_coords1,d_coords2, mode) :
    """ input : list of lists which contains the coords of every atoms
    output : distance matrix
    """
    nbres1 = len(d_coords1["reslist"])
    nbres2 = len(d_coords2["reslist"])
    distmat = numpy.zeros((nbres1,nbres2))
    
    for i in range(nbres1) :
        d_coordi = d_coords1[d_coords1["reslist"][i]]
        for j in range(nbres2) :
            d_coordj = d_coords2[d_coords2["reslist"][j]]
            dij = computeDist_dico(d_coordi, d_coordj, mode)
            distmat[i][j], distmat[j][i] = dij, dij

    return distmat


def getResid(contactsList, dPDB) :

    resPairs = []
    
    for pair in contactsList :
        resPairs.append([dPDB["reslist"][pair[0]], dPDB["reslist"][pair[1]]])        

    return resPairs

def getResid2(contactsList, dPDB1,dPDB2) :

    resPairs = []
    
    for pair in contactsList :
        resPairs.append([dPDB1["reslist"][pair[0]], dPDB2["reslist"][pair[1]]])        

    return resPairs



#####################################################################################
#
#                       MAIN
#
#####################################################################################




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
    chain1 = sys.argv[sys.argv.index("-chain1")+1]
except:                                        
    chain1 = "A"
try:
    chain2 = sys.argv[sys.argv.index("-chain2")+1]
except:                                        
    chain2 = "B"
try:      
    mode = sys.argv[sys.argv.index("-mode")+1] # if mode = 'atom', compute the distance between all the atoms 
except:                                        # of the two residues and return the smallest distance 
    mode = "atom"                              # if mode = 'center', compute the distance between the two centers
                                               # of mass of the two residues and return it

try:
    outplot1 = sys.argv[sys.argv.index("-outplot1")+1]

except:
    outplot1 = "distance_matrix1.eps"
try:
    outplot2 = sys.argv[sys.argv.index("-outplot2")+1]

except:
    outplot2 = "distance_matrix2.eps"

try:
    outname1 = sys.argv[sys.argv.index("-opairs1")+1]

except:
    outname1 = "contactsPairs1.txt"
try:
    outname2 = sys.argv[sys.argv.index("-opairs2")+1]

except:
    outname2 = "contactsPairs2.txt"

try:
    seuil = float(sys.argv[sys.argv.index("-seuil")+1])

except:
    seuil = 5.0


try:
    atomlist = sys.argv[sys.argv.index("-atomtype")+1]
    atomtype = []
    atomtype = atomlist.split("[")[1].split("]")[0].split(",")

except:
    atomtype = ["all"]





# Computes distances between every residues, stores them in distmat (array) and returns a contact matrix in eps format
#=======================================================================================

# parses the pdb file
dPDB = parsePDBMultiChains(infile)


# computes the distances
matdist1 = computeContactMatrix(dPDB[chain1], mode)
matdist2 = computeContactMatrix2(dPDB[chain1],dPDB[chain2], mode)
# plots and stores the distances
pcolor(matdist1)
print matdist1[1]
print len(matdist1[1])
savefig('contactMatrix_dicopg1.eps',figsize=(4,4),dpi=100)
pcolor(matdist2)
print matdist2[1]
print len(matdist2[1])
savefig('contactMatrix_dicopg2.eps',figsize=(4,4),dpi=100)
# extracts the contacts
contacts1 = extractContactResidues(matdist1, seuil)
contacts2 = extractContactResidues2(matdist2, seuil)
resPairs1 = getResid(contacts1, dPDB[chain1])
resPairs2 = getResid2(contacts2, dPDB[chain1],dPDB[chain2])

# writes the pairs in contact
fout1 = open(outname1, "w")

for pairs in resPairs1 :
    fout1.write("%s\t%s\n"%(pairs[0],pairs[1]))

fout1.close()
print outname2
print resPairs2
fout2 = open(outname2, "w")

for pairs in resPairs2 :
    fout2.write("%s\t%s\n"%(pairs[0],pairs[1]))

fout2.close()
