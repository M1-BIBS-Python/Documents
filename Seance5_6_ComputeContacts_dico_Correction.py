#!/usr/bin/env python

import random, math, numpy, string, sys
from pylab import *
import structureTools_Correction

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
            dij = structureTools_Correction.computeDist_dico(d_coordi, d_coordj, mode)
            distmat[i][j], distmat[j][i] = dij, dij
            j +=1

    return distmat            


def getResid(contactsList, dPDB) :

    resPairs = []
    
    for pair in contactsList :
        resPairs.append([dPDB["reslist"][pair[0]], dPDB["reslist"][pair[1]]])        

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
    mode = sys.argv[sys.argv.index("-mode")+1] # if mode = 'atom', compute the distance between all the atoms 
except:                                        # of the two residues and return the smallest distance 
    mode = "atom"                              # if mode = 'center', compute the distance between the two centers
                                               # of mass of the two residues and return it

try:
    outplot = sys.argv[sys.argv.index("-oplot")+1]

except:
    outplot = "distance_matrix.eps"

try:
    outname = sys.argv[sys.argv.index("-opairs")+1]

except:
    outname = "contactsPairs.txt"

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
dPDB = structureTools_Correction.parsePDBMultiChains(infile)
print "nb res ", len(dPDB["reslist"])

# computes the distances
print "computing the distances according to the", mode
matdist = computeContactMatrix(dPDB, mode)

# plots and stores the distances
pcolor(matdist)
print matdist[1]
print len(matdist[1])
savefig('contactMatrix_dicopg.eps',figsize=(4,4),dpi=100)

# extracts the contacts
contacts = structureTools_Correction.extractContactResidues(matdist, seuil)
resPairs = getResid(contacts, dPDB)

# writes the pairs in contact
fout = open(outname, "w")

for pairs in resPairs :
    fout.write("%s\t%s\n"%(pairs[0],pairs[1]))

fout.close()

