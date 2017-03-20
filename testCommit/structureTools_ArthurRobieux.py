#!/usr/bin/env  python
# -*- coding:utf8 -*-

"""
Auteur : Arthur ROBIEUX
E-mail : arthur.robieux@gmail.com

Objectif : Fonction permettant de parser une protéine au format .pdb !

Contenu fichier .pdb :
1-6 : name / 7-11 : number ID / 12-16 : Nom Atome / 17 : altLoc / 18 - 20 : Residu name / ID : 22
23-26 : resSeq / 31-38 : x / 39-46: y / 47-54 : z /

Hiérarchie de ce qu'on extrait : ID - resSeq - NomAtom - Coordonnés + n°id
"""

import math, string


def parsePDBMultiChains(infile):

    # lecture du fichier PDB
    f = open(infile, "r")
    lines = f.readlines()
    f.close()

    # var init
    chaine = True
    firstline = True
    prevres = None
    dPDB = {}
    dPDB["reslist"] = []
    dPDB["chains"] = []

    # parcoure le PDB
    for line in lines:
        if line[0:4] == "ATOM":
            chain = line[21]
            if not chain in dPDB["chains"]:
                dPDB["chains"].append(chain)
                dPDB[chain] = {}
                dPDB[chain]["reslist"] = []
            curres = "%s" % (line[22:26]).strip()
            if not curres in dPDB[chain]["reslist"]:
                dPDB[chain]["reslist"].append(curres)
                dPDB[chain][curres] = {}
                dPDB[chain][curres]["resname"] = string.strip(line[17:20])
                dPDB[chain][curres]["atomlist"] = []
            atomtype = string.strip(line[12:16])
            dPDB[chain][curres]["atomlist"].append(atomtype)
            dPDB[chain][curres][atomtype] = {}
            # print "cures ", curres
            # print dPDB[chain][curres]

            dPDB[chain][curres][atomtype]["x"] = float(line[30:38])
            dPDB[chain][curres][atomtype]["y"] = float(line[38:46])
            dPDB[chain][curres][atomtype]["z"] = float(line[46:54])
            dPDB[chain][curres][atomtype]["id"] = line[6:11].strip()

    return dPDB
