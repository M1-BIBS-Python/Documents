#!/usr/bin/env python
#-*- coding : utf8 -*-

"""
Author: Severine Liegeois
Contact: 
Date: 06/03/2017
Description: 
- A function parsing a PDB file into a dictionary
"""

from math import sqrt
import string
import re

def ParsingPDB (pdbFile):
		
	infile = open(pdbFile, "r")
	lines = infile.readlines() 						# cree une liste dont chaque element est une ligne du fichier

	molecule = {}									# dictionnaire le plus externe
	chainList = []
	rList = []
	
	cptAlt = False

	for line in lines:
		if line[:4:] == 'ATOM':						# si la ligne commence par 'ATOM'
			
			if cptAlt == False:
				alt = line[17]
				cptAlt = True
				
			if line[17] == alt:
				chaine = line[21]
				
				if chaine not in chainList:
					chainList.append(chaine)
					molecule[chaine] = {}
					resList=[]
				
				curres = line[22:26].strip()
			
				if curres not in resList :
					resList.append(curres)
					molecule[chaine][curres] = {}
					atomList=[]
					rList=[]
				
				atom = line[13:16].strip()
				if atom not in atomList:
					atomList.append(atom)
					molecule[chaine][curres][atom] = {}
				
				if line[17:20] not in rList:
					rList.append(line[17:20].strip())
				
					
				molecule[chaine][curres][atom]['x'] = float(line[30:38])
				molecule[chaine][curres][atom]['y'] = float(line[38:46])
				molecule[chaine][curres][atom]['z'] = float(line[46:54])
				molecule[chaine][curres][atom]['id'] = line[6:11].strip()
			
				molecule[chaine][curres]['resname'] = rList
				molecule[chaine]['reslist'] = resList
				molecule[chaine][curres]['atomlist'] = atomList
			
	infile.close()
	return(molecule)
	

def ContactResidus(dico):
	
	atomicMass = dict()
	atomicMass = {"C": float(12.0107), "N": float(14.0067), 
				  "O": float(15.9994), "S": float(32.065),
				  "OH": float(17.0073), "NH": float(15.0146)}
	
	x = 0
	y = 0
	z = 0
	massTot = 0
	CM = {}
	
	for chain in dico.keys():
		for res in dico[chain].keys():
			if res != 'reslist':
				for atom in dico[chain][res].keys():
					alpha = dico[chain][res]['CA']
					if atom != 'resname' and atom != 'atomlist':
						for key in atomicMass.keys():
							if re.match(key, atom):
								mass = atomicMass[key]
						x += mass*dico[chain][res][atom]['x']
						y += mass*dico[chain][res][atom]['y']
						z += mass*dico[chain][res][atom]['z']
						massTot += mass

	CM['x'] = x/massTot
	CM['y'] = y/massTot
	CM['z'] = z/massTot
	return CM


			


if __name__ == '__main__':
	dico = ParsingPDB("../Data/arginine.pdb")
	print ContactResidus(dico)
	

	
			
			
		
