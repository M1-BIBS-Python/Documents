#!/usr/bin/env python
#-*- coding : utf8 -*-

"""
Author: Severine Liegeois
Contact: 
Date: 06/03/2017
Description: 
- A function parsing a PDB file into a dictionary
"""

def ParsingPDB (pdbFile):
	
	print("Ce programme vous permet de transformer un fichier .pdb au format ATOM en un dictionnaire manipulable par Python.")
	
	infile = open(pdbFile, "r")
	lines = infile.readlines() 						# cree une liste dont chaque element est une ligne du fichier

	molecule = {}									# dictionnaire le plus externe
	chainList = []
	rList = []
	

	for line in lines:
		if line[:4:] == 'ATOM':						# si la ligne commence par 'ATOM'
			
			chaine = line[21]
			if chaine not in chainList:
				chainList.append(chaine)
				molecule[chaine] = {}
				resList=[]
			
			res = line[23:26]
			if res not in resList :
				resList.append(res)
				molecule[chaine][res] = {}
				atomList=[]
				rList=[]
				
			atom = line[13:16]
			if atom not in atomList:
				atomList.append(atom)
				molecule[chaine][res][atom] = {}
				
			if line[17:20] not in rList:
				rList.append(line[17:20])
				
					
			molecule[chaine][res][atom]['x'] = line[31:38]
			molecule[chaine][res][atom]['y'] = line[39:46]
			molecule[chaine][res][atom]['z'] = line[47:54]
			molecule[chaine][res][atom]['id'] = line[7:11]
			
			molecule[chaine][res]['resname'] = rList
			molecule[chaine]['reslist'] = resList
			molecule[chaine][res]['atomlist'] = atomList
			
	infile.close()
	return(molecule)
	
print(ParsingPDB("1EJH.pdb"))


			



	
			
			
		
