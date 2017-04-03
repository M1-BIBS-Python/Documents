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
import matplotlib.pyplot as plt
import numpy as np

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
				alt = line[16]
				cptAlt = True
				
			if line[16] == alt:
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
					molecule[chaine][curres]['resname'] = line[17:20].strip()
				
				atom = line[13:16].strip()
				if atom not in atomList:
					atomList.append(atom)
					molecule[chaine][curres][atom] = {}
				
				
				
				molecule[chaine]['reslist'] = resList
				molecule[chaine][curres]['atomlist'] = atomList
					
				molecule[chaine][curres][atom]['x'] = float(line[30:38])
				molecule[chaine][curres][atom]['y'] = float(line[38:46])
				molecule[chaine][curres][atom]['z'] = float(line[46:54])
				molecule[chaine][curres][atom]['id'] = line[6:11].strip()
			
				
			
	infile.close()
	return(molecule)
	

def CentreMasse(dico):
	
	atomicMass = dict()
	atomicMass = {"C": float(12.0107), "N": float(14.0067), 
				  "O": float(15.9994), "S": float(32.065),
				  "OH": float(17.0073), "NH": float(15.0146)}
	
	x = 0
	y = 0
	z = 0
	massTot = 0
	CM_res = {}
	for atom in dico.keys():
		if atom != 'resname' and atom != 'atomlist':
			for key in atomicMass.keys():
				if re.match(key, atom):
					mass = atomicMass[key]
			x += mass*dico[atom]['x']
			y += mass*dico[atom]['y']
			z += mass*dico[atom]['z']
			massTot += mass
	# coordonnees du centre de masse du residu:
	CM_res['x'] = x/massTot
	CM_res['y'] = y/massTot
	CM_res['z'] = z/massTot

	return CM_res

def CentreGravite(dico):
	x = 0
	y = 0
	z = 0
	nbAtomes = 0
	CG_res = {}
	for atom in dico.keys():
		if atom != 'resname' and atom != 'atomlist' and atom != 'temp':
			x += dico[atom]['x']
			y += dico[atom]['y']
			z += dico[atom]['z']
			nbAtomes += 1
	# coordonnees du centre de gravite du residu:
	CG_res['x'] = x/nbAtomes
	CG_res['y'] = y/nbAtomes
	CG_res['z'] = z/nbAtomes
	
	return CG_res

# distance entre 2 points dans l'espace:
def Distance(res1, res2):
	return sqrt((res1['x']-res2['x'])*(res1['x']-res2['x'])+(res1['y']-res2['y'])*(res1['y']-res2['y'])+(res1['z']-res2['z'])*(res1['z']-res2['z']))

# retourne la matrice des distances 2 a 2 des residus:
def MatriceDistances(dico):
	seuil = 4
	for chain in dico.keys():
		distMatrix = []
		enContact = {}
		for res1 in dico[chain].keys():
			if res1 != 'reslist':
				r1 = CentreGravite(dico[chain][res1])
				listr1 = [] # liste qui contiendra les distances entre r1 et tous les autres residus
				for res2 in dico[chain].keys():
					if res2 != 'reslist':
						r2 = CentreGravite(dico[chain][res2])
						listr1.append(Distance(r1, r2))
						if (Distance(r1,r2) <= seuil and Distance(r1,r2) != 0 and (res2, res1) not in enContact.keys()):
							enContact[res1,res2] = Distance(r1, r2)
				distMatrix.append(listr1)
				
		
		data = np.array(distMatrix).reshape(len(listr1), len(listr1))
		heatmap = plt.pcolor(data, cmap="RdBu_r")
		
		
		plt.show() # montre la heatmap dans une fenetre
		#plt.savefig("/home/tp-home007/sliegeo/Documents/Python/heatmap.png")
		return enContact


# extrait les residus de l'interface d'un complexe proteique
def InterfaceComplexe(dico, seuil):
	# calcul des centres de gravite de chaque residu du complexe
	for chain in dico.keys():
		for res in dico[chain].keys():
			if res != 'reslist':
				dico[chain][res]['centre'] = CentreGravite(dico[chain][res])
	
	# calcul de la distance entre les residus de la chaine 1 et ceux de la chaine 2
	distance = {}
	nbCharges = 0
	nbPolaires = 0
	nbHydrophobes = 0
	nbContactC = 0
	nbContactP = 0
	nbContactHH = 0
	for chain1 in dico.keys():
		for chain2 in dico.keys():
			if chain1 != chain2:
				for res1 in dico[chain1].keys():
					for res2 in dico[chain2].keys():
						if (res1 != 'reslist' and res2 != 'reslist'):
							d = Distance(dico[chain1][res1]['centre'], dico[chain2][res2]['centre'])
							
							if d <= seuil:
								dico[chain1][res1]['temp'] = 1
								dico[chain2][res2]['temp'] = 1
								distance[(chain1, res1), (chain2, res2)] = d
								if ChargeResidu(res1):
									nbCharges += 1
								if ChargeResidu(res2):
									nbCharges += 1
								if (ChargeResidu(res1) and ChargeResidu(res2)):
									nbContactC += 1
									
								if PolariteResidu(res1):
									nbPolaires += 1
								if PolariteResidu(res2):
									nbPolaires += 1
								if (PolariteResidu(res1) and PolariteResidu(res2)):
									nbContactP += 1
								
								if HydrophobiciteResidu(res1):
									nbHydrophobes += 1
								if HydrophobiciteResidu(res2):
									nbHydrophobes += 1
								if ((HydrophobiciteResidu(res1) and not(HydrophobiciteResidu(res2))) or (not(HydrophobiciteResidu(res1)) and HydrophobiciteResidu(res2))):
									nbContactHH += 1
								
							else:
								dico[chain1][res1]['temp'] = 0
								dico[chain2][res2]['temp'] = 0
	return dico

						
"""
	print "Nombre residus charges : " + str(nbCharges)
	print "Nombre residus polaires : " + str(nbPolaires)
	print "Nombre residus hydrophobes : " + str(nbHydrophobes)
	print "Nombre de contacts polaire-polaire : " + str(nbContactP)
	print "Nombre de contacts charge-charge : " + str(nbContactC)
	print "Nombre de contacts hydrophobe-hydrophile : " + str(nbContactHH)		
	"""	
	
def ChargeResidu(res):
	resCharges = ["GLU", "HIS", "LYS", "ASP", "ARG"]
	if res.upper() in resCharges:
		return True
	else:
		return False
		
def PolariteResidu(res):
	resPolaires = ["ASP", "GLU", "HIS", "LYS", "ASN", "GLN", "ARG", "SER", "TYR"]
	if res.upper() in resPolaires:
		return True
	else:
		return False

def HydrophobiciteResidu(res):
	resHydrophobes = ["ALA", "PHE", "GLY", "ILE", "LEU", "MET", "PRO", "VAl"]
	if res.upper() in resHydrophobes:
		return True
	else:
		return False



if __name__ == '__main__':
	dico = ParsingPDB("/home/tp-home007/sliegeo/Documents/Python/2w18.pdb")
	#print dico
	"""for chain in dico.keys():
		for res in dico[chain].keys():
			if res != 'reslist':
				for key in dico[chain][res].keys():
					print key"""
	#print "Centre de masse : " + str(CentreMasse(dico))
	#print "Centre de gravite : "+ str(CentreGravite(dico))
	#print MatriceDistances(dico)
	print InterfaceComplexe(dico,5)
	
			
			
		
