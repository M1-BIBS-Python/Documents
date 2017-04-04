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
import sys
import os

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
	
#####################################
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

#######################################
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

########################################
# extrait les residus de l'interface d'un complexe proteique
def InterfaceComplexe(dico, seuil):
	# calcul des centres de gravite de chaque residu du complexe
	for chain in dico.keys():
		for res in dico[chain].keys():
			if res != 'reslist':
				dico[chain][res]['centre'] = CentreGravite(dico[chain][res])
	
	# calcul de la distance entre les residus de la chaine 1 et ceux de la chaine 2
	distance = dict()
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
								if ChargeResidu(dico[chain1][res1]['resname']):
									nbCharges += 1
								if ChargeResidu(dico[chain1][res2]['resname']):
									nbCharges += 1
								if (ChargeResidu(dico[chain1][res1]['resname']) and ChargeResidu(dico[chain1][res2]['resname'])):
									nbContactC += 1
									
								if PolariteResidu(dico[chain1][res1]['resname']):
									nbPolaires += 1
								if PolariteResidu(dico[chain1][res2]['resname']):
									nbPolaires += 1
								if (PolariteResidu(dico[chain1][res1]['resname']) and PolariteResidu(dico[chain1][res2]['resname'])):
									nbContactP += 1
								
								if HydrophobiciteResidu(dico[chain1][res1]['resname']):
									nbHydrophobes += 1
								if HydrophobiciteResidu(dico[chain1][res2]['resname']):
									nbHydrophobes += 1
								if ((HydrophobiciteResidu(dico[chain1][res1]['resname']) and not(HydrophobiciteResidu(dico[chain1][res2]['resname']))) or (not(HydrophobiciteResidu(dico[chain1][res1]['resname'])) and HydrophobiciteResidu(dico[chain1][res2]['resname']))):
									nbContactHH += 1
								
							else:
								dico[chain1][res1]['temp'] = 0
								dico[chain2][res2]['temp'] = 0	
	
	"""		
	print "Nombre residus charges : " + str(nbCharges)
	print "Nombre residus polaires : " + str(nbPolaires)
	print "Nombre residus hydrophobes : " + str(nbHydrophobes)
	print "Nombre de contacts polaire-polaire : " + str(nbContactP)
	print "Nombre de contacts charge-charge : " + str(nbContactC)
	print "Nombre de contacts hydrophobe-hydrophile : " + str(nbContactHH)"""	
	return dico
	
	
	
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

###########################################
# Script test partie 1 : calcul du RMSD entre les 3 atomes de deux structures
def RMSD_part1(dico1, dico2):
	nb_pairs = 0
	somme = 0
	for chain in dico1.keys():
		for res in dico1[chain].keys():
			if (res != 'reslist' and res != 'reslist'):
				for atom in dico1[chain][res].keys():
					if (atom != 'atomlist'  and atom != 'resname'):
						d = Distance(dico1[chain][res][atom], dico2[chain][res][atom])
						somme += d*d
						nb_pairs += 1
	
	rmsd = sqrt(somme/nb_pairs)
	return rmsd

# Partie 2 : calcul du RMSD entre deux proteines entieres
def RMSD_part2(file1, file2, calc_atom):
	dico1 = ParsingPDB(file1)
	dico2 = ParsingPDB(file2)
	
	nb_pairs = 0
	nb_res = 0
	somme = 0
	for chain in dico1.keys():
		for res in dico1[chain].keys():
			if (res != 'reslist' and res != 'reslist'):
				nb_res += 1
				for atom in dico1[chain][res].keys():
					if re.match(calc_atom, atom):
						d = Distance(dico1[chain][res][atom], dico2[chain][res][atom])
						somme += d*d
						nb_pairs += 1
	
	rmsd = sqrt(somme/nb_pairs)
	return rmsd # pour la partie 3
	"""
	print "Nombre d'atomes pour le calcul du RMSD : " + str(nb_pairs)
	print "\nRMSD = " + str(rmsd)
	print "\nNombre de residus de la proteine : " + str(nb_res)
	print "\n"
	"""
	

def usage():
	print "Obligatory arguments : -pdb1 <first pdb file> -pdb2 <second pdb file>\n"
	print "Optional arguments (atom used in the RMSD calculation) : -atom <name of atom>\n"

# Partie 3 : calcule automatiquement le RMSD entre toutes les paires de structures stockees dans le repertoire "structures"
def RMSD_part3(dossier, atom):
	outfile = open("rmsd.txt", "w")
	for file1 in os.listdir(dossier):
		if not file1.startswith('.'):
			for file2 in os.listdir(dossier):
				if not file2.startswith('.'):
					if file1 != file2:
						rmsd = RMSD_part2(os.path.join(dossier ,file1), os.path.join(dossier,file2), atom)
						outfile.write(str(file1)+"\t"+ str(file2)+"\t"+str(rmsd)+"\n")

	outfile.close()

#def usage():
	usage = """ 
	Computes the RMSD between all the pairs of protein structures (in pdb format) contained in a directory.
	Input = a directory containing files of superposed protein structures
	Output = a text file with 3 columns: respectively the names of the 2 pdb files used for the calculation
			 and the corresponding RMSD
	
	Obligatory argument : -dir <path of the directory>
						  -> absolute or relative path
						  
	Optional argument : -atom <name of atom>
						-> the name of the atom used in the RMSD calculation
						-> default : "CA"
	"""
	#print usage


# Partie 4 : modifier le script de la partie 2 de facon a calculer le RMSD independamment pour chaque residu de la proteine
def RMSD_part4(file1, file2, calc_atom):
	dico1 = ParsingPDB(file1)
	dico2 = ParsingPDB(file2)
	
	RMSD = []
	list_res = []
	
	outfile = open("rmsd_part4.txt", "w")
	outfile.write("Residue\tRMSD\n")
		
	for chain in dico1.keys():
		for res in dico1[chain].keys():
			nb_pairs = 0
			somme = 0
			if (res != 'reslist' and res != 'reslist'):
				for atom in dico1[chain][res].keys():
					if re.match(calc_atom, atom):
						d = Distance(dico1[chain][res][atom], dico2[chain][res][atom])
						somme += d*d
						nb_pairs += 1
				rmsd = sqrt(somme/nb_pairs)
				outfile.write(str(res)+"\t"+str(rmsd)+"\n")
				RMSD.append(rmsd)
				list_res.append(res)
				
	outfile.close()
	x = np.array(list_res)
	y = np.array(RMSD)
	plt.plot(x, y)
	plt.show()
	











if __name__ == '__main__':
	#dico = ParsingPDB("/home/tp-home007/sliegeo/Documents/Python/2w18.pdb")
	#print dico
	#print "Centre de masse : " + str(CentreMasse(dico))
	#print "Centre de gravite : "+ str(CentreGravite(dico))
	#print MatriceDistances(dico)
	#print InterfaceComplexe(dico,5)
	#dico_rouge = ParsingPDB("/home/tp-home007/sliegeo/Documents/Python/Seance8_RMSD_RayonG/Seance8_data/rouge.pdb")
	#dico_bleu = ParsingPDB("/home/tp-home007/sliegeo/Documents/Python/Seance8_RMSD_RayonG/Seance8_data/bleu.pdb")
	#print RMSD_part1(dico_rouge, dico_bleu)
	
	
	try:
		infile1 = sys.argv[sys.argv.index("-pdb1")+1]
		infile2 = sys.argv[sys.argv.index("-pdb2")+1]
		print "pdb files to treat:", infile1, infile2
	except:    
		usage()
		print "ERROR: please enter the names of the pdb files"
		sys.exit()
	try:
		atom = sys.argv[sys.argv.index("-atom")+1]
	except:
		atom = "CA"
		
	
	#RMSD_part2(infile1, infile2, atom)
	"""
	try:
		input_dir = sys.argv[sys.argv.index("-dir")+1]
	except:
		usage()
		print "ERROR: please enter the path of the directory containing the pdb files"
		sys.exit()
	try:
		atom = sys.argv[sys.argv.index("-atom")+1]
	except:
		atom = "CA"
		"""
	#RMSD_part3(input_dir, atom)
	
	RMSD_part4(infile1, infile2, atom)
