#!/usr/bin/env python
#-*- coding : utf8 -*-

'''
Author: Marine Duhamel
Contact: marine.c.duhamel@gmail.com
Date: 27/02/2017
Usage: Distance map
Command : python TP5.py fichier.pdb > results.txt
'''







import sys #In order to use exit() and argv
import string
from math import sqrt

def PDB_Parsing(pdb_file):
	
	
	#Final dictionnary creation
	dddd_protein = {}
	#Intermediate list creation
	#l_chain stores the protein chain names
	l_chains = []
	#l_reslist contains the residus IDs
	l_reslist=[]
	#l_atoms contains the atom names in a single residu
	l_atoms=[]
	
	#File opening checking
	try:
		#In order to read the file line after line
		IN = open(pdb_file, "r")
		#Each line is contained in a list
		#Chaque ligne est contenu dans une liste
		line = IN.readlines()
	except : 
		print("Le fichier n'a ppu etre charge correcement. Verifier qu'il existe bien et relancez votre programme.")
		sys.exit(0)
	
	#For each element in the list, the programm find the lines from the file.txt which correspond to a residu
	for element in line:
		if element[:4:]=='ATOM':
			#If a chain appears for the first time, a new key named as and an associated sub-dictionnary it is created
			#And the residu list is reinitiated in order to be filled again with residu number associated to the new chain
			#And the chain name is added to the chain list
			if element[21] not in dddd_protein.keys():
				l_reslist=[]
				dddd_protein[element[21]] = {}
				l_chains.append(element[21])
			#if the residu is not already present in the residu sub-dictionnary keys, a new sub-dictionnary associated to this residu is created
			#And the atom list is reinitalised in order to be filled again with the new residu atoms associated
			if element[23:26] not in dddd_protein[element[21]].keys():
				res = element[23:26].strip()
				dddd_protein[element[21]][res] = {}	
				l_atoms=[]
				#If the residu is not present in the residu list, it is added
				if res not in l_reslist:
					l_reslist.append(res)
			#For each atom absent of the atom sub-dictionnary keys,
			#the concerned atom is added to the atom list
			#the residu name is stored in the residu sub-dictionnary, associated to the resname key
			#a new sub-dictionnary associated to this atom is created
			#the coordinates and atom's ID are added to this new sub-dictionnary
			#float() check that the strings contain float in order to avoid a mistake
			if element[13:16] not in dddd_protein[element[21]][element[23:26]].keys():
				l_atoms.append(element[13:16])
				dddd_protein[element[21]][res]['resname'] = element[17:20]
				dddd_protein[element[21]][res][element[13:16]] = {}
				dddd_protein[element[21]][res][element[13:16]]['x'] = element[32:38]
				dddd_protein[element[21]][res][element[13:16]]['y'] = element[40:46]
				dddd_protein[element[21]][res][element[13:16]]['z'] = element[48:54]
				dddd_protein[element[21]][res][element[13:16]]['id'] = element[23:26]
			
			#the atom list associated to the residu sub-dictionnary is updated with the intermediate atom list
			dddd_protein[element[21]][res]['atomlist'] = l_atoms
			#the residu list associated to the chain sub-dictionnary is updated with the intermediate residu list
			dddd_protein[element[21]]['reslist'] = l_reslist
			#the chain list associated to the dictionnary is updated with the intermediate chain list
			dddd_protein['chains'] = l_chains
		

	IN.close()
	
	return dddd_protein





'''
La suite ne fonctionne pas :
Traceback (most recent call last):
  File "TP5.py", line 135, in <module>
    print distance(dddd_PDB)
  File "TP5.py", line 106, in distance
    l_distance[res1] = []
TypeError: list indices must be integers, not str
'''




def distance(dddd_PDB):
	
	for res1 in dddd_PDB['A'].keys():
		if (res1 != 'reslist'):
			l_distance = []
			l_distance[res1] = []
			l_distance_all = []
			for res2 in dddd_PDB['A'].keys():
				if (res1 != res2 and res2 != 'reslist'): 
					for atom1 in dddd_PDB['A'][res1].keys():
						if (atom1 != 'resname' and atom1 != 'atomlist'):
							for atom2 in dddd_PDB['A'][res2].keys():
								if (atom2 != 'resname' and atom2 != 'atomlist'):
									x = (dddd_PDB['A'][res1][atom1]['x'] - dddd_PDB['A'][res2][atom2]['x']) **2
									y = (dddd_PDB['A'][res1][atom1]['y'] - dddd_PDB['A'][res2][atom2]['y']) **2
									z = (dddd_PDB['A'][res1][atom1]['z'] - dddd_PDB['A'][res2][atom2]['z']) **2
									l_distance_all.append(sqrt(x + y + z))
							d_distance[res1].append(min(l_distance_all))
			
	return d_distance
						
	
'''
def contact(liste_dist, valeur):
	
	for res1 in liste_dist:
		for res2 in list.dist[res1]:
			if dico_dist[res1][res2] <= valeur:
			...
'''	
					

if __name__ == '__main__':
	fichier = sys.argv[1]
	dddd_PDB = PDB_Parsing(fichier)
	print distance(dddd_PDB)
	


