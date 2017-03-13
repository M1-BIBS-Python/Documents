#!/usr/bin/env python
#-*- coding : utf8 -*-

'''
Author: Marine Duhamel
Contact: marine.c.duhamel@gmail.com
Date: 27/02/2017
Usage: PDB file parsing
Command : python TP3.py fichier.pdb
'''



import sys #In order to use exit() and argv


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
				dddd_protein[element[21]][element[23:26]] = {}	
				l_atoms=[]
				#If the residu is not present in the residu list, it is added
				if element[23:26] not in l_reslist:
					l_reslist.append(element[23:26])
			#For each atom absent of the atom sub-dictionnary keys,
			#the concerned atom is added to the atom list
			#the residu name is stored in the residu sub-dictionnary, associated to the resname key
			#a new sub-dictionnary associated to this atom is created
			#the coordinates and atom's ID are added to this new sub-dictionnary
			#float() check that the strings contain float in order to avoid a mistake
			if element[13:16] not in dddd_protein[element[21]][element[23:26]].keys():
				l_atoms.append(element[13:16])
				dddd_protein[element[21]][element[23:26]]['resname'] = element[17:20]
				dddd_protein[element[21]][element[23:26]][element[13:16]] = {}
				dddd_protein[element[21]][element[23:26]][element[13:16]]['x'] = element[32:38]
				dddd_protein[element[21]][element[23:26]][element[13:16]]['y'] = element[40:46]
				dddd_protein[element[21]][element[23:26]][element[13:16]]['z'] = element[58:64]
				dddd_protein[element[21]][element[23:26]][element[13:16]]['id'] = element[23:26]
			
			#the atom list associated to the residu sub-dictionnary is updated with the intermediate atom list
			dddd_protein[element[21]][element[23:26]]['atomlist'] = l_atoms
			#the residu list associated to the chain sub-dictionnary is updated with the intermediate residu list
			dddd_protein[element[21]]['reslist'] = l_reslist
			#the chain list associated to the dictionnary is updated with the intermediate chain list
			dddd_protein['chains'] = l_chains
		

	IN.close()
	
	return dddd_protein




fichier = sys.argv[1]


print PDB_Parsing(fichier)

