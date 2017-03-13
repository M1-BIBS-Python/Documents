#!/usr/bin/env python
#-*- coding : utf3 -*-
import sys

def lirePDB(a):


	#### Dictionnaires ####
	residu = dict()
	chaine = dict()
	atome = dict()
	####	Fin dico   ####
	conf=""

	with open(a, "r") as fichier:
		fic = fichier.readlines()
		for l in fic:
			
			if (conf=="" and l[0:4]=="ATOM"):
				conf=l[16]
				
			if (l[0:4]=="ATOM" and l[16]==conf): #Scan des chaines d'interet
				##### Recuperation des donnees ###
				res=l[24:31].strip()
				nomAtome = l[13:16].strip()
				chName=l[20:23].strip()
				
				info={'ID': l[9:12].strip(),
						'x'	: l[32:39].strip(),
						'y' : l[40:47].strip(),
						'z' : l[48:55].strip()
					}
				######## Fin recuperation ########
				
				 # Chaine non repertoriee donc on l'ajoute
				if (chName not in chaine.keys()):
					chaine[chName] = {}
				
				## Residu non encore repertorie
				if(res not in chaine[chName].keys()):
					chaine[chName][res]={}
				
				chaine[chName][res][nomAtome]=info

	return chaine


print lirePDB(sys.argv[1])


