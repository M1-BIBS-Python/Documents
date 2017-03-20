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
				
				info={'ID': (float)(l[9:12].strip()),
						'x'	: (float)(l[32:39].strip()),
						'y' : (float)(l[40:47].strip()),
						'z' : (float)(l[48:55].strip())
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


def centreMasse(b):
	
	for i in b.keys():# On parcours la chaine
		for j in b[i].keys():
			
			cdm = {'x':0,'y':0,'z':0}
			c=0
			for k in b[i][j].keys():
				
				cdm['x']+=b[i][j][k]['x']
				cdm['y']+=b[i][j][k]['y']
				cdm['z']+=b[i][j][k]['z']
				c+=1
				
			cdm['x']/=c
			cdm['y']/=c
			cdm['z']/=c
			
			b[i][j]["cdm"]=cdm
			
	return b
	
	
def dessiner(c):
	
	
	

if __name__ == '__main__':
	monDico = dict()
	monDico = lirePDB(sys.argv[1])
	centreMasse(monDico);

