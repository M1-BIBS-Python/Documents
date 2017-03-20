#!/usr/bin/env python
#-*- coding : utf8 -*-
#Auteur : Kevin De Azevedo
#Date : 13/03/2017

import string
import sys

def parserPDB(infile):
	
	fichier = open(infile,"r")
	lines = fichier.readlines() #Chaque element est une ligne du fichier

	d_proteine= {}
	d_proteine["nchaine"]=[]

	for i in lines:
		if i[0:4]=="ATOM": #On verifie qu'on est bien sur la ligne d'un atome
			
			nchaine=i[21:22] #On enregistre le numero de chaine
			#On verifie que la chaine n'est pas deja enregistree dans le dictionnaire pour ne pas ecraser l'information
			if not nchaine in d_proteine["nchaine"]:
				d_proteine["nchaine"].append(nchaine)
				d_proteine[nchaine]={}
				d_proteine[nchaine]["position"]=[]
			
			#Pour chaque chaine, on enregistre les positions
			position=i[23:26].strip()
			if not position in d_proteine[nchaine]["position"]:
				d_proteine[nchaine]["position"].append(position)
				d_proteine[nchaine][position]={}
				d_proteine[nchaine][position]["resname"]=line[17:20].strip()
				d_proteine[nchaine][position]["atome"]=[]
			
			#Pour chaque position on enregistre les atomes
			atome=i[13:16].strip()
			if not atome in d_proteine[nchaine][position]["atome"]:
				d_proteine[nchaine][position]["atome"].append(atome)
				d_proteine[nchaine][position][atome]={}
			
			#Pour chaque atome, on enregistre ses coordonnees et son ID
			d_proteine[nchaine][position][atome]["x"]=i[32:38]
			d_proteine[nchaine][position][atome]["y"]=i[40:46]
			d_proteine[nchaine][position][atome]["z"]=i[48:54]
			d_proteine[nchaine][position][atome]["ID"]=i[9:11]
		
	return d_proteine
	fichier.close()

#Autres fonctions



#Test

parserPDB("1EJH.pdb")

if __name__ == " __main__ ":
	import json
	print ("\n Identifiants des chaines de la proteine 1EJH :")
	print (parserPDB("1EJH.pdb")["nchaine"])
	print ("\n Identifiants des residus de la chaine A de la proteine 1EJH ")
	print (parserPDB("1EJH.pdb")["A"]["position"])

