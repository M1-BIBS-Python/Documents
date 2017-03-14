#!/usr/bin/env python
#-*- coding : utf8 -*-
#Auteur : Kevin De Azevedo
#Date : 13/03/2017

import string
import sys

def parserPDB(infile):
	
	fichier = open("infile","r")
	lines = fichier.readlines() #Chaque element est une ligne du fichier
	fichier.close()

	d_proteine= dict()

	for i in lines:
		if i[0:4]=="ATOM": #On verifie qu'on est bien sur la ligne d'un atome
			
			nchaine=i[21:22] #On enregistre le numero de chaine
			#On verifie que la chaine n'est pas deja enregistree dans le dictionnaire pour ne pas ecraser l'information
			if d_proteine.has_key(nchaine)==False:
				d_proteine[nchaine]={}
			
			#Pour chaque chaine, on enregistre les positions
			position=i[23:26]
			if d_proteine[nchaine].has_key(position)==False:
				d_proteine[nchaine][position]={}
			
			#Pour chaque position on enregistre les atomes
			atome=i[13:16]
			if d_proteine[nchaine][position].has_key(atome)==False:
				d_proteine[nchaine][position][atome]={}
			
			#Pour chaque atome, on enregistre ses coordonnees et son ID
			d_proteine[nchaine][position][atome]["x"]=i[32:38]
			d_proteine[nchaine][position][atome]["y"]=i[40:46]
			d_proteine[nchaine][position][atome]["z"]=i[48:54]
			d_proteine[nchaine][position][atome]["ID"]=i[9:11]
			
	return d_proteine

#Autres fonctions


#Test

if __name__ == " __main__ ":
	import json
	print ("\n Identifiants des chaines de la proteine 1 EJH :")
	print (parserPDB("1EJH.pdb")[nchaine])
	print ("\n Identifiants des residus de la chaine A de la proteine 1EJH ")
	print (parserPDB("1EJH.pdb")["A"][position])
