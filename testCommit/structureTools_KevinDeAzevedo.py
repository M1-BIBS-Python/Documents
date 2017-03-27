#!/usr/bin/env python
#-*- coding : utf8 -*-
#Auteur : Kevin De Azevedo
#Date : 13/03/2017

import string
import sys
import math

def parserPDB(infile):
	
	fichier = open(infile,"r")
	lines = fichier.readlines() #Chaque element est une ligne du fichier

	d_proteine= {}
	d_proteine["nchaine"]=[]

	for i in lines:
		if i[0:4]=="ATOM": #On verifie qu'on est bien sur la ligne d'un atome
			
			nchaine=i[21] #On enregistre le numero de chaine
			#On verifie que la chaine n'est pas deja enregistree dans le dictionnaire pour ne pas ecraser l'information
			if not nchaine in d_proteine["nchaine"]:
				d_proteine["nchaine"].append(nchaine)
				d_proteine[nchaine]={}
				d_proteine[nchaine]["position"]=[]
			
			#Pour chaque chaine, on enregistre les positions
			position=i[22:26].strip()
			if not position in d_proteine[nchaine]["position"]:
				d_proteine[nchaine]["position"].append(position)
				d_proteine[nchaine][position]={}
				d_proteine[nchaine][position]["residu"]=i[17:20].strip()
				d_proteine[nchaine][position]["atome"]=[]
			
			#Pour chaque position on enregistre les atomes
			atome=i[12:16].strip()
			if not atome in d_proteine[nchaine][position]["atome"]:
				d_proteine[nchaine][position]["atome"].append(atome)
				d_proteine[nchaine][position][atome]={}
			
			#Pour chaque atome, on enregistre ses coordonnees et son ID
			d_proteine[nchaine][position][atome]["x"]=i[30:38]
			d_proteine[nchaine][position][atome]["y"]=i[38:46]
			d_proteine[nchaine][position][atome]["z"]=i[46:54]
			d_proteine[nchaine][position][atome]["ID"]=i[6:11]
	fichier.close()	
	return d_proteine
	

#Autres fonctions

def CalculeDistance(dico):
	
	distance=[]
	for nchaine in dico["nchaine"]:
			for position1 in dico[nchaine]["position"]:
				for position2 in dico[nchaine]["position"]:
					if position1 != position2:
						for atome in dico[nchaine][position2]["atome"]:
							if atome=="CA":
								x1=float(dico[nchaine][position1][atome]["x"])
								x2=float(dico[nchaine][position2][atome]["x"])
								y1=float(dico[nchaine][position1][atome]["y"])
								y2=float(dico[nchaine][position2][atome]["y"])
								z1=float(dico[nchaine][position1][atome]["z"])
								z2=float(dico[nchaine][position2][atome]["z"])
								d=math.sqrt(math.pow((x1-x2),2)+math.pow((y1-y2),2)+math.pow((z1-z2),2))
								distance.append(d)
	return distance
	
def Interface(infile):
	fichier = open(infile,"r")
		

#Test

import json
print ("\n Identifiants des chaines de la proteine 1EJH :")
print (parserPDB("/home/kazevedo/Documents/M1BIBS/S2/Python/GitRepo/Documents/Data/1EJH.pdb")["nchaine"])
print ("\n Identifiants des residus de la chaine A de la proteine 1EJH ")
print (parserPDB("/home/kazevedo/Documents/M1BIBS/S2/Python/GitRepo/Documents/Data/1EJH.pdb")["A"]["position"])
dicoProt = parserPDB("/home/kazevedo/Documents/M1BIBS/S2/Python/GitRepo/Documents/Data/1EJH.pdb")
#print dicoProt
print (CarteContact(dicoProt))

