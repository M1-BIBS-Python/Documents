#!/usr/bin/env python
#-*- coding : utf8 -*-
#Auteur : Kevin De Azevedo
#Date : 13/03/2017

import string
import sys
import math
import numpy

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
			d_proteine[nchaine][position][atome]["x"]=float(i[30:38])
			d_proteine[nchaine][position][atome]["y"]=float(i[38:46])
			d_proteine[nchaine][position][atome]["z"]=float(i[46:54])
			d_proteine[nchaine][position][atome]["ID"]=float(i[6:11])
	fichier.close()	
	return d_proteine
	

#Autres fonctions

def CentreMasse(dico):
	Masse = dict()
	Masse = {"C": float(12.0107), "N": float(14.0067), 
			 "O": float(15.9994), "S": float(32.065),
			 "OH": float(17.0073), "NH": float(15.0146)}
	
	CM={}
	x=0
	y=0
	z=0
	massetotale=0
	
	for atome in dico["atome"]:
		for key in Masse.keys():
			if atome==key:
				x+=Masse[key]*dico[atome]["x"]
				y+=Masse[key]*dico[atome]["y"]
				z+=Masse[key]*dico[atome]["z"]
				massetotale+=Masse[key]
						
	CM["x"]=x/massetotale
	CM["y"]=y/massetotale
	CM["z"]=z/massetotale
	
	return CM
				
	

def CalculeDistance(dico):
	
	distance=[]
	for nchaine in dico["nchaine"]:
			for position1 in dico[nchaine]["position"]:
				p1=CentreMasse(dico[nchaine][position1])
				for position2 in dico[nchaine]["position"]:
					p2=CentreMasse(dico[nchaine][position2])
					if position1 != position2:
						x1=p1["x"]
						x2=p2["x"]
						y1=p1["y"]
						y2=p2["y"]
						z1=p1["z"]
						z2=p2["z"]
						d=math.sqrt(math.pow((x1-x2),2)+math.pow((y1-y2),2)+math.pow((z1-z2),2))
						distance.append(d)
	return distance
		

#Test

import json
print ("\n Identifiants des chaines de la proteine 1EJH :")
print (parserPDB("/home/kazevedo/Documents/M1BIBS/S2/Python/GitRepo/Documents/Data/1EJH.pdb")["nchaine"])
print ("\n Identifiants des residus de la chaine A de la proteine 1EJH ")
print (parserPDB("/home/kazevedo/Documents/M1BIBS/S2/Python/GitRepo/Documents/Data/1EJH.pdb")["A"]["position"])
dicoTest = CalculeDistance(parserPDB("/home/kazevedo/Documents/M1BIBS/S2/Python/GitRepo/Documents/Data/1EJH.pdb"))
print dicoTest

