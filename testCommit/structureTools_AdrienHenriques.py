#!/usr/bin/env python
#-*- coding:utf8 -*-

""" Programme du TP5 """
from math import sqrt
#### PARCER PDB #######
from parcing import ParcerPDB

# Utilisation du parcer (main)


	
	#### Fait la moyenne de chaques x , y et z pour tout les atomes a chaques positions
def PointPosition(resultats):
	position=dict()
	for chain in resultats.keys():
		for residu in resultats[chain].keys():
			position[residu]={}
			compteur=0
			sommeX=0
			sommeY=0
			sommeZ=0
			for atome in resultats[chain][residu].keys():
				sommeX =sommeX+float(resultats[chain][residu][atome]["x"])
				sommeY =sommeY+float(resultats[chain][residu][atome]["y"])
				sommeZ =sommeZ+float(resultats[chain][residu][atome]["z"])
				compteur+=1
			position[residu]["x"]=sommeX/compteur
			position[residu]["y"]=sommeY/compteur
			position[residu]["z"]=sommeZ/compteur
			
	return(position)
			
	
			
	### Calcule la distance entre 2 positions donnée	
#def CalculeDistance(position):
#	distance=dict()
#	for residu in position.keys():
#		for resi in position.keys():
			#if (residu not in distance.keys()) and (resi not in distance[residu].keys()):
#				distance[residu][resi]={}
#				distance[residu][resi]=sqrt((position[residu]["x"]-position[resi]["x"])^2-(position[residu]["y"]-position[resi]["y"])^2-(position[residu]["z"]-position[resi]["z"])^2)
#		print position
	
#	return()
	



fichier = raw_input ("Saississez le nom de votre fichier avec le format (ex: arginine.pdb):")
ParcerPDB(fichier)

dddd_result = ParcerPDB(fichier)
position=PointPosition(dddd_result)
#CalculeDistance(position)


for residu in position:
	for resi in position:
		distance=sqrt((position[residu]["x"]-position[resi]["x"])^2-(position[residu]["y"]-position[resi]["y"])^2-(position[residu]["z"]-position[resi]["z"])^2)
		print distance







