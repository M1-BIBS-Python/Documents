#! /usr/bin/env python
# -*- coding: utf-8 -*-
import math

def CDM(dico1,dico2):
	dico_mass=dict()
	dico_mass['C']=12
	dico_mass['H']=1
	dico_mass['O']=16
	dico_mass['N']=14
	
	x1=0
	y1=0
	z1=0
	residu1=dico1.keys()
	residu2=dico2.keys()
	for keys1 in dico1[residu1[0]]["atomlist"]:
		x1=x1+float(dico1[residu1[0]][keys1]["x"])
		y1=y1+float(dico1[residu1[0]][keys1]["y"])
		z1=z1+float(dico1[residu1[0]][keys1]["z"])
		
	x1=x1/cpt
	y1=y1/cpt
	z1=z1/cpt
	cpt=0
	x2=0
	y2=0
	z2=0
	for keys2 in dico2[residu2[0]]["atomlist"]:
		x2=x2+float(dico2[residu2[0]][keys2]["x"])
		y2=y2+float(dico2[residu2[0]][keys2]["y"])
		z2=z2+float(dico2[residu2[0]][keys2]["z"])
		cpt=cpt+1
	x2=x2/cpt
	y2=y2/cpt
	z2=z2/cpt
	dist=math.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
	return dist
	
	
	


if __name__ = "main":
	
	
	
	#exécutables...
	
