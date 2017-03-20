#! /usr/bin/env python
# -*- coding: utf-8 -*-
import math



def dico(lines):
	dico_prot = dict()

	for line in lines:
		if line[0:4]=="ATOM":
			chain = line[21]
	
			if chain not in dico_prot.keys():
				dico_prot[chain] = dict()
				dico_prot[chain]["reslist"]=[]
	
			res = line[22:26]
			if res not in dico_prot[chain].keys():
				dico_prot[chain][res] = dict()
				dico_prot[chain]["reslist"].append(res)
		
			residu = line[17:20]
			if residu not in dico_prot[chain][res].keys():
				dico_prot[chain][res][residu] = dict()
				dico_prot[chain][res][residu]["atomlist"]=[]
	
			atom = line[13:16]
			if atom not in dico_prot[chain][res][residu].keys():
				dico_prot[chain][res][residu][atom] = dict()
				dico_prot[chain][res][residu]["atomlist"].append(atom)
		
			dico_prot[chain][res][residu][atom]["x"]=line[31:38]
			dico_prot[chain][res][residu][atom]["y"]=line[39:46]
			dico_prot[chain][res][residu][atom]["z"]=line[47:54]

	return dico_prot

	
	
	
	




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
	
