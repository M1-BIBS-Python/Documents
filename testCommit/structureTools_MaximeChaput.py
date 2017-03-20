#!/usr/bin/env python
#-*- coding : utf8 -*-
#By: maxime chaput droit:  cc by
#lecture d'un fichier pdb

#myfile= open("arginine.pdb" , "r")
import pylab as pl
import numpy as np

def parser(FILE): 								#lit un fichier pdb et creer un dictionnaire
	myfile= open(FILE, "r")
	lines= myfile.readlines()
	dddd_proteine= dict()
	for line in lines:
		if line[0:4]=="ATOM":
#			dddd_proteine["chains"]=[]
			chain=line[21]									#chain
			if chain not in dddd_proteine.keys():
				dddd_proteine[chain]={}
				dddd_proteine[chain]["reslist"]=[]
#				dddd_proteine["chains"].append(chain)
				location_indicateur=line[16]						#permet d'eviter d'utilise des les location alternative
			resids=line[22:27]
			resids = resids.strip()
			if resids not in dddd_proteine[chain].keys():
				dddd_proteine[chain]["reslist"].append(resids)
				dddd_proteine[chain][resids]={}						#position du a
				dddd_proteine[chain][resids]["atomlist"]=[]
				resname=line[17:19]
				dddd_proteine[chain][resids]["resname"]=resname
			if line[17]==location_indicateur:
				atome=line[13:16]
				dddd_proteine[chain][resids]["atomlist"].append(atome.strip())
				dddd_proteine[chain][resids][atome]={}
				dddd_proteine[chain][resids][atome]["x"]=line[31:38]
				dddd_proteine[chain][resids][atome]["y"]=line[39:46]
				dddd_proteine[chain][resids][atome]["z"]=line[47:54]
				dddd_proteine[chain][resids][atome]["id"]=line[9:11]
	#print dddd_proteine
	return dddd_proteine


def remplitmatrice(d_prot,tab):
	for resa in d_prot["A"]["reslist"]:
		for resb in d_prot["A"]["reslist"]:
			print("resa="+str(resa)+"resb="+str(resb))
			print("d_prot[A][resa]="+str(d_prot["A"][resa]))
			tab[i][j]=dist(d_prot["A"][resa],d_prot["A"][resb])
	return tab

def dist (d_res1, d_res2):
	print("d_res1="+str(d_res1))
	xa=float(d_res1['CA']["x"])
	ya=float(d_res1['CA']["y"])
	za=float(d_res1['CA']["z"])
	xb=float(d_res2['CA']["x"])
	yb=float(d_res2['CA']["y"])
	zb=float(d_res2['CA']["z"])
	return sqrt((xa- xb)^2+(ya-yb)^2+(za-zb)^2)


if __name__ == "__main__":


	d_dico= dict()
	d_dico=parser("../Data/1EJH.pdb" )
	#d_dico=parser("../Data/arginine.pdb")
	#print d_dico
	tailleMat=len(d_dico["A"]["reslist"])
	print tailleMat
	print d_dico.keys
	tab=np.zeros((tailleMat,tailleMat))
	print tab
	tab = remplitmatrice(d_dico,tab)
	print tab
	pl.pcolor(tab)
	pl.colorbar()
	pl.show()
	

	


		
	#compteur=0
	#for i in lines:  
	#	print i,'\n'
	
