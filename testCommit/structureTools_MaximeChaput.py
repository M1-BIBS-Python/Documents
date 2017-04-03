#!/usr/bin/env python
#-*- coding : utf8 -*-
#By: maxime chaput droit:  cc by
#lecture d'un fichier pdb

#myfile= open("arginine.pdb" , "r")
import pylab as pl
import numpy as np
import math

def parser(FILE): 								#lit un fichier pdb et creer un dictionnaire
	myfile= open(FILE, "r")
	lines= myfile.readlines()
	dddd_proteine= dict()
	dddd_proteine["chains"]=[]
	for line in lines:
		if line[0:4]=="ATOM":
			chain=line[21]									#chain
			if chain not in dddd_proteine.keys():
				dddd_proteine[chain]={}
				dddd_proteine[chain]["reslist"]=[]
				dddd_proteine["chains"].append(chain)
				location_indicateur=line[16]						#permet d'eviter d'utilise des les location alternative
			resids=line[22:27]
			resids = resids.strip()
			if resids not in dddd_proteine[chain].keys():
				dddd_proteine[chain]["reslist"].append(resids)
				dddd_proteine[chain][resids]={}						#position du a
				dddd_proteine[chain][resids]["atomlist"]=[]
				resname=line[17:19]
				dddd_proteine[chain][resids]["resname"]=resname
			if line[16]==location_indicateur:
				atome=line[13:16]
				atome=atome.strip()
				dddd_proteine[chain][resids]["atomlist"].append(atome)
				dddd_proteine[chain][resids][atome]={}
				dddd_proteine[chain][resids][atome]["x"]=line[31:38]
				dddd_proteine[chain][resids][atome]["y"]=line[39:46]
				dddd_proteine[chain][resids][atome]["z"]=line[47:54]
				dddd_proteine[chain][resids][atome]["id"]=line[9:11]
				dddd_proteine[chain][resids][atome]["fact_temperature"]=line[61:66]
				dddd_proteine[chain][resids][atome]["temp"]=line[60:65]
				dddd_proteine[chain][resids][atome]["occupation"]=line[60:65]
				dddd_proteine[chain][resids][atome]["atome_signe"]=line[77:78]
	#print dddd_proteine
	return dddd_proteine


def remplitmatrice_chain(d_prot,tab,chains):
	i=0
	for resa in d_prot[chains]["reslist"]:
		j=0
		for resb in d_prot[chains]["reslist"]:
			j=j+1
			tab[i][j]=dist(d_prot[chains][resa],d_prot[chains][resb])
		i=i+1
	return tab

def remplitmatrice_all(d_prot,tab):
	i = (-1)

	for chains in d_prot["chains"]:
		for resa in d_prot[chains]["reslist"]:
			i += 1
			j = 0
			for chains2 in d_prot["chains"]:
				for resb in d_prot[chains2]["reslist"]:
					#~ print "J =",j
					tab[i][j]=dist(d_prot[chains][resa],d_prot[chains2][resb])
					j=j+1
				
	return tab





def dist (d_res1, d_res2):
	#~ print("d_res1="+str(d_res1))
	xa=float(d_res1['CA']["x"])
	ya=float(d_res1['CA']["y"])
	za=float(d_res1['CA']["z"])
	xb=float(d_res2['CA']["x"])
	yb=float(d_res2['CA']["y"])
	zb=float(d_res2['CA']["z"])
	#~ print ("xa"+str(xa))
	#~ print (xb)
	#~ print (xa-xb)
	return math.sqrt((xa- xb)**2+(ya-yb)**2+(za-zb)**2)

def ecrire_pdb (nom,d_final):
	nom+=".txt"
	f= open(nom,"w")
	for i in range(1,9):
		for j in range(9):
			f.write(" ")
		f.write(str(i))
	f.write("\n")
	for i in range(8):
		for j in range(1,11):
			f.write(str(j%10))
	for res in d_final["A"]["reslist"]:
		for atome_bis in d_final["A"][res]["atomlist"]:
			f.write("\nATOM")
			for i in range(6-len(res)):
				f.write(" ")
			f.write(d_final["A"][res][atome_bis]["id"])
			f.write("  ")
			f.write(atome_bis)
			
			f.write(res)
		
	f.close()







if __name__ == "__main__":


	d_dico= dict()
	d_dico=parser("../Data/1EJH.pdb" )
	#~ d_dico=parser("../Data/arginine.pdb")
	#~ print d_dico
	tailleMat=0
	for chains in d_dico["chains"]:
		tailleMat=tailleMat+len(d_dico[chains]["reslist"])
	print tailleMat
	#~ print d_dico.keys
	tab=np.zeros((tailleMat,tailleMat))
	#~ print tab
	#~ tab = remplitmatrice_all(d_dico,tab)
	tab = remplitmatrice_chain(d_dico,tab,"A")
	print tab
	pl.pcolor(tab)
	pl.colorbar()
	pl.show()
	#~ ecrire_pdb("testo",d_dico)


		
	#compteur=0
	#for i in lines:  
	#	print i,'\n'
	
