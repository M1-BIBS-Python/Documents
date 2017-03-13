#!/usr/bin/env python
#-*- coding : utf8 -*-
#By: maxime chaput droit:  cc by
#lecture d'un fichier pdb

#myfile= open("arginine.pdb" , "r")
import pylab

def parser(FILE): 
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


def dist(**d_prot):
	res=len(d_prot["A"]["reslist"])
	print res
	distance=[]
	for i in range(res):
		for j in range(res):
			a=int(d_prot["A"]["reslist"][j])
			b=int(d_prot["A"]["reslist"][i])
			distance.append(abs(a- b))
	return distance

def graph(*resultat):
	a=0
	while a<= len(resultat):
		for i in range(180):
			if resultat[a]<=20:
				print "1",
			elif resultat[a]<=20:
				print "2",
			else :
				print "3",
			a=a+1
		print "\n"

	



if __name__ == "__main__":


	d_dico= dict()
	d_dico=parser("../Data/1EJH.pdb" )
	#d_dico=parser("../Data/arginine.pdb")
	#print d_dico
	tab=[]
	tab=dist(**d_dico)
	#graph(*tab)
	#print dddd_proteine




		
	#compteur=0
	#for i in lines:  
	#	print i,'\n'
	
