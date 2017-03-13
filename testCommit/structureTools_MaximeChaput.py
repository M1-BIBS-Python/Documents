#!/usr/bin/env python
#-*- coding : utf8 -*-
#By: maxime chaput droit:  cc by
#lecture d'un fichier pdb

#myfile= open("arginine.pdb" , "r")












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
			if line[16]==location_indicateur:
				atome=line[13:16]
				dddd_proteine[chain][resids]["atomlist"].append(atome.strip())
				dddd_proteine[chain][resids][atome]={}
				dddd_proteine[chain][resids][atome]["x"]=line[31:38]
				dddd_proteine[chain][resids][atome]["y"]=line[39:46]
				dddd_proteine[chain][resids][atome]["z"]=line[47:54]
				dddd_proteine[chain][resids][atome]["id"]=line[9:11]
	print dddd_proteine
	return dddd_proteine





	



if __name__ == "__main__":



	#parser("/home/tp-home001/mchaput/Documents/python/test2.pdb" )
	parser("/home/tp-home001/mchaput/Documents/python/arginine.pdb" )
	#print dddd_proteine




		
	#compteur=0
	#for i in lines:  
	#	print i,'\n'
	
