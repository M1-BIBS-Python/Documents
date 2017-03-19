#!/usr/bin/env python

#Author: GHARBI Houssem

#Contact: houssem-eddine.gharbi@u-psud.fr
#Date: 20/02/2017
#Description: Parceur d'un fichier PDB


import math, string


def parsePDBMultiChains(infile) :

    # lecture du fichier PDB 
    f = open(infile, "r")
    lines = f.readlines()
    f.close()


    # var init
    
    dPDB = {} #initiation d'un dic
   
    dPDB["listChains"] = []  #dPDB est un dic qui a comme cle listChain qui prend une liste comme valeur


    # parcoure le PDB   
    for line in lines :
        if line[0:4] == "ATOM" :
			
            chain = line[21] #chain contient la lettre de la chaine Ex:A
            
            if not chain in dPDB["listChains"] :
				
                dPDB["listChains"].append(chain) #rajouter A a la liste
                
                dPDB[chain] = {} # creer une cle chain dans dPDB (Ex:A) qui prend un dic comme valeur
                
                dPDB[chain]["reslist"] = [] # ce dic a une cle reslist qui prend une liste comme valeur
                
            curres = "%s"%(line[22:26]).strip() #curres contient le num du residu
            
            if not curres in dPDB[chain]["reslist"] :
				
                dPDB[chain]["reslist"].append(curres) #ajouter curres a la liste
                
                dPDB[chain][curres] = {} #on creer un dic dans dPDB[chain] pour la cle curres Ex: -3
                
                dPDB[chain][curres]["resname"] = string.strip(line[17:20]) #ce dernier dic a pour cle resname
                
                dPDB[chain][curres]["atomlist"] = [] #ce dernier dic a pour cle atomelist qui prend pour valeur un liste
                
            atomtype = string.strip(line[12:16])
            
            dPDB[chain][curres]["atomlist"].append(atomtype) #ajouter l'atome a la liste
            
            dPDB[chain][curres][atomtype] = {} #creer un dic dans dPDB[chain][curres] pour la cle atomtype Ex:C
            
 
		# ce dernier dic a pour cle x, y, z, et id 
            dPDB[chain][curres][atomtype]["x"] = float(line[30:38]) 
            dPDB[chain][curres][atomtype]["y"] = float(line[38:46])
            dPDB[chain][curres][atomtype]["z"] = float(line[46:54])
            dPDB[chain][curres][atomtype]["id"] = line[6:11].strip()

    return dPDB
    
    
#test:
res=parsePDBMultiChains("arginine.pdb")
print res
 
   
