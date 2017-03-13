#!/usr/bin/env python
#-*- coding : utf8 -*-


"""
Author: Hugo ARBES
Contact: hugo.arbes@u-psud.fr
Date: 13/03/2017
Description: Script containing many useful functions for parsing and processing PDF files (i.e. 3D structure of proteins)
Licence (que je garde car elle me semble parfaite): DSSL (http://dssl.flyounet.net/)
"""



def Parcer_PDB(entree) :
	"""
	Fonction permettant de creer un dictionnaire definissant une proteine en prenant un fichier PDB en parametre

	Hugo ARBES le 06/03/2017 : Parser de fichier PDB
	"""
	
	try :								# Test d'ouverture du fichier pdb d'entree
		pdb = open(entree,"r")
		lines = pdb.readlines()
		pdb.close()
	
	except :
		print "Le fichier ne s'est pas ouvert correctement"
	
	
	dico = {}
	dico["Chaine"] = []
	
	for line in lines:
		
		if(lines[i][0:4] == "ATOM"):
			
			if((lines[i][16] == "A") or (lines[i][16] == " ")):
				
				
				chain = lines[i][21:22]
				
				if not Chain in dico["Chaine"]:
					dico["Chaine"].append(Chain)
					dico[chain]= {}
					dico[chain]["Residus"] = []
						
				residu = lines[i][24:26]
				
				if not Residu in dico[chain]["Residus"]:
					dico[chain]["Residus"].append(residu)
					dico[chain][residu]= {}
					dico[chain][residu]["Atome"] = []
				
				atom = lines[i][13:15]
				
				if not Atom in dico[chain][residu]["Atome"]:
					dico[chain][residu]["Atome"].append(Atom)
					dico[chain][residu][atom] = {}
				
				dico[chain][residu][atom]["Name"] = lines[i][13:15]
				dico[chain][residu][atom]["x"] = lines[i][31:38]
				dico[chain][residu][atom]["y"] = lines[i][39:46]
				dico[chain][residu][atom]["z"] = lines[i][47:54]
	
	
	return dico


# Autres fonctions de TP






if __name__ == '__main__':
	
	import json
	print("Donnees pour l'arginine : \n"+json.dumps(Parcer_PDB("1EJH.pdb"), indent = 4))



