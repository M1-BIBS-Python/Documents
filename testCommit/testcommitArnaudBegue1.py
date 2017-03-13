#!/usr/bin/env python
#-*- coding : utf8 -*-

import sys

	#test

def lire_pdb(fichier):

	f= open(fichier, "r")

	data=f.readlines()

	# creation des dictionnaires
	chain={} 

	# creation flag pour identifier la conformation

	flag=0

	for myline in data:
		

		if (myline[0:4] == "ATOM" and flag == 0):
			
			conf=myline[16]
			
			flag=1
			
			
		if (myline[0:4] == "ATOM" and myline[16] == conf):
			
			#si la chaine n'est pas connue on l'a cree
			
			mychain=myline[20:23].strip() #nom de la chaine			
			if( mychain not in chain.keys()): # si la chaine n'existe pas 			
				chain[mychain]={}
				
			
			myresidu=myline[24:31].strip() # le numero de residu
			if( myresidu not in chain[mychain].keys()): # si le residu n'existe pas 					  			
				chain[mychain][myresidu] = {}
			

			nomatome=myline[13:16].strip() #nom de l'atome	
			info={
				"X":myline[32:39].strip(),
				"y":myline[40:47].strip(),
				"z":myline[48:55].strip(),
				"id":myline[9:12].strip()
			}

			chain[mychain][myresidu][nomatome] = info	
			
			f.close()
	
	## Affichage

	affichage= str(chain).replace("},","},\n")
	affichage= str(affichage).replace("},","},\n")		
	affichage=str(affichage).replace("},","},\n")		

	print affichage	
	
		
	return chain ;

lire_pdb(sys.argv[1])

