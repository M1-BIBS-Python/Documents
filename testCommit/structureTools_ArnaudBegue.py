#!/usr/bin/env python
#-*- coding : utf8 -*-

import sys
import math

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
				"X":(float)(myline[32:39].strip()),
				"Y":(float)(myline[40:47].strip()),
				"Z":(float)(myline[48:55].strip()),
				"id":myline[9:12].strip()
			}

			chain[mychain][myresidu][nomatome] = info	
			
			f.close()
	
		
	return chain ;

def affichage(chain):
	## Affichage

	affichage= str(chain).replace("},","},\n")
	affichage= str(affichage).replace("},","},\n")		
	affichage=str(affichage).replace("},","},\n")		

	print affichage	


def cbm(chain):

	for i in chain.keys():
		for j in chain[i].keys():
			cbm={'X':0 , 'Y':0, 'Z':0}
			div=0
			
			for k in chain[i][j].keys():
				
				cbm['X']+=chain[i][j][k]['X']
				cbm['Y']+=chain[i][j][k]['Y']
				cbm['Z']+=chain[i][j][k]['Z']
				div=div+1
			
			cbm['X']/=div	
			cbm["Y"]/=div
			cbm["Z"]/=div
			chain[i][j]["cbm"]=cbm
			
	return chain
	
def distance(chain):
	distance=dict()
	for i in chain.keys():
		for j in chain [i].keys():
			for k in chain [i].keys():
				if(j != k and str(i+":"+k+"&"+j) not in distance) :
					# correspond aux coordonnee du cbm d'un residu de la chaine	aux carres
					x1=chain[i][j]['cbm']['X']
					y1=chain[i][j]['cbm']['Y']
					z1=chain[i][j]['cbm']['Z']
					# correspond aux coordonnee du cbm d'un residu de la chaine	aux carres
					x2=chain[i][k]['cbm']['X']
					y2=chain[i][k]['cbm']['Y']
					z2=chain[i][k]['cbm']['Z']
					
					dist=math.sqrt(pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2))
					distance[str(i+":"+j+"&"+k)]=dist
	
	return distance
	


		
		
	
#~ def graphe(chain):

	


	
if __name__=='__main__':
	a=lire_pdb(sys.argv[1])
	
	b=cbm(a)
	
	c=distance(b)
	affichage(c)
	
