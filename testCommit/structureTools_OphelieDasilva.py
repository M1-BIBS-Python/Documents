#!/ usr/bin/env python
#coding: utf-8

"""
Author: Ophélie Da Silva
Contact: 
Date: 06/03/2017
Description: 
- Fichier contenant l'ensemble des fonctions nécessaires à blablabla
"""

# Importations
import sys
from math import sqrt
import numpy as np, math, matplotlib.pylab as plt

# Function definition
def parserPDB(input_f):
	''' 
	Fonction qui transforme un fichier pdb en dictionnaire directement 
	utilisable via python 
	'''
	
	# Verification ouverture fichier
	try:
		IN = open(input_f)
		lines = IN.readlines()
		IN.close()
	except:
		print("Erreur ouverture fichier")
		sys.exit(0) # Arret execution programme
    
	lchaine = list()
	lresidus = list()
	latome = list()
	lres = list()
	d_proteine = dict()

	for element in lines :
		if element[:4:] == "ATOM" :
			chaine = element[21]
			if chaine not in lchaine :
				lchaine.append(chaine)
				d_proteine[chaine] = dict()
				lresidus = []
		
			residus = element[24:26]
			if residus not in lresidus :
				lresidus.append(residus)
				d_proteine[chaine][residus] = dict()
				latome=[]
				lres = []
		
			d_proteine[chaine]["reslist"] = lresidus
		
			atome = element[13:16:]
			if atome not in latome :
				latome.append(atome)
		
			if element[17:20] not in lres : 
				lres.append(element[17:20])

			d_proteine[chaine][residus][atome] = dict()
			d_proteine[chaine][residus][atome]["x"] = float(element[31:38:])
			d_proteine[chaine][residus][atome]["y"] = float(element[39:46:])
			d_proteine[chaine][residus][atome]["z"] = float(element[47:54:])
			d_proteine[chaine][residus][atome]["id"] = element[9:11:]
			d_proteine[chaine]["resname"] = lres
			d_proteine[chaine][residus]["atomelist"] = latome
	
	return d_proteine








# Demonstration des fonctions locales
if __name__ == '__main__':
	d_prot = parserPDB("/home/tp-home001/odasilv/Bureau/python/GitRepo/Documents/Data/5iby.pdb")
#~ print d_prot
	daa_min = list()
	
	for i in d_prot.keys() : # pour chaque chaine A, B,C...
		l_aa = d_prot[i].keys()
		l_aa.remove('reslist')
		l_aa.remove('resname')
		
		lx_moy = list()
		ly_moy = list()
		lz_moy = list()
		
		datomeall = list()
			
		for j in l_aa : # pour chaque aa
			
			l_atome = d_prot[i][j].keys()
			l_atome.remove('atomelist')

			lx = list()
			ly = list()
			lz = list()
			nb_element = 0
						
			for k in l_atome : # pour chaque atome
				lx.append(d_prot[i][j][k]["x"])
				ly.append(d_prot[i][j][k]["y"])
				lz.append(d_prot[i][j][k]["z"])
				
				nb_element += 1
			
			# barycentre des acides aminés
			lx_moy.append(sum(lx) / nb_element)
			ly_moy.append(sum(ly) / nb_element)
			lz_moy.append(sum(lz) / nb_element)
		
		#~ # distance entre les atomes
		#~ datome = list()
		#~ daa = list()

		#~ for i in range(len(lx) - 1) :
			#~ for j in range(i+1,len(lx)) :
				#~ distance = sqrt((lx[i]-lx[j])**2 + (ly[i]-ly[j])**2 + (lz[i]-lz[j])**2)
				#~ datome.append(distance)
		#~ datomeall += datome
		
	#~ # distance entre les residus
	#~ if (len(lx_moy) != 1) :
		#~ for i in range(len(lx_moy) - 1) :
			#~ for j in range(i+1,len(lx_moy)) :
				#~ distance = sqrt((lx_moy[i]-lx_moy[j])**2 + (ly_moy[i]-ly_moy[j])**2 + (lz_moy[i]-lz_moy[j])**2)
				#~ daa.append(distance)

	if (len(lx_moy) != 1) :
		m_distance = np.zeros((len(lx_moy),len(lx_moy)))
		for i in range(len(lx_moy) - 1) :
			for j in range(i+1,len(lx_moy)) :
				d = sqrt((lx_moy[i]-lx_moy[j])**2 + (ly_moy[i]-ly_moy[j])**2 + (lz_moy[i]-lz_moy[j])**2)
				m_distance[i,j] = d
				m_distance[j,i] = d
			m_distance[i][i] = 0
print m_distance

#~ # Ensemble des distances entres les aa
#~ print daa
#~ print len(daa)
#~ print len(lx_moy)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_aspect('equal')
plt.imshow(m_distance, interpolation='nearest', cmap=plt.cm.coolwarm)
plt.colorbar()
plt.show()


