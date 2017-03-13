#!/usr/bin/python
#-*- coding : utf8 -*-
# @author : alexandra.benamar@yahoo.fr
# date: 06/03/2017

import os                      # gestion de fichiers et de dossiers
import sys                     # gestion des erreurs et des arguments
import math

def parser_pdb():

    """ Cette fonction a pour but de parser un fichier de type pdb afin
    d'en recuperer les informations sur les differents atomes qui composent
    les acides amines.

    Usage : ./tp3.py <fichier.pdb> """

    # Le nom du fichier est passe en argument
    if len(sys.argv) != 2:
        sys.exit("ERREUR : il faut exactement un argument.")


    # Test d'ouverture du fichier
    try:
        f = open(sys.argv[1],'r')
    except:
        print "Erreur, le fichier n'a pas pu s'ouvrir"
        sys.exit(1)

    lignes  = f.readlines()
    nombre_lignes=len(lignes)


    # Test fichier vide
    if nombre_lignes==0:
        print "Erreur, le fichier est vide"
        sys.exit(1)

    # Initialisation d'un dictionnaire
    dict = {}

    # Test ATOM dans le fichier
    compteur=0

    for i in range(nombre_lignes):

        # Pour toutes les lignes qui commencent par ATOM (celles qui ont des atomes)
        if (lignes[i][0:4] == "ATOM"):

			compteur+=1

			nom_aa=lignes[i][21:22]
			if dict.has_key(nom_aa)==False:
				dict[nom_aa] = {}

        # On enregistre la position de l'atome
			position=lignes[i][23:26]
			if dict[nom_aa].has_key(position)==False:
				dict[nom_aa][position] = {}

        # On enregistre le nom de l'atome
			nom_atome=lignes[i][13:16]
			if dict[nom_aa][position].has_key(nom_atome)==False:
            # On cree un dictionnaire avec toutes les informations qui nous interessent sur l'atome
				dict[nom_aa][position][nom_atome] = {}

				dict[nom_aa][position][nom_atome]['position'] = lignes[i][9:13]
				dict[nom_aa][position][nom_atome]['x'] = lignes[i][31:38]
				dict[nom_aa][position][nom_atome]['y'] = lignes[i][40:48]
				dict[nom_aa][position][nom_atome]['z'] = lignes[i][48:56]

    # Test presence d'ATOM
    if compteur==0:
        print "Le fichier ne contient pas d'ATOM"
        sys.exit(1)

    # Fermeture du fichier
    f.close()

    print compteur,

    # Affichage
    #for i in dict.keys():
        #print dict[i]

    return dict


parser_pdb()












#Cette fonction va nous servir a calculer la distance la plus courte entre les atomes de deux acides amines
def Minimum_distance(dict):
	#on prend en entree un dictionnaire 
	#on souhaite un autre dictionnaire en sortie
	
	#distance = sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2)
	
	#on boucle sur les 
	compteur = 0
	print dict.keys()
	for cle in dict.keys():
		print cle
		#ici on boucle sur le nom de l'acide amine
		compteur = compteur + 1
		for cle2 in dict[cle].keys():
			print cle2
			#maintenant on boucle sur la position, c'est l'acide amine que l'on compare aux autres
			for cle3 in dict[cle].keys():
				print cle3
				
				if cle2 != cle3:
					#on test que nous ne comparions pas l'acide amine avec lui meme
					for cle4 in dict[cle][cle3].keys():
						print cle4
						#on boucle sur les differents atomes des acides amines
						
						x1 = float(dict[cle][cle2][cle4]['x'])
						x2 = float(dict[cle][cle3][cle4]['x'])
						
						y1 = float(dict[cle][cle2][cle4]['y'])
						y2 = float(dict[cle][cle3][cle4]['y'])
						
						z1 = float(dict[cle][cle2][cle4]['z'])
						z2 = float(dict[cle][cle3][cle4]['z'])
						
						distance = math.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
						print distance
						#print compteur
	return distance
				


dict = parser_pdb()

Minimum_distance(dict)

