#! /usr/bin/env python
# -*- coding: utf-8 -*-

def getvoisins(a,b):
    """but : retourner les coord des cases voisines de la case de coord a,b
    input : 2 entiers correspondant aux coord phi et theta de la case dont on veut les voisins
    output : liste de liste : liste de coord [phi,theta] de chacune des cases voisines
    """ 
    return [[a+5,b+5],[a-5,b-5],[a+5,b],[a-5,b],[a,b+5],[a,b-5],[a+5,b-5],[a-5,b+5]]

#créer des classes
def createClass(listE, bestscore, nbcl) :
    """but : un dictionnaire permettant de classer les énergies en nombre de classes que les utilisateurs souhaitaient
        input : liste d'energie ordonnees de facon croissante, meilleur score, nb de classes
        output : dico 
    """
    classe={}
    compteur=0
    enerseuil=0
    while len(classe) != nbcl:
        compteur=compteur+1
        classe[compteur]=[]
        enerseuil=(nbcl-compteur)*(bestscore/nbcl)
        for energie in listE:
            if energie <= enerseuil:
                classe[compteur].append(energie)
        #print "classe[compteur]:",classe
    dico_etoclass={}# dictionnaire d'énergie à la classe
    for key in classe:
        for elem in classe[key]:
            if not elem in dico_etoclass:
                dico_etoclass[elem]=key
            
    return dico_etoclass
