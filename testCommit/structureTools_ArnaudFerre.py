#!/ usr/bin/env python
#coding: utf-8

"""
Author: Arnaud Ferré
Contact: arnaud.ferre@u-psud.fr
Date: 20/02/2017
Description: Script containing many useful functions for parsing and processing PDF files (i.e. 3D structure of proteins)
Licence: DSSL (http://dssl.flyounet.net/)
"""


import string 
import sys # Pour accéder à exit()


def parsePDBMultiChains(infile) :
    """
        Cette fonction permet de charger un fichier PDB (Protein Data Bank) au format ATOM.
        Puis de parser son contenu (structure 3D d'une molécule) pour le stocker dans une variable Python.
        
        Paramètre(s) :
            - infile : emplacement du fichier à charger et parser
        
        Valeur renvoyée :
            - dddd_PDB : dictionnaire contenant 
        
        Pour plus d'informations sur les fichiers PDB et en particulier le format ATOM, veuillez vous reporter à :
        http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM 
    """

    ### On vérifie que l'ouverture du fichier se passe correctement :
    try:
        f = open(infile, "r")
        lines = f.readlines()
        f.close()
    except:
        print("Le fichier n'a pu être chargé correctement. Vérifiez que le fichier existe bien et relancez votre programme.")
        sys.exit(0) ### Stoppe simplement l'exécution du programme.

        
    dddd_PDB = {}
    dddd_PDB["chains"] = []
    
    cmptAltLoc = False
    for line in lines :
        if line[0:4] == "ATOM" :
            
            if cmptAltLoc == False:
                altLoc = line[17]
                cmptAltLoc = True
            
            if line[17] == altLoc:
            
                chain = line[21]

                if not chain in dddd_PDB["chains"] :
                    dddd_PDB["chains"].append(chain)
                    dddd_PDB[chain] = {}
                    dddd_PDB[chain]["reslist"] = []

                curres = line[22:26].strip() 

                if not curres in dddd_PDB[chain]["reslist"] :
                    dddd_PDB[chain]["reslist"].append(curres)
                    dddd_PDB[chain][curres] = {}
                    dddd_PDB[chain][curres]["resname"] = line[17:20].strip()
                    dddd_PDB[chain][curres]["atomlist"] = []

                atomtype = line[12:16].strip()
                dddd_PDB[chain][curres]["atomlist"].append(atomtype)
                dddd_PDB[chain][curres][atomtype] = {}

                # On pourrait aussi vérifier que les chaînes de caractères continnent bien des float pour éviter une erreur.
                dddd_PDB[chain][curres][atomtype]["x"] = float(line[30:38])
                dddd_PDB[chain][curres][atomtype]["y"] = float(line[38:46])
                dddd_PDB[chain][curres][atomtype]["z"] = float(line[46:54])
                
                dddd_PDB[chain][curres][atomtype]["id"] = line[6:11].strip()

    return dddd_PDB


###
# Ici, vous écrirez vos prochaines fonctions
###


### Ici, vous pourrez testez vos fonctions :
if __name__ == "__main__":
    
    # Pour afficher une structure de façon un peu plus esthétique :
    import json
    print("Données pur l'arginine : \n"+json.dumps(parsePDBMultiChains("arginine.pdb"), indent = 4))
    
    print("\nIdentifiants des chaînes de la protéine 1EJH :")
    print(parsePDBMultiChains("1EJH.pdb")["chains"])
    print("\nIdentifiants des résidus de la chaîne A de la protéine 1EJH :")
    print(parsePDBMultiChains("1EJH.pdb")["A"]["reslist"])
    
	