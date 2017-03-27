#!/usr/bin/env python
#-*- coding : utf8 -*-

"""
Contact: 
Date: 06/03/2017
Description: 
- A function parsing a PDB file into a dictionary
"""
from math import sqrt

def ParsingPDB (pdbFile):

    print("Ce programme vous permet de transformer un fichier .pdb au format ATOM en un dictionnaire manipulable par Python.")

    infile = open(pdbFile, "r")
    lines = infile.readlines() 						# cree une liste dont chaque element est une ligne du fichier

    molecule = {}									# dictionnaire le plus externe
    chainList = []
    rList = []


    for line in lines:
        if line[:4:] == 'ATOM':						# si la ligne commence par 'ATOM'
        
            chaine = line[21]
            if chaine not in chainList:
                chainList.append(chaine)
                molecule[chaine] = {}
                resList=[]
        
            res = line[23:26]
            if res not in resList :
                resList.append(res)
                molecule[chaine][res] = {}
                atomList=[]
                rList=[]
            
            atom = line[13:16]
            if atom not in atomList:
                atomList.append(atom)
                molecule[chaine][res][atom] = {}
            
            if line[17:20] not in rList:
                rList.append(line[17:20])
    
            molecule[chaine][res][atom]['x'] = line[31:38]
            molecule[chaine][res][atom]['y'] = line[39:46]
            molecule[chaine][res][atom]['z'] = line[47:54]
            molecule[chaine][res][atom]['id'] = line[7:11]
        
            molecule[chaine][res]['resname'] = rList
            molecule[chaine]['reslist'] = resList
            molecule[chaine][res]['atomlist'] = atomList
        
    infile.close()
    return(molecule)

    
def Distance_res(dico1,dico2,mode):
    #Calcul la distance entre 2 résidus
    #Dico1 correspondant au premier residus
    #Dico2 correspondant au deixième residus
    #Mode : les differents modes de calculs de distance
    
    #Mode de calcul 1 : distance la plus courte
    if mode==1:
        min=200 #Distance la plus courte
        for atom1 in dico1["atomlist"]: #Calcul de la distance entre les 2 residus
            coor1=[dico1[atom1]["x"],dico1[atom1]["y"],dico1[atom1]["z"]]
            for atom2 in dico2["atomlist"]:
                coor2=[dico2[atom2]["x"],dico2[atom2]["y"],dico2[atom2]["z"]]
                distance=Distance(coor1[0],coor1[1],coor1[2],coor2[0],coor2[1],coor[2])
                print(distance)
                if distance<min:
                    min=distance
    
    
    #Mode de calcul 2 : entre les centres de masse
    
    elif mode==2:
        #Calcul de la distance du centre vers les deux residus
        #Moyenne de chaque distance
        distance=centre_masse(coor1[0],coor1[1],coor1[2],coor2[0],coor2[1],coor[2],0,0,0)
        
    return(distance)

def centre_masse(x1,y1,z1,x2,y2,z2,x_centre,y_centre,z_centre):
    #Permet de calculer la distance entre le centre et deux residus fournis en argument
    #x_centre,y_centre et z_centre sont les coordonnées du centre
    d1=list()
    d2=list()
    d1=Distance(x1,x_centre,y1,y_centre,z1,z_centre)
    d2=Distance(x2,x_centre,y2,y_centre,z2,z_centre)
    
    moyenne=(d1+d2)/2
    return(moyenne)

def Distance(x1,y1,z1,x2,y2,z2):#Calcule la distance entre deux points
    return(sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)


def Ecrituretxt():#ecriture dans fichier txt
    fichier_txt=open("fichier.txt","w")
    
    liste_residus
    nbr_residus
    residus_polaires
    residus_hydrophobes
    
    fichier_txt.write(liste_residus)
    fichier_txt.write(nbr_resiuds)
    fichier_txt.write(residus_polaires)
    fichier_txt.write(residus_hydrophobes)
    fichier_txt.close()

def EcritureDansPdb():#ecriture dans le fichier pdb dans lequel les residus de linterface sont etiquetes
    fichier_pdb=open("fichier.pdb","w")
     fichier_pdb.close()



if __name__ == '__main__':
    
    res1=ParsingPDB("/Users/mathildebertrand/Desktop/M1-S2/python/tp3/Documents/arginine.pdb")
    #Dico1 correspondant au premier residus
    
    res2=ParsingPDB("/Users/mathildebertrand/Desktop/M1-S2/python/tp3/Documents/arginine.pdb")
    #Dico2 correspondant au deuxième residus

    print(Distance_res(res1,res2,1))
    print(Distance_res(res1,res2,2))
 
   
    




