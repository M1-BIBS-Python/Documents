#! /usr/bin/env python
# -*- coding: utf-8 -*-
import math



def dico(lines):
    dico_prot = dict()

    for line in lines:
        if line[0:4]=="ATOM":
            chain = line[21]
    
            if chain not in dico_prot.keys():
                dico_prot[chain] = dict()
                dico_prot[chain]["reslist"]=[]
    
            res = line[22:26]
            if res not in dico_prot[chain].keys():
                dico_prot[chain][res] = dict()
                dico_prot[chain]["reslist"].append(res)
        
            residu = line[17:20]
            if residu not in dico_prot[chain][res].keys():
                dico_prot[chain][res][residu] = dict()
                dico_prot[chain][res][residu]["atomlist"]=[]
    
            atom = line[13:16]
            if atom not in dico_prot[chain][res][residu].keys():
                dico_prot[chain][res][residu][atom] = dict()
                dico_prot[chain][res][residu]["atomlist"].append(atom)
        
            dico_prot[chain][res][residu][atom]["x"]=line[31:38]
            dico_prot[chain][res][residu][atom]["y"]=line[39:46]
            dico_prot[chain][res][residu][atom]["z"]=line[47:54]

    return dico_prot
    


def calculdist(dico1,dico2):
    dist=[]
    residu1=dico1.keys()
    residu2=dico2.keys()
    for keys1 in dico1[residu1[0]]["atomlist"]:
        for keys2 in dico2[residu2[0]]["atomlist"]:
            dist.append(math.sqrt((float(dico2[residu2[0]][keys2]["x"])-float(dico1[residu1[0]][keys1]["x"]))**2+(float(dico2[residu2[0]][keys2]["y"])-float(dico1[residu1[0]][keys1]["y"]))**2+(float(dico2[residu2[0]][keys2]["z"])-float(dico1[residu1[0]][keys1]["z"]))**2))
    return min(dist)    
    
    

def CDM(dico1,dico2):
    cpt=0
    x1=0
    y1=0
    z1=0
    residu1=dico1.keys()
    residu2=dico2.keys()
    for keys1 in dico1[residu1[0]]["atomlist"]:
        x1=x1+float(dico1[residu1[0]][keys1]["x"])
        y1=y1+float(dico1[residu1[0]][keys1]["y"])
        z1=z1+float(dico1[residu1[0]][keys1]["z"])
        cpt=cpt+1
    x1=x1/cpt
    y1=y1/cpt
    z1=z1/cpt   
    cpt=0
    x2=0
    y2=0
    z2=0
    for keys2 in dico2[residu2[0]]["atomlist"]:
        x2=x2+float(dico2[residu2[0]][keys2]["x"])
        y2=y2+float(dico2[residu2[0]][keys2]["y"])
        z2=z2+float(dico2[residu2[0]][keys2]["z"])
        cpt=cpt+1
    x2=x2/cpt
    y2=y2/cpt
    z2=z2/cpt
    dist=math.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
    return dist
    
    
    


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import numpy as np

    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f", dest="f",help ="Veuillez fournir le chemin vers votre fichier pdb")
    parser.add_option("-d", dest="d",help ="DistanceCourte / CentreDeMasse ?")
    parser.add_option("-s", dest="s",help ="Seuil de distance")
    parser.add_option("-o", dest="o",help ="Le chemin et nom du fichier de sortie")
    parser.print_help()
    (options, args)=parser.parse_args()
    f = options.f
    d = options.d
    s = options.s
    o = options.o

    filin=open (f,"r")
    filout=open(o,"w")
    lines=filin.readlines()

######################################################
#   creer dictionnaire a partir du fichier pdb
######################################################

    dic={}
    dic=dico(lines)
    
######################################################
#   dictionnaire des proprietes des acides amines
######################################################

    dprot=dict()
    dpropriete=dict()
    dprot["charge"]=[['LYS'],['ARG'],['HIS'],['ASP'],['GLU']]
    dprot["polaire"]=[['SER'],['THR'],['CYS'],['MET'],['ASN'],['GLN']]
    dprot["hydrophobe"]=[['ALA'],['ILE'],['LEU'],['VAL'],['PHE'],['TRP'],['TYR'],['GLY'],['PRO']]
    
    for keys in dprot:
        for elements in dprot[keys]:
            dpropriete["%s"%elements]=keys

######################################################
#   matrice de distance
######################################################

    residulist=[]
    for chain in dic:
        for i in dic[chain]["reslist"]:
            residulist.append(dic[chain][i].keys())
          
    matrix=[]
    matrix=np.zeros((len(residulist),len(residulist)))
    res_interface=[]
    hydrophobes_hydrophiles=0
    polaires_polaires=0
    charges_charges=0
    res=[]
    for chain in dic:
        for i in range(len(dic[chain]["reslist"])):
            for j in range(len(dic[chain]["reslist"])):
                if d=="DistanceCourte":
                    matrix[i,j]=calculdist(dic[chain][dic[chain]["reslist"][i]],dic[chain][dic[chain]["reslist"][j]])
                    if matrix[[i,j]]<=s: #stocker les residus de l'interface dans une liste
                        if dic[chain]["reslist"][i] not in res:
                            res_interface.append(dic[chain][dic[chain]["reslist"][i]].keys()) 
                            res.append(dic[chain]["reslist"][i])
                        elif dic[chain]["reslist"][j] not in res:
                            res_interface.append(dic[chain][dic[chain]["reslist"][j]].keys())
                            res.append(dic[chain]["reslist"][j])
                         
                        if dpropriete["%s"%dic[chain][dic[chain]["reslist"][i]].keys()]=="hydrophobe" and dpropriete["%s"%dic[chain][dic[chain]["reslist"][j]].keys()]!="hydrophobe": #on peut considerer que les residus non hydrophobes sont hydrophiles?
                            hydrophobes_hydrophiles=hydrophobes_hydrophiles+1
                        elif dpropriete["%s"%dic[chain][dic[chain]["reslist"][i]].keys()]!="hydrophobe" and dpropriete["%s"%dic[chain][dic[chain]["reslist"][j]].keys()]=="hydrophobe":
                            hydrophobes_hydrophiles=hydrophobes_hydrophiles+1
                        elif dpropriete["%s"%dic[chain][dic[chain]["reslist"][i]].keys()]=="polaire" and dpropriete["%s"%dic[chain][dic[chain]["reslist"][j]].keys()]=="polaire":
                            polaires_polaires=polaires_polaires+1
                        elif dpropriete["%s"%dic[chain][dic[chain]["reslist"][i]].keys()]=="charge" and dpropriete["%s"%dic[chain][dic[chain]["reslist"][j]].keys()]=="charge":
                            charges_charges=charges_charges+1

                elif d=="CentreDeMasse":
                    matrix[i,j]=CDM(dic[chain][dic[chain]["reslist"][i]],dic[chain][dic[chain]["reslist"][j]])
                    if matrix[[i,j]]<=s: #stocker les residus de l'interface dans une liste
                        if dic[chain]["reslist"][i] not in res:
                            res_interface.append(dic[chain][dic[chain]["reslist"][i]].keys()) 
                            res.append(dic[chain]["reslist"][i])
                        elif dic[chain]["reslist"][j] not in res:
                            res_interface.append(dic[chain][dic[chain]["reslist"][j]].keys())
                            res.append(dic[chain]["reslist"][j])
                            
                        if dpropriete["%s"%dic[chain][dic[chain]["reslist"][i]].keys()]=="hydrophobe" and dpropriete["%s"%dic[chain][dic[chain]["reslist"][j]].keys()]!="hydrophobe": #on peut considerer que les residus non hydrophobes sont hydrophiles?
                            hydrophobes_hydrophiles=hydrophobes_hydrophiles+1
                        elif dpropriete["%s"%dic[chain][dic[chain]["reslist"][i]].keys()]!="hydrophobe" and dpropriete["%s"%dic[chain][dic[chain]["reslist"][j]].keys()]=="hydrophobe":
                            hydrophobes_hydrophiles=hydrophobes_hydrophiles+1
                        elif dpropriete["%s"%dic[chain][dic[chain]["reslist"][i]].keys()]=="polaire" and dpropriete["%s"%dic[chain][dic[chain]["reslist"][j]].keys()]=="polaire":
                            polaires_polaires=polaires_polaires+1
                        elif dpropriete["%s"%dic[chain][dic[chain]["reslist"][i]].keys()]=="charge" and dpropriete["%s"%dic[chain][dic[chain]["reslist"][j]].keys()]=="charge":
                            charges_charges=charges_charges+1
    
    plt.pcolor(matrix, cmap='gist_rainbow')
    plt.colorbar()
    plt.show()
    
###################################################################################
#	ecriture dans un fichier 
###################################################################################
    filout.write("residus à l'interface"+"\t"+"reslist de ce résidu"+"\n")
    for i in range(len(res)):
        filout.write("%s"%res_interface[i]+"\t\t\t"+"%s"%res[i]+"\n")
    filout.write("\n\n"+"nombre de contact hydrophobes-hydrophiles:"+"%d"%hydrophobes_hydrophiles)
    filout.write("\n"+"nombre de contact polaires-polaires:"+"%d"%polaires_polaires)
    filout.write("\n"+"nombre de contact charges-charges:"+"%d"%charges_charges)
    
    
