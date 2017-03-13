#!/usr/bin/env python
#-*-coding:utf8 -*-

"""
Licence : Rechercher ce qu'il faut noter avec le CC
Auteur : Julien Roziere 
Date : 06/02/2017
Description : Programme pour débuter python.
"""

#Partie 2

print 'Hello World!' 
print 1
print 2.3
print[4,5,6]

print 'Hello World!' + str(2.0) #str permet de convertir le 2.0 en chaîne de caractères


for i in range (20):
	print 'a',   #affiche a et b 20fois
	print 'b',
	

for i in range (20):
	print 'a',
print 'b',		#affiche 20 fois a et une fois b


#Partie 3

a=15
print a , type (a)

b=2.34
print b , type(b) #la fonction type donne le type de la variable 

c="Hello"
print c , type (c)

la=list() #liste est un type 
lb=[]
lc=[1,2.0,"trois"]
print la,lb,lc,type(la)

#Travail sur des listes

l1=[1,2,3]
print len (l1)
l1.append(5) #rajoute 5 à la fin de la liste
print l1
l1[3]=4 #modifie la position 3 de la liste en 4
print l1


l2=list(l1) #Si j'avais fait l2=l1 alors dans ce cas ce serait les adresses qui seraient pointées et donc dans le cas où je modifie l1 je risque de modifier l2
l3 = range(10) #tous les entiers de 0 à 10 (le 10 est exclu)
l4 = range(2,7) #tous les entiers de 2 à 7 avec le 7 exclu
l5 = range (3,40,7) #comporte les valeurs de 3 à 40 de 7 en 7
l6 = 14+15
print l2,l3,l4,l5,l6


a=[1,2,3]
b=a
b.append(4) #cela indique bien que lorsque je fais b=a c'est en fait l'adresse qui est désigné car si je modifie le b ou le a (b dans ce cas) l'autre est également modifié
print a


c=[13]
d=c
print c,d
c[0]=14
print c,d

#Partie 4

for i in range(10): #affiche les entiers de 0 à 10 (10 exclu)
	print i,
	

liste=["pierre","feuille","ciseaux"]
for w in liste:
	print w,	


#Exercice entrainement

ta="ta"
ti="ti"

for i in range(100):
	print "ta",

for i in range(50):
	print "ti",


for i in range(100):
	print ta+ti,


abc='abcdefghijklmnopqrstuvwxyz'

for i in range(0,len (abc),2):
	print abc[i],


"""for i in range(len(abc),0,-1):
	print abc[i], #à revoir, il faudrait inverser la liste """



#4
#voir la correction

#5

#voir le poly

#6

a = [8,4,12,3,7,2]

for i in range(0,len (a)-1):
	if a[i+1]<=a[i]:
		mini=a[i+1]
print mini	
	
#Fonction pour le tri
n = len(a)

for i in range(n) :                  
	k = i
	for j in range(i+1,n) :
		if a[k]>a[j] :  k = j
		a[k],a[i] = a[i],a[k]
	
print a		

#Exercice pour la redondance

lb=[4,10,4,8,3,19,1,7,7,19,11,3,21,7]

lc=list()

for b in lb:
	if b not in lc:
		lc.append(b)
print lc		

#7

import random
ADN=range(100)

for i in range(100):
	ADN[i]=random.randint(0,3)
	if ADN[i]==0: ADN[i]='A'
	elif ADN[i]==1: ADN[i]='T'
	elif ADN[i]==2: ADN[i]='C'
	elif ADN[i]==3: ADN[i]='G'
print ADN

ADNc=range(100)

for i in range(100):
	if ADN[i]=='A': ADNc[i]='T'
	elif ADN[i]=='T': ADNc[i]='A' 
	elif ADN[i]=='G': ADNc[i]='C'
	elif ADN[i]=='C': ADNc[i]='G'
print ADNc	

nA=0
nT=0
nG=0
nC=0

for i in range(100):
	if ADN[i]=='A': 
		nA=nA+1
	elif ADN[i]=='T': 
		nT=nT+1
	elif ADN[i]=='G': 
		nG=nG+1
	elif ADN[i]=='C': 
		nC=nC+1
print nA,nT,nC,nG

#Partie 5

#obtention d'un 6
de=0
i=0
while de<6:
	de=random.randint(1,6)
	i+=1
print i	

#deux 6
de=0
i=0
nb6=0
while nb6<2:
	de=random.randint(1,6)
	i+=1
	if de==6: nb6+=1
	
print i	

#deux 6 consécutifs
de=0
i=0
consecutif=0
while consecutif<2:
	de=random.randint(1,6)
	i+=1
	if de==6: consecutif+=1
	else: consecutif=0
	
print i	


#Partie 6

#voir énoncé et correction
