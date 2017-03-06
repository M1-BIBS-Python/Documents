#!/usr/bin/env python
#-*- coding : utf8 -*-

'''
Author: Marine Duhamel
Contact: marine.c.duhamel@gmail.com
Date: 20/02/2017
'''

seq="ATCGTTGCAAGCTATGCCAATCCCATGGGCTTAGGACTACTGATCGACGT"


def Nb_occ(ma_sequence):
	'''Ceci est une fonction qui prend en argument une liste contenant une sequence nucleotidique
	et retourne un dictionnaire qui contient la sequence, et le nombre d'occurence des mots de 2 lettres'''
	d_occ={}
	'''Integration de la sequence au dictionnaire'''
	d_occ["complete"]=ma_sequence
	'''chaque mot de 2 lettres present dans la sequence est identifie et compte le nombre d'occurences'''
	for i in range(len(ma_sequence)-1):
		mot2=seq[i]+seq[i+1]
		if mot2 in d_occ.keys():
			d_occ[mot2] +=1
		else:
			d_occ[mot2] = 1
	return d_occ

print Nb_occ(seq)
