#!/ usr/bin/env python
#coding: utf-8

"""
Author: Arnaud Ferré
Contact: arnaud.ferre@u-psud.fr
Date: 09/02/2017
Description: 
- A clean example of a Python script architecture for functions definition
- Used to test commit command
Licence: DSSL
"""

# Importations
import math
# Autres syntaxes/utilisations possibles :
# from math import sqrt
# import math as mat 


# Function definition
def hypotenuse(a, b):
    """
    Calculate the size of the hypotenuse of a triangle with the 
    2 other sides of size a and b. 
    """
    
    result = math.sqrt(a**2 + b**2)
    
    return result


# Test/Demo of the local functions:
if __name__ == '__main__':
    
    print("hyp(1,1)="+str(hypotenuse(1,1)))
	print("hyp(2,2)="+str(hypotenuse(2,2)))
	print("hyp(1,2)="+str(hypotenuse(1,2)))
	