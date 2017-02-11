#!/ usr/bin/env python
#coding: utf-8

"""
Author: Arnaud Ferré
Contact: arnaud.ferre@u-psud.fr
Date: 09/02/2017
Description: 
- A clean example of a Python script architecture
- Used to test commit command
Licence: DSSL
"""

# Importations
import math


# Function definition
def hypotenuse(a, b):
    """
    Calculate the size of the hypotenuse of a triangle with the 
    2 other sides of size a and b. 
    """
    
    result = math.sqrt(a**2 + b**2)
    
    return result


# Is executed only if it is this script which is executed by the interpretor 
# (not if it is contained in an imported file) 
# Useful to do unit test in every script!
if __name__ == '__main__':
    
    print(hypotenuse(1,1))

	