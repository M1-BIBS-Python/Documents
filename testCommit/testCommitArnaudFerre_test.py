#!/ usr/bin/env python
#coding: utf-8

"""
Author: Arnaud Ferré
Contact: arnaud.ferre@u-psud.fr
Date: 09/02/2017
Description: 
- A clean example of a Python script architecture for test
- Used to test commit command
Licence: DSSL
"""


# External importations:
import math
# Importations of fonctions to test:
from testCommitArnaudFerre_fonctions.py import *


# tests:
if __name__ == '__main__':
    
	# Tests for hypotenuse fonction:
	assert(hypotenuse(0,0) == 0)
    assert(hypotenuse(1,1) == math.sqrt(2))
	assert(hypotenuse(2,2) == math.sqrt(8))
	assert(hypotenuse(3,4) == 5)
	