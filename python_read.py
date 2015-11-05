#!/usr/bin/env python

import astropy
import math
#astropy.test()

from math import *

#print dir(math)

import re

#program to make plots in python

print ''
print ''
print ''
print '**************************'
print '*                        *'
print '*  Starting Python code  *'
print '*                        *'
print '**************************'
print ''
print ''

print sqrt(25)


"""def shut_down(s):
    if s=="yes":
        return "Shutting down"
    elif s=="no":
        return "Shutdown aborted"
    else:
        return "Sorry"
    """

#letters = ['a', 'b', 'c']
#letters.append('d')

#print "Hens", 25 + 30 / 6
#print 10+5

#a = "string"
#b = 10
#print "S=%s Number=%d" % (a, b)
#print "Number=%f" % b
#type

Datadir = "./data/"

"""with open(Datadir+"liwhite_2009.txt", "r") as file:
    content = file.readlines()

file = open(Datadir+"liwhite_2009.txt")
for line in file:
       cmass,mass2,mass3,cphi,cerror  = line.strip().split()
       print mass3
"""

#10**2 - 10^2

# variable.lower() - call method lower of variable
#srt() turn into string


with open(Datadir+"liwhite_2009.txt") as file:
    for line in file:  #Line is a string       
        numbers_str = line.split() #split the string on whitespace, return a list of numbers (as strings)       
        if len(line.split())>0: 
            values = [float(x) for x in numbers_str]  #convert numbers to floats
      

cmass,mass2,mass3,cphi,cerror=values

print mass3
