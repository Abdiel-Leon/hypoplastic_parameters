#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 17:04:45 2023

@author: abdiel
"""
import numpy
import matplotlib.pyplot as plt

class Hypoplastic_Factors:
    
    def __init__(self, *args):
        self.phi_c = args[0]
        self.hs = args[1]
        self.ncoeff = args[2]
        self.ei0 = args[3]
        self.ec0 = args[4]
        self.ed0 = args[5]
        self.beta = args[6]
        self.alpha = args[7]
        
    
    def Compute_fb(self):
        n = 1000
        fb = numpy.zeros(n)
        p = numpy.zeros(n)
        
        # pressure loop
        for i in range(n): 
            p[i] = 100*i
            a =numpy.sqrt(3.0) * (3.- numpy.sin(self.phi_c))/(2.0*numpy.sqrt(2.0)*numpy.sin(self.phi_c))       
            fb[i] =self.hs/self.ncoeff *  pow(self.ei0/self.ec0, self.beta) * pow(3.*p[i]/self.hs,(1.0-self.ncoeff))\
            *pow( ( 3.0 + a*a - numpy.sqrt(3.0)*a * pow((self.ei0 - self.ed0)/(self.ec0 - self.ed0), self.alpha)),-1) 
            
        return p,fb
    
    
def main():
 #args = (30*3.1416/180,1e9,0.40,1,0.85,0.54,1,0.12)   
 args = (20*3.1416/180,1e9,0.40,1,0.85,0.54,1,0.12)   
 
 pressure = Hypoplastic_Factors(*args)  
 pf, fbf = pressure.Compute_fb()   
 plt.plot(pf, fbf)
 
 
main() 