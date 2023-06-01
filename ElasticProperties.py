#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 09:20:06 2023

@author: marco
"""

import numpy as np
import sys
import pandas as pd
import os

def get_elastic_tensor(filename):
    ''' Reads the elastic tensor from the OUTCAR. 
    Args:
        filename : the name of the vasp OUTCAR
    Returns:
        elastic_tensor : 6x6 tensor of the elastic moduli
    '''
    f = open(filename,"r")
    lines = f.readlines()
    f.close()
    copy = False
    elastic_tensor = []
    for line in lines:
        inp = line.split()
        if inp == []:
            continue 
        if len(inp) < 4 or len(inp) > 7:
            continue
        if len(inp) == 4 and inp[0] == 'TOTAL':
            copy = True
        if copy:
            if len(inp) == 7 and len(inp[0]) == 2:
                elastic_tensor.append(inp[1:])
    return np.asarray(elastic_tensor).astype(np.float)




def VoigtMat():
    a = np.asarray([[0, 5, 4], [5, 1, 3], [4, 3, 2]])
    return a
def SVoigtCoeff(p,q):
    return 1/(np.ceil((p+1.)/3.)*np.ceil((q+1.)/3.))
def Smat(SVoigt):
    VM = VoigtMat()
    SM = np.zeros(shape=(3,3,3,3))
    for i in 0,1,2:
        for j in 0,1,2:
            for k in 0,1,2:
                for l in 0,1,2:
                    SM[i,j,k,l] = SVoigtCoeff(VM[i, j], VM[k, l]) \
                    *SVoigt[VM[i, j], VM[k, l]]
    return SM
def dirVector(theta,phi):
    a =  [np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)]
    return a
    
def dirVector2(theta,phi,chi):
    a = [ np.cos(theta)*np.cos(phi)*np.cos(chi) - np.sin(phi)*np.sin(chi), \
    np.cos(theta)*np.sin(phi)*np.cos(chi) + \
    np.cos(phi)*np.sin(chi), - np.sin(theta)*np.cos(chi)]
    return a
def YoungModulus(theta,phi,Sij):
    a = dirVector(theta, phi)
    denom = 0
    SM = Smat(Sij)
    for i in 0,1,2:
        for j in 0,1,2:
            for k in 0,1,2:
                for l in 0,1,2:
                    denom = denom + a[i]*a[j]*a[k]*a[l]*SM[i,j,k,l]
    mod = 1/denom
    return mod

def ShearModulus(theta,phi,chi,Sij):
    a = dirVector(theta, phi)
    b = dirVector2(theta,phi,chi)
    denom = 0
    SM = Smat(Sij)
    for i in 0,1,2:
        for j in 0,1,2:
            for k in 0,1,2:
                for l in 0,1,2:
                    denom = denom + a[i]*b[j]*a[k]*b[l]*SM[i,j,k,l]
    mod = 1/(4*denom)
    return mod
def linearCompressibility(theta, phi, Sij):
    a = dirVector(theta, phi)
    comp = 0
    SM = Smat(Sij)
    for i in 0,1,2:
        for j in 0,1,2:
            for k in 0,1,2:
                    comp = comp + a[i]*a[j]*SM[i,j,k,k]
    return comp
path0='/Users/marco/Dropbox/Work/WORKING_DIR_UNIMELB/GraphiteSiBatteries/Data'
# 
# os.path.

# sys.exit(-1)

path=path0+'/Anodes2ndPhases/AnodesElastic/SiLiC2/'

elastic_tensor = get_elastic_tensor(path+'OUTCAR')

Cij = elastic_tensor/10


Sij = np.linalg.inv(Cij)


df = pd.DataFrame(Cij)

df.to_csv(path+'Cij.csv', index=False)


Kv = ((Cij[0,0] + Cij[1,1] + Cij[2,2]) + 2 * (Cij[0,1] + Cij[1,2] + Cij[2,0])) / 9
Kv

Kr = 1/((Sij[0,0] + Sij[1,1] + Sij[2,2]) + 2 * (Sij[0,1] + Sij[1,2] + Sij[2,0])) 
Kr
     

Gv = (4 * (Cij[0,0] + Cij[1,1] + Cij[2,2]) - 4 * (Cij[0,1] + Cij[1,2] + Cij[2,0]) + 3 * (Cij[3,3] + Cij[4,4] + Cij[5,5]))/15
Gv

Gr = 15 / (4 * (Sij[0,0] + Sij[1,1] + Sij[2,2]) - 4 * (Sij[0,1] + Sij[1,2] + Sij[2,0]) + 3 * (Sij[3,3] + Sij[4,4] + Sij[5,5]))
Gr

Kvrh = (Kv + Kr)/2
Kvrh

Gvrh = (Gv + Gr)/2
Gvrh

mu = (3 * Kvrh - 2 * Gvrh) / (6 * Kvrh + 2 * Gvrh )
mu

# data={'Voigt bulk modulus': Kv, 'Reuss bulk modulus ': Kr, 'Voigt shear modulus': Gv, 'Reuss shear modulus':Gr, 'Voigt-Reuss-Hill bulk modulus': Kvrh, 'Voigt-Reuss-Hill shear modulus': Gvrh, 'Isotropic Poisson ratio': mu}
data={ Kv,Kr, Gv, Gr,  Kvrh,  Gvrh,  mu}



df2 = pd.DataFrame(data,  index=['Voigt bulk modulus', 'Reuss bulk modulus ', 'Voigt shear modulus', 'Reuss shear modulus', 'Voigt-Reuss-Hill bulk modulus', 'Voigt-Reuss-Hill shear modulus', 'Isotropic Poisson ratio'])
df2.to_csv(path+'output.csv')

print(path,df)

print(df2)

sys.exit(-1)

############################################

import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pylab
from pylab import rcParams
# import cubehelix

rcParams['figure.figsize'] = 10, 8
matplotlib.rc('xtick', labelsize=20) 
matplotlib.rc('ytick', labelsize=20)
cx1 = plt.cm.get_cmap('plasma')
cx2 = plt.cm.get_cmap('plasma_r')



a = YoungModulus(0.,0.4475,Sij)


a = np.arange(0,2*np.pi,np.pi/60)
YoungsDirectional = []
for i in a:
    for j in a:
        YoungsDirectional.append(YoungModulus(i,j,Sij))
print(max(YoungsDirectional))

x=y=np.arange(0,2*np.pi,np.pi/10)
X, Y = np.meshgrid(x, y)
point_names = ['0','$\pi/3$','$2\pi$/3', '$\pi$','$4\pi$/3','$5\pi$/3','$2\pi$']
Z = YoungModulus(X,Y,Sij)
im = plt.imshow(Z,cmap=cx2,extent=(0,6,0,6),origin='lower')
cset = plt.contour(Z,np.arange(0,50,10),linewidths=2,cmap=cx1,extent=(0,6,0,6),origin='lower')
pylab.colorbar(im)
#plt.clabel(cset, inline=True, fmt='%1.1f', fontsize=20)
plt.xlabel('$\phi$', fontsize=30)
plt.ylabel('$\Theta$', fontsize=30)
plt.xticks(np.arange(0,6),point_names)
plt.yticks(np.arange(0,6),point_names)
plt.savefig(path+'YoungsModulus.eps')
plt.show()

x=y=np.arange(0,2*np.pi,np.pi/10)
X, Y = np.meshgrid(x, y)
Z = linearCompressibility(X,Y,Sij)
print(np.amin(Z))
im = plt.imshow(Z,cmap=cx2,extent=(0,6,0,6),origin='lower')
cset = plt.contour(Z,np.arange(np.amin(Z),np.amax(Z),0.002),linewidths=2,cmap=cx1,extent=(0,6,0,6),origin='lower')
pylab.colorbar(im)
plt.xlabel('$\phi$', fontsize=30)
plt.ylabel('$\Theta$', fontsize=30)
plt.xticks(np.arange(0,6),point_names)
plt.yticks(np.arange(0,6),point_names)
plt.savefig(path+'LinearCompress.eps')
#plt.clabel(cset, inline=True, fmt='%1.1f', fontsize=20)
plt.show()


x=y=np.arange(0,2*np.pi,np.pi/10)
i = 0
fig, axs = plt.subplots(2,2)
for ax in axs.ravel():
    X, Y = np.meshgrid(x, y)
    Z = ShearModulus(X,Y,i*np.pi/2,Sij)
    cs = ax.contourf(X, Y, Z, cmap=cx2, origin='lower')
    im = plt.imshow(Z,cmap=cx2,extent=(0,6,0,6),origin='lower')
    cset = ax.contour(Z,np.arange(np.amin(Z),np.amax(Z),10),linewidths=2,cmap=cx2,extent=(0,6,0,6),origin='lower')
    #fig.colorbar(cs, ax=ax, shrink=0.9)
    #ax.set_title("extend = %s" % extend)
    ax.locator_params(nbins=4)
    i = i + 1

#print np.amin(Z)

#pylab.colorbar(im)
#plt.clabel(cset, inline=True, fmt='%1.1f', fontsize=20)
plt.show()


