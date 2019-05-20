# Plot for Ariana Email Model C
# Only for degeneracies:
import matplotlib.pyplot as plt
import numpy as np
import os

fig,ax1= plt.subplots(nrows=1,ncols=1,figsize=(12,5))

data = np.loadtxt('./old_mc.dat')
x = data[:,0]/24
y = data[:,1]
ax1.scatter(x,y,s=20,label="old_C++_MC result")

data = np.loadtxt('./new_mc.dat')
x = data[:,0]/24
y = data[:,1]
ax1.scatter(x,y,s=20,marker='+',label="new_Julia_WL result")


data = np.loadtxt('./ariana_result.dat')
x = data[:,0]/24
y = data[:,1]
ax1.scatter(x,y,s=20,marker='x',label="Ariana's result")

ax1.legend(loc='best')
ax1.set_title("Energy Spectrum for Model B")
ax1.set_xlabel("Energy")
ax1.set_ylabel("Degeneracies")
ax1.set_ylim(0,max(y))
ax1.set_xlim(0,1)

plt.savefig('fig00.png')