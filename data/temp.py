import matplotlib.pyplot as plt
import numpy as np

fig,ax1 = plt.subplots(nrows=1,ncols=1,figsize=(12,5))

plt.rcParams['figure.dpi'] = 600 
data = np.loadtxt('./temp.dat')
x = [i for i in range(0,25)]
y = data
# x = data[:,0]
# y = data[:,1]
ax1.scatter(x,y,label='temp')

ax1.legend(loc='best')
ax1.set_title("temp_plot")
ax1.set_xlabel(r"$\beta$")
ax1.set_ylabel("lnZ")
ax1.set_xlim(0,)
ax1.set_ylim(0,)

plt.savefig("temp.png")
#plt.show()