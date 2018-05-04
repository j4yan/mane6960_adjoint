
import sys
import matplotlib.pyplot as plt
import numpy

font = {'family': 'serif',
        'color': 'black',
        'weight': 'bold',
        'size': 30,}

#  linestyle = "ko-"
#  linestyle = "ko-"
labels = ["p=1", "p=2", "p=3", "p=4"]
varnames = [r'$\phi_{\rho}$', r'$\phi_{\rho u}$', r'$\phi_e$']
f = []
data = []
plot = []
ivar = 1
dir = sys.argv[1]
fnames = ["p1.dat", "p2.dat", "p3.dat", "p4.dat"]
#  for i in range(len(sys.argv)-1):
#      f.append(sys.argv[i+1])
for p, fname in enumerate(fnames):
    f.append(dir+"/"+fname)
    data.append(numpy.loadtxt(f[-1], comments = "\"", delimiter=" "))
    #  plot.append(plt.plot(data[-1][:,0], data[-1][:,1], linestyle, linewidth=2., label=label, markersize=10))
    plot.append(plt.plot(data[-1][:,0], data[-1][:,ivar], linewidth=2., label=labels[p], markersize=10))

plt.legend(loc='best', handlelength=4, fontsize=16, ncol=1)
plt.gcf().subplots_adjust(left=0.20)
plt.gcf().subplots_adjust(bottom=0.15)
plt.tick_params(axis='both', labelsize=16)
#  plt.xlim(0.0, 1.2)
#  plt.ylim(0.0, 9.0)
plt.xlabel('x', fontdict=font)
plt.ylabel(varnames[ivar-1], fontdict=font)
plt.grid()
plt.show()
