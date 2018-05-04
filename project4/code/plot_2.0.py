
import sys
import matplotlib.pyplot as plt
import numpy

font = {'family': 'serif',
        'color': 'black',
        'weight': 'bold',
        'size': 20,}

#  linestyle = "ko-"
#  linestyle = "ko-"
labels = ["t=0", "t=0.25", "t=0.5", "t=2.0"]
f = []
data = []
plot = []
ivar1 = 0
ivar2 = 3

#  n_data = int(sys.argv[1])
#  for i in range(n_data):
    #  f.append("psi_" + str(i) + ".dat")
for i in range(len(sys.argv)-1):
    f.append(sys.argv[i+1])
    data.append(numpy.loadtxt(f[-1], comments = "\"", delimiter=" "))
    #  plot.append(plt.plot(data[-1][:,0], data[-1][:,1], linestyle, linewidth=2., label=label, markersize=10))
    plot.append(plt.plot(data[-1][:,ivar1], data[-1][:,1], linewidth=2., label=r'$\psi_{\rho}$', markersize=10))
    plot.append(plt.plot(data[-1][:,ivar1], data[-1][:,2], linewidth=2., label=r'$\psi_{\rho u}$', markersize=10))
    plot.append(plt.plot(data[-1][:,ivar1], data[-1][:,3], linewidth=2., label=r'$\psi_{E}$', markersize=10))



#  plt.legend(loc='best', handlelength=1, fontsize=16)
plt.gcf().subplots_adjust(left=0.25)
plt.gcf().subplots_adjust(bottom=0.15)
plt.tick_params(axis='both', length=3, labelsize=16)
#  plt.xlim(0.0, 1.2)
#  plt.ylim(0.0, 9.0)
plt.xlabel('x', fontdict=font)
#  plt.ylabel(r'$\psi_{\rho u}$', fontdict=font)
plt.ylabel(r'$\psi$', fontdict=font)
#  plt.xscale('log')
#  plt.yscale('log')
plt.legend(loc='best', handlelength=1, fontsize=20)
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker
import matplotlib.ticker as mticker
ax=plt.gca()
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
g = lambda x,pos : "${}$".format(f._formatSciNotation('%1.10e' % x))
plt.gca().yaxis.set_major_formatter(mticker.FuncFormatter(g))
plt.grid()
plt.show()
