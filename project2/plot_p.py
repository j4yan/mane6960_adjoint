import sys
import matplotlib.pyplot as plt
#  import numpy
import numpy as np

font = {'family': 'serif',
        'color': 'black',
        'weight': 'bold',
        'size': 28,}

plt.figure(figsize=(6,8))

p = [0, 1, 2, 3]
plt_t = [0]*4
x = [0]*4
y = [0]*4

if len(sys.argv) < 2:
    print("Usage: directory")
    sys.exit()

directory = sys.argv[1]
print (directory)

fn = []
data = []
plots = []
linestyles = []
labels = []
fn.append(directory + '/j1j2_p1.dat')
fn.append(directory + '/j1j2_p2.dat')
fn.append(directory + '/j1j2_p3.dat')
fn.append(directory + '/j1j2_p4.dat')
linestyles.append('ko-')
linestyles.append('ks-')
linestyles.append('k^-')
linestyles.append('k<-')
labels.append('p=1')
labels.append('p=2')
labels.append('p=3')
labels.append('p=4')

i = 0
ivar = 1
for f in fn:
    data.append(np.loadtxt(fn[i]))
    data[i][:,1] = abs(data[i][:,1])
    numMeshes = len(data[i])
    plots.append(plt.plot(data[i][:,0], data[i][:,ivar], linestyles[i], markerfacecolor='none', markeredgewidth=1.5, linewidth=1.5, label=labels[i], markersize=10))

    p_avg = 0.0
    for l in range(numMeshes-1):
        p_rate =  (np.log(data[i][l+1, ivar]) -np.log(data[i][l, ivar]))/(np.log(data[i][l+1, 0]) -np.log(data[i][l, 0]))
        print("p = ", i+1, ", p_rate = ", p_rate)

    p_avg =  (np.log(data[i][0, ivar]) -np.log(data[i][-1, ivar]))/(np.log(data[i][0, 0]) -np.log(data[i][-1, 0]))
    print("p_avg = ", p_avg)
    i += 1
    print()


numData_arr = [5, 5, 1, 1]
for i in range(4):
    slope = (i+1) * 2 + 0
    numOfPoints = 1
    #  numData = len(data[i]) - 0
    numData = numData_arr[i]
    if i > 1:
        numData = 2
    w1 = 1.1
    w2 = 0.01
    x[0] = w2*data[i][numData-1,0] + w1*data[i][numData-numOfPoints-1, 0]
    x[1] = x[0]
    x[2] = w1*data[i][numData-1,0] + w2*data[i][numData-numOfPoints-1, 0]
    x[3] = x[0]
    y[0] = data[i][numData-numOfPoints-1, ivar]
    y[1] = np.power(x[2], slope) / np.power(x[0], slope) * y[0]
    y[2] = y[1]
    y[3] = y[0]


    #  d = min(data[i][numData-1, 1], data[i][numData-1, 1])
    #  x[0] = w1*data[i][numData-1,0] + w2*data[i][numData-numOfPoints-1, 0]
    #  #  y[0] = data[i][numData-1, 2]
    #  y[0] = d
    #  x[1] = w2*data[i][numData-1,0] + w1*data[i][numData-numOfPoints-1, 0]
    #  y[1] = y[0]
    #  x[2] = x[1]
    #  y[2] = np.power(x[2], i+2)/np.power(x[0], i+2)*y[0]
    #  x[3] = x[0]
    #  y[3] = y[0]

    plt_t[i] = plt.plot(x, y, "k")

plt.annotate('2:1', xy=(0.15, 0.68), xycoords="axes fraction", fontsize=22)
plt.annotate('3:1', xy=(0.15, 0.26), xycoords="axes fraction", fontsize=22)
plt.annotate('4:1', xy=(0.825, 0.435), xycoords="axes fraction", fontsize=22)
plt.annotate('5:1', xy=(0.825, 0.25), xycoords="axes fraction", fontsize=22)

#  plt.legend(loc=4, handlelength=3, fontsize=20)
plt.gcf().subplots_adjust(left=0.225, right=0.95, bottom=0.1, top=0.95)
plt.tick_params(axis='both', labelsize=18)
numData = len(data[1])
plt.xlim(0.8*data[1][numData-1,0], 1.2*data[1][0,0])
plt.ylim(1e-13, 1.e-1)
plt.xlabel('h', fontdict=font)
plt.ylabel('$J_1$ error', fontdict=font)
plt.locator_params(axis='x', numticks=6)
plt.xscale('log')
plt.yscale('log')
plt.grid()

fout = sys.argv[1] + ".png"
plt.savefig(fout)
#  plt.show()


