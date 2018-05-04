import sys
import matplotlib.pyplot as plt
#  import numpy
import numpy as np

font = {'family': 'serif',
        'color': 'black',
        'weight': 'bold',
        'size': 28,}

plt.figure(figsize=(8,8))

p = [0, 1, 2, 3]
plt_t = [0]*4
x = [0]*4
y1 = [0]*4
y2 = [0]*4
y3 = [0]*4

if len(sys.argv) < 2:
    print("Usage: directory")
    sys.exit()

directory = sys.argv[1]
print(directory)

fn = []
data = []
plots = []
linestyles = []
labels = []
fn.append(directory + '/p1_J1.dat')
#  fn.append(directory + '/p2_J1.dat')
#  fn.append(directory + '/p3_J2.dat')
#  fn.append(directory + '/j1j2_p4.dat')
linestyles.append('k-')
linestyles.append('ko')
linestyles.append('k--')
linestyles.append('k^-')
linestyles.append('k<-')
labels.append('functional error')
labels.append('error estimate')
labels.append('corrected func. error')
labels.append('p=4')

i = 0
ivar = 1
for f in fn:
    data.append(np.loadtxt(fn[i]))
    data[i][:,1] = abs(data[i][:,1])
    numMeshes = len(data[i])
    plots.append(plt.plot(data[i][:,0], data[i][:,ivar], linestyles[0], markerfacecolor='none', markeredgewidth=1.5, linewidth=1.5, label=labels[0], markersize=10))
    plots.append(plt.plot(data[i][:,0], data[i][:,ivar+1], linestyles[1], markerfacecolor='none', markeredgewidth=1.5, linewidth=1.5, label=labels[1], markersize=10))
    plots.append(plt.plot(data[i][:,0], data[i][:,ivar+2], linestyles[2], markerfacecolor='none', markeredgewidth=1.5, linewidth=1.5, label=labels[2], markersize=10))

    p_avg = 0.0
    for l in range(numMeshes-1):
        p_rate =  (np.log(data[i][l+1, ivar]) -np.log(data[i][l, ivar]))/(np.log(data[i][l+1, 0]) -np.log(data[i][l, 0]))
        print("p = ", i+1, ", p_rate = ", p_rate)
        p_rate =  (np.log(data[i][l+1, ivar+1]) -np.log(data[i][l, ivar+1]))/(np.log(data[i][l+1, 0]) -np.log(data[i][l, 0]))
        print("p = ", i+1, ", p_rate = ", p_rate)
        p_rate =  (np.log(data[i][l+1, ivar+2]) -np.log(data[i][l, ivar+2]))/(np.log(data[i][l+1, 0]) -np.log(data[i][l, 0]))
        print("p = ", i+1, ", p_rate = ", p_rate)

    p_avg =  (np.log(data[i][0, ivar]) -np.log(data[i][-1, ivar]))/(np.log(data[i][0, 0]) -np.log(data[i][-1, 0]))
    print("p_avg = ", p_avg)
    i += 1
    print()


#  numData_arr = [4, 4, 4, 1]
slope = [2, 2, 4]
for i in range(1):
    #  slope = (i+1) * 2 + 0
    numOfPoints = 1
    numData = len(data[0]) - 0
    #  numData = numData_arr[i]
    w1 = 1.1
    w2 = 0.01
    x[0] = w2*data[0][numData-1,0] + w1*data[0][numData-numOfPoints-1, 0]
    x[1] = x[0]
    x[2] = w1*data[0][numData-1,0] + w2*data[0][numData-numOfPoints-1, 0]
    x[3] = x[0]
    y1[0] = data[0][numData-numOfPoints-1, ivar]
    y1[1] = np.power(x[2], slope[0]) / np.power(x[0], slope[0]) * y1[0]
    y1[2] = y1[1]
    y1[3] = y1[0]
    y2[0] = data[0][numData-numOfPoints-1, ivar+1]
    y2[1] = np.power(x[2], slope[1]) / np.power(x[0], slope[1]) * y2[0]
    y2[2] = y2[1]
    y2[3] = y2[0]
    y3[0] = data[0][numData-numOfPoints-1, ivar+2]
    y3[1] = np.power(x[2], slope[2]) / np.power(x[0], slope[2]) * y3[0]
    y3[2] = y3[1]
    y3[3] = y3[0]


    plt_t[i] = plt.plot(x, y1, "k")
    plt_t[i] = plt.plot(x, y2, "k")
    plt_t[i] = plt.plot(x, y3, "k")

plt.annotate('2:1', xy=(0.25, 0.575), xycoords="axes fraction", fontsize=22)
plt.annotate('4:1', xy=(0.25, 0.10), xycoords="axes fraction", fontsize=22)
#  plt.annotate('5:1', xy=(0.825, 0.25), xycoords="axes fraction", fontsize=22)

plt.legend(loc='best', handlelength=1, fontsize=20)
plt.gcf().subplots_adjust(left=0.175, right=0.95, bottom=0.1, top=0.95)
plt.tick_params(axis='both', labelsize=18)
numData = len(data[0])
plt.xlim(0.8*data[0][numData-1,0], 1.2*data[0][0,0])
plt.ylim(1e-10, 1.e-1)
plt.xlabel('h', fontdict=font)
plt.ylabel('$J_1$ error', fontdict=font)
plt.locator_params(axis='x', numticks=6)
plt.xscale('log')
plt.yscale('log')
plt.grid()

fout = sys.argv[1] + ".png"
#  plt.savefig(fout)
plt.show()


