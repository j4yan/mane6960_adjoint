
import sys
import matplotlib.pyplot as plt
import numpy

font = {'family': 'serif',
        'color': 'black',
        'weight': 'bold',
        'size': 30,}

#  linestyle = "ko-"
#  linestyle = "ko-"

ivar = 1
dir = sys.argv[1]
degrees = [1, 2, 3, 4]

for p in degrees:
    plt.gcf().clear()
    f = []
    data = []
    plot = []
    fnames = []
    labels = []

    numelems = [20, 40, 80]
    for n in numelems:
        fnames.append('p' + str(p) + '_' + str(n) + ".dat")
        labels.append("numelem=" + str(n))
    
    for n, fname in enumerate(fnames):
        f.append(dir+"/"+fname)
        print(f[-1])
        data.append(numpy.loadtxt(f[-1], comments = "\"", delimiter=" "))
        plot.append(plt.plot(data[-1][:,0], data[-1][:,ivar], linewidth=2., label=labels[n]))

    plt.legend(loc='best', handlelength=4, fontsize=16, ncol=1)
    plt.gcf().subplots_adjust(left=0.20)
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.tick_params(axis='both', labelsize=16)
    #  plt.xlim(0.0, 1.2)
    #  plt.ylim(-0.75, 0.75)
    plt.ylim(-0.06, 0.08)
    plt.xlabel('x', fontdict=font)
    plt.ylabel('$DJ_2 / DA$', fontdict=font)
    plt.grid()
    #  plt.show()
    fout = 'DJ2DA_p' + str(p) + '.png'
    plt.savefig(fout)
