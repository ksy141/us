import numpy as np
import glob
import pandas as pd
kj_to_kcal = 0.239006

files_tmp = glob.glob('pmf[1-4].dat')
files     = [x.split('.')[0] for x in files_tmp]
pmf = {}

for file in files:
    a = pd.read_csv(file + '.dat', delim_whitespace=True, header=None, usecols=[0,1], comment='#')
    pmf[file] = {}
    pmf[file]['x'] = a[0].values
    pmf[file]['y'] = a[1].values

### X post-processing
convert = False
for file in files:
    if (pmf[file]['x'] > 0).all():
        pmf[file]['x'] *= -1
        convert = True

### Y post-processing
for file in files:
    if convert:
        pmf[file]['y'] -= pmf[file]['y'][-1]
    else:
        pmf[file]['y'] -= pmf[file]['y'][0]

for file in files:
    pmf[file]['y'] *= kj_to_kcal


### Calculate average

pmf['pmf'] = {}
pmf['pmf']['x'] = pmf['pmf1']['x']

ys = []
for file in files:
    ys.append(np.array(pmf[file]['y']))
ys = np.array(ys)

pmf['pmf']['y'] = np.average(ys, axis=0)
pmf['pmf']['std'] = np.std(ys, axis=0)


### Write a file
ofile = open('pmf', 'w')
for i in range(len(pmf['pmf']['x'])):
    x = pmf['pmf']['x']
    y = pmf['pmf']['y']
    std = pmf['pmf']['std']
    ofile.write('{: 6.3f} {: 6.3f} {: 6.3f} {: 6.3f} {: 6.3f}\n'.format(x[i], y[i], std[i], y[i]-std[i], y[i]+std[i]))
ofile.close()


### Calculate mean force
a     = np.loadtxt('pmf')
diff  = (a[2:,1]-a[:-2,1]) / (a[2:,0]-a[:-2,0])
ofile = open('mf', 'w')
for x, y in zip(a[1:-1,0], diff):
    ofile.write('{: 6.3f} {: 6.3f}\n'.format(x, y))
ofile.close()




