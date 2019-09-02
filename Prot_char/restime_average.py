###Residence Time Plotting###

import numpy as np
import matplotlib.pyplot as plt
from mpltools import style
from mpltools import layout
import matplotlib as mpl

a = open('outputA.txt','r')
restimea=[]
for x in a:
    x = [int(i) for i in x.split()]
    if x:
        restimea.append(np.mean(x))
    elif not x:
        restimea.append(0)

a.close()

b = open('outputB.txt','r')
restimeb=[]
for x in b:
    x = [int(i) for i in x.split()]
    if x:
        restimeb.append(np.mean(x))
    elif not x:
        restimeb.append(0)

b.close()

combineab = np.maximum(restimea,restimeb)

c = open('outputC.txt','r')
restimec=[]
for x in c:
    x = [int(i) for i in x.split()]
    if x:
        restimec.append(np.mean(x))
    elif not x:
        restimec.append(0)

c.close()

combineabc = np.maximum(combineab,restimec)

d = open('outputD.txt','r')
restimed=[]
for x in d:
    x = [int(i) for i in x.split()]
    if x:
        restimed.append(np.mean(x))
    elif not x:
        restimed.append(0)
d.close()

combineabcd = np.maximum(combineabc,restimed)

e = open('outputE.txt','r')
restimee=[]
for x in e:
    x = [int(i) for i in x.split()]
    if x:
        restimee.append(np.mean(x))
    elif not x:
        restimee.append(0)
e.close()

combineabcde = np.maximum(combineabcd,restimee)

f = open('outputF.txt','r')
restimef=[]
for x in f:
    x = [int(i) for i in x.split()]
    if x:
        restimef.append(np.mean(x))
    elif not x:
        restimef.append(0)
f.close()

combineabcdef = np.maximum(combineabcde,restimef)

g = open('outputG.txt','r')
restimeg=[]
for x in g:
    x = [int(i) for i in x.split()]
    if x:
        restimeg.append(np.mean(x))
    elif not x:
        restimeg.append(0)
g.close()

combineabcdefg = np.maximum(combineabcdef,restimeg)

h = open('outputH.txt','r')
restimeh=[]
for x in h:
    x = [int(i) for i in x.split()]
    if x:
        restimeh.append(np.mean(x))
    elif not x:
        restimeh.append(0)
h.close()

combine = np.maximum(combineabcdefg,restimeh)

combine=np.nan_to_num(combine)
resid=[]
for i in range(1,437):
    resid.append(i+2)

time=[x*20 / 1000  for x in combine]

def plot(x, y, axis):
    axis.plot(x,y, marker='o', color='black', linestyle='None')
    axis.set_xlim([0, 436])
    axis.xaxis.set_ticks(np.arange(0,436,50))
    plt.rcParams.update({'font.size': 12})

fig = plt.gcf()
fig.set_size_inches(6,6)
axes1= fig.add_subplot(111)

mydict={'res5':[resid,time,axes1]}

for item in mydict:
    plot(mydict[item][0], mydict[item][1], mydict[item][2])


axes1.set_ylabel('Residence Time (ns)', fontsize=14)
axes1.set_xlabel('Residue Number', fontsize=14)
axes1.tick_params(axis='x', pad=6)
axes1.tick_params(axis='both', which='major', labelsize=12)
axes1.tick_params(axis='both', which='minor', labelsize=12)

fig.savefig('res.pdf', format='pdf', bbox_inches='tight',dpi=500)
