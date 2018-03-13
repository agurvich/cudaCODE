import matplotlib.pyplot as plt
import numpy as np 
from all_utils import nameAxes
from distinct_colours import get_distinct

input_values = np.genfromtxt('sample_input.txt')
input_values = input_values.T[:-1].T # strip the last row which has a hanging delimiter

method = 'rk4'
vals = np.genfromtxt('output.txt',unpack=1)

tss = vals[0]

nelements = np.sum(tss==tss[0])
tss = tss.reshape(-1,nelements).T
ts = tss[0] # same for everyone

nYs = len(vals[1:]) #number of ys
ysss = []
for i,col in enumerate(vals[1:]):
    ysss+=[col.reshape(-1,nelements).T]

ysss = np.array(ysss) # 3 neqn  x 4096 elements x 200 time steps

ics = np.genfromtxt('sample_input.txt',unpack=1)
y0s,g0s=np.array([ics[:len(ics)/2],ics[len(ics)/2:]])
#ics = ics.reshape(ics.shape[0],2,-1) # 4096 elements x ys/gs x 3 neqn
#ics = np.transpose(ics,axes=[0,2,1]) # 4096 elements x 3 neqn x ys/gs 

#plt.plot(ts,ysss[2][0],lw=3)

def eqn1(t,y0,g0):
    return y0 + g0*t

def eqn2(t,y0,g0):
    return -1./g0*np.sin(g0*t)+y0

def eqn3(t,y00,y10,y20,g00,g10,g20):
    return (0.5*g00*t**2 + (y10+y00)*t + (1.0/g10**2)*np.cos(g10*t) )*g20 + y20

colors = get_distinct(nYs)
maxnum = 16
fig,axs = plt.subplots(3,1,sharex=True)
## calculate y0
ax = axs[0]
theo_y1ss =[]
for element in ics.T[:maxnum]:
    y0,g0 = element[0],element[len(element)/2 ]
    theo_y1ss+=[eqn1(ts,y0,g0)]

ax.plot(ts,ysss[0][0],c=colors[0],lw=3)
ax.plot(ts,theo_y1ss[0],'k--',lw=3)
nameAxes(ax,None,0,'y0',supertitle = method)

ax = axs[1]
theo_y2ss = []
## caclulate y1
for element in ics.T[:maxnum]:
    y0,g0 = element[1],element[1+len(element)/2 ]
    theo_y2ss+=[eqn2(ts,y0,g0)]

ax.plot(ts,ysss[1][0],c=colors[1],lw=3)
ax.plot(ts,theo_y2ss[0],'k--',lw=3)
nameAxes(ax,None,0,'y1',supertitle = method)

ax = axs[2]
## calculate y2
theo_y3ss = []
for element in ics.T[:16]:
    y3s = eqn3(ts,*element)
    theo_y3ss+=[y3s]

ax.plot(ts,ysss[-1][0],c=colors[2],lw=3)
ax.plot(ts,theo_y3ss[0],'k--',lw=3)
nameAxes(ax,None,'t (sec)','y2',supertitle = method)
plt.subplots_adjust(hspace=0)
fig.set_size_inches(16,9)

plt.savefig('curve_overlay')

"""
#rts,rys0,rys1 = np.genfromtxt('results.txt',unpack=1)

fn = lambda ts: np.exp(3./2 * ts**2)
fn = lambda ts: -np.sin(ts)

plt.plot(ts,ys1,'r',lw=3,label='Riemann')
plt.plot(rts,rys1,'g:',lw=3,label='RK4')
plt.plot(ts,fn(ts),'k--',lw=3)

plt.plot(ts,[0]*len(ts))
plt.plot([0,0],[-2,1])

nameAxes(plt.gca(),None,'t',None,supertitle = 'dt = %.1e'%(ts[1]-ts[0]),logflag=(1,0),make_legend=1)

plt.savefig('compare')
"""
