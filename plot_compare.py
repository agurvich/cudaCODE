import matplotlib.pyplot as plt
import numpy as np 
from all_utils import nameAxes
from distinct_colours import get_distinct

input_values = np.genfromtxt('sample_input.txt')
input_values = input_values.T[:-1].T # strip the last row which has a hanging delimiter

method = 'rk4'
dt = 0.05

vals = np.genfromtxt('%s_%.2f_output.txt'%(method,dt),unpack=1)

tss = vals[0]

nelements = np.sum(tss==tss[0])
tss = tss.reshape(-1,nelements).T
ts = tss[0] # same for everyone
assert dt == ts[1]-ts[0]

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
maxnum = nelements

## calculate y0
theo_y1ss =[]
for element in ics.T[:maxnum]:
    y0,g0 = element[0],element[len(element)/2 ]
    theo_y1ss+=[eqn1(ts,y0,g0)]

theo_y2ss = []
## caclulate y1
for element in ics.T[:maxnum]:
    y0,g0 = element[1],element[1+len(element)/2 ]
    theo_y2ss+=[eqn2(ts,y0,g0)]

## calculate y2
theo_y3ss = []
for element in ics.T[:maxnum]:
    y3s = eqn3(ts,*element)
    theo_y3ss+=[y3s]

def plotCurveOverlay(index):
    fig,axs = plt.subplots(3,1,sharex=True)
    ax = axs[0]
    ax.plot(ts,ysss[0][index],c=colors[0],lw=3)
    ax.plot(ts,theo_y1ss[index],'k--',lw=3)
    nameAxes(ax,None,0,'y0',supertitle = method)

    ax = axs[1]
    ax.plot(ts,ysss[1][index],c=colors[1],lw=3)
    ax.plot(ts,theo_y2ss[index],'k--',lw=3)
    nameAxes(ax,None,0,'y1',supertitle = method)

    ax = axs[2]
    ax.plot(ts,ysss[-1][index],c=colors[2],lw=3)
    ax.plot(ts,theo_y3ss[index],'k--',lw=3)
    nameAxes(ax,None,'t (sec)','y2',supertitle = method)
    plt.subplots_adjust(hspace=0)
    fig.set_size_inches(16,9)

    plt.savefig('plots/%s_%.2f_curve_overlay'%(method,dt))

for cindex in xrange(16):
    plotCurveOverlay(cindex)

fig,axs = plt.subplots(3,1,sharex=True)
####### calculate residuals over time

## calculate y0
ax = axs[0]
ax.plot(ts,theo_y1ss[0]-ysss[0][0],c=colors[0],lw=3)
ax.plot(ts,0*ts,'k--',lw=3)
nameAxes(ax,None,0,'yt0-y0',supertitle = method)

ax = axs[1]
ax.plot(ts,theo_y2ss[0]-ysss[1][0],c=colors[1],lw=3)
ax.plot(ts,0*ts,'k--',lw=3)
nameAxes(ax,None,0,'yt1-y1',supertitle = method)

ax = axs[2]
ax.plot(ts,theo_y3ss[0]-ysss[-1][0],c=colors[2],lw=3)
ax.plot(ts,0*ts,'k--',lw=3)
nameAxes(ax,None,'t (sec)','yt2-y2',supertitle = method)
plt.subplots_adjust(hspace=0)
fig.set_size_inches(16,9)

plt.savefig('plots/%s_%.2f_residuals'%(method,dt))

## calculate percent error
plt.figure()
ax = plt.gca()

mega_residuals0 =[(theo_y1ss[i]-ysss[0][i])/theo_y1ss[i] for i in xrange(maxnum)] 
mega_residuals1 =[(theo_y2ss[i]-ysss[1][i])/theo_y2ss[i] for i in xrange(maxnum)] 
mega_residuals2 =[(theo_y3ss[i]-ysss[2][i])/theo_y3ss[i] for i in xrange(maxnum)] 

ax.plot(ts,0*ts,'k--',lw=3)

ax.plot(ts,mega_residuals0[0],lw=3,c=colors[0])
ax.plot(ts,mega_residuals1[0],lw=3,c=colors[1])
ax.plot(ts,mega_residuals2[0],lw=3,c=colors[2])


nameAxes(ax,None,'t (sec)','|(yt-y)/yt|',subtitle = method)
fig.set_size_inches(8,4.5)
plt.savefig('plots/%s_%.2f_percent_errors'%(method,dt))

fig,axs = plt.subplots(4,4,sharey=True,sharex=True)

axs = axs.flatten()

for i in xrange(16):
    ax = axs[i]
    ax.plot(ts,np.abs(mega_residuals0[i]),lw=3,c=colors[0])
    ax.plot(ts,np.abs(mega_residuals1[i]),lw=3,c=colors[1])
    ax.plot(ts,np.abs(mega_residuals2[i]),lw=3,c=colors[2])

ax.set_ylim(0,0.05)
plt.subplots_adjust(hspace=0,wspace=0)

fig.set_size_inches(24,24)
plt.savefig('plots/%s_%.2f_mega_percent_errors'%(method,dt))

## avg percent error
plt.figure()
ax = plt.gca()

ax.plot(ts,0*ts,'k--',lw=3)
ax.plot(ts,np.mean(np.abs(mega_residuals0),axis=0),lw=3,c=colors[0])
ax.plot(ts,np.mean(np.abs(mega_residuals1),axis=0),lw=3,c=colors[1])
ax.plot(ts,np.mean(np.abs(mega_residuals2),axis=0),lw=3,c=colors[2])

nameAxes(ax,None,'t (sec)',r'$\langle$|(yt-y)/yt|$\rangle$',subtitle = method,supertitle='%d systems'%maxnum)
fig.set_size_inches(8,4.5)
plt.savefig('plots/%s_%.2f_avg_percent_errors'%(method,dt))
