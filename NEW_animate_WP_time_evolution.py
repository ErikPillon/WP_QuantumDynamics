import pandas
import numpy as np
import numpy.linalg as LA
import h5py


def Potential_tully(x, A=0.01, B=1.6, C=0.005, D=1.0):
    """
    One dimensional potential used for the implementation of Tully's algorithm
    """
    V = np.zeros((2,2))
    if x>0:
        V[0,0] = A*(1-np.exp(-B*x))
    else:
        V[0,0] = -A*(1-np.exp(B*x))
    V[1,1] = -V[0,0]
    V[0,1] = C*np.exp(-D*x**2)
    V[1,0] = V[0,1]
    return V

evals = np.zeros((2,2048))
xvect = np.linspace(-30,30,2048)

for count in range(2048):
    e_values1, e_vectors1 = LA.eigh(Potential_tully(xvect[count]))
    evals[0,count] = e_values1[0] 
    evals[1,count] = e_values1[1] 

############## Input Reading for Animation ######################
# read data
dataframe = pandas.read_csv("fort.66",header=None, delim_whitespace=True)
dataset = dataframe.values

ngrid=int(dataset[0,0])
nstep=int(dataset[0,1])

print(ngrid, nstep)
print()

xgrid=[0.0]*ngrid
time=[0.0]*nstep
PsiPsi1=np.zeros((nstep,ngrid))
PsiPsi2=np.zeros((nstep,ngrid))

#sort dataset
for i in range(1,nstep+1):
    ii=(i-1)*ngrid+i
    time[i-1]=float(dataset[ii,0])
    #print(i-1,"ii",ii,"time",time[i-1])

    count=0
    for j in range(1,ngrid+1):
        jj=ii+j
        if dataset[jj,1] != "time":
            count=count+1
            PsiPsi1[i-1,j-1]=float(dataset[jj,1]) #1st state
            PsiPsi2[i-1,j-1]=float(dataset[jj,2])  #2nd state
            #print(time[i-1],count,dataset[jj,0],PsiPsi[i-1,j-1])

    if i == 1:
        for k in range(0,ngrid):
            kk=k+2
            xgrid[k]=float(dataset[kk,0])
           #print(i,time[i-1],k,xgrid[k])

#######################  Animation ##################################

maxval=np.amax(PsiPsi1)
maxval2=np.amax(PsiPsi2)
print(maxval, "maxval")

PsiPsi2=PsiPsi2*maxval/maxval2
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
plt.style.use('seaborn-pastel')



fig = plt.figure()
fig.suptitle('Vertically stacked subplots')

ax1 = plt.subplots(2,1,1)
ax1 = plt.axes(xlim=(-20, 20), ylim=(-0.012, 0.015))
line, = ax1.plot([], [], lw=3)
line2, = ax1.plot([], [], lw=3)
line3, = ax1.plot([], [], lw=3)
line4, = ax1.plot([], [], lw=3)
ax2 = plt.subplots(2,2,3)
ax3 = plt.subplots(2,2,4)


def init():
    line.set_data([], [])
    return line,

def animate(i):
    x = xgrid
    y1 = PsiPsi1[i][:]+evals[0,:]
    y2 = PsiPsi2[i][:]+evals[1,:]
    y3 = evals[0,:]
    y4 = evals[1,:]
    line.set_data(x, y1)
    line2.set_data(x, y2)
    line3.set_data(x, y3)
    line4.set_data(x, y4)
    return line,line2

anim = FuncAnimation(fig, animate, init_func=init,frames=nstep, interval=100, blit=True)


anim.save('WP_evolve.gif', writer='imagemagick')

# anim.save('WP_evolve.mp4', fps=30, extra_args=['-vcodec', 'libx264'])


print("Wavepacket animation is sucessfully generated.")
print("Open the file: WP_evolve.gif")
