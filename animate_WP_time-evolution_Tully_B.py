import pandas
import numpy as np
import numpy.linalg as LA
import h5py
import os

from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
plt.style.use('seaborn-pastel')
            
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

def Dual_avoided_crossing(x):
    """
    One dimensional potential used for the implementation of Tully's algorithm
    This corresponds to Model B of Tully's model.
    """
    A, B, C, D = 0.1, 0.28, 0.015, 0.06
    E_0 = 0.05
    V = np.zeros((2,2))
    V[0,0] = 0.0
    V[1,1] = -A*np.exp(-B*x**2)+E_0
    V[0,1] = V[1,0] = C*np.exp(-D*x**2)
    return V

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


evals = np.zeros((2,2048))
xvect = np.linspace(-30,30,2048)

for count in range(2048):
    e_values1, e_vectors1 = LA.eigh(Dual_avoided_crossing(xvect[count]))
    evals[0,count] = e_values1[0] 
    evals[1,count] = e_values1[1] 


current = os.getcwd()

with os.scandir(current+"/data/Dual_Avoided_Crossing/") as it:
    for entry in it:
        if entry.name.endswith(".data") and entry.is_file():
            print(entry.name)
            # print(entry.name, entry.path)
            current_k = float((entry.name).split(".")[0])
            current_k += float((entry.name).split(".")[1])/10

            ############## Input Reading for Animation ######################
            # read data
            dataframe = pandas.read_csv(current+"/data/Dual_Avoided_Crossing/"+entry.name,header=None, delim_whitespace=True)
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
                        PsiPsi1[i-1,j-1]=float(dataset[jj,3]) #1st state
                        PsiPsi2[i-1,j-1]=float(dataset[jj,4])  #2nd state
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
            
            
            fig = plt.figure()
            fig.suptitle('Exact wavepacket dynamics \n $\\hbar k=$ %f'%(current_k))
            fig.supxlabel('Nuclear Coordinates')
            fig.supylabel('Energy')
            
            maximum = maxval2*1.2+0.01
            ax = plt.axes(xlim=(-20, 20), ylim=(-0.055, 0.07))
            line, = ax.plot([], [], lw=3)
            line2, = ax.plot([], [], lw=3)
            line3, = ax.plot([], [], lw=3)
            line4, = ax.plot([], [], lw=3)

            # produce the final sequence
            anim = FuncAnimation(fig, animate, init_func=init,frames=nstep, interval=100, blit=True)
            # save the final gif 
            anim.save('WP_evolve_Tully_B_%f.gif'%(current_k), writer='imagemagick')

# anim.save('WP_evolve.mp4', fps=30, extra_args=['-vcodec', 'libx264'])


print("Wavepacket animation is sucessfully generated.")
print("Open the file: WP_evolve_Tully_B_[k].gif")
