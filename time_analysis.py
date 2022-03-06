import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
################## Input Reading for Analysis ######################
# read data
dataframe = pd.read_csv("fort.66",header=None, delim_whitespace=True)
dataset = dataframe.values

ngrid=int(dataset[0,0])
nstep=int(dataset[0,1])

print(ngrid, nstep)
print()

xgrid=[0.0]*ngrid
time=[0.0]*nstep
PsiPsi1=np.zeros((nstep,ngrid))
PsiPsi2=np.zeros((nstep,ngrid))

Psi1 = np.zeros(nstep)
Psi2 = np.zeros(nstep)

#sort dataset

for i in range(1,nstep+1):
    ii=(i-1)*ngrid+i
    time[i-1]=float(dataset[ii,0])
    #print(i-1,"ii",ii,"time",time[i-1])

    Psi1[i-1] = sum([float(dataset[k,3]) for k in range(ii+1,ii+ngrid+1)])
    Psi2[i-1] = sum([float(dataset[k,4]) for k in range(ii+1,ii+ngrid+1)])
   # count=0
   # for j in range(1,ngrid+1):
   #     jj=ii+j
   #     if dataset[jj,1] != "time":
   #         count=count+1
   #         PsiPsi1[i-1,j-1]=float(dataset[jj,1])
   #         PsiPsi2[i-1,j-1]=float(dataset[jj,2])
   #         #print(time[i-1],count,dataset[jj,0],PsiPsi[i-1,j-1])

   # if i == 1:
   #     for k in range(0,ngrid):
   #         kk=k+2
   #         xgrid[k]=float(dataset[kk,0])
   #        #print(i,time[i-1],k,xgrid[k])
####################################################################
print(time)
print("last 5 steps of Psi1: ", Psi1[-5:-1])
print("last 5 steps of Psi2: ", Psi2[-5:-1])
plt.figure()
plt.plot(time,Psi1)
plt.plot(time,Psi2)
plt.show()


