import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
################## Input Reading for Analysis ######################
current = os.getcwd()

trasm = []
trasm2 = []
k_range = []
with os.scandir(current+"/data/Masses/d10/") as it:
    for entry in it:
        if entry.name.endswith(".data") and entry.is_file():
            print(entry.name)
            # print(entry.name, entry.path)
            current_k = float((entry.name).split(".")[0])
            current_k += float((entry.name).split(".")[1])/10

            k_range.append(current_k)

            # read data
            dataframe = pd.read_csv(current+"/data/Masses/d10/"+entry.name,header=None, delim_whitespace=True)
            dataset = dataframe.values
            print(dataframe.shape[0])
            ngrid=int(dataset[0,0])
            nstep=int(dataset[0,1])

            print(ngrid, nstep)
            print()

            xgrid=[0.0]*ngrid
            time=[0.0]*nstep
            
            Psi1 = np.zeros(nstep)
            Psi2 = np.zeros(nstep)
            
            #sort dataset
            for i in range(1,nstep+1):
                ii=(i-1)*ngrid+i
#                print(i)
#                print(ii)
#                print(dataset[ii,:])
                time[i-1]=float(dataset[ii,0])
                #print(i-1,"ii",ii,"time",time[i-1])
            
                Psi1[i-1] = sum([float(dataset[k,3]) for k in range(ii+1,ii+ngrid+1)])
                Psi2[i-1] = sum([float(dataset[k,4]) for k in range(ii+1,ii+ngrid+1)])

            trasm.append(Psi2[-1])
            trasm2.append(Psi1[-1])
            plt.figure()
            plt.title("Exact WP dynamics (A) (red. mass)  \n k="+str(current_k))
            plt.plot(Psi1)
            plt.plot(Psi2)
            plt.xlabel('time (atomic units)')
            plt.ylabel('Population distribution')
            plt.legend(['groundstate', 'first ES'])
            plt.savefig(str(current_k)+".png")
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



######################### Plot data ################################
plt.figure()
plt.plot(k_range,trasm,"+")
plt.plot(k_range,trasm2,"*")
plt.title("Transition probability (A): Exact Quantum Dynamics")
plt.grid()
plt.xlabel("Initial momentum")
plt.ylabel("Exited state population")
plt.savefig("trasmission_EQD.png")
# plt.show()
####################################################################

######################### Write data into external file ############
f = open("analysis_new.out", "w")
f.write("Below some information the data of the exact dynamics analysis\n")
f.write("System: Simple Avoided Crossing, with increas mass\n")
# f.write("Coupling coefficient: C=0.003")
f.write("The following data refer to the transmission from the lower state to the upper state.")
f.write("\n")
# np.savetxt(f,k_range)
f.write("k_range= {} \n".format(k_range))
f.write("trasm= {} \n".format(trasm))
# f.write("trasm", trasm)
f.write("\n")
f.write("We used {} time steps".format(nstep)) 
f.close()
####################################################################
