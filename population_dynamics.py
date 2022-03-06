import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

################## Input Reading for Analysis ######################
# read data
dataframe = pd.read_csv("data/Dual_Avoided_Crossing/40.0.data",header=None, delim_whitespace=True)
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
####################################################################

col_to_read = ["C2"]

dataframe_ehr = pd.read_csv("data/Dual_Avoided_Crossing/ehrenfest/ehr.40.0.csv", usecols=col_to_read, delimiter=',')
dataframe_ehr_minus = pd.read_csv("data/Dual_Avoided_Crossing/ehrenfest/ext_ehr_minus.40.0.csv", usecols=col_to_read,  delimiter=',')
dataframe_ehr_plus = pd.read_csv("data/Dual_Avoided_Crossing/ehrenfest/ext_ehr_plus.40.0.csv", usecols=col_to_read,  delimiter=',')

time2 = [ 0.1*x for x in range(0,len(dataframe_ehr))]
time3 = [ 0.03*x for x in range(0,len(dataframe_ehr_minus))]
time4 = [ 0.03*x for x in range(0,len(dataframe_ehr_plus))]

plt.figure()
plt.title('Population dynamics \n hk=40 ')
plt.xlabel('time (a.u.)')

plt.plot(time,Psi2)
plt.plot(time2,dataframe_ehr)
plt.plot(time3,dataframe_ehr_minus)
plt.plot(time4,dataframe_ehr_plus)
plt.legend(['ED','Ehr','Ehr min','Ehr plus'])
plt.savefig('pop_dynamics_40k.png', dpi=300)
plt.show()


######################### Write data into external file ############
f = open("pdyn_analysis.out", "w")
f.write("Below some information the data of the exact population dynamics analysis\n")
f.write("System: Dual Avoided Crossing\n")
# f.write("Coupling coefficient: C=0.003")
f.write("The following data refer to the transmission from the lower state to the upper state.")
f.write("\n")
# np.savetxt(f,k_range)
f.write("k= 40")
f.write("time = {}\n".format(time))
f.write("Psi1= {} \n".format(Psi1))
# f.write("trasm", trasm)
f.write("\n")
f.write("We used {} time steps".format(nstep)) 
f.close()
####################################################################
