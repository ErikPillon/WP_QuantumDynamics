import matplotlib.pyplot as plt
import numpy as np

k = 12
Y = np.loadtxt("fort.11")

kinetic = Y[:,0]
potential = Y[:,1]
total = kinetic+potential

plt.plot(kinetic)
plt.plot(potential)
plt.plot(total)
plt.xlabel("time (a.u.)")
plt.ylabel('Energy (a.u.)')
plt.title("Total Energy of the exact wavepacket dynamics \n k={}".format(k))
plt.legend(['kinetic','potential', 'total'])
plt.savefig("exact_kinetic_{}k".format(k), dpi=300)
plt.show()
