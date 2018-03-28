import numpy as np
Q = np.loadtxt("Q.dat")
#print(Q)
#print("The Order Parameter is: ")
print(max(np.linalg.eigh(Q)[0]))
