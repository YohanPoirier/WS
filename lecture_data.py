import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


fichier = open("Wave_elevation_ 1_PYW.dat","r")
lignes = fichier.readlines()
lignes = lignes[2:]


Nt = len(lignes)

data = np.zeros((4,Nt))

for i in range(Nt):
    data[:,i] = list(map(float,lignes[i].split()))

for i in range(3):
    plt.plot(data[0], data[1+i], color = "b")
    
    
fichier = open("Wave_elevation_ 1_PYW2.dat","r")
lignes = fichier.readlines()
lignes = lignes[2:]


Nt = len(lignes)

data = np.zeros((4,Nt))

for i in range(Nt):
    data[:,i] = list(map(float,lignes[i].split()))

for i in range(3):
    plt.plot(data[0], data[1+i], color = "r")
    
    
plt.show()