import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt



def plot_elevation(ax, fichier, ylim) :
    lignes = fichier.readlines()
    lignes = lignes[2:]
    
    print(lignes[0])
    data = np.zeros((4,len(lignes)))
    
    for i in range(len(lignes)):
        data[:,i] = list(map(float,lignes[i].split()))

    ax.plot(data[0], data[1], label = "Incidente")
    ax.plot(data[0], data[2], label = "Perturbation")
    ax.plot(data[0], data[3], label = "Totale")
    
    ax.set_ylim(ylim)
    ax.legend()

    fichier.close()
    
def plot_elevation_pert(ax, fichier, ylim):
    lignes = fichier.readlines()
    lignes = lignes[2:]
    
    data = np.zeros(len(lignes))
    
    for i in range(len(lignes)):
        data[:,i] = list(map(float,lignes[i].split()))
    
    ax.plot(data[0], data[2], label = "Perturbation")
    
    ax.set_ylim(ylim)
    ax.legend()

    fichier.close()

    
    
# fichier = open("Output_NL/Wave_elevation_ 1.dat","r")
# lignes = fichier.readlines()
# lignes = lignes[2:]
# 
# 
# Nt = len(lignes)
# 
# data = np.zeros((4,Nt))
# 
# for i in range(Nt):
#     data[:,i] = list(map(float,lignes[i].split()))
# 
# ax2.plot(data[0], data[1], label = "Incidente")
# ax2.plot(data[0], data[2], label = "Perturbation")
# ax2.plot(data[0], data[3], label = "Totale")
# 
# 

#ax2.set_ylim([-0.03,0.03])
    
N_ax = 2

f, ax = plt.subplots(N_ax, 1, sharex=True)
ylim = [-0.03,0.03]




fichier = open("Output_boucle/Wave_elevation_ 1.dat","r")

plot_elevation(ax[0], fichier, ylim)

fichier = open("Output_boucle/Wave_elevation_ 2.dat","r")

plot_elevation(ax[1], fichier, ylim)

plt.show()