import numpy as np


def save_mesh(L_P, L_T, filename):
    
    
    N_p = L_P.shape[0]
    N_f = L_T.shape[0]
    
    
    fichier = open(filename,"w")
    
    fichier.write("# vtk DataFile Version 2.0\n")
    fichier.write("maillage\n")
    fichier.write("ASCII\n")
    fichier.write("\n")
    fichier.write("DATASET UNSTRUCTURED_GRID\n")
    fichier.write("\n")
    
    fichier.write("POINTS {} float\n".format(N_p))
    
    for i in range(N_p):
        fichier.write("{} {} {}\n".format(L_P[i][0], L_P[i][1], L_P[i][2]))
        
    fichier.write("CELLS {} {}\n".format(N_f, 4*N_f))
    
    for i in range(N_f):
        fichier.write("3 {} {} {}\n".format(L_T[i][0], L_T[i][1], L_T[i][2]))
        
    fichier.write("CELL_TYPES {}\n".format(N_f))
    
    for i in range(N_f):
        fichier.write("5\n")
        
    fichier.close()
        

def add_mesh(L_values, name_values, filename):
    
    N_p = L_values.shape[0]

    fichier = open(filename,"a")
    
    fichier.write("POINT_DATA {}\n".format(N_p))
    fichier.write("Scalars {} float\n".format(name_values))
    fichier.write("LOOKUP_TABLE default\n")
    
    for i in range(N_p):
        fichier.write("{}\n".format(L_values[i]))
        
    fichier.close()