from __future__ import division
from sys import path
from WSC_mod import *


import VTK

import cas_test as ct
import numpy as np
import time

import ecriture_data as ec


def calcul(lineairefs, lineairebody, nt, io, filestate_in = "", filestate_out = "") :

    # Reading the input files
    if filestate_in == "" :
        api_wsc.api_execution('ws.in','test.geom')
    else :
        api_wsc.api_execution('ws.in','test.geom', True, filestate_in)
        
    if lineairefs :
        parameters.lineairefs = -1
    else :
        parameters.lineairefs = 0
        
    if lineairebody :
        parameters.lineairebody = -1
    else :
        parameters.lineairebody = 0
        
    parameters.nt = nt
    
    # Creating the geometry
    api_wsc.api_geometry()

    
    # Creating the mesh
    api_wsc.api_mesh()
    
    N_f, N_n, N_body = api_wsc.import_mesh_dim()
    print("bla : ", N_f, N_n, N_body,"---------------------------------------")
    L_P, L_ds, L_T, L_G, L_N, L_Rmax, L_type = api_wsc.import_mesh(N_f, N_n)



    # Calling all the subroutines which are used before the beginning of the temporal loop.
    api_wsc.pre_temporal_loop()
          

    t1 = time.time()
    

    api_wsc.ti = parameters.t0
    api_wsc.api_write_wp(io,1)
    

    # Temporal loop
    for jt in range(1,parameters.nt+1): # +1 to reach jt = nt
    
    
        
        N_f, N_n, N_body = api_wsc.import_mesh_dim()
        
        

        
        print("bla : ", N_f, N_n, N_body,"---------------------------------------")
        
        t2 = time.time()
        print()
        print("pas de temps : {} ({}s) ".format(jt, t2-t1))
        print("Temps : ", api_wsc.ti)
        

        
        # RK4 loop
        for jk in range(1,4+1): # +1 to reach jk = 4
        
        
            print(parameters.deformmesh)
            print("*******Debut", jk,jt)
            input()
    
            # Computation of the time step of the RK4 step and the current time
            api_wsc.api_time_step_current_time(jt,jk)

            # Influence coefficient parameter
            api_wsc.api_rci_manager(jk)
    
            # Lissage
            if((jk == 1) and (jt !=1)):
                api_wsc.api_lissage(jt)
    
            
            #----------------------------------------------------------------------------------
            #           Computation of the time differentiation of the state vector
            #----------------------------------------------------------------------------------
    
            # Computing the incident flow
            api_wsc.api_incident()
    
            # Plots
            if(jk == 1):
                api_wsc.api_plots()
                
            
    
            # New mesh if necessary
            if ((jk == 1) and (jt != 1) and parameters.deformmesh == -1):
                # New mesh
                api_wsc.api_regeneration_mesh()
                

                # Deallocation and new allocation
                api_wsc.new_allocation(jt)
    
                # Zeroing the quantities for the new points
                api_wsc.api_zeroing_new_points()
            elif(jk != 1):
                print("Pif")
                api_wsc.boolremesh = False # 0 = False
    
            
            print("*******regene", jk,jt)
            input()

            
            
            # Reinitialisation de Mesh0 et Ecoulement0 avant la resolution de la passe RK1
            if (jk == 1):
                api_wsc.initialization_mesh0_ecoulement0()
    
            # Body condition
            api_wsc.api_bodycondition()
    
            # Calcul des vitesses des noeuds du maillage (updating the velocity of the nodes)
    
            if (parameters.deformmesh == -1): # -1 = True
                api_wsc.api_meshvel()
            
    
    
            
            #Resolution du probleme surfacique
            api_wsc.api_solbvp(False)
    
            
            # Computation of the gradient on the floater (surface and normal gradient)
            api_wsc.api_gradient()
    
            # Computation of DPhiDt and DEtaDt and spatial differentiations
            api_wsc.api_derive()
    
            # Hydrodynamic loads
            if (parameters.freebodies == -1): # -1 = True
    
                #                     FreeBodyMotion
    
                # Allocation of A, B and Sol
                api_wsc.allocation_freebodymotion()
    
                # Computation of q'%perturbation
                api_wsc.api_accconvec()
    
                # Computation of the inertia matrix of the body
                api_wsc.api_inertia()
    
                # Initialisation of the linear system
                # api_wsc.api_systlin_freebody()
    
                #                    SystLin_FreeBody
    
                # Computation of M, CK, CT, Q and Th
                api_wsc.computation_m_ck_ct_q_th()
    
                # Building of A
                api_wsc.api_building_a()
    
                # Building of B
                api_wsc.api_building_b()
    
                # Initialization of the solution
                api_wsc.api_initialization_solution()
    
                # Preconditioning.
                api_wsc.api_preconditioning()
    
                # Deallocation.
                api_wsc.api_deallocating()
    
                #                  End of SystLin_freeBody
    
                # Resolution of the linear system
                api_wsc.solving_systlin()
    
                # Distribution of the solutions
                api_wsc.api_postsl_freebody()
    
                # PlotBVPPhit and deallocation
                api_wsc.api_plotbvpphit_deallocation()
    
                #                  End of FreeBodyMotion
    
            else: # not(freebodies) = Still floater
                if (parameters.finitedifference_forcecalcul == -1): # -1 = True
                    if (jk == 1): # Only for the first RK4 step
                        api_wsc.api_finite_differences(jt)
                else: # not(FiniteDifference_ForceCalcul)
                    api_wsc.forcebodymotion_solbvp()
    
            # Writting of the hydrodynamic loads in the output files
            if (jk == 1):
                api_wsc.api_plotforces()
    
            # Creating of the internal state vector for this RK step
            api_wsc.api_rk_manager(jk)
    
            #----------------------------------------------------------------------------------
            #               Preparation of the state vector for the next RK4 step
            #----------------------------------------------------------------------------------
    
            if (jk != 4): # No updating for the 4th RK4 step.
    
                # Computation of the time step of the RK4 step and the current time
                api_wsc.api_time_step_current_time(jt,jk+1)
    
                # Position of the floater
                api_wsc.api_bodymotion()
    
                # Mise a jour de la vitesse du corps au temps ti
                api_wsc.api_bodyvelocity()
    
                # Phi_p, Eta_p
                api_wsc.api_free_surface_rk_step(jk)
    
                # Mise a jour de la position des noeuds du maillage a partir de la vitesse des noeuds au temps ti (deformation)
                if (parameters.deformmesh == -1): # -1 = True
    
                    api_wsc.api_remesh()
                    


    
        # Time-stepping
        api_wsc.api_time_stepping()
   
        # Writting the time info in the command window
        api_wsc.api_writting_time_info(jt)
        
        api_wsc.api_write_wp(io,1)


    
    if filestate_out != "":
        api_wsc.api_write_state(filestate_out,jt)

    # Closing of the output files
    api_wsc.api_closing()
    # Deallocating
    api_wsc.deallocating_temporal_loop()

    t2 = time.time()
    
    print("Temps total : {}".format(t2-t1))
    
    
## Main

api_wsc.api_open_file_debug()



N_iterations = 10

N_ordis = 2

io = 7988

api_wsc.api_execution('ws.in','test.geom')

pas_ordi = parameters.nt // N_ordis # A affiner

api_wsc.api_parareal_init('ws.in','test.geom', N_iterations, N_ordis)



# Initialisation

for i_ordi in range(N_ordis):
    
    print()
    print()
    print("----Ordinateur n°{}----------".format(i_ordi))
    print()
    print()

    filename = "output_WP/WP_G_{}_{}.dat".format(0, i_ordi)
    api_wsc.api_open_file_wp(filename, io)
        
        
    if i_ordi == 0 :
        filestate_in = ""
    else:
        filestate_in = "lambda_{}_{}.dat".format(0, i_ordi)
    
    filestate_out = "lambda_{}_{}.dat".format(0, i_ordi+1)
    
    #calcul(True, True, pas_ordi, io, filestate_in, filestate_out) 
    
    ## A SUPPRIMER
    calcul(False, False, pas_ordi, io, filestate_in, filestate_out)
    
    ## 
    api_wsc.api_parareal_save_g(0, i_ordi)

    api_wsc.api_close_file_wp(io)
    
    

print()
print("Fin de l'init !!! ")
print()


## A SUPPRIMER


i_iter = 0
i_ordi = 0
    
print()
print()
print("----Ordinateur n°{}----------".format(i_ordi))
print()
print()

    
filename = "output_WP/WP_F_{}_{}.dat".format(i_iter, i_ordi)
api_wsc.api_open_file_wp(filename, io)

if i_ordi == 0 :
    filestate_in = ""
else:
    filestate_in = "lambda_{}_{}.dat".format(i_iter, i_ordi)

filestate_out = ""

print("****Calcul")
input()
calcul(False, False, pas_ordi, io, filestate_in, filestate_out) 



# api_wsc.api_parareal_save_f(i_ordi)
# 
# api_wsc.api_close_file_wp(io)
#         
##

# # Iterations
# 
# for i_iter in range(N_iterations) :
#     
#     # Calcul fin (en parallele)
#     
#     for i_ordi in range(N_ordis) :
#         
#         
#         print()
#         print()
#         print("----Ordinateur n°{}----------".format(i_ordi))
#         print()
#         print()
#     
#             
#         filename = "output_WP/WP_F_{}_{}.dat".format(i_iter, i_ordi)
#         api_wsc.api_open_file_wp(filename, io)
#         
#         if i_ordi == 0 :
#             filestate_in = ""
#         else:
#             filestate_in = "lambda_{}_{}.dat".format(i_iter, i_ordi)
#     
#         filestate_out = ""
# 
#         print("****Calcul")
#         calcul(False, False, pas_ordi, io, filestate_in, filestate_out) 
# 
#         api_wsc.api_parareal_save_f(i_ordi)
#         
#         api_wsc.api_close_file_wp(io)
#         
#         
#     # Calcul grossier (en sequentiel)
#     
#     for i_ordi in range(N_ordis):
#         
#         filename = "output_WP/WP_G_{}_{}.dat".format(i_iter+1, i_ordi)
#         api_wsc.api_open_file_wp(filename, io)
#         
#         if i_ordi == 0 :
#             filestate_in = ""
#         else:
#             filestate_in = "lambda_{}_{}.dat".format(i_iter+1, i_ordi)
#         
#         filestate_out = ""
# 
#         calcul(True, True, pas_ordi, io, filestate_in, filestate_out) 
# 
#         api_wsc.api_parareal_save_g(i_iter+1, i_ordi)
# 
#         filestate_out = "lambda_{}_{}.dat".format(i_iter + 1, i_ordi+1)
#         api_wsc.api_parareal_calcul_lambda(i_iter, i_ordi, 1, filestate_out) #jt = 1 (arbitraire, a verifier)
#         
#         api_wsc.api_close_file_wp(io)
# 
#     input()
    
api_wsc.api_close_file_debug()
    
# for lo in range(3) :
#     
#     nt = 5
#     filestate = ""
#     
#     if lo > 0 :
#         filestate = "State_{}.dat".format(nt)
#     
#     calcul(True, True, nt, filestate) 
#     
# 
# 
