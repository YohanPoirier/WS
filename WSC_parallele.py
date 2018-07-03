from __future__ import division
from sys import path
from WSC_mod import *


import VTK

import cas_test as ct
import numpy as np
import time

import ecriture_data as ec


for lo in range(1) :
    
    
    if lo == 0 :
        # Reading the input files
        api_wsc.api_execution('ws.in','test.geom')
        parameters.lineairefs = 0
        parameters.lineairebody = 0
    else :
        # Reading the input files
        api_wsc.api_execution('ws.in','test.geom', True, "State_{}.dat".format(parameters.nt))
        parameters.lineairefs = 0
        parameters.lineairebody = 0
    
    
    # Creating the geometry
    api_wsc.api_geometry()
    
    
    # Creating the mesh
    api_wsc.api_mesh()
    
    
    # Calling all the subroutines which are used before the beginning of the temporal loop.
    api_wsc.pre_temporal_loop()
          
    t1 = time.time()

    # Temporal loop
    for jt in range(1,parameters.nt+1): # +1 to reach jt = nt
        
        
        t2 = time.time()
        print()
        print("pas de temps : {} ({}s) ".format(jt, t2-t1))
        print("Temps : ", api_wsc.ti)
        
            
        # RK4 loop
        for jk in range(1,4+1): # +1 to reach jk = 4
    
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
                api_wsc.boolremesh = False # 0 = False
    
    
    
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
        
        
        if jt == parameters.nt :
            filename = "State_{}.dat".format(jt)
            print("Temps", api_wsc.ti)
            api_wsc.api_write_state(filename,jt)

    
    # Closing of the output files
    api_wsc.api_closing()

    # Deallocating
    api_wsc.deallocating_temporal_loop()

    t2 = time.time()
    
    print("Temps total : {}".format(t2-t1))


