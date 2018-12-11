from __future__ import division
from sys import path
from WSC_mod import *
from sys import exit


import VTK

import cas_test as ct
import numpy as np
import time

import ecriture_data as ec


def calcul(grossier, nt, io_WP, io_F, filestate_in = "", filestate_out = "") :


    global aff

    
        
    t1 = time.time()
    ta = time.time()
    
    
    # Reading the input files
    if filestate_in == "" :
        api_wsc.api_execution('ws.in','test.geom')
    else :
        api_wsc.api_execution('ws.in','test.geom', True, filestate_in)
        
      
      
      
    # A SUPPRIMER
    change_size(dx1)


    tb = time.time()
        

    
    if grossier :
        parameters.lineairefs = -1
        parameters.lineairebody = -1
        parameters.is_bs_gs = 0
        parameters.bool_coarse = -1
    else :
        parameters.lineairefs = 0
        parameters.lineairebody = 0
        parameters.is_bs_gs = -1
        parameters.bool_coarse = 0
        

    parameters.nt = nt
    
    
    # Creating the geometry
    api_wsc.api_geometry()
    
    
    
    

    if filestate_in == "" :
    # Creating the mesh
        api_wsc.api_mesh()

    else :
    # Importing the mesh
        api_wsc.api_import_mesh_ecoulement()
        
        if not(grossier):
            api_wsc.api_geominit()
            

  
    N_f, N_n, N_body = api_wsc.import_mesh_dim()
    api_wsc.print_file("N_ noeuds = {} , N_facettes = {}, N_body = {}".format(N_n, N_f, N_body), 1111)
    
 
    tc = time.time()



    # Calling all the subroutines which are used before the beginning of the temporal loop.
    api_wsc.pre_temporal_loop_parareal()
    

    td = time.time()
    #api_wsc.save_mesh_ci()

    api_wsc.ti = parameters.t0


    

    if (grossier):
        # Initialization of mesh_ref
        if (parameters.bool_coarse_init == 0) :
            api_wsc.api_init_lineaire()
            parameters.bool_coarse_init = -1
        # Interpolation to mesh_ref
        else:
            api_wsc.api_interpolation_mesh_ref()


    api_wsc.api_parareal_write_wp(io_WP,1)



    te = time.time()
    
    
    ich = 0  #  asuppr
    
    
    # Temporal loop
    for jt in range(1,parameters.nt+1): # +1 to reach jt = nt
    

        api_wsc.print_file("Iteration n {}".format(jt), 1111)

        t2 = time.time()
        
        
        # Peut etre a supprimer
        if jt == 1 :
            api_wsc.forcedremesh = -1
        else :
            api_wsc.forcedremesh = 0
        
     


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
            #if ((jk == 1) and (jt != 1) and parameters.deformmesh == -1):
            if ((jk == 1) and parameters.deformmesh == -1):
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
                
                if not(grossier):
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
                
                # if i_ordi == 1 :
                #     api_wsc.aff = 1
                # else:
                #     api_wsc.aff = 0
                    
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
                    
                    
            if (jk == 1) :
                
                api_wsc.api_parareal_write_force(io_F,2)
                
     
        # Time-stepping
        api_wsc.api_time_stepping()
   
        # Writting the time info in the command window
        api_wsc.api_writting_time_info(jt)
        
        
        


        api_wsc.api_parareal_write_wp(io_WP,1)
        

    tf = time.time()
    
    
    
    

      
    if filestate_out != "":
        
        api_wsc.api_write_state(filestate_out,jt)
        
        #api_wsc.api_write_state_with_interpolation(filestate_out,jt)

    # Closing of the output files
    api_wsc.api_closing()
    
    # Deallocating
    api_wsc.deallocating_temporal_loop()
    t2 = time.time()
    # 
    
    
    print("")
    print("Temps")
    print(tb-ta,tc-tb,td-tc,te-td,tf-te)
    print("")
    
    
    print("Temps total : {}".format(t2-t1))

 
def change_size(dx1):
    
    
    K = dx1/0.4
    
    
    parameters.dx1 = K*0.4
    parameters.dx2 = K*0.03
    parameters.dx3 = K*0.03
    
    parameters.dx2domain = K*0.05
    
    # parameters.hrefx = K*0.6
    # parameters.hrefy = K*0.6
    # 
    # parameters.hrefxfs = K*0.05
    # parameters.hrefyfs = K*0.05
    # parameters.dsi = K*0.002
    
 
    
    
    
    
    
    
## Main


L_N_ordis = [1, 2, 4 , 10, 20, 40, 100, 200, 400]

L_N_ordis = [10]

L_dx1 = [0.4]

api_wsc.api_open_file_debug()
    
aff = False


for dx1 in L_dx1 :
    for N_ordis in L_N_ordis :
    
    
        api_wsc.print_file("", 1111)
        api_wsc.print_file("----------------------------------------------------------------", 1111)
        api_wsc.print_file("Nombre d'ordinateurs {}".format(N_ordis), 1111)
        api_wsc.print_file("dx1 {}".format(dx1), 1111)
        api_wsc.print_file("----------------------------------------------------------------", 1111)
        api_wsc.print_file("", 1111)
        
        
        
        
        N_iterations = 8
        #N_ordis = 1
        
        io_WP = 7988
        io_F = 4524
        

        
        api_wsc.api_parareal_init('ws.in','test.geom', N_iterations, N_ordis)
        

        parameters.affich = 0
        
        parameters.bool_coarse_init = 0
        
        pas_ordi = parameters.nt // N_ordis # A affiner
        
        
        
        ta = time.time()
        
        t_g = 0
        t_f = 0
        
        
        # Initialisation
        
        for i_ordi in range(N_ordis):
            
        
            api_wsc.print_file("", 1111)
            api_wsc.print_file("----------------------------------------------------------------", 1111)
            api_wsc.print_file("Initialisation, Ordinateur n {}".format( i_ordi + 1), 1111)
            api_wsc.print_file("----------------------------------------------------------------", 1111)
            api_wsc.print_file("", 1111)
        
        
            filename = "output_WP/WP_G_{}_{}.dat".format(0, i_ordi)
            api_wsc.api_parareal_open_file_wp(filename, io_WP)
            
            filename = "output_WP/force_G_{}_{}.dat".format(0, i_ordi)
            api_wsc.api_parareal_open_file_force(filename, io_F)
                
                
            if i_ordi == 0 :
                filestate_in = ""
            else:
                filestate_in = "output_WP/lambda_{}_{}.dat".format(0, i_ordi)
            
            filestate_out = "output_WP/lambda_{}_{}.dat".format(0, i_ordi+1)
        
        
            calcul(True, pas_ordi, io_WP, io_F, filestate_in, filestate_out) 
        

            api_wsc.api_parareal_save_g(i_ordi)
            

            api_wsc.api_parareal_close_file(io_WP)
            api_wsc.api_parareal_close_file(io_F)
            
            
            

        
        api_wsc.api_parareal_next_g(N_ordis) 
            
        api_wsc.print_file("", 1111)
        api_wsc.print_file("----------------------------------------------------------------", 1111)
        api_wsc.print_file("Fin de l initialisation", 1111)
        api_wsc.print_file("----------------------------------------------------------------", 1111)
        api_wsc.print_file("", 1111)
        

        tb = time.time()
        
        api_wsc.print_file("Chrono init : {}".format(tb-ta), 1111)
        

        # Iterations
        i_iter = 0
        converge = False
            
        while not converge :
            
            
            # Calcul fin (en parallele)
            
            
            ta = time.time()
        
            for i_ordi in range(N_ordis) :
                
        
        
                if i_ordi == 8 :
                    aff = True
                else :
                    aff = False
                    
                    
                api_wsc.print_file("", 1111)
                api_wsc.print_file("----------------------------------------------------------------", 1111)
                api_wsc.print_file("Fin : Iteration n {}, Ordinateur n {}".format(i_iter+1, i_ordi + 1), 1111)
                api_wsc.print_file("----------------------------------------------------------------", 1111)
                api_wsc.print_file("", 1111)
        
                    
                filename = "output_WP/WP_F_{}_{}.dat".format(i_iter, i_ordi)
                api_wsc.api_parareal_open_file_wp(filename, io_WP)
                
                filename = "output_WP/force_F_{}_{}.dat".format(i_iter, i_ordi)
                api_wsc.api_parareal_open_file_force(filename, io_F)
                    
        
        
                if i_ordi == 0 :
                    filestate_in = ""
                else:
                    filestate_in = "output_WP/lambda_{}_{}.dat".format(i_iter, i_ordi)
            
                filestate_out = ""
        
                calcul(False, pas_ordi, io_WP, io_F, filestate_in, filestate_out) 
        
                api_wsc.api_parareal_save_f(i_ordi)
                
                api_wsc.api_parareal_close_file(io_WP)
                api_wsc.api_parareal_close_file(io_F)
                
            tb = time.time()


            t_f += tb-ta
                
            # Calcul grossier (en sequentiel)
            
            
            ta = time.time()
            
        
            converge = True
            
            for i_ordi in range(N_ordis):
                
                
                api_wsc.print_file("", 1111)
                api_wsc.print_file("----------------------------------------------------------------", 1111)
                api_wsc.print_file("Grossier : Iteration n {}, Ordinateur n {}".format(i_iter+1, i_ordi + 1), 1111)
                api_wsc.print_file("----------------------------------------------------------------", 1111)
                api_wsc.print_file("", 1111)
                
                filename = "output_WP/WP_G_{}_{}.dat".format(i_iter+1, i_ordi)
                api_wsc.api_parareal_open_file_wp(filename, io_WP)
                
                filename = "output_WP/force_G_{}_{}.dat".format(i_iter+1, i_ordi)
                api_wsc.api_parareal_open_file_force(filename, io_F)
                
                
                
                if i_ordi == 0 :
                    filestate_in = ""
                else:
                    filestate_in = "output_WP/lambda_{}_{}.dat".format(i_iter+1, i_ordi)
                
                filestate_out = ""
        
                t111 = time.time()
                
                calcul(True, pas_ordi, io_WP, io_F, filestate_in, filestate_out) 
        
                api_wsc.api_parareal_save_g(i_ordi)
        
                filestate_out = "output_WP/lambda_{}_{}.dat".format(i_iter + 1, i_ordi+1)
                
                t112 = time.time()

                erreur = api_wsc.api_parareal_calcul_lambda(i_iter, i_ordi, 1, filestate_out) #jt = 1 (arbitraire, a verifier)
                
                
                t113 = time.time()
                
                
                print("Temps lambda : ", t112-t111, t113-t112)
                
                if erreur > 0.001 :
                    
                    converge = False
                    
                    
                    
                api_wsc.print_file("Erreur : {}".format(erreur), 1111)
                
                api_wsc.api_parareal_close_file(io_WP)
                api_wsc.api_parareal_close_file(io_F)
                
                
            api_wsc.api_parareal_next_g(N_ordis) 
                
        
        
            tb = time.time()
                
            
            t_g += tb-ta
            
            
            i_iter += 1
        
            api_wsc.print_file("Chrono : {}, {}".format(t_f, t_g), 1111)
        
        
        
        api_wsc.api_parareal_close()
            
    
        
api_wsc.api_close_file_debug()
        
