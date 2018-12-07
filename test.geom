# GENERAL PARAMETERS ------------------------------------------------------------------------------
    Body is present                             	:   is_body         :   T
    Number of bodies                             	:   NBodies         :   1
    Index of body                               	:   Int_body        :   2
    Resolution point different from COG                 :   iFixPoint       :   F
    Position of the resolution point                    :   FixPointPos     :   [0. 0. 0.]
    Body Forward Velocity				:   Forward_Velocity:	0.
# DOMAIN GEOMETRY ---------------------------------------------------------------------------------
    Type of the domain (1: , 2:, 3: , 4:)               :   idtype          :   2
    Depth (m)                                           :   Ldom_3          :   1.
    Size discretization for the domain (Mesh_type = 2)  :   dx2Domain       :   0.05
# Rectangular (if idtype = 1 or 3)    
    Length (m)                                          :   Ldom_1          :   22.2
    Large (m)                                           :   Ldom_2          :   22.2
    Large (m)                                           :   Ldom_4          :   0.
# Cylinder (if idtype = 2 or 4)
    Radius (m)                                          :   Ldom_5          :   3.0
# BODY GEOMETRY -----------------------------------------------------------------------------------
    Body free of move                           	:   free_body       :   T
    Number of discrete point on the sphere      	:   Nphisphere      :   25
    Size discretization for the floater	(Mesh_type=2)   :   dx2		    :   0.03
    Mesh file has to be read (T) or not (F)		:   is_Mesh_file    :   F
    File for inertia                                    :   Mesh_file       :   Sphere_WSC_clipped.dat
    Type of geometry (1:cube, 2:cylinder, 5:axisym)     :   igtype          :   2
    Initial position (m)                                :   Position_1      :   [0. 0. 0.0]
    Initial angular position (rad)                      :   Position_2      :   [0. 0. 0.]
    Cog				                         :   PositionG       :   [0. 0. 0.0]
    Mass of the system (kg)                             :   mass_body       :   64.4
    Inertia file has to be read (1) or not (0)		:   is_Inertia_File :   T
    File for inertia                                    :   file_inertia    :   Inertia_Body_1.dat
    Inertia on FixPoint                                 :   is_BodyA        :   F
# Dimension of the cylinder (if igeom = 2)   
    Length (m)                                          :   Lgeom_1         :   1.
    Radius (m)                                          :   Lgeom_2         :   0.2
    Radius (m)                                          :   Lgeom_3         :   0.
# Dimension of the cube (if igeom = 1)
    Length (m)                                          :   Lgeom_1         :   0.
# Definition of the AxiSym (if igeom = 5)
    File                                                :   file_axisym     :   wavestar_axisym.dat
# Definition of displacement (if free_body = F)    
    Constante veocity (m/s)                             :   Vcst            :   0.
    Amplitude (m)                               	:   Acorps          :   0.1
    Pulsation (rad/s)                           	:   wcorps          :   3.2
    Phase (rad)                           		:   Phicorps        :   0.
# PTO Model    
    Stiffness                                   	:   Raideur         :   0
    Spring length                               	:   LRessort        :   0
    Position                                    	:   PRessort        :   [0. 0 0]
    Damping coefficient                         	:   B_PTO           :   0
    Viscous force, Morison coefficient          	:   Cd_Morison      :   0
    Degree of Freedom                           	:   DOF             :   [F F F F F F]
#--------------------------------------------------------------------------------------------------