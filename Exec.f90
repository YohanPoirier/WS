module Exec_mod
use Houle_RF
use Constantes
use FonctionsCommunes
use iso_c_binding

implicit none

interface read_param
    module procedure read_param_i, read_param_rp, read_param_l, read_param_ld, read_param_id, read_param_rpd, read_param_c
end interface

interface read_param_only
    module procedure read_param_i_only, read_param_rp_only, read_param_l_only, read_param_ld_only, read_param_id_only, read_param_rpd_only, read_param_c_only
end interface

contains

subroutine Execution(fileparam,filegeom,InputData,get_State)
        
    character (len=50),intent(in)       :: fileparam,filegeom   ! Input files.
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(inout) :: InputData            ! Input data.
    logical,intent(in)                  :: get_State            ! True if the state input file is present, false otherwise.
    
    character(len=8)                    :: date_exec            ! Current date.
    character(len=10)                   :: time_exec            ! Current time.
    integer                             :: ierror               ! Error flag.

    ! This subroutine reads the input files and creates the structure InputData.
    
    call date_and_time(date_exec,time_exec) 
    
    print*,''
    print*,'------------------------- Prog. Weak Scatterer --------------------'
    print*,''
    write(*,*) 'Project Name :'
    write(*,*) 'User         :'
    write(*,*) 'Date         :',date_exec(7:8),'/',date_exec(5:6),'/',date_exec(1:4),' ',&
    &                           time_exec(1:2),':',time_exec(3:4),':',time_exec(5:6)
    write(*,*) 'exec.   dir. :'
    write(*,*) 'input dir.   :'
    write(*,*) 'output dir.  :'
    write(*,*) ''
    
    ! Reading NBodies
    call read_NBodies(filegeom)
    
    ! Reading nWP and the number of waves in case of Airy waves.
    call read_nWP(fileparam)
    
    ! Creating and initializing of InputData
    call New_InputData(InputData,Nbodies,nWP)
    call Initialization_InputData(InputData)
    
    ! Reading of the input file *.in.
    call read_input(fileparam,InputData,ierror)
            
    ! Reading of the input file *.geom.
    call read_geom(filegeom,InputData)
    
    ! Extra parameters.
    call Define_Extra_Data()
    
    ! Opening the output file Intersection_curves.dat
    call Opening_Intersection(get_State)
    
    ! Error
    ierror = 0
    if(ierror/=0)then
        ierror = 100
        goto 9999
    endif
    ierr = 0
    
    9999 continue
    if(ierror/=0)then
        write(*,99),ierror
    endif
99  format('** error #',i3,' : Main')
    
end subroutine Execution

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! READ INPUT
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_input(fileparam,InputData,ierror)
    
    character(len=50),intent(in)        :: fileparam                ! Input file *.in.
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(inout) :: InputData                ! Input data.
    integer,intent(inout)               :: ierror                   ! Error flag.
    
    character(len=120)                  :: line                     ! Line.
    real(rp)                            :: file_version             !
    real(rp),dimension(3)               :: tab_dx                   ! Structure to get the data back.
    real(rp),dimension(2)               :: tab_dbl, tab_href, tab_T ! Structures to get the data back.
    real(rp),dimension(2)               :: tab_hrefFS               ! Structure to get the data back.
    integer                             :: nc                       ! Loop parameter.
    
    ! This subroutine reads the input file *.in.
    
    ierror = 0
    
    ! Opening
    open(unit=ioparam,file=fileparam,iostat=ierror)
    if(ierror/=0) goto 9999
    
    read(ioparam,'(a)',end=1000,err=9999) line
    ! PROJECT
    call read_param(ioparam,file_version)    
    call read_param(ioparam,filename)
    
    read(ioparam,'(a)',end=1000,err=9999) line
    ! INCIDENT WAVE
    call read_param(ioparam,Htype)
    if(Htype.eq.2)then
        print*,"Water depth for the RF wave (input parameter: prof) must match with the water depth given when running waverf.exe."
    end if
    call read_param(ioparam,profondeur)    
    read(ioparam,'(a)',end=1000,err=9999) line
    call read_param(ioparam,NHoule)
    call read_param(ioparam,DirectionAiry,NHoule)
    call read_param(ioparam,w,NHoule)
    call read_param(ioparam,Aphi,NHoule)
    call read_param(ioparam,PhiWave,NHoule)
    print*,"The modification of the phase of the incident wave is only for the Airy waves."
    read(ioparam,'(a)',end=1000,err=9999) line
    call read_param(ioparam,filecoeff)
    call read_param(ioparam,DirectionRF)
    call read_param(ioparam,PhiWaveRF)
    
    ! Extraction des coefficients de houle de Rienecker et Fenton.
    if (Htype.eq.2) then
        call extract_WaveRF(HouleRF)
    end if
    
    read(ioparam,'(a)',end=1000,err=9999) line
    ! MESH DEFINITION
    call read_param(ioparam,cuve_ferme)
    call read_param(ioparam,symmetry)
    call read_param(ioparam,bottom_sym)
    call read_param(ioparam,LAbs)
    LAbs = max(LAbs,Epsilon2)
    call read_param(ioparam,Mesh_type)
    if(Mesh_type.eq.2 .and. not(cuve_ferme))then
        print*,"Open domain not available with Mesh_type = 2."
    end if
    if(Mesh_type.eq.3)then
        print*,"The weakscatterer code must be called in using Python to use Meshmagick."
    end if
    if(Mesh_type.ne.1 .or. Mesh_type.ne.2 .or. Mesh_type.ne.3)then
        print*,"This value of Mesh_type is not available."
    end if
    read(ioparam,'(a)',end=1000,err=9999) line    
    call read_param(ioparam,decalage)
    call read_param(ioparam,Nth)
    call read_param(ioparam,Nr)
    read(ioparam,'(a)',end=1000,err=9999) line
    call read_param(ioparam,tab_dx,2)
    call read_param(ioparam,d_bl)
    call read_param(ioparam,tab_href,2)
    call read_param(ioparam,tab_hrefFS,2)
    call read_param(ioparam,dsi)
    call read_param(ioparam,NormOverDref)
    call read_param(ioparam,idref)
    if(idref.eq.0 .or. idref.eq.1)then
        print*,"Size of the panels found by a boundary problem solved, adapted in case of several floaters."
    else if(idref.eq.2)then
        print*,"Size of the panels found by a linear approximation, adapted in case of a single floater."
        if(Nbodies.gt.1)then
            print*,"idref = 2 would be better in case of several floaters"
            pause
        end if
    end if
    call read_param(ioparam,RemeshFS)
    call read_param(ioparam,ForcedRemeshFS)
    if(ForcedRemeshFS .and. not(RemeshFS))then
        print*,"If ForcedRemeshFS = T, RemeshFS must also be T."
        pause
    end if
        
    dx1 = tab_dx(1) ; dx3 = tab_dx(2)
    dx4 = dx1
    hrefx = tab_href(1) ; hrefy = tab_href(2)
    hrefxFS = tab_hrefFS(1) ; hrefyFS = tab_hrefFS(2)
        
    ! NUMERICAL PARAMETERS
    read(ioparam,'(a)',end=1000,err=9999) line
    call read_param(ioparam,CCI)
    call read_param(ioparam,is_BS)
    call read_param(ioparam,SplineOrder)
    if(SplineOrder.eq.0)then
        print*,"The interpolation algorithm uses thin plate b-splines."
    else if(SplineOrder.eq.3)then
        print*,"The interpolation algorithm uses polynomial b-splines of order 3."
    else if(SplineOrder.eq.4)then
        print*,"The interpolation algorithm uses polynomial b-splines of order 4."
    else
        print*,"This order of b-splines is not available. Please change it. SplineOrder = ",SplineOrder
        call exit()
    end if
    call read_param(ioparam,DPNormals)
    if(DPNormals.lt.0._RP .or. DPNormals.gt.1._RP)then
        print*,"DPNormals must be between 0 and 1."
        call exit()
	end if
    call read_param(ioparam,lineaireFS)
    call read_param(ioparam,lineaireBody)
    call read_param(ioparam,RK_fige)
    call read_param(ioparam,is_latching)
    call read_param(ioparam,t0)
    call read_param(ioparam,dt)
    call read_param(ioparam,nt)
    call read_param(ioparam,tab_T,2)
    call read_param(ioparam,Tout)
    call read_param(ioparam,nout)
    call read_param(ioparam,Solv)
    call read_param(ioparam,TolGMRES)
    call read_param(ioparam,nliss_input)
    nliss = nliss_input
    call read_param(ioparam,hydrostatique)
    call read_param(ioparam,NThreads)
    
#ifdef _OPENMP
    if(NThreads.le.0)then
        print*,"NThreads must be higher or equal to 1."
        call exit()
    else if(NThreads.eq.1)then
        print*,"One thread will be used: sequential computation with OpenMP."
    else
        print*,NThreads," threads will be used: parallel computation with OpenMP."
    end if
#else
    print*,"One thread will be used: sequential computation."
#endif
    
    
    T1 = tab_T(1) ; T2 = tab_T(2)    
    
    ! WAVE PROBES
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line ! nWP
    do nc = 1,nWP
        call read_param(ioparam,InputData%PositionWP(1:3,nc),3)
    end do
    
    ! PHYSICAL CONSTANTS
    read(ioparam,'(a)',end=1000,err=9999) line
    call read_param(ioparam,ro)
    call read_param(ioparam,g)
    
    ! OUTPUT PARAMETERS 
    read(ioparam,'(a)',end=1000,err=9999) line
    call read_param(ioparam,iprint)
    call read_param(ioparam,idebug)
    call read_param(ioparam,iwhoule)
    call read_param(ioparam,iwenergy)
    call read_param(ioparam,iwmbody)
    call read_param(ioparam,iwmonitor)
    call read_param(ioparam,icheck)
    call read_param(ioparam,iwevol)
    call read_param(ioparam,iwevolRemesh)
    call read_param(ioparam,iwCartGrid)

    1000 continue

    ! Closing
    close(ioparam)
    
    9999 continue
    if(ierror/=0)then
        write(*,99) ierror
    endif
    99 format('** error #',i3,' : pb. read input.')

end subroutine read_input

subroutine read_NBodies(filegeom)
    
    character(len=50),intent(in)    :: filegeom ! Input file *.geom.
    
    character(len=50)               :: line     ! Line
    integer                         :: ierror   ! Error flag
    
    ! This subroutine reads NBodies in the input file *.geom.
    
    ierror = 0._rp
    
    ! Opening
    open(unit=iogeom,file=filegeom)
       
    read(iogeom,'(a)',end=1000,err=9999) line
    read(iogeom,'(a)',end=1000,err=9999) line ! is_body
    call read_param(iogeom,NBodies)
    
    1000 continue

    ! Closing
    close(iogeom)
    
    9999 continue
    if(ierror/=0)then
        write(*,99) ierror
    endif
    99 format('** error #',i3,' : pb. read Nbodies.')
    
end subroutine read_NBodies

subroutine read_nWP(fileparam)

    character(len=50),intent(in)    :: fileparam    ! Input file *.in.
    
    character(len=50)               :: line         ! Line
    integer                         :: ierror       ! Error flag

    ! This subroutine reads nWP in the input file *.ws.
    
    ierror = 0
    
    ! Opening
    open(unit=ioparam,file=fileparam,iostat=ierror)
    if(ierror/=0) goto 9999
    
    read(ioparam,'(a)',end=1000,err=9999) line
    ! PROJECT
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    
    read(ioparam,'(a)',end=1000,err=9999) line
    ! INCIDENT WAVE
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line 
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    
    read(ioparam,'(a)',end=1000,err=9999) line
    ! MESH DEFINITION
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line  
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
        
    ! NUMERICAL PARAMETERS
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
    read(ioparam,'(a)',end=1000,err=9999) line
        
    ! WAVE PROBES
    read(ioparam,'(a)',end=1000,err=9999) line
    call read_param(ioparam,nWP)
    
    1000 continue

    ! Closing
    close(ioparam)
    
9999 continue
    if(ierror/=0)then
        write(*,99) ierror
    endif
99 format('** error #',i3,' : pb. read nWP.')

end subroutine read_nWP

subroutine read_geom(filegeom,InputData)
    
    character(len=50),intent(in)        :: filegeom         ! Input file *.geom.
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(inout) :: InputData        ! Input data.
    
    character(len=50)                   :: line             ! Line.
    integer                             :: j,nc             ! Loop parameters.
    integer                             :: ierror           ! Error flag.
    
    ! This subroutine reads the input file *.geom.
    
    ierror = 0._rp
    
    ! Opening
    open(unit=iogeom,file=filegeom)
    
    ! General parameters
    read(iogeom,'(a)',end=1000,err=9999) line
    call read_param(iogeom,is_body)
    read(iogeom,'(a)',end=1000,err=9999) line !NBodies
    call read_param(iogeom,Int_body)
    call read_param(iogeom,iFixPoint)
    call read_param(iogeom,FixPointPos,3)
    call read_param(iogeom,Forward_Velocity)
    
    if(abs(Forward_Velocity).gt.Epsilon)then
        print*,"The forward speed is only along the x-axis."
        pause
    end if
    
    ! DOMAIN GEOMETRY
    read(iogeom,'(a)',end=1000,err=9999) line
    call read_param(iogeom,idtype)
    call read_param(iogeom,Ldom(3))
    if(dx1.gt.Ldom(3))then
        print*,"Panel size on the tank is bigger than the tank depth!"
        pause
    end if
    call read_param(iogeom,dx2Domain)
    read(iogeom,'(a)',end=1000,err=9999) line
    call read_param(iogeom,Ldom(1))
    call read_param(iogeom,Ldom(2))
    call read_param(iogeom,Ldom(4))
    read(iogeom,'(a)',end=1000,err=9999) line
    call read_param(iogeom,Ldom(5))
    
    ! Consistency between the length of the absorbing beach and the radius of the cylindrical domain.
    if((Mesh_type.eq.2 .or. Mesh_type.eq.3) .and. (idtype.eq.2 .or. idtype.eq.4) .and. LAbs.gt.Ldom(5))then
        print*,"The length of the absorbing beach is greated then the radius of the cylindrical domain."
        print*,"Please, modify that."
        pause
    end if
    
    ! Modification to do with a rectangular tank.
    if(idtype.ne.2 .and. idtype.ne.4)then
        print*,"create_dref_grid shoud be adapted to rectangular tank!" ! Ldom(5) only works with circular tank. Lgeom(2) only works with circular floater.
        print*,"Absorbing beach may be not well defined for a rectangular tank!"
        print*,"hrefxFS and hrefyFS are not defined for a rectangular tank (square free surface)!" ! Right now the free surface is defined as the first disc of a cylinder. If the FS is a square, hrefxFS and hrefyFS should be used with the first face of a cube: cf MeshModule/Mesh2d. Same thing when writing the Cartesian grid.
        pause
    end if
    
    ! is_body
    if(is_body .and. NBodies.eq.0)then
        print*,"is_body = ",is_body," and NBodies = ",NBodies,". These parameters should match."
        call exit()
    else if(not(is_body) .and. NBodies.ge.1)then
        print*,"is_body = ",is_body," and NBodies = ",NBodies,". Consequently the bodies are above the free surface."
    end if
    
    ! Water depth
    if(abs(profondeur - Ldom(3)).gt.Epsilon)then
        print*,"Water depth = ",profondeur
        print*,"Ldom(3) = ",Ldom(3)
        print*,"Water depth and Ldom_3 are not equal. Be sure of your input parameters."
    end if
    
    ! BODY GEOMETRY
    do nc = 1,NBodies
        read(iogeom,'(a)',end=1000,err=9999) line
        call read_param(iogeom,InputData%free_body(nc))
        call read_param(iogeom,InputData%Nphisphere(nc))
        call read_param(iogeom,InputData%dx2(nc))
        call read_param(iogeom,InputData%is_Mesh_File(nc))
        call read_param(iogeom,InputData%Mesh_file(nc))
        call read_param(iogeom,InputData%igtype(nc))
        if(InputData%igtype(nc) == 1 .and. Symmetry)then ! Cube
            print*,""
            print*,"Cube is not available with symmetry."
            print*,""
            call exit()
        end if
        if(InputData%igtype(nc) == 5 .and. idref.eq.2)then
            print*,"Lgeom_2 needs to be filled in case of idref = 2 (linear approximation for the panel size)."
        end if
        call read_param(iogeom,InputData%Position(1:3,1,nc),3)
        call read_param(iogeom,InputData%Position(1:3,2,nc),3) ! Degrees
        
        ! Deg -> Rad
        InputData%Position(1:3,2,nc) = (PI/180._RP)*InputData%Position(1:3,2,nc)
        
        call read_param(iogeom,InputData%PositionG(1:3,nc),3) ! In the body frame WITHOUT ROTATION.
        
        ! Vector: Position of the center of the mesh - Center of gravity at t = 0 s. Useful to update the frames of the geometries the subroutine Remesh.
        InputData%GOt0(1:3,nc) = InputData%Position(1:3,1,nc) - InputData%PositionG(1:3,nc)
                
        call read_param(iogeom,InputData%Mass(nc))
        if(InputData%Mass(nc).lt.0._RP)then
            print*,"Mass of the floater ",nc," is its displacement."
        end if
        call read_param(iogeom,InputData%is_Inertia_File(nc))
        call read_param(iogeom,InputData%file_inertia(nc))
        call read_param(iogeom,InputData%is_IBodyA(nc))
        read(iogeom,'(a)',end=1000,err=9999) line
        call read_param(iogeom,InputData%Lgeom(1,nc))
        call read_param(iogeom,InputData%Lgeom(2,nc))
        call read_param(iogeom,InputData%Lgeom(3,nc))
        read(iogeom,'(a)',end=1000,err=9999) line
        if(InputData%igtype(nc)==1)then
            call read_param(iogeom,InputData%Lgeom(1,nc))
        else    
            read(iogeom,'(a)',end=1000,err=9999) line
        endif
        read(iogeom,'(a)',end=1000,err=9999) line
        call read_param(iogeom,InputData%file_axisym(nc))
        read(iogeom,'(a)',end=1000,err=9999) line
        call read_param(iogeom,InputData%Vcst(nc))
        call read_param(iogeom,InputData%ACorps(nc))
        call read_param(iogeom,InputData%wCorps(nc))
        call read_param(iogeom,InputData%PhiCorps(nc))
        read(iogeom,'(a)',end=1000,err=9999) line
        call read_param(iogeom,InputData%Raideur(nc))
        call read_param(iogeom,InputData%Lressort(nc))
        call read_param(iogeom,InputData%Pressort(1:3,nc),3)
        call read_param(iogeom,InputData%MuPTO(nc))
        call read_param(iogeom,InputData%Cd_Morison(nc))
        call read_param(iogeom,InputData%dll_dir(1:6,nc),6)
        
        ! Supplementaru parameters and verifications.
        
        ! Approximation
        InputData%MovingBody(nc) = InputData%free_body(nc) .or. (not(InputData%free_body(nc)).and.abs(InputData%ACorps(nc)).gt.Epsilon)
        if (not(lineaireFS) .and. lineaireBody .and. InputData%MovingBody(nc)) print*, "Attention, l'approximation choisie est SL non linéaire, corps linéaire et corps mobile"      
        
        ! LambdACorps
        if(InputData%WCorps(nc).gt.epsilon)then
            InputData%LambdACorps(nc) = 2._RP*Pi*g/(InputData%WCorps(nc)**2)
        endif
        
        ! Number of degrees of freedom.
        InputData%Ndll(nc) = 0
        do j = 1,6
            if(InputData%dll_dir(j,nc))then
                InputData%Ndll(nc) = InputData%Ndll(nc) + 1
            end if
        end do
                
        if(Symmetry)then
            if(InputData%dll_dir(2,nc) .or.InputData%dll_dir(4,nc).or.InputData%dll_dir(6,nc))then
                print*, 'Error : Symmetry and degrees of freedom do not match. Only surge, heave and pitch motions are available.'
                call exit()
            end if
        end if
        
        ! Forced motion.
        if(not(InputData%free_body(nc)))then
            if(InputData%Ndll(nc).gt.1)then
                print*,"If a forced motion is wanted, only one dof can be chosen. Please verify the dof of the body ",nc
                call exit()
            end if
            
            if(abs(InputData%ACorps(nc)).gt.Epsilon .and. abs(InputData%Vcst(nc)).gt.Epsilon)then
                print*,"Constant velocity and oscillating motion are both activated. Please choose only one kind of forced motion for body ",nc
                call exit()
            end if
            
            if(abs(InputData%ACorps(nc)).gt.Epsilon)then
                print*,"Harmonic motion of the body ",nc
                print*,"x(t) = Acorps*sin(Wcorps*t + Phicorps)"
            end if
            
        end if
        
        ! Inertia.
        if(not(InputData%is_Inertia_File(nc)))then
            print*,"The inertia for the body ",nc," is computed by the code."
        else
            if(Symmetry)then
                print*,"The mass matrix of the body ",nc," is given in input for the WHOLE body, not for the half one."
            else
                print*,"The mass matrix of the body ",nc," is given in input."
            end if
        end if
        
    end do
    
    read(iogeom,'(a)',end=1000,err=9999) line
    
    ! -- SUPPLEMENTARY PARAMETERS -------------------------------------
    
    !do nc = 1,NBodies
    !    if (not(InputData%MovingBody(nc))) then
    !        lineaireBody = .true.
    !        print*, "Corps fixe --> forcage de l'approximation linéaire sur le corps (pas de déplacement du maillage)"
    !    end if
    !end do
    
    FreeBodies = .false.
    do nc = 1,NBodies
        if(InputData%free_body(nc))then
            FreeBodies = .true. ! If at least one body is free to move, FreeBodies = .true.
        end if
    end do
    
    1000 continue
    
    9999 continue
    if(ierror/=0)then
        write(99,*) ierror
    endif
    
    ! Closing
    close(iogeom)

    99  format('** error #',i3,' : p. read file parameters.')
        
end subroutine read_geom

subroutine Define_Extra_Data()
    
    real(rp)    :: val  ! Real.
    integer     :: j    ! Loop parameters.
    
    ! This subroutine initializes extra data which are neither in *.in nor in *.geom.
    
    ! Approximation
    if (Htype.eq.0) then
        if (not(lineaireFS)) print*, "Houle incidente nulle --> forcage de l'approximation linéaire sur la surface libre"
        lineaireFS = .true.
    end if
    
    ! Valisations solveur et dérivées locales.
    t_spline = .false.
    t_solver = .false.
    
    ! Wave.
    if(Htype.eq.1)then
        if(NHoule.gt.10)then
            print*,"The maximum number of Airy waves is 10."
            call exit()
        else
            do j = 1,NHoule
                if(abs(Aphi(j)).lt.Epsilon)then
                    print*,"If the amplitude of the wave is null, please use Htype = 0."
                    call exit()
                end if
            end do
        end if
    end if
    
    ! Détermination de la longueur d'onde si l'on est en profondeur finie.
    if(Htype.eq.1)then
        val = profondeur*maxval(w(:))**2/(4._rp*Pi**2*g)
        infinite_depth = val > 0.8
        if(infinite_depth)then
            print*,"Infinite depth approximation."
        end if
        if(infinite_depth)then
            konde = w**2/g
        else
            call find_konde(konde,w,profondeur,NHoule)
        endif
        lambda(1:NHoule) = 2._rp*Pi/konde(1:NHoule)
    endif
    
    NMAX = NHMAX
    Nhoulemx = NHMAX
    
    if(idebug.ge.1)then
        iwmesh1 = 1
    endif
    
    ! Ordre of the splines.
    ! Nordre(0) = Thin plate.
    ! Nordre(1) = Not available.
    ! Nordre(2) = Not available.
    ! Nordre(3) = Polynomial ordre 3.
    ! Nordre(4) = Polynomial ordre 4.
    Nordre(0:4) = [4,0,0,7,11] 
    
    ! SUPPLEMENTARY PARAMETERS
    Ra = 0._rp
    
end subroutine

subroutine read_param_i(iunit,val)
    
    integer,intent(in)      :: iunit    ! Index of the file.
    integer,intent(inout)   :: val      ! Integer.
    
    character(len=120)      :: line     ! Line.
    integer                 :: ind      ! Index of the ":"
    
    ! This subroutine reads an integer written after a ":".
    
    ! Reading of the line.
    read(iunit,'(a)') line
    
    ! Index of ":".
    ind = scan(line,':',.true.)
    
    ! Reading of the integer after the ":".
    read(line(ind+1:len(line)),*) val
    
end subroutine read_param_i

subroutine read_param_i_only(iunit,val)
    
    integer,intent(in)      :: iunit    ! Index of the file.
    integer,intent(inout)   :: val      ! Integer.
    
    character(len=120)      :: line     ! Line.
        
    ! This subroutine reads an integer only (no ":").
    
    ! Reading of the line.
    read(iunit,'(a)') line
        
    ! Reading of the integer after the ":".
    read(line(1:len(line)),*) val
    
end subroutine read_param_i_only

subroutine read_param_rp(iunit,val)
    
    integer,intent(in)      :: iunit! Index of the file.
    real(rp),intent(inout)  :: val  ! Real.
    
    character(len=120)      :: line ! Line.
    integer                 :: ind  ! Index of the ":"
    
    ! This subroutine reads a real written after a ":".
    
    ! Reading of the line.
    read(iunit,'(a)') line    
    
    ! Index of ":".
    ind = scan(line,':',.true.)
    
    ! Reading of the real after the ":".
    read(line(ind+1:len(line)),*) val
    
end subroutine read_param_rp

subroutine read_param_rp_only(iunit,val)
    
    integer,intent(in)      :: iunit! Index of the file.
    real(rp),intent(inout)  :: val  ! Real.
    
    character(len=120)      :: line ! Line.
        
    ! This subroutine reads a real only (no ":").
    
    ! Reading of the line.
    read(iunit,'(a)') line    
    
    ! Reading of the real.
    read(line(1:len(line)),*) val
    
end subroutine read_param_rp_only

subroutine read_param_l(iunit,val)
    
    integer,intent(in)      :: iunit    ! Index of the file.
    logical,intent(inout)   :: val      ! Boolean.
    
    character(len=120)      :: line     ! Line.
    integer                 :: ind      ! Index of the ":"
    
    ! This subroutine reads a Boolean written after a ":".
    
    ! Reading of the line.
    read(iunit,'(a)') line
    
    ! Index of ":".
    ind = scan(line,':',.true.)
    
    ! Reading of the Boolean after the ":".
    read(line(ind+1:len(line)),*) val
    
end subroutine read_param_l   

subroutine read_param_l_only(iunit,val)
    
    integer,intent(in)      :: iunit    ! Index of the file.
    logical,intent(inout)   :: val      ! Boolean.
    
    character(len=120)      :: line     ! Line.
        
    ! This subroutine reads a Boolean only (no ":").
    
    ! Reading of the line.
    read(iunit,'(a)') line
    
    ! Reading of the Boolean after the ":".
    read(line(1:len(line)),*) val
    
end subroutine read_param_l_only

subroutine read_param_ld(iunit,val,ndim)
    
    integer,intent(in)                      :: iunit            ! Index of the file.
    integer,intent(in)                      :: ndim             ! Dimension of the vector of Booleans.
    logical,dimension(ndim),intent(inout)   :: val              ! Vector of Booleans.
    
    character(len=120)                      :: line             ! Line.
    integer                                 :: ind, ind1, ind2  ! Indexes of the ":", "[" and "]".
    character(len=20)                       :: strval           ! String.
    
    ! This subroutine reads a vector of Booleans written after a ":" and between "[" and "]".
    
    ! Reading of the line.
    read(iunit,'(a)') line
    
    ! Index of ":".
    ind = scan(line,':',.true.)
    
    ! Reading after ":".
    strval = line(ind+1:len(line))
    
    ! Index of "[".
    ind1 = scan(strval,'[')
    
    ! Index of "]".
    ind2 = scan(strval,']')
    
    ! Reading of the vector of Booleans between "[" and "]".
    strval = strval(ind1+1:ind2-1)
    read(strval,*) val(1:ndim)
    
end subroutine read_param_ld  

subroutine read_param_ld_only(iunit,val,ndim)
    
    integer,intent(in)                      :: iunit            ! Index of the file.
    integer,intent(in)                      :: ndim             ! Dimension of the vector of Booleans.
    logical,dimension(ndim),intent(inout)   :: val              ! Vector of Booleans.
    
    character(len=120)                      :: line             ! Line.
    character(len=20)                       :: strval           ! String.
    
    ! This subroutine reads a vector of Booleans only (no ":"; "[" or "]").
    
    ! Reading of the line.
    read(iunit,'(a)') line
    
    ! Reading of the vector of Booleans between "[" and "]".
    strval = line(1:len(line))
    read(strval,*) val(1:ndim)
    
end subroutine read_param_ld_only

subroutine read_param_rpd(iunit,val,ndim)
    
    integer,intent(in)                      :: iunit            ! Index of the file.
    integer,intent(in)                      :: ndim             ! Dimension of the vector of reals.
    real(rp),dimension(ndim),intent(inout)  :: val              ! Vector of reals.
    
    character(len=120)                      :: line             ! Line.
    integer                                 :: ind, ind1, ind2  ! Indexes of the ":", "[" and "]".
    character(len=50)                       :: strval           ! String.
    
    ! This subroutine reads a vector of reals written after a ":" and between "[" and "]".
    
    ! Reading of the line.
    read(iunit,'(a)') line
    
    ! Index of ":".
    ind = scan(line,':',.true.)
    
    ! Reading after ":".
    strval = line(ind+1:len(line))
    
    ! Index of "[".
    ind1 = scan(strval,'[') 
    
    ! Index of "]".
    ind2 = scan(strval,']')
    
    ! Reading of the vector of reals between "[" and "]".
    strval = strval(ind1+1:ind2-1)
    read(strval,*) val(1:ndim)
    
end subroutine read_param_rpd

subroutine read_param_rpd_only(iunit,val,ndim)
    
    integer,intent(in)                      :: iunit    ! Index of the file.
    integer,intent(in)                      :: ndim     ! Dimension of the vector of reals.
    real(rp),dimension(ndim),intent(inout)  :: val      ! Vector of reals.
    
    character(len=300)                      :: line     ! Line.
    character(len=300)                      :: strval   ! String.
        
    ! This subroutine reads a vector of reals only (no ":"; "[" or "]").
    
    ! Reading of the line.
    read(iunit,'(a)') line
    
    ! Reading of the vector of reals
    strval = line(1:len(line))
    read(strval,*) val(1:ndim)
        
end subroutine read_param_rpd_only

subroutine read_param_id(iunit,val,ndim)
    
    integer,intent(in)                      :: iunit            ! Index of the file.
    integer,intent(in)                      :: ndim             ! Dimension of the vector of integers.
    integer,dimension(ndim),intent(inout)   :: val              ! Vector of integers.
    
    character(len=300)                      :: line             ! Line.
    integer                                 :: ind, ind1, ind2  ! Indexes of the ":", "[" and "]".
    character(len=300)                      :: strval           ! String.
    
    ! This subroutine reads a vector of integers written after a ":" and between "[" and "]".
    
    ! Reading of the line.
    read(iunit,'(a)') line
    
    ! Index of ":".
    ind = scan(line,':',.true.) 
    
    ! Reading after ":".
    strval = line(ind+1:len(line))
    
    ! Index of "[".
    ind1 = scan(strval,'[')
    
    ! Index of "]".
    ind2 = scan(strval,']')
    
    ! Reading of the vector of integers between "[" and "]".
    strval = strval(ind1+1:ind2-1)
    read(strval,*) val(1:ndim)
    
end subroutine read_param_id

subroutine read_param_id_only(iunit,val,ndim)
    
    integer,intent(in)                      :: iunit            ! Index of the file.
    integer,intent(in)                      :: ndim             ! Dimension of the vector of integers.
    integer,dimension(ndim),intent(inout)   :: val              ! Vector of integers.
    
    character(len=300)                      :: line             ! Line.
    character(len=300)                      :: strval           ! String.
    
    ! This subroutine reads a vector of integers only (no ":"; "[" or "]").
    
    ! Reading of the line.
    read(iunit,'(a)') line
    
    ! Reading of the vector of integers.
    strval = line(1:len(line))
    read(strval,*) val(1:ndim)
    
end subroutine read_param_id_only

subroutine read_param_c(iunit,val)
    
    integer,intent(in)              :: iunit    ! Index of the file.
    character(len=50),intent(inout) :: val      ! String.
    
    character(len=120)              :: line     ! Line.
    integer                         :: ind      ! Index of the ":".
    
    ! This subroutine reads a string written after a ":".
    
    ! Reading of the line.
    read(iunit,'(a)') line
    
    ! Index of ":".
    ind = scan(line,':',.true.)
    
    ! Reading of the vstring after the ":".
    val = trim(adjustl(line(ind+1:len(line))))
    
end subroutine read_param_c

subroutine read_param_c_only(iunit,val)
    
    integer,intent(in)              :: iunit    ! Index of the file.
    character(len=50),intent(inout) :: val      ! String.
    
    character(len=120)              :: line     ! Line.
        
    ! This subroutine reads a string only (no ":").
    
    ! Reading of the line.
    read(iunit,'(a)') line
    
    ! Reading of the vstring after the ":".
    val = trim(adjustl(line(1:len(line))))
    
end subroutine read_param_c_only

subroutine extract_WaveRF(HouleRF)

    !f2py integer*1, dimension(1000)    :: HouleRF
    type(THouleRF), intent(out)         :: HouleRF      ! RF wave.
    
    logical                             :: ex           ! Parameter.
    integer                             :: j, k, ios    ! Loop parameters.
    
    ! This subroutine extracts the RF wave from its input file.
    
    ! Test de présence de fichier de coefficients de Houle.
    inquire(file=filecoeff,exist=ex) 
    
    if(ex) then 
        
        ! Opening
        open(unit=13, file=filecoeff, access='sequential', form='formatted', iostat = ios)
        if (ios/=0) stop "Erreur à l'ouverture du fichier de houle"
        
        ! Reading.
        read (unit=13, fmt='(F15.6, E17.7,3E15.6)') HouleRF%Landa, HouleRF%Hcc, HouleRF%k, HouleRF%T, HouleRF%C
        HouleRF%Landa = HouleRF%Landa*profondeur
        HouleRF%Hcc = HouleRF%Hcc*profondeur
        HouleRF%k = HouleRF%k/profondeur
        HouleRF%T = HouleRF%T/sqrt(g/profondeur)
        HouleRF%C = HouleRF%C*sqrt(g*profondeur)
        read (unit=13, fmt='(2E15.7, I)') HouleRF%Cp, HouleRF%Cc, HouleRF%N
        HouleRF%Cp = HouleRF%Cp*sqrt(g*profondeur)
        HouleRF%Cc = HouleRF%Cc*sqrt(g*profondeur)
        do j=1,HouleRF%N+1
            read(unit=13, fmt='(I,2E15.7)') k, HouleRF%An(k), HouleRF%Bn(k)
            HouleRF%An(k) = HouleRF%An(k)*profondeur
            HouleRF%Bn(k) = HouleRF%Bn(k)*profondeur*sqrt(g*profondeur)
        end do
        HouleRF%Bn(0) = HouleRF%Bn(0)/profondeur
        HouleRF%d = profondeur
        
        ! Closing.
        close (unit=13)
        print*, 'Houle Rienecker and Fenton'
        print*, 'Period = ', HouleRF%T
        print*, 'Wave length = ', HouleRF%Landa
        print*, 'Amplitude = ', HouleRF%Hcc*0.5_RP
        print*, 'Depth = ', profondeur
    else
        print*, 'Pas de fichier de houle RF correspondant'
    end if
    
    ! HouleRF%Index is fixed equal to 1.
    
end subroutine Extract_WaveRF

subroutine read_Inertia(file_inertia,Mesh,NumBody,ierror)
    
    character(len=50),intent(in)    :: file_inertia ! Inertia file.
    !f2py integer*1, dimension(1000):: Mesh
    type(TMaillage),intent(inout)   :: Mesh         ! Mesh.
    Integer,intent(in)              :: NumBody      ! Body number.
    integer,intent(inout)           :: ierror       ! Error flag.
        
    integer                         :: i            ! Loop parameters.
    
    ! This subroutine reads the total mass matrix in the inertia files.
    
    ierror = 0
    
    ! Opening.
    open(unit=ioinertia,file=file_inertia,iostat=ierror)
    
    9999 continue
    if(ierror/=0)then
        close(ioinertia)
        print*,"Problem in reading the inertia file of the body: ",NumBody
        return
    endif
    
    ! Filling.
    do i=1,6
        read(ioinertia,*) Mesh%Body(NumBody)%IBody(i,1:6)
    end do
    
    ! Closing.
    close(ioinertia)
    
end subroutine read_Inertia

subroutine Opening_Intersection(get_State)
    
    logical,intent(in)  :: get_State    ! True if the state input file is present, false otherwise.
    
    integer             :: ios          ! Parameter.
    
    ! This subroutine opens the output file Intersection_curves.dat.
    
    ! Opening
    if(not(get_State))then
        open(unit=ioIntersection,file='Intersection_curves.dat')
        write(ioIntersection,*) 'Title = "Intersection curves between the floaters and the free surface"'
        write(ioIntersection,*) 'VARIABLES = "X","Y","Z","Nedge","E1","E2","E3","E4","E5"'
    else
        open(unit=ioIntersection,file='Intersection_curves.dat',Access = 'append',Status='old', iostat=ios)
        if (ios/=0) stop "Error at the opening of the intersection output file."
    end if
    
end subroutine Opening_Intersection


end module Exec_mod
