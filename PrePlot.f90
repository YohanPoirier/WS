module PrePlot
use Constantes
use Structuresdonnees
use FonctionsCommunes
use Incident_mod
use Houle_RF
use Airy
use GeomStruct
use Exec_mod
use Spline
implicit none

contains

! ------------------------------------------------------------
! PrePlots : ouverture et préparation des fichier de sortie
! ------------------------------------------------------------

subroutine PrePlots(InputData,NBody,get_State)
    
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(in)    :: InputData    ! Input data.
    integer,intent(in)                  :: NBody        ! Mesh%NBody.
    logical,intent(in)                  :: get_State    ! True if the state input file is present, false otherwise.
    
    integer                             :: ios,nc       ! Parameters.
    character(len=200)                  :: fileMBody    ! Temporary file name for the motion.
    character(len=200)                  :: fileLoads    ! Temporary file name for the loads.
    character(len=200)                  :: fileWP       ! Temporary file name for the wave elevation.
    
    ! This subroutine opens the output files.
    
    ! Output.dat
    if(not(get_State))then
        open(unit=33,file='output_'//filename, iostat=ios)
        if (ios/=0) stop "Erreur à l'ouverture du fichier de sortie principal"
        write(33,fmt='(50a)') 'Title = "Mesh de la cuve"'
        write(33,fmt='(150a)') 'VARIABLES = "X","Y","Z","VX","VY","VZ","Eta_Perturbation","Phi_Perturbation","GPhi_Perturbation_x","GPhi_Perturbation_y","GPhi_Perturbation_z","Eta_Incident","Phi_incident","GPhi0x","GPhi0y","GPhi0z","DPhiDt","DPhiDn","GEtapx","GEtapy","GEtapz","DEtaDt","DDPhiDnDt","DGradPhiSqDn","Pression"' !,"DPhiDt","DPhiDn","GEtax","GEtay","GEtaz","DEtaDt","DEta0Dt","DPhi0Dt","GEta0x","GEta0y","GEta0z","GPhi0x","GPhi0y","GPhi0z"'
    else
        open(unit=33,file='output_'//filename,Access = 'append',Status='old', iostat=ios)
        if (ios/=0) stop "Erreur à l'ouverture du fichier de sortie principal"
    end if
    
    ! Not used.
    if(iwmatv)then
        open(unit=102,file='mesh_displacement_'//filename,iostat=ios)
        if(ios/=0) stop "Erreur à l'ouverture du fichier de sortie principal"
        write(102,'(50a)') 'Title = "Mesh_displacement"'
        write(102,'(50a)') 'VARIABLES = "X","Y","Z","Eta0","GPhi0(1)","GPhi0(2)","GPhi0(3)","V_mesh(1)","V_mesh(2)","V_mesh(3)","N(1)","N(2)","N(3)","dPhi0dn(1)","dPhi0dn(2)","dPhi0dn(3)"'
    endif
    
    ! Not used.
    if(iwmesh)then
        open(unit=iomesh,file='mesh_'//filename,iostat=ios)
        if(ios/=0) stop "Erreur a l'ouverture du fichier de sortie Maillage"
        write(iomesh,fmt='(50a)') 'Title = "Remaillage de la cuve"'
        write(iomesh,fmt='(130a)') 'VARIABLES = "X","Y","Z","Nx","Ny","Nz","IX","IY","IZ","JX","JY","JZ","typeFrontiere","Npanneau","Damping","NVoisin","Ndouble","Vx","Vy","Vz","DL1","DL2","DL3","DL4","DL5","Ordre"' !,"Ndouble","double"'
    endif
    
    if(iwmetric)then
        open(unit=iometric,file='metric_'//filename,iostat=ios)
        if(ios/=0) stop "Erreur à l'ouverture du fichier de metric"
        write(iometric,fmt='(50a)') 'Title = "Metric associes a la qualite du maillage"'
        write(iometric,fmt='(130a)') 'VARIABLES = "X","Y","Z","FSIZE","FSHAPE"'
    endif

    if(iwenergy)then
        if(not(get_State))then
            open(unit=ioenergy,file='Energy_'//filename,iostat=ios)
            if(ios/=0) stop "Erreur à l'ouverture du fichier energy"
            write(ioenergy,fmt='(50a)') 'Title = "Energy"'
            write(ioenergy,fmt='(50a)') 'VARIABLES = "t","E_TOTAL","EK_SL","EP_SL","EK_SM","EK_Body","VGx",VGy","VGz","Gx","Gy","Gz"'
        else
            open(unit=ioenergy,file='Energy_'//filename,Access = 'append',Status='old',iostat=ios)
            if(ios/=0) stop "Erreur à l'ouverture du fichier energy"
        end if
    endif

    if(iwphi)then
        open(unit=iophi,file='PhiCorps_'//filename,iostat=ios)
        if(ios/=0) stop "Erreur à l'ouverture du fichier PhiCorps"
        write(iophi,fmt='(50a)') 'Title = "PhiCorps"'
        write(iophi,fmt='(50a)') 'VARIABLES = "X",Y","Z","PhiRef","PhiSol"'
    endif

    if(iwBVP)then
        open(unit=ioBVP,file='BVP_'//filename, iostat=ios)
        if (ios/=0) stop "Erreur à l'ouverture du fichier de gradient"
        write(ioBVP,fmt='(50a)') 'Title= "Sortie BVP"'
        write(ioBVP,fmt='(50a)') 'VARIABLES = "X","Y","Z","Phi_p","Phi_i","DPhiDn_p","DPhiDn_i"'
    endif
    
    if(iwTest_BVP)then
        open(unit=ioTest_BVP,file='Test_BVP_'//filename, iostat=ios)
        if (ios/=0) stop "Erreur à l'ouverture du fichier de gradient"
        write(ioTest_BVP,fmt='(50a)') 'Title= "Sortie BVP"'
        write(ioTest_BVP,fmt='(50a)') 'VARIABLES = "X","Y","Z","Phi_p","Phi_i","DPhiDn_p","DPhiDn_i","Error"'
    endif
    
    if(iwDerive)then
        open(unit=ioBVP,file='Derive_'//filename, iostat=ios)
        if (ios/=0) stop "Erreur à l'ouverture du fichier de calcul des dérivées"
        if(not(get_State))then
            write(ioBVP,fmt='(50a)') 'Title= "Sortie BVP"'
            write(ioBVP,fmt='(50a)') 'VARIABLES = "X","Y","Z","Phi_p","DPhiDt_p","Phi_i","DPhiDt_i"'
        end if
    endif

    if(iwpress3D)then
        open(unit=iopress3D,file='Pression3D_'//filename,iostat=ios)
        if(ios/=0) stop "Erreur à l'ouverture du fichier Pression3D."
        if(not(get_State))then
            write(iopress3D,fmt='(50a)') 'Title = "Pression3D"'
            write(iopress3D,fmt='(50a)') 'VARIABLES = "X","Y","Z","P","P_t1","P_t2","P_t3"'
        end if
    endif
    
    if(iwmbody)then
        do nc = Int_Body,NBody
            write(fileMBody,'("MBody_",i2,"_",50a)'),nc-1,trim(filename)
            if(not(get_State))then
                open(unit=iombody+nc,file=fileMBody,iostat=ios)
                if(ios/=0) stop "Erreur lors de l'ouverture du fichier MBody"
                write(iombody+nc,fmt='(50a)') 'Title = "Mouvement du Corps"'
                write(iombody+nc,fmt='(50a)') 'VARIABLES = "t","X","Y","Z","RX","RY","RZ","VX","VY","VZ","VRX","VRY","VRZ","AX","AY","AZ","ARX","ARY","ARZ"'
            else
                open(unit=iombody+nc,file=fileMBody,Access = 'append',Status='old')
            end if
        end do
    endif
    
    if(iwphit)then
        open(unit=iophit,file='BVPPhit_'//filename,iostat=ios)
        if(ios/=0) stop "Erreur à l'ouverture du fichier BVPPhit"
        if(not(get_State))then
            write(iophit,fmt='(50a)') 'Title = "Termes du probleme aux limites sur Phit"'
            write(iophit,fmt='(50a)') 'VARIABLES = "X","Y","Z","DDPhi0DnDt","DDPhiDnDt","DGradPhiSqDn","AccN"'
        end if
    endif

    if(iwpress .and. is_body)then
        if(not(InputData%free_body(1)) .and. FiniteDifference_ForceCalcul)then
            if(not(get_State))then
                print*,"InputData%free_body(1) was used. No MB in the part of the subroutine."
                open(unit=iopress,file='Pression_P_DF_'//filename,iostat=ios)
                if(ios/=0) stop "Erreur à l'ouverture du fichier Pression_P_DF"
                write(iopress,fmt='(50a)') 'Title = "Pression"'
                write(iopress,fmt='(50a)') 'VARIABLES = "t","FX","FY","FZ","FX0","FY0","FZ0","FX1","FY1","FZ1","FX2","FY2","FZ2","Mx","My","Mz"'
            else
                open(unit=iopress,file='Pression_P_DF_'//filename,Access = 'append',Status='old',iostat=ios)
                if(ios/=0) stop "Erreur à l'ouverture du fichier Pression_P_DF"
            end if
        else   
            if(not(get_State))then
                open(unit=iopress,file='Pression_P_BVP_'//filename,iostat=ios)
                if(ios/=0) stop "Erreur à l'ouverture du fichier Pression_P_BVP"
                write(iopress,fmt='(50a)') 'Title = "Pression"'
                write(iopress,fmt='(50a)') 'VARIABLES = "t","FX","FY","FZ","FX0","FY0","FZ0","FX1","FY1","FZ1","FX2","FY2","FZ2","Mx","My","Mz"'
            else
                open(unit=iopress,file='Pression_P_BVP_'//filename,Access = 'append',Status='old',iostat=ios)
                if(ios/=0) stop "Erreur à l'ouverture du fichier Pression_P_BVP"
            end if
        end if
    end if

    if(iwdphidt)then
        if (FiniteDifference_ForceCalcul) then
            open(unit=iodphidt,file='DPhiDt_DF_'//filename,iostat=ios)
            if(ios/=0) stop "Erreur à l'ouverture du fichier DPhiDt_DF"
            if(not(get_State))then
                write(iodphidt,fmt='(50a)') 'Title = "DPhiDt sur le Corps"'
                write(iodphidt,fmt='(50a)') 'VARIABLES = "X","Y","Z","Phi","DpPhiDt","DPhiDt"'
            end if
        else
            open(unit=iodphidt2,file='DPhiDt_BVP_'//filename,iostat=ios)
            if(ios/=0) stop "Erreur à l'ouverture du fichier DPhiDt_BVP"
            if(not(get_State))then
                write(iodphidt2,fmt='(50a)') 'Title = "DPhiDt sur le Corps"'
                write(iodphidt2,fmt='(50a)') 'VARIABLES = "X","Y","Z","Phi","DpPhiDt","DPhiDt"'
            end if
        end if
    end if

    if(iwts)then
        open(unit=iots, file='Ts.dat')
        if(not(get_State))then
            write(iots,fmt='(50a)') 'Title = "Ts"'
            write(iots,fmt='(50a)') 'VARIABLES = "t","Ts(1)","CTdSM"'
        end if
        close(iots)
    endif

    if(iwsize)then
        open(unit=iosize,file='size_mesh.log')
        if(not(get_State))then
            write(iosize,*) 't fshape fsize_min fsize_max'
        end if
    endif

    if(iwmonitor)then
        if(not(get_State))then
            open(unit=iomonitor,file='monitoring.log')
            write(iomonitor,'(50a)') 'Title = "Monitoring"'
            write(iomonitor,'(100a)') 'VARIABLES = "t","Body mesh quality","FS mesh quality","Crossing the free surface","Forced remeshing","FS remeshing asked","Remeshing asked","Remeshing worked"'
        else
            open(unit=iomonitor,file='monitoring.log',Access = 'append',Status='old')
        end if
    end if

    if(iwmass)then
        if(not(get_State))then
            open(unit=iomass,file='mass_monitoring.dat')
            write(iomass,'(50a)') 'Title = "Mass monitoring"'
            write(iomass,'(50a)') 'VARIABLES = "t","mass"'
        else
            open(unit=iomass,file='mass_monitoring.dat',Access = 'append',Status='old')
        end if
    end if

    if(iwwave)then
        do nc = 1,nWP
            write(fileWP,'("Wave_elevation_",i2,".dat")'),nc
            if(not(get_State))then
                open(unit=iowave+nc,file=fileWP,iostat=ios)
                if(ios/=0) stop "Erreur lors de l'ouverture du fichier Wave_elevation"
                write(iowave+nc,'(50a)') 'Title = "Wave elevation"'
                write(iowave+nc,'(50a)') 'VARIABLES = "t","Eta_Incident","Eta_Perturbation","Eta"'
            else
                open(unit=iowave+nc,file=fileWP,Access = 'append',Status='old',iostat=ios)
                if(ios/=0) stop "Erreur lors de l'ouverture du fichier Wave_elevation"
            end if
        end do
    end if
        
    if(iwLoads)then
        do nc = Int_Body,NBody
            write(fileLoads,'("Forces_",i2,"_",50a,".dat")'),nc-1,trim(filename)
            if(not(get_State))then
                open(unit=ioLoads+nc,file=fileLoads,iostat=ios)
                if(ios/=0) stop "Erreur lors de l'ouverture du fichier Forces"
                write(ioLoads+nc,fmt='(50a)') 'Title = "Loads"'
                write(ioLoads+nc,'(100a)') 'VARIABLES = "t","FWX","FWY","FWZ","MWX","MWY","MWZ","FhsX","FhsY","FhsZ","MhsX","MhsY","MhsZ","FhydroX","FhydroY","FhydroZ","MhydroX","MhydroY","MhydroZ","FPTOX","FPTOY","FPTOZ","MPTOX","MPTOY","MPTOZ","FMorisonX","FMorisonY","FMorisonZ","MMorisonX","MMorisonY","MMorisonZ","FTotalX","FTotalY","FTotalZ","MTotalX","MTotalY","MTotalZ"'
            else
                open(unit=ioLoads+nc,file=fileLoads,Access = 'append',Status='old',iostat=ios)
                if(ios/=0) stop "Erreur lors de l'ouverture du fichier Forces"
            end if
        end do
    end if
    
end subroutine PrePlots

! -----------------------------------------------------------------
! close_output : fermeture des fichiers de sortie
! -----------------------------------------------------------------

subroutine close_output(time,NBody)
    
    real(rp),intent(in) :: time     ! Current time
    integer,intent(in)  :: NBody    ! Mesh%NBody
    
    integer             :: ios,nc   ! Parameters
    
    ! This subroutine closes all the output files.
    
    close(33)
    close(101)
    close(11)
    close(58)

    if(iwmatv)         close(102)
    if(iwmesh)         close(iomesh)
    if(iwmetric)       close(iometric)
    if(iwpress)        close(iopress)
    if(iwenergy)       close(ioenergy)
    if(iwphi)          close(iophi)
    if(iwBVP)          close(ioBVP)
    if(iwTest_BVP)     close(ioTest_BVP)
    if(iwDerive)       close(ioDerive)
    if(iwpress3D)      close(iopress3D)
    if(iwphit)         close(iophit)
    if(iwsize)         close(iosize)
    if(iwmonitor)      close(iomonitor)
    if(iwmass)         close(iomass)
    if(iwIntersection) close(ioIntersection)
        
    
    if(iwwave)then
        do nc = 1,nWP
            close(iowave+nc)
        end do
    end if
    
    
    
    if(iwmbody)then
        do nc = Int_Body,NBody
            close(iombody+nc)
        end do
    end if
    
    
    
    if(iwLoads)then
        do nc = Int_Body,NBody
            close(ioLoads+nc)
        end do
    end if
    
    
    
    open(unit=17,file='Temps_'//filename, iostat=ios)
    if (ios/=0) stop "Erreur à l'ouverture du fichier Temps"
    write(17,fmt='(a,e,a)')'Temps total du calcul : ', time, ' secondes'
    close(17)

end subroutine close_output

subroutine Plots(t, Mesh, Ecoulement,InputData)
    
    real(rp)                            :: t            ! Current time.
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage)                     :: Mesh         ! Mesh.
    !f2py integer*1, dimension(1000)    :: Ecoulement
    type(TECoulement)                   :: Ecoulement   ! Flow parameters.
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(in)    :: InputData    ! Input data.
        
    ! This subroutine updates the output files.
    
    ! Mesh.
    if(iwhoule)   call PlotPropHoule(t, Mesh, Ecoulement)
    
    ! Motion, velocity and acceleration of the floaters.
    if(iwmbody)   call PlotMBody(t, Mesh)
    
    ! Wave elevations.
    if(iwwave)    call PlotWaveElevation(t,InputData,Mesh,Ecoulement)
    
end subroutine Plots

! ---------------------------------------------------------------------------
! PlotEcoul (non utilise)
! ---------------------------------------------------------------------------
subroutine PlotEcoul(Mesh, Ecoulement, t, filen)
! Paramètres :
!f2py integer*1, dimension(1000) :: Mesh
type(TMaillage), intent(in) :: Mesh
!f2py integer*1, dimension(1000) :: Ecoulement
type(TEcoulement), intent(in) :: Ecoulement
real(rp), intent(in) :: t
character(len=50), intent(in), optional :: filen
! Variables :
integer :: j, ios, iotemp
character(len=10) :: num
! Begin
if (present(filen)) then
    iotemp = 3
    open(unit=iotemp,file=filen, iostat=ios)
    if (ios/=0) stop "Erreur à l'ouverture du fichier de gradient"
    write(iotemp,fmt='(50a)') 'Title= "Sortie Dérivées SL"'
    write(iotemp,fmt='(50a)') 'VARIABLES = "X","Y","Z","Phi_p","DPhiDt_p","Phi_i","DPhiDt_i"'
    write(iotemp,fmt='(a,i,a,i,a)') 'Zone N =', Mesh%Nsys, ', E=', Mesh%Nfsys, ' , ET=TRIANGLE, F=FEPOINT'
else
    iotemp = ioBVP
    write( num, '( f0.4 )' ) t
    write(iotemp,fmt='(a,i,a,i,a)')  'Zone T = "'//trim(num)//'seconds", N =', Mesh%Nsys, ', E=', Mesh%Nfsys, ' , ET=TRIANGLE, F=FEPOINT, STRANDID = 1, SOLUTIONTIME = '//trim(num)
end if 
do j=1,Mesh%Nsys
    write(iotemp,'(7E)') Mesh%Tnoeud(j)%Pnoeud, Ecoulement%Phi(j)%perturbation, Ecoulement%DPhiDt(j)%perturbation, Ecoulement%Phi(j)%incident, Ecoulement%DPhiDt(j)%incident
end do
do j=1,Mesh%Nfsys
    write(iotemp,'(4I)') Mesh%Tfacette(j)%Tnoeud
end do
if (present(filen)) close(unit=iotemp)
! End
end subroutine PlotEcoul

subroutine PlotMaill(filename, Mesh, t, iout)
    
    character(len=50), intent(in)   :: filename         ! File name.
    !f2py integer*1, dimension(1000):: Mesh
    type(TMaillage), intent(in)     :: Mesh             ! Mesh.
    real(rp),intent(in), optional   :: t                ! Current time.
    integer, optional               :: iout             ! Output parameter.
    
    integer                         :: j, ios, iunit,k  ! Loop parameters.
    character(len=10)               :: num              ! Solution time.
    integer                         :: Nv               ! Number of 1st order neighbours.
        
    ! This subroutine writes the Tecplot file of the mesh.
    
    iunit = ioMaillage
    if(present(iout)) iunit = iout

    if(.not.present(t))then
        open(unit=ioMaillage,file=filename, iostat=ios)
        if (ios/=0) stop "PlotMaill: error at the opening."
        write(ioMaillage,fmt='(50a)') 'Title = "Remaillage de la cuve"'
        write(ioMaillage,fmt='(200a)') 'VARIABLES = "X","Y","Z","Nx","Ny","Nz","IX","IY","IZ","JX","JY","JZ","typeFrontiere","Npanneau","Damping","NVoisin 1","NVoisin 1,2","Ndouble","DL1","DL2","DL3","DL4","DL5","PC","MOB","Nfacette","Aire"' !,"Ndouble","double"'
        write(ioMaillage,fmt='(a,i,a,i,a)') 'Zone N =', Mesh%Nnoeud, ', E=', Mesh%Nfacette, ' , ET=TRIANGLE, F=FEPOINT' 

        do j = 1,Mesh%Nnoeud
            ! Number of 1st order neighbours.
            Nv = 0
            k = 2
            do while (Mesh%Tnoeud(j)%Tvoisin(k,2).eq.1) ! Loop over the 1st order neighbours.
                Nv = Nv + 1
                k = k + 1
            end do
            ! Mesh%Tnoeud(j)%NVoisin(2) = number of all neighbours (1st order and 2nd order).
                write(ioMaillage,'(12E,2I,E,3I,5E,3I,E)') Mesh%Tnoeud(j)%Pnoeud, Mesh%Tnoeud(j)%Normale,Mesh%Tnoeud(j)%Plocal(1:3,1,1),Mesh%Tnoeud(j)%Plocal(1:3,2,1),&
            &   Mesh%Tnoeud(j)%TypeNoeud, Mesh%Tnoeud(j)%Npanneau, Mesh%Tnoeud(j)%Damping, Nv,Mesh%Tnoeud(j)%NVoisin(2), Mesh%Tnoeud(j)%Ndouble, &
            &   Mesh%Tnoeud(j)%DLocal(1:5),merge(1,0,Mesh%Tnoeud(j)%control_point),merge(1,0,Mesh%Tnoeud(j)%mobility),Mesh%Tnoeud(j)%Nfacette,Mesh%Tnoeud(j)%Aire
        end do
        do j = 1,Mesh%Nfacette
            write(ioMaillage,'(4I)') Mesh%Tfacette(j)%Tnoeud
        end do 
        close(unit=3)

    else

        Write( num, '( f0.4 )' ) t
        write(iunit,fmt='(a,i,a,i,a)') 'Zone T =  "'//trim(num)//'seconds", N =', Mesh%Nnoeud, ', E=', Mesh%Nfacette, ' , ET=TRIANGLE, F=FEPOINT, STRANDID = 1, SOLUTIONTIME = '//trim(num)
        do j=1,Mesh%Nnoeud
                write(iunit,'(12E,2I,E,2I,3E,5E,I)') Mesh%Tnoeud(j)%Pnoeud, Mesh%Tnoeud(j)%Normale,Mesh%Tnoeud(j)%Plocal(1:3,1,1),Mesh%Tnoeud(j)%Plocal(1:3,2,1),&
            &   Mesh%Tnoeud(j)%TypeNoeud, Mesh%Tnoeud(j)%Npanneau, Mesh%Tnoeud(j)%Damping, Mesh%Tnoeud(j)%NVoisin(2), Mesh%Tnoeud(j)%Ndouble , Mesh%Tnoeud(j)%Velocity, &
            &   Mesh%Tnoeud(j)%DLocal(1:5),Mesh%Tnoeud(j)%Ordre 
        end do
        do j=1,Mesh%Nfacette
            write(iunit,'(4I)') Mesh%Tfacette(j)%Tnoeud
        end do 
        
    endif

end subroutine PlotMaill

subroutine ExportMesh(Mesh, BodyOnly)
    
    !f2py integer*1, dimension(1000):: Mesh
    type(TMaillage),intent(in)      :: Mesh                                         ! Mesh.
    logical, intent(in), optional   :: BodyOnly                                     ! = True: exporting only the mesh of the bodies, total mesh otherwise.
    
    integer                         :: Nnoeud, Nfacette, Noeud1, Facette1, ios,j,nc ! Parameters.
    logical                         :: isBodyOnly                                   ! = True: exporting only the mesh of the bodies, total mesh otherwise.
    character(len=50)               :: fileMesh                                     ! Mesh file.
    
    ! This subroutine writes the mesh of the bodies in the output file *.dat.
    
    isBodyOnly = .false.
    if (present(BodyOnly)) isBodyOnly = BodyOnly 
    
    if (isBodyOnly) then
        
        ! Mesh of the bodies only.
        
        do nc = Int_Body,Mesh%NBody
            Nnoeud = Mesh%Body(nc)%IndBody(3) - Mesh%Body(nc)%IndBody(1) + 1
            Nfacette = Mesh%Body(nc)%IndBody(4) - Mesh%Body(nc)%IndBody(2) + 1
            Noeud1 = Mesh%Body(nc)%IndBody(1)
            Facette1 = Mesh%Body(nc)%IndBody(2)
        
            write(fileMesh,'("Mesh_Body_",i2,".dat")'),nc-1
            
            ! Opening
            open(unit=ioMeshBody+nc,file=fileMesh, iostat=ios)
            if (ios/=0)then
                print*,"ExportMesh: error at the opening for the body ",nc-1
                stop
            end if
            
            ! Writing
            write(ioMeshBody+nc,fmt='(50a)') 'Title= "Mesh"'
            write(ioMeshBody+nc,fmt='(150a)') 'VARIABLES = "X","Y","Z"'
            write(ioMeshBody+nc,fmt='(a,i,a,i,a)') 'Zone N =', Nnoeud, ', E=', Nfacette, ' , ET=TRIANGLE, F=FEPOINT' 

            do j=Noeud1,Noeud1+Nnoeud-1
                write(ioMeshBody+nc,'(I,3E)') j-(Noeud1-1), Mesh%Tnoeud(j)%Pnoeud
            end do
            do j=Facette1,Facette1+Nfacette-1
                write(ioMeshBody+nc,'(4I)') [Mesh%Tfacette(j)%TNoeud,Mesh%Tfacette(j)%TNoeud(3)] - (Noeud1-1)
            end do 
            
            ! Closing
            close(unit=ioMeshBody+nc)
            
        end do
        
    else
        
        ! Total mesh.
        
        Nnoeud = Mesh%Nnoeud
        Nfacette = Mesh%Nfacette
        Noeud1 = 1
        Facette1 = 1
        
        fileMesh = "Mesh_Total.dat"
        
        ! Opening
        open(unit=ioMeshBody,file=fileMesh, iostat=ios)
        if (ios/=0) stop "ExportMesh: error at the opening."
        
        ! Writing
        write(ioMeshBody,fmt='(50a)') 'Title= "Mesh"'
        write(ioMeshBody,fmt='(150a)') 'VARIABLES = "X","Y","Z"'
        write(ioMeshBody,fmt='(a,i,a,i,a)') 'Zone N =', Nnoeud, ', E=', Nfacette, ' , ET=TRIANGLE, F=FEPOINT' 

        do j=Noeud1,Noeud1+Nnoeud-1
            write(ioMeshBody,'(I,3E)') j-(Noeud1-1), Mesh%Tnoeud(j)%Pnoeud
        end do
        do j=Facette1,Facette1+Nfacette-1
            write(ioMeshBody,'(4I)') [Mesh%Tfacette(j)%TNoeud,Mesh%Tfacette(j)%TNoeud(3)] - (Noeud1-1)
        end do 
        
        ! Closing
        close(unit=ioMeshBody)
                
    end if
    
end subroutine ExportMesh

subroutine ExportMesh_STL(Mesh,fgeom_vect,BodyOnly)
    
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage),intent(in)          :: Mesh                                             ! Mesh.
    !f2py integer*1, dimension(1000)    :: fgeom_vect
    type(type_GeomVect),intent(in)      :: fgeom_vect                                       ! Geometry of the floaters.
    logical, intent(in), optional       :: BodyOnly                                         ! = True: exporting only the mesh of the bodies, total mesh otherwise.
    
    integer                             :: Nnoeud, Nfacette, Noeud1, Facette1, ios, j, k,nc ! Parameters.
    character(len=50)                   :: fileMesh                                         ! Mesh file.
    logical                             :: isBodyOnly                                       ! = True: exporting only the mesh of the bodies, total mesh otherwise.
    
    ! This subroutine writes the mesh of the bodies in the output file *.stl.
    
    isBodyOnly = .false.
    if (present(BodyOnly)) isBodyOnly = BodyOnly

    if (isBodyOnly) then
        
        ! Mesh of the bodies only.
        
        do nc = Int_Body,Mesh%NBody
            if(fgeom_vect%Active(nc-1))then
                Nnoeud = Mesh%Body(nc)%IndBody(3) - Mesh%Body(nc)%IndBody(1) + 1
                Nfacette = Mesh%Body(nc)%IndBody(4) - Mesh%Body(nc)%IndBody(2) + 1
                Noeud1 = Mesh%Body(nc)%IndBody(1)
                Facette1 = Mesh%Body(nc)%IndBody(2)
            
                write(fileMesh,'("Mesh_Body_",i2,".stl")'),nc-1
            
                ! Opening
                open(unit=ioMeshBodSTL+nc,file=fileMesh, iostat=ios)
                if (ios/=0)then
                    print*,"ExportMesh_STL: error at the opening for the body ",nc-1
                    stop
                end if
            
                ! Writing
                write(ioMeshBodSTL+nc,fmt='(50a)') 'solid name'

                do j=Facette1,Facette1+Nfacette-1
                    write(ioMeshBodSTL+nc,fmt=90) Mesh%Tfacette(j)%Normale
                    write(ioMeshBodSTL+nc,fmt='(50a)') '    outer loop'
                    do k=1,3
                        write(ioMeshBodSTL+nc,fmt=91) Mesh%Tnoeud(Mesh%Tfacette(j)%TNoeud(k))%Pnoeud
                    end do
                    write(ioMeshBodSTL+nc,fmt='(50a)') '    endloop'
                    write(ioMeshBodSTL+nc,fmt='(50a)') 'endfacet'
                end do
                write(ioMeshBodSTL+nc,fmt='(50a)') 'endsolid name'
                90 format('facet normal ',3F)
                91 format('     vertex ',3F)
   
                ! Closing
                close(unit=ioMeshBodSTL)
            end if
        end do
            
    else
        
        ! Total mesh.
        
        Nnoeud = Mesh%Nnoeud
        Nfacette = Mesh%Nfacette
        Noeud1 = 1
        Facette1 = 1
        
        fileMesh = "Mesh.stl"
        
        ! Opening
        open(unit=ioMeshBodSTL,file=fileMesh, iostat=ios)
        if (ios/=0) stop "ExportMesh_STL: error at the opening."
        
        ! Writing
        write(ioMeshBodSTL,fmt='(50a)') 'solid name'

        do j=Facette1,Facette1+Nfacette-1
            write(ioMeshBodSTL,fmt=92) Mesh%Tfacette(j)%Normale
            write(ioMeshBodSTL,fmt='(50a)') '    outer loop'
            do k=1,3
                write(ioMeshBodSTL,fmt=93) Mesh%Tnoeud(Mesh%Tfacette(j)%TNoeud(k))%Pnoeud
            end do
            write(ioMeshBodSTL,fmt='(50a)') '    endloop'
            write(ioMeshBodSTL,fmt='(50a)') 'endfacet'
        end do
        write(ioMeshBodSTL,fmt='(50a)') 'endsolid name'
        92 format('facet normal ',3F)
        93 format('     vertex ',3F)
   
        ! Closing
        close(unit=ioMeshBodSTL)
        
    end if
    
end subroutine ExportMesh_STL

subroutine PlotBVP(Mesh, Ecoulement, t, TestorNot, filen)
    
    !f2py integer*1, dimension(1000)            :: Mesh
    type(TMaillage), intent(in)                 :: Mesh                 ! Mesh.
    !f2py integer*1, dimension(1000)            :: Ecoulement
    type(TEcoulement), intent(in)               :: Ecoulement           ! Flow parameters.
    real(rp), intent(in)                        :: t                    ! Current time.
    logical,intent(in)                          :: TestorNot            !  = True: BVP from Test_Solveur, = false, only BVP.
    character(len=50), intent(in), optional     :: filen                ! Output file name.
    
    integer                                     :: j, nc, ios, iotemp   ! Loop parameters.
    real(rp)                                    :: MaxE                 ! Maximum of Phi or DPhiDn.
    real(rp), dimension(Mesh%Nsys)              :: Error                ! Relative errors.
    character(len=10)                           :: num                  !
    
    ! This subroutine writes the Phi and DPhiDn for the mesh and if it is the test the relative errors (Eq 4.21 of LL) in the case of a Rankine source located above the free surface.
    
    if (present(filen)) then
        iotemp = 3
        open(unit=iotemp,file=filen, iostat=ios)
        if (ios/=0) stop "Erreur à l'ouverture du fichier de gradient"
        write(iotemp,fmt='(50a)') 'Title= "Sortie BVP"'
        write(iotemp,fmt='(50a)') 'VARIABLES = "X","Y","Z","Phi_p","Phi_i","DPhiDn_p","DPhiDn_i","Error"'
        write(iotemp,fmt='(a,i,a,i,a)') 'Zone N =', Mesh%Nsys, ', E=', Mesh%Nfsys, ' , ET=TRIANGLE, F=FEPOINT'
    else
        if(TestorNot)then
            iotemp = ioTest_BVP ! Test_BVP
        else
            iotemp = ioBVP ! Only BVP.
        end if
        write( num, '( f0.4 )' ) t
        write(iotemp,fmt='(a,i,a,i,a)')  'Zone T = "'//trim(num)//'seconds", N =', Mesh%Nsys, ', E=', Mesh%Nfsys, ' , ET=TRIANGLE, F=FEPOINT, STRANDID = 1, SOLUTIONTIME = '//trim(num)
    end if
    
    ! Error Calculation
    ! FS
    MaxE = 0._RP
    do j=Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)
        MaxE = max(MaxE,abs(Ecoulement%DPhiDn(j)%Incident))
    end do
    do j=Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)
        if (abs(MaxE).lt.Epsilon) MaxE = 1._RP
        Error(j) = 100._RP*abs(Ecoulement%DPhiDn(j)%Perturbation - Ecoulement%DPhiDn(j)%Incident)/abs(MaxE)
    end do
    
    ! Bodies
    if (cuve_ferme) then
        do nc=1,Mesh%NBody
            MaxE = 0._RP
            do j=Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                MaxE = max(MaxE,abs(Ecoulement%Phi(j)%Incident))
            end do
            if (abs(MaxE).lt.Epsilon) MaxE = 1._RP
            do j=Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                Error(j) = 100._RP*abs(Ecoulement%Phi(j)%Perturbation - Ecoulement%Phi(j)%Incident)/abs(MaxE)
            end do
        end do
    else
        do nc=1,Mesh%NBody
            MaxE = 0._RP
            do j=Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                MaxE = max(MaxE,abs(Ecoulement%Phi(j)%Incident))
            end do
            if (abs(MaxE).lt.Epsilon) MaxE = 1._RP
            do j=Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                Error(j) = 100._RP*abs(Ecoulement%Phi(j)%Perturbation - Ecoulement%Phi(j)%Incident)/abs(MaxE)
            end do
        end do
    end if
    
    ! Plots
    do j=1,Mesh%Nsys
        write(iotemp,'(8E)') Mesh%Tnoeud(j)%Pnoeud, Ecoulement%Phi(j)%Perturbation, Ecoulement%Phi(j)%Incident, Ecoulement%DPhiDn(j)%Perturbation, Ecoulement%DPhiDn(j)%Incident, Error(j)
    end do
    do j=1,Mesh%Nfsys
        write(iotemp,'(4I)') Mesh%Tfacette(j)%Tnoeud
    end do 
    if (present(filen)) close(unit=iotemp)
    
end subroutine PlotBVP

subroutine PlotDerive(filen, Mesh, Ecoulement)
! Paramètres :
character(len=50), intent(in) :: filen
!f2py integer*1, dimension(1000) :: Mesh
type(TMaillage), intent(in) :: Mesh
!f2py integer*1, dimension(1000) :: Ecoulement
type(TEcoulement), intent(in) :: Ecoulement
! Variables :
!character(len=50) :: fileout
integer :: j, ios
!fileout='Remaillage_'//filename
open(unit=3,file=filen, iostat=ios)
if (ios/=0) stop "Erreur à l'ouverture du fichier de gradient"
write(3,fmt='(50a)') 'Title= "Sortie BVP"'
write(3,fmt='(50a)') 'VARIABLES = "X","Y","Z","Eta","Phi","DEtaDt","DPhiDt"' !,"Ndouble","double"'
write(3,fmt='(a,i,a,i,a)') 'Zone N =', Mesh%Nsys, ', E=', Mesh%Nfsys, ' , ET=TRIANGLE, F=FEPOINT' 
do j=1,Mesh%Nsys
    write(3,'(7E)') Mesh%Tnoeud(j)%Pnoeud, Ecoulement%Eta(j)%perturbation, Ecoulement%Phi(j)%perturbation, Ecoulement%DEtaDt(j)%perturbation, Ecoulement%DPhiDt(j)%perturbation
end do
!write(3,'(3E,1I)') 0, 0, 0, 0
do j=1,Mesh%Nfsys
    write(3,'(4I)') Mesh%Tfacette(j)%Tnoeud
end do 
close(unit=3)
end subroutine PlotDerive

!subroutine PlotPression(t, Mesh, FpreN, FpreN0, FpreN1, FpreN2, FpreF, FpreF0, FpreF1, FpreF2)
!! Parameters
!real(rp), intent(in) :: t
!!f2py integer*1, dimension(1000) :: Mesh
!type(TMaillage), intent(in) :: Mesh
!real(rp), dimension(:,:), intent(in) :: FpreN, FpreN0, FpreN1, FpreN2
!real(rp), dimension(:,:), intent(in), optional :: FpreF, FpreF0, FpreF1, FpreF2
!! Locals
!integer :: j
!real(rp) :: moitie
!real(rp), dimension(8) :: Fpression
!
!moitie = 1._RP
!if (symmetry) moitie = 2._RP
!Fpression = 0._RP
!if(is_body)then
!    do j=Mesh%Body(Int_body)%IndBody(1),Mesh%Body(Int_body)%IndBody(3)
!        Fpression(1) = Fpression(1) + FpreN(3,j)*moitie
!        Fpression(2) = Fpression(2) + FpreN0(3,j)*moitie
!        Fpression(3) = Fpression(3) + FpreN1(3,j)*moitie
!        Fpression(4) = Fpression(4) + FpreN2(3,j)*moitie
!    end do
!    if (present(FpreF)) then
!        do j=Mesh%Body(Int_body)%IndBody(2),Mesh%Body(Int_body)%IndBody(4)
!            Fpression(5) = Fpression(5) + FpreF(3,j)*moitie
!            Fpression(6) = Fpression(6) + FpreF0(3,j)*moitie
!            Fpression(7) = Fpression(7) + FpreF1(3,j)*moitie
!            Fpression(8) = Fpression(8) + FpreF2(3,j)*moitie
!        end do
!    end if
!endif
!if (present(FpreF)) then
!    write(58, '(9E)')  t, Fpression
!else
!    write(58, '(5E)')  t, Fpression(1:4)
!end if
!! End
!end subroutine PlotPression

subroutine PlotPress(t, Mesh)
implicit none
real(rp), intent(in) :: t
!f2py integer*1, dimension(1000) :: Mesh
type(TMaillage), intent(in) :: Mesh
!   Begin
write(iopress,'(13E)') t, Mesh%Body(Int_Body)%FBody(1:12) !, Mesh%Body(Int_Body)%Moment(1:6)
! End
end subroutine PlotPress

subroutine PlotPropHoule(t, Mesh, Ecoulement)
    use Incident_mod
    use Houle_RF
    use Airy
    
    !f2py integer*1, dimension(1000)            :: Mesh
    type(TMaillage), intent(in)                 :: Mesh         ! Mesh.
    !f2py integer*1, dimension(1000)            :: Ecoulement
    type(TEcoulement), intent(inout), optional  :: Ecoulement   ! Flow parameters.
    real(rp), intent(in)                        :: t            ! Current time.
    
    integer                                     :: j            ! Loop parameter.
    character(len=20)                           :: num          ! Solution time.
    real(rp)                                    :: pression     ! Pression from Bernoulli equation.
    real(rp)                                    :: TimeRamp     ! Ramp function.
    real(rp), dimension(3)                      :: GPhip, GPhi0 ! Gradients.
    real(rp)                                    :: DPhiDt, GPhi2! DPhiDt and square gradient.
    
        
    ! This subroutine updates the output file output.dat
    
    write( num, '( f0.4 )' ) t
    
    write(unit=33,fmt='(a,i,a,i,a)') 'Zone T = "'//trim(num)//'seconds", N =', Mesh%Nnoeud, ', E=', Mesh%Nfacette,&
    &  ' , ET=TRIANGLE, F=FEPOINT, STRANDID = 1, SOLUTIONTIME = '//trim(num)
    
    ! Time ramp.
    call Ramp_Management(t,TimeRamp)
    
    ! Nodes.
    do j = 1,Mesh%Nnoeud
        
        if(Mesh%Tnoeud(j)%typeNoeud /= 0) Ecoulement%Eta(j)%Perturbation = 0._rp
        
        GPhip = Ecoulement%GPhi(:,j)%perturbation
        GPhi0 = TimeRamp*Ecoulement%GPhi(:,j)%incident
        DPhiDt = TimeRamp*Ecoulement%DPhiDt(j)%incident + Ecoulement%DPhiDt(j)%perturbation
        GPhi2 = dot_product(GPhi0+GPhip,GPhi0+GPhip) - dot_product(GPhip,GPhip)
        
        pression = -ro*(DPhiDt + 0.5_rp*GPhi2 + g*Mesh%Tnoeud(j)%Pnoeud(3))
        
        write(33,'(27E)') Mesh%Tnoeud(j)%Pnoeud(1:2), Mesh%Tnoeud(j)%Pnoeud(3) + 0*Ecoulement%Eta(j)%perturbation,Mesh%Tnoeud(j)%Velocity(1:3),&
                    & Ecoulement%Eta(j)%perturbation, Ecoulement%Phi(j)%perturbation, Ecoulement%GPhi(:,j)%perturbation,&
                    & Ecoulement%Eta(j)%incident , Ecoulement%Phi(j)%incident, Ecoulement%GPhi(:,j)%incident &
                    &, Ecoulement%DPhiDt(j)%perturbation, Ecoulement%DPhiDn(j)%perturbation, Ecoulement%GEta(:,j)%perturbation, Ecoulement%DEtaDt(j)%perturbation,&
                    & Ecoulement%DDPhiDnDt(j)%perturbation, Ecoulement%DGradPhiSqDn(j)%perturbation,pression
          
    end do
    
    ! Panels.
    do j=1,Mesh%Nfacette
        write(33,'(3I)') Mesh%Tfacette(j)%Tnoeud
    end do
    
end subroutine PlotPropHoule

subroutine PlotPression3D(t, Mesh, Ecoulement)
    real(rp), intent(in) :: t
    !f2py integer*1, dimension(1000) :: Mesh
    type(TMaillage), intent(in) :: Mesh
    !f2py integer*1, dimension(1000) :: Ecoulement
    type(TEcoulement), intent(in) :: Ecoulement
!   local
    integer :: j
    integer :: Ind_sl0, Ind_sl1, Ind_body0, Ind_body1
    real(rp) :: GPhi2, DPhiDt
    real(rp), dimension(3) :: M, GPhi0, GPhip
    real(rp), dimension(:), allocatable :: press
    real(rp), dimension(:,:), allocatable :: TPress
    character(len=10) :: num

    allocate(press(Mesh%Nnoeud))
    allocate(TPress(3,Mesh%Nnoeud))

    Ind_sl0 = Mesh%FS%IndFS(1)
    Ind_sl1 = Mesh%FS%IndFS(3)
    Ind_body0 = Mesh%Body(Int_body)%IndBody(1)
    Ind_body1 = Mesh%Body(Int_Body)%IndBody(3)

    do j=Ind_sl0,Ind_sl1    
        M = Mesh%Tnoeud(j)%Pnoeud
        GPhi0 = Ecoulement%GPhi(1:3,j)%incident
        GPhip = Ecoulement%GPhi(1:3,j)%perturbation
        DPhiDt = Ecoulement%DPhiDt(j)%incident + Ecoulement%DPhiDt(j)%perturbation
        DPhiDt = DPhiDt - dot_product(Mesh%Tnoeud(j)%Velocity,GPhip)
        GPhi2 = dot_product(GPhi0+GPhip,GPhi0+GPhip)
        press(j) = -ro*(DPhiDt + 0.5*GPhi2 + g*M(3))
        TPress(1:3,j) = [-ro*DPhiDt, -0.5*ro*GPhi2, -ro*g*M(3)]
    enddo

    do j=Ind_body0,Ind_body1
        M = Mesh%Tnoeud(j)%Pnoeud
        GPhi0 = Ecoulement%GPhi(1:3,j)%incident
        GPhip = Ecoulement%GPhi(1:3,j)%perturbation
        DPhiDt = Ecoulement%DPhiDt(j)%incident + Ecoulement%DPhiDt(j)%perturbation
        GPhi2 = dot_product(GPhi0+GPhip,GPhi0+GPhip)
        press(j) = -ro*(DPhiDt + 0.5*GPhi2 + g*M(3))
        TPress(1:3,j) = [-ro*DPhiDt, -ro*0.5*GPhi2, -ro*g*M(3)]
    enddo

    write(num, '(f0.4)') t
!
    write(iopress3D,fmt='(a,i,a,i,a)') 'Zone T = "'//trim(num)//'seconds", N =', Mesh%Nnoeud, ', E=', Mesh%Nfacette,&
  &  ' , ET=TRIANGLE, F=FEPOINT, STRANDID = 1, SOLUTIONTIME = '//trim(num)

    do j=1,Mesh%Nnoeud
        write(iopress3D,'(7e)') Mesh%Tnoeud(j)%Pnoeud(1:3),press(j),TPress(1:3,j)
    enddo
!
    do j=1,Mesh%Nfacette
        write(iopress3D,'(3I)') Mesh%Tfacette(j)%TNoeud
    end do

    if(allocated(press)) deallocate(press) 
    if(allocated(TPress)) deallocate(TPress)


end subroutine PlotPression3D

subroutine PlotDPhiDt(Mesh, Ecoulement, t, iotruc)
    
    !f2py integer*1, dimension(1000) :: Mesh
    type(TMaillage),intent(in) :: Mesh
    !f2py integer*1, dimension(1000) :: Ecoulement
    type(TEcoulement),intent(in) :: Ecoulement
    real(rp),intent(in)  :: t
    integer :: iotruc
!   Locals
    integer :: j
    character(len=10) :: num

    write(num, '(f0.4)') t
    write(iotruc,fmt='(a,i,a,i,a)') 'Zone T = "'//trim(num)//'seconds", N =', Mesh%Body(Int_body)%IndBody(3)-Mesh%Body(Int_body)%IndBody(1)+1,&
    &', E=', Mesh%Body(Int_body)%IndBody(4)-Mesh%Body(Int_body)%IndBody(2)+1,  ' , ET=TRIANGLE, F=FEPOINT, STRANDID = 1, SOLUTIONTIME = '//trim(num)

    do j=Mesh%Body(Int_body)%IndBody(1),Mesh%Body(Int_body)%IndBody(3)
        write(iotruc,'(5E)') Mesh%Tnoeud(j)%Pnoeud(1:3), Ecoulement%Phi(j)%Perturbation, Ecoulement%DpPhiDt(j)%Perturbation, Ecoulement%DPhiDt(j)%Perturbation
    enddo

    do j=Mesh%Body(Int_body)%IndBody(2),Mesh%Body(Int_body)%IndBody(4)
        write(iotruc,'(3I)') Mesh%Tfacette(j)%Tnoeud-Mesh%Body(Int_body)%IndBody(1)+1
    enddo
end subroutine PlotDPhiDt

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PlotMBody
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine PlotMBody(t,Mesh)
    implicit none
    real(rp),intent(in)                 :: t    ! Current time
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage),intent(in)          :: Mesh ! Mesh
    
    integer                             :: nc   ! Loop parameter
    
    ! This subroutine writes the data in the output file "output_filename.dat"
    
    do nc = Int_Body,Mesh%NBody
        write(iombody+nc,'(19f16.8)') t,Mesh%Body(nc)%CSolv(1:6),Mesh%Body(nc)%VBody(1:6),Mesh%Body(nc)%ABody(1:6)
    end do
    
end subroutine PlotMBody

subroutine PlotBVPPhit(t,Mesh,Ecoulement)
    implicit none
    real(rp),intent(in) :: t
    !f2py integer*1, dimension(1000) :: Mesh
    type(TMaillage),intent(in) :: Mesh
    !f2py integer*1, dimension(1000) :: Ecoulement
    type(TEcoulement),intent(in) :: Ecoulement
!   local
    integer :: j
    integer :: Ind_body0, Ind_body1, Ind_face0, Ind_face1
    integer :: Nb, Nf
    real(rp) :: AccN
    character(len=50) :: num
    

    Ind_body0 = Mesh%Body(Int_Body)%IndBody(1)
    Ind_body1 = Mesh%Body(Int_Body)%IndBody(3)
    Nb = Ind_body1 - Ind_body0 +1

    Ind_face0 = Mesh%Body(Int_Body)%IndBody(2)
    Ind_face1 = Mesh%Body(Int_Body)%IndBody(4)
    Nf = Ind_face1 - Ind_face0 + 1

    write(num, '(f0.4)') t

    write(iophit,fmt='(a,i,a,i,a)') 'Zone T = "'//trim(num)//'seconds", N =', Nb, ', E=', Nf,&
  &  ' , ET=TRIANGLE, F=FEPOINT, STRANDID = 1, SOLUTIONTIME = '//trim(num)

    do j=Ind_body0,Ind_body1
        AccN = dot_product(Mesh%Body(Int_Body)%ABody(1:3),Mesh%Tnoeud(j)%Normale(1:3))
        write(iophit,'(7e)') Mesh%Tnoeud(j)%Pnoeud(1:3), Ecoulement%DDPhiDnDt(j)%incident, Ecoulement%DDPhiDnDt(j)%perturbation,&
&                           Ecoulement%DGradPhiSqDn(j)%perturbation, AccN
    enddo

    do j=Ind_face0,Ind_face1
        write(iophit,'(3i)') Mesh%Tfacette(j)%Tnoeud(1:3)-Ind_body0+1
    enddo

end subroutine PlotBVPPhit

subroutine PlotMonitor(t,MeshQualityBodies,MeshQualityFS,CrossingFS,ForcedRemesh,boolRemeshFS,boolRemesh,ierror)
    
    real(rp),intent(in) :: t                                    ! Current time.
    logical,intent(in)  :: MeshQualityBodies                    ! = false: no remeshing of the bodies, = true: remeshing of the bodies.
    logical,intent(in)  :: MeshQualityFS                        ! = false: no remeshing of the free surface, = true: remeshing of the free surface.
    logical,intent(in)  :: CrossingFS                           ! A body crossed the free surface (True) or not (False).
    logical,intent(in)  :: ForcedRemesh                         ! Force the remeshing if necessary (crossing the free surface).
    logical,intent(in)  :: boolRemesh                           ! Boolean to know if the mesh of the bodies or/and the free surface was regenerated (True) or not (False).
    logical,intent(in)  :: boolRemeshFS                         ! = true: remeshing the free surface, = false: no remeshing.
    integer,intent(in)  :: ierror                               ! Error flag to know if the remeshing worked (=0) or not (!=0).
    
    integer             :: int1,int2,int3,int4,int5,int6,int7   ! Integers.
    
    ! This subroutine writes the information in case of remeshing the output file manitoring.dat.
    
    ! Body mesh quality.
    if(MeshQualityBodies)then
        int1 = 1
    else
        int1 = 0
    endif
    
    ! Free surface mesh quality.
    if(MeshQualityFS)then
        int2 = 1
    else
        int2 = 0
    endif
    
    ! Crossing the free surface.
    if(CrossingFS)then
        int3 = 1
    else    
        int3 = 0
    endif
    
    ! Forced remeshing.
    if(ForcedRemesh)then
        int4 = 1
    else    
        int4 = 0
    endif
    
    ! FS Remeshing asked.
    if(boolRemesh)then
        int5 = 1
    else    
        int5 = 0
    endif
    
    ! Remeshing asked.
    if(boolRemesh)then
        int6 = 1
    else    
        int6 = 0
    endif
    
    ! Remeshing worked.
    if(ierror/=0) then
        int7 = 0
    else
        int7 = 1
    end if
    
    write(iomonitor,'(f8.4,7i4)') t, int1, int2, int3, int4, int5, int6, int7

end subroutine  PlotMonitor


subroutine PlotEcoulement(Mesh, Ecoulement)
    
    !f2py integer*1, dimension(1000) :: Mesh
    type(TMaillage),intent(in) :: Mesh
    !f2py integer*1, dimension(1000) :: Ecoulement
    type(TEcoulement),intent(in) :: Ecoulement
!   Locals
    integer :: j, iotruc, ios
    character(len=10) :: num

    iotruc = 3226
    open(unit=iotruc,file='Ecoulement_'//filename,iostat=ios)
    if(ios/=0) stop "Erreur à l'ouverture du fichier Ecoulement"
    write(iotruc,fmt='(50a)') 'Title = "Composantes incidentes sur le Corps"'
    write(iotruc,fmt='(50a)') 'VARIABLES = "X","Y","Z","Phi","DDPhiDzDt","DGPhiDzx","DGPhiDzy","DGPhiDzz"'
    !,"DPhiDn","GPhix","GPhiy","GPhiz","GPhi2u","GPhi2v","GPhi2uu","GPhi2vv","GradGrad(1,1)","GradGrad(2,2)","GradGrad(3,3)","DGPhiDz1","DGPhiDz2","DGPhiDz3","DPhiDt","DDPhiDnDt"'
    
    write(num, '(f0.4)') 32._RP
    !write(iotruc,fmt='(a,i,a,i,a)') 'Zone T = "'//trim(num)//'seconds", N =', Mesh%Body(Int_body)%IndBody(3)-Mesh%Body(Int_body)%IndBody(1)+1,&
    !&', E=', Mesh%Body(Int_body)%IndBody(4)-Mesh%Body(Int_body)%IndBody(2)+1,  ' , ET=TRIANGLE, F=FEPOINT, STRANDID = 1, SOLUTIONTIME = '//trim(num)
    !write(iotruc,fmt='(a,i,a,i,a)') 'Zone T = "'//trim(num)//'seconds", N =', Mesh%Nnoeud,', E=', Mesh%Nfacette,  ' , ET=TRIANGLE, F=FEPOINT, STRANDID = 1, SOLUTIONTIME = '//trim(num)
    write(iotruc,fmt='(a,i,a,i,a)') 'Zone N =', Mesh%Body(Int_body)%IndBody(3)-Mesh%Body(Int_body)%IndBody(1)+1, ', E=', Mesh%Body(Int_body)%IndBody(4)-Mesh%Body(Int_body)%IndBody(2)+1, ' , ET=TRIANGLE, F=FEPOINT' 
    do j=Mesh%Body(Int_body)%IndBody(1),Mesh%Body(Int_body)%IndBody(3)
    !do j=1,Mesh%Nnoeud
        write(iotruc,'(25E)') Mesh%Tnoeud(j)%Pnoeud(1:3), Ecoulement%Phi(j)%Incident, Ecoulement%DDPhiDnDt(j)%Incident, Ecoulement%DGPhiDz(:,j)%incident
        !, Ecoulement%DPhiDn(j)%Incident&
        !&, Ecoulement%GPhi(1:3,j)%Incident, Ecoulement%GPhi2(1:4,j)%Incident&
        !&, Ecoulement%GradGrad(1,1,j)%Incident, Ecoulement%GradGrad(2,2,j)%Incident, Ecoulement%GradGrad(3,3,j)%Incident&
        !&, Ecoulement%DGPhiDz(:,j)%Incident, Ecoulement%DPhiDt(j)%Incident, Ecoulement%DDPhiDnDt(j)%Incident
    end do
    do j=Mesh%Body(Int_body)%IndBody(2),Mesh%Body(Int_body)%IndBody(4)
    !do j=1,Mesh%Nfacette
        write(iotruc,'(3I)') Mesh%Tfacette(j)%Tnoeud-Mesh%Body(Int_body)%IndBody(1)+1
    end do
    close(iotruc)
end subroutine PlotEcoulement

subroutine PlotWaveElevation(t,InputData,Mesh,Ecoulement)
    
    real(rp),intent(in)                 :: t                ! Current time.
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(in)    :: InputData        ! Input data.
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage)                     :: Mesh             ! Mesh.
    !f2py integer*1, dimension(1000)    :: Ecoulement
    type(TECoulement)                   :: Ecoulement       ! Flow parameters.
        
    real(rp)                            :: Eta_Incident     ! Incident wave elevation.
    real(rp)                            :: Eta_Pertubation  ! Perturbation wave elevation.
    integer                             :: nc,k             ! Loop parameters.
    real(rp)                            :: dist,PjPk        ! Distances between two points.
    integer                             :: Na, Nvoisin      ! Spline ordre and number of neighbours.
    integer                             :: ClosestPoints    ! Closest point of the mesh for the position of a wave probe.
    real(rp),allocatable                :: Pvoisin(:,:)     ! Position of the neighbours.
    real(rp),allocatable                :: A(:,:), B(:)     ! Matrix A and vector B for the B-spline problem.
    integer                             :: info             !
    character(len=1)                    :: trans            !
    integer                             :: lda, nrhs        ! Size of the matrix A and number of rhs.
    integer, allocatable                :: ipiv(:)          !
        
    ! This subroutine writes the wave elevation output files.
    
    ! Computation of the wave elevations
    do nc = 1,nWP
        
        ! Eta_Incident.
        select case(Htype)
        case(0) ! Still water
            Eta_Incident = 0._RP
        case(1) ! Airy waves
            call Computation_EtaAiry(InputData%PositionWP(1:3,nc),t,Eta_Incident)
        case(2) ! RF
            call Computation_EtaRF(InputData%PositionWP(1:3,nc), t,Eta_Incident)
        end select
        
        ! Eta_Perturbation from B-spline interpolation of Ecoulement%Eta_p
        Eta_Pertubation = 0._RP
        ClosestPoints = 0
        dist = 999._RP
        do k = Mesh%FS%IndFS(1), Mesh%FS%IndFS(3)
            PjPk = norm2(InputData%PositionWP(1:3,nc)-Mesh%Tnoeud(k)%Pnoeud)
            if(PjPk .le. dist)then
                dist = PjPk
                ClosestPoints = k
            end if
        end do
        
        if(is_BS)then
            
            ! Spline ordre at the closest point of the wave probe.
            Na = Nordre(Mesh%Tnoeud(ClosestPoints)%Ordre)
            
            ! Number of neighbours of the closest point.
            NVoisin = Mesh%Tnoeud(ClosestPoints)%NVoisin(2)-1 ! -1 because the point itself is counted.
            
            ! Allocations
            allocate(Pvoisin(3,NVoisin+1), B(NVoisin+Na), A(NVoisin+Na,NVoisin+Na))
                                    
            ! Pvoisin and B
            do k=1,NVoisin+1 
                
                ! Neighbours of the points.
                Pvoisin(1:3,k) = Mesh%Tnoeud(abs(Mesh%Tnoeud(ClosestPoints)%TVoisin(k,1)))%Pnoeud ! k = 1 is Mesh%Tnoeud(ClosestPoints) itself.
                
                if(Mesh%Tnoeud(abs(Mesh%Tnoeud(ClosestPoints)%TVoisin(k,1)))%NPanneau.ne.0)then
                    print*,"PlotWaveElevation: Point is not on the FS:"
                    pause
                end if
                
                if (Mesh%Tnoeud(ClosestPoints)%TVoisin(k,1).lt.0) Pvoisin(2,k) = - Pvoisin(2,k)
                B(k) = Ecoulement%Eta(abs(Mesh%Tnoeud(ClosestPoints)%TVoisin(k,1)))%perturbation
                
            end do
            B(Nvoisin+2:Nvoisin+Na) = 0._RP
            
            ! A
            call SplineMatrix(Nvoisin, Mesh%Tnoeud(ClosestPoints)%Ordre, Pvoisin, A)
            
            ! Inversion of the linear system
            lda = Nvoisin + Na
            allocate(ipiv(lda))
            ipiv = 0
            trans = 'n'
            nrhs = 1 ! Number of rhs.
            call dgetrf(lda,lda,A,lda,ipiv,info)
            if (info.ne.0) then
                print*,'WARNING Interpolation: info =', info, lda, Mesh%Tnoeud(ClosestPoints)%NVoisin(2), Na
                print*,Mesh%Tnoeud(ClosestPoints)%TVoisin(1:Mesh%Tnoeud(ClosestPoints)%NVoisin(2),1)
                pause
            end if
            call dgetrs(trans,lda,nrhs,A,lda,ipiv,B,lda,info) ! The B-spline coefficients are in B.
            deallocate(ipiv)
            
            ! Eta_p
            call SplineF(Nvoisin, Mesh%Tnoeud(ClosestPoints)%Ordre, B(1:Nvoisin+Na), InputData%PositionWP(1:3,nc), Pvoisin, Eta_Pertubation) ! InputData%PositionWP(1:3,nc) is where SplineF is computed.
            
            ! Deallocations
            if(allocated(Pvoisin)) deallocate(Pvoisin)
            if(allocated(B)) deallocate(B)
            if(allocated(A)) deallocate(A)
            
        else
            print*,"PlotWaveElevation: interpolation algorithm only works with B-splines."
            call exit()
        end if
        
        ! Writing.
        write(iowave+nc,'(4E)') t,Eta_Incident,Eta_Pertubation,Eta_Incident + Eta_Pertubation
        
    end do
    
end subroutine PlotWaveElevation

subroutine Write_Inertia(Mesh,NumBody)

    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage),intent(in)          :: Mesh         ! Mesh.
    integer,intent(in)                  :: NumBody      ! Number of the body.
    
    integer                             :: j            ! Loop parameter.
    integer                             :: ios          ! Output flag.
    character(len=50)                   :: fileInertia  ! Inertia file.
    
    ! This subroutine writes the inertia matrix in the ouput file Inertia_Body_Computed_by_WSC.dat
    
        
    write(fileInertia,'("Inertia_Body_",i2,"_Computed_by_WSC.dat")'),NumBody-1
        
    ! Opening.
    open(unit=ioWriteInertia+NumBody,file=fileInertia, iostat=ios)
    if (ios/=0)then 
        print*,"Write_Inertia: error at the opening for the body ",NumBody-1
        stop
    end if
        
    ! Writing.
    write(ioWriteInertia+NumBody,fmt='(50a)') 'Title = "Inertia"'
        
    do j=1,6
        write(ioWriteInertia+NumBody,'(3E,3E,3E,3E,3E,3E)') Mesh%Body(NumBody)%IBody(j,1:6)
    end do
    
    ! Closing.
    close(unit=ioWriteInertia+NumBody)
                
    
end subroutine Write_Inertia

subroutine Write_State(Mesh,Ecoulement,t,jt,Starting_time,jFiltering)
    
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage),intent(in)          :: Mesh         ! Mesh.
    !f2py integer*1, dimension(1000)    :: Ecoulement
    type(TEcoulement),intent(in)        :: Ecoulement   ! Flow parameters.
    real(rp),intent(in)                 :: t            ! Current time.
    integer,intent(in)                  :: jt           ! Time loop parameter.
    real(rp),intent(in)                 :: Starting_time! Starting time saved.
    integer,intent(in)                  :: jFiltering   ! Counter to wait for 200 time steps before using the original parameters in case of a body crossing the free surface.
    
    integer                             :: nc,j,p       ! Loop parameters.
    integer                             :: ios          ! Output flag.
    
    ! This subroutine writes the state vectors (time, position of the floaters, velocity of the floaters, position of the nodes, perturbation velocity potential on the free surface, pertubation wave elevation on the free surface) in State.txt.
    
    if(iwState)then
        
        ! Opening.
        ! Two back-up files is safer.
        if(modulo(jt,2).lt.Epsilon)then ! jt is even.
            open(unit=ioState,file="State_Even.dat",iostat=ios)
        else ! jt is odd.
            open(unit=ioState,file="State_Odd.dat",iostat=ios)
        end if
        if(ios/=0) stop "Write_State: error in the opening of State.txt"
        write(ioState,'(50a)') 'Title = "State vector"'
        
        ! Starting time.
        write(ioState,'(E)') Starting_time
        
        ! Current time.
        write(ioState,'(E)') t
        
        ! Current jt.
        write(ioState,'(I)') jt
                
        ! jFiltering.
        write(ioState,'(I)') jFiltering
        
        ! Nnoeud
        write(ioState,'(I)') Mesh%Nnoeud
        
        ! Active.
        do nc = 1,Mesh%NBody
            write(ioState,'(L)') Mesh%Body(nc)%Active
        end do
        
        ! MBody.
        do nc = Int_Body,Mesh%NBody
            write(ioState,'(6E)') Mesh%Body(nc)%MBody
        end do
        
        ! CSolv.
        do nc = Int_Body,Mesh%NBody
            write(ioState,'(6E)') Mesh%Body(nc)%CSolv
        end do
        
        ! GBody.
        do nc = Int_Body,Mesh%NBody
            write(ioState,'(3E)') Mesh%Body(nc)%GBody
        end do
        
        ! VBody.
        do nc = Int_Body,Mesh%NBody
            write(ioState,'(6E)') Mesh%Body(nc)%VBody
        end do
        
        ! ----------------------------------------------
        !                   Free surface 
        ! ----------------------------------------------
        
        ! IndFS.
        write(ioState,'(4I)') Mesh%FS%IndFS(1),Mesh%FS%IndFS(2),Mesh%FS%IndFS(3),Mesh%FS%IndFS(4)
        
        ! PNoeud, NVoisin and TVoisin on the free surface.
        do j = Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)
            write(ioState,'(3E)') Mesh%Tnoeud(j)%Pnoeud
            write(ioState,'(2I)') Mesh%Tnoeud(j)%NVoisin
            do p = 1,Mesh%Tnoeud(j)%NVoisin(2)
                write(ioState,'(2I)') Mesh%Tnoeud(j)%TVoisin(p,1),Mesh%Tnoeud(j)%TVoisin(p,2)
            end do
        end do
        
        ! Phi_p, Eta_p.
        do j = Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)
            write(ioState,'(2E)') Ecoulement%Phi(j)%perturbation,Ecoulement%Eta(j)%perturbation
        end do
        
        ! ----------------------------------------------
        !              Bodies (tank included) 
        ! ----------------------------------------------
        
        do nc = 1,Mesh%NBody ! Tank included.
            
            if(Mesh%Body(nc)%Active)then
            
                ! IndBody
                write(ioState,'(4I)') Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(2),Mesh%Body(nc)%IndBody(3),Mesh%Body(nc)%IndBody(4)
            
                ! PNoeud on the body mesh.
                do j = Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                    write(ioState,'(3E)') Mesh%Tnoeud(j)%Pnoeud
                end do
            
            end if
            
        end do
        
        ! Closing.
        close(ioState)
        
    end if
    
end subroutine Write_State


subroutine Write_State2(Mesh,Ecoulement,t,jt,Starting_time,jFiltering, filename)
    
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage),intent(in)          :: Mesh         ! Mesh.
    !f2py integer*1, dimension(1000)    :: Ecoulement
    character(len=50),intent(in)        :: filename    ! State input file.
    type(TEcoulement),intent(in)        :: Ecoulement   ! Flow parameters.
    real(rp),intent(in)                 :: t            ! Current time.
    integer,intent(in)                  :: jt           ! Time loop parameter.
    real(rp),intent(in)                 :: Starting_time! Starting time saved.
    integer,intent(in)                  :: jFiltering   ! Counter to wait for 200 time steps before using the original parameters in case of a body crossing the free surface.
    
    integer                             :: nc,j,p       ! Loop parameters.
    integer                             :: ios          ! Output flag.
    
    ! This subroutine writes the state vectors (time, position of the floaters, velocity of the floaters, position of the nodes, perturbation velocity potential on the free surface, pertubation wave elevation on the free surface) in State.txt.
    

    ! Opening.
    open(unit=ioState,file=filename,iostat=ios)

    if(ios/=0) stop "Write_State: error in the opening of State.txt"
    write(ioState,'(50a)') 'Title = "State vector"'
        
    ! Starting time.
    write(ioState,'(E)') Starting_time
        
    ! Current time.
    write(ioState,'(E)') t
        
    ! Current jt.
    write(ioState,'(I)') jt
                
    ! jFiltering.
    write(ioState,'(I)') jFiltering
        
    ! Nnoeud
    write(ioState,'(I)') Mesh%Nnoeud
        
    ! Active.
    do nc = 1,Mesh%NBody
        write(ioState,'(L)') Mesh%Body(nc)%Active
    end do
        
    ! MBody.
    do nc = Int_Body,Mesh%NBody
        write(ioState,'(6E)') Mesh%Body(nc)%MBody
    end do
        
    ! CSolv.
    do nc = Int_Body,Mesh%NBody
        write(ioState,'(6E)') Mesh%Body(nc)%CSolv
    end do
        
    ! GBody.
    do nc = Int_Body,Mesh%NBody
        write(ioState,'(3E)') Mesh%Body(nc)%GBody
    end do
        
    ! VBody.
    do nc = Int_Body,Mesh%NBody
        write(ioState,'(6E)') Mesh%Body(nc)%VBody
    end do
        
    ! ----------------------------------------------
    !                   Free surface 
    ! ----------------------------------------------
        
    ! IndFS.
    write(ioState,'(4I)') Mesh%FS%IndFS(1),Mesh%FS%IndFS(2),Mesh%FS%IndFS(3),Mesh%FS%IndFS(4)
        
    ! PNoeud, NVoisin and TVoisin on the free surface.
    do j = Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)
        write(ioState,'(3E)') Mesh%Tnoeud(j)%Pnoeud
        write(ioState,'(2I)') Mesh%Tnoeud(j)%NVoisin
        do p = 1,Mesh%Tnoeud(j)%NVoisin(2)
            write(ioState,'(2I)') Mesh%Tnoeud(j)%TVoisin(p,1),Mesh%Tnoeud(j)%TVoisin(p,2)
        end do
    end do
        
    ! Phi_p, Eta_p.
    do j = Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)
        write(ioState,'(2E)') Ecoulement%Phi(j)%perturbation,Ecoulement%Eta(j)%perturbation
    end do
        
    ! ----------------------------------------------
    !              Bodies (tank included) 
    ! ----------------------------------------------
        
    do nc = 1,Mesh%NBody ! Tank included.
            
        if(Mesh%Body(nc)%Active)then
            
            ! IndBody
            write(ioState,'(4I)') Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(2),Mesh%Body(nc)%IndBody(3),Mesh%Body(nc)%IndBody(4)
            
            ! PNoeud on the body mesh.
            do j = Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                write(ioState,'(3E)') Mesh%Tnoeud(j)%Pnoeud
            end do
            
        end if
            
    end do
        
    ! Closing.
    close(ioState)

    
end subroutine Write_State2

subroutine read_State_InputData(fileState,InputData,t,jt0)
    
    character(len=50),intent(in)        :: fileState    ! State input file.
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct)               :: InputData    ! Input data.
    real(rp),intent(inout)              :: t            ! New stating time (t0).
    integer,intent(inout)               :: jt0          ! Time loop parameter.
    
    integer                             :: nc           ! Loop parameter.
    real(rp),dimension(6)               :: Vect         ! Vector of size 6.
    integer                             :: ios          ! Output flag.
    character(len=50)                   :: line         ! Line.
    
    ! This subroutine reads the state input file in order to update InputData. 
    
    if(iwState)then
        
        ! Opening.
        open(unit=ioState,file=fileState,iostat=ios)
        if(ios/=0) stop "read_State_InputData: error in the opening of state input file."
                
        ! Title.
        read(ioState,'(a)',end=1000,err=9999) line
        
        ! Starting time.
        read(ioState,'(a)',end=1000,err=9999) line
        
        ! Current time.
        call read_param_only(ioState,t)
        
        ! jt0 (current jt).
        call read_param_only(ioState,jt0)
                
        ! jFiltering.
        read(ioState,'(a)',end=1000,err=9999) line
        
        ! Nnoeud.
        read(ioState,'(a)',end=1000,err=9999) line
        
        ! Active.
        do nc = 1,(NBodies+1) ! NBodies + 1 = Mesh%NBody.
            read(ioState,'(a)',end=1000,err=9999) line
        end do
        
        ! MBody.
        do nc = 1,NBodies
            read(ioState,'(a)',end=1000,err=9999) line
        end do
        
        ! CSolv.
        do nc = 1,NBodies
            call read_param_only(ioState,Vect,6)
            InputData%Position(1:3,1,nc) = Vect(1:3) ! Linear position.
            InputData%Position(1:3,2,nc) = Vect(4:6) ! Angular position in radians.
        end do
        
        ! GBody.
        do nc = 1,NBodies
            call read_param_only(ioState,InputData%PositionG(1:3,nc),3)
        end do
        
        ! Closing.
        close(ioState)
        
    end if
        
    1000 continue
     
    9999 continue
    
end subroutine read_State_InputData

subroutine read_State(fileState,Mesh,Ecoulement,Starting_time,jFiltering)
    
    character(len=50),intent(in)        :: fileState    ! State input file.
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage),intent(inout)       :: Mesh         ! Mesh.
    !f2py integer*1, dimension(1000)    :: Ecoulement
    type(TEcoulement),intent(inout)     :: Ecoulement   ! Flow parameters.
    real(rp),intent(inout)              :: Starting_time! Starting time saved.
    integer,intent(inout)               :: jFiltering   ! Counter to wait for 200 time steps before using the original parameters in case of a body crossing the free surface.
        
    integer                             :: nc,j,p       ! Loop parameters.
    character(len=50)                   :: line         ! Line.
    integer                             :: ios          ! Output flag.
    integer                             :: Nnoeud       ! Number of nodes. 
    integer,dimension(4)                :: Ind          ! Vector to store the node/panel indexes.
    real(rp),dimension(2)               :: Phi_Eta      ! Vector to store Phi_p and Eta_p.
    
    ! This subroutine reads the state input file. 
    
    if(iwState)then
        
        ! Opening.
        open(unit=ioState,file=fileState,iostat=ios)
        if(ios/=0) stop "read_State: error in the opening of state input file."
                
        ! Title.
        read(ioState,'(a)',end=1000,err=9999) line
        
        ! Starting time.
        call read_param_only(ioState,Starting_time)
        
        ! Current time.
        read(ioState,'(a)',end=1000,err=9999) line
        
        ! jt0 (current jt).
        read(ioState,'(a)',end=1000,err=9999) line
        
        ! jFiltering.
        call read_param_only(ioState,jFiltering)
                
        ! Nnoeud.
        call read_param_only(ioState,Nnoeud)
          
        ! New TMaillage.
        call NewMaillage(Mesh,Nnoeud,NBodies+1) ! +1 for the tank.
        Mesh%Nnoeud = Nnoeud
        
        ! NBody
        Mesh%NBody = NBodies+1 ! 1 for the tank.
        
        ! Active.
        do nc = 1,Mesh%NBody
            call read_param_only(ioState,Mesh%Body(nc)%Active)
        end do
        
        ! MBody.
        do nc = Int_Body,Mesh%NBody
            call read_param_only(ioState,Mesh%Body(nc)%MBody,6)
        end do
        
        ! CSolv.
        do nc = Int_Body,Mesh%NBody
            call read_param_only(ioState,Mesh%Body(nc)%CSolv,6)
        end do
        
        ! GBody.
        do nc = Int_Body,Mesh%NBody
            call read_param_only(ioState,Mesh%Body(nc)%GBody,3)
        end do
        
        ! VBody.
        do nc = Int_Body,Mesh%NBody
            call read_param_only(ioState,Mesh%Body(nc)%VBody,6)
        end do
        
        ! ----------------------------------------------
        !                   Free surface 
        ! ----------------------------------------------
        
        ! IndFS.
        Ind = 0._RP
        call read_param_only(ioState,Ind,4)
        do j = 1,4
            Mesh%FS%IndFS(j) = Ind(j)
        end do
        
        ! PNoeud,NVoisin and TVoisin on the free surface.
        do j = Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)
            call read_param_only(ioState,Mesh%Tnoeud(j)%Pnoeud,3)
            call read_param_only(ioState,Mesh%Tnoeud(j)%NVoisin,2)
            do p = 1,Mesh%Tnoeud(j)%NVoisin(2)
                call read_param_only(ioState,Mesh%Tnoeud(j)%TVoisin(p,1:2),2)
            end do
        end do
        
        ! New TEcoulement.
        call NewEcoulement(Ecoulement,Mesh%Nnoeud)
        call IniEcoulement(Ecoulement,Mesh%Nnoeud, 0._RP) ! %incident and %perturbation are 0.
        
        ! Phi_p, Eta_p.
        Phi_Eta = 0._RP
        do j = Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)
            call read_param_only(ioState,Phi_Eta,2)
            Ecoulement%Phi(j)%perturbation = Phi_Eta(1)
            Ecoulement%Eta(j)%perturbation = Phi_Eta(2)
        end do
        
        ! ----------------------------------------------
        !              Bodies (tank included) 
        ! ----------------------------------------------
        
        do nc = 1,Mesh%NBody ! Tank included.
            
            ! IndBody
            Ind = 0._RP
            call read_param_only(ioState,Ind,4)
            do j = 1,4
                Mesh%Body(nc)%IndBody(j) = Ind(j)
            end do
            
            ! PNoeud on the body mesh.
            do j = Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                call read_param_only(ioState,Mesh%Tnoeud(j)%Pnoeud,3)
            end do
            
        end do
        
        ! Closing.
        close(ioState)
        
    end if
    
    print*,"State input file read."
    
    1000 continue
     
    9999 continue
    
end subroutine read_State

subroutine Extract_Mesh_Body(filemesh,NewMesh,fgeom_vect,InputData)
    
    character(len=150),intent(in)       :: filemesh     ! Tecplot format input file.
    !f2py integer*1, dimension(1000)    :: NewMesh
    type(TMaillage),intent(inout)       :: NewMesh      ! Body Mesh only (no Mesh%FS).
    !f2py integer*1, dimension(1000)    :: fgeom_vect
    type(type_GeomVect),intent(inout)   :: fgeom_vect   ! Geometry of the floaters (not the domain).
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(inout) :: InputData    ! Input data.
    
    integer                             :: ios          ! Opening file parameter.
    character(len=50)                   :: line         ! Line.
    integer                             :: N,E          ! Number of nodes and panels.
    integer                             :: idx          ! Index to search N and E.
    integer                             :: j,k          ! Loop parameters.
    integer,dimension(4)                :: Tnoeud       ! Table of nodes for each panel.
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage)                     :: Mesh         ! Body Mesh only (no Mesh%FS).
        
    ! This subroutine extracts a body mesh from a Tecplot format input file.
    
    ! Constants used.
    NBodies = 1
    idebug = 0
    is_body = .true.
    is_immerged = .false.
    
    ! Opening mesh file.
    open(ioTmpBodyMesh,file=filemesh, iostat=ios)
    if (ios/=0) stop "Extract_Mesh_Body: Error in the opening of the body mesh file."
    
    ! Headlines.
    read(ioTmpBodyMesh,'(a)',end=1000,err=9999) line
    read(ioTmpBodyMesh,'(a)',end=1000,err=9999) line
    read(ioTmpBodyMesh,'(a)',end=1000,err=9999) line
    
    ! Number of nodes and panels.
    read(ioTmpBodyMesh,'(a)',end=1000,err=9999) line
    idx = index(line, 'N=      ') + len('N=      ')
    read(line(idx:idx+3), '(I4)') N
    idx = index(line, 'E=      ') + len('E=      ')
    read(line(idx:idx+3), '(I4)') E
    
    ! New mesh structure.
    call NewMaillage(Mesh,PointMax,2) ! 1 = tank (useless), 2 = body.
    Int_Body = 2
    Mesh%NBody = 2
    Mesh%Nnoeud = N
    Mesh%Nfacette = E
    Mesh%Nsys = N
    Mesh%FS%IndFS(2) = 1 ! 1st FS panel.
    Mesh%FS%IndFS(4) = -1 ! last FS panel. -1 not to go in the for loops.
    Mesh%FS%IndFS(1) = 1 ! 1st FS node.
    Mesh%FS%IndFS(3) = -1 ! Last FS node. -1 not to go in the for loops.
    Mesh%Body(1)%IndBody(2) = 1 ! 1st tank panel.
    Mesh%Body(1)%IndBody(4) = -1 ! Last tank panel. -1 not to go in the for loops.
    Mesh%Body(1)%IndBody(1) = 1 ! 1st tank node.
    Mesh%Body(1)%IndBody(3) = -1 ! Last tank node. -1 not to go in the for loops.
    Mesh%Body(2)%IndBody(2) = 1 ! 1st body panel.
    Mesh%Body(2)%IndBody(4) = E ! Last body panel.
    Mesh%Body(2)%IndBody(1) = 1 ! 1st body node.
    Mesh%Body(2)%IndBody(3) = N ! Last body node.
    
    Mesh%Body(2)%GBody(1:3) = [0._RP,0._RP,-1._RP]
    print*,"Extract_Mesh_Body: GBody hard coded."
    
    ! Nodes.
    do j = 1,N ! Not Mesh%Nnoeud because it is updated in the loop due to the twin nodes.
        read (ioTmpBodyMesh,*) Mesh%Tnoeud(j)%Pnoeud(1:3)
        Mesh%Tnoeud(j)%typeNoeud = 1
        Mesh%Tnoeud(j)%Npanneau = Int_Body
    end do
    
    ! Panels.
    print*,"Extract_Mesh_Body INFO:"
    k = 0
    do j = 1,Mesh%Nfacette
        read (ioTmpBodyMesh,*) Tnoeud(1:4)
        if(Tnoeud(4).eq.Tnoeud(1))then
            Mesh%Tfacette(j)%Tnoeud(1:3) = Tnoeud(1:3)
        else
            print*,"Panel ",j," is not square."
            k = k + 1
        end if
    end do
    print*,k," quadrangular panels."
    
    ! New TMaillage structure.
    call NewMaillage(NewMesh,Mesh%Nnoeud,2) ! 1 = tank (useless), 2 = body.
    call CopyMaillage(NewMesh,Mesh)
    call DelMaillage(Mesh)
    
    ! Allocating of fgeom_vect
    call init_GeomVect(fgeom_vect)
    fgeom_vect%Active(1) = .true. ! 1st floater (not the tank) is active.
    
    close(ioTmpBodyMesh)
    
1000 continue    
     
9999 continue     

end subroutine Extract_Mesh_Body

end module PrePlot