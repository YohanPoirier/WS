module GeomDef
use GeomStruct
use GeomFonct
use GeomGen
use GeomAxiSym
use SolvNum
use GeomWigley
use Parameters

implicit none

! --------------- GeomDef --------------
 
contains
    
subroutine Generation_Geometry(fgeom_vect,fdomaine,nface,tab2,n_tab2,rep0,InputData,ierror,n_tab)
    
    !f2py integer*1, dimension(1000)                :: fdomaine
    type(type_geom),intent(out)                     :: fdomaine     ! Geometry of the domain.
    !f2py integer*1, dimension(1000)                :: fgeom_vect
    type(type_GeomVect),intent(out)                 :: fgeom_vect   ! Geometry of the floaters.
    integer,intent(out)                             :: nface        ! Number of faces in both the floater and the domain.
    !f2py integer*1, dimension(1000)                :: tab2
    type(chaine_point_pt),dimension(100),intent(out):: tab2         ! Table of the pointers toward the intersection points.
    integer,intent(out)                             :: n_tab2,n_tab ! Number of intersection curves and lines.
    !f2py integer*1, dimension(1000)                :: rep0
    type(repere3d), intent(out)                     :: rep0         ! Inertial frame.
    !f2py integer*1, dimension(1000)                :: InputData
    type(InputDataStruct),intent(inout)             :: InputData    ! Input data.
    integer                                         :: ierror       ! Error flag.
    
    real(rp),dimension(3)                           :: O            ! Origine of the inertial frame.
    !f2py integer*1, dimension(1000)                :: VX,VY,VZ
    type(vector)                                    :: VX,VY,VZ     ! Unit vector of the basis of the inertial frame.
    integer                                         :: nline        ! Number of lines in both the floater and the doamin.
    integer                                         :: nface_old    ! Number of faces for the previous body (in order to compute nface_vect.
    integer                                         :: jj           ! Loop parameter.
    logical                                         :: AboveFS      ! All parts of the geometry are above the free surface (True) or not (False).
        
    ! This subroutine generates the geometry of both the floater and the domain and the intersection line.
    

    
    if(Mesh_type.eq.2)then ! Method of Camille: floater immerged partially.
        
        ! Allocating of fgeom_vect
        call init_GeomVect(fgeom_vect)
                
        ! Creating of the inertial frame
        O   = [0._rp,0._rp,0._rp]
        VX%coord  = [1._rp,0._rp,0._rp] ! Previously only VX
        VY%coord  = [0._rp,1._rp,0._rp] ! Previously only VY
        VZ%coord  = [0._rp,0._rp,1._rp] ! Previously only VZ
        call assign_repere(1,O,VX%coord,VY%coord,VZ%coord,0._rp,0._rp,rep0)
        
        ! Initalization
        nface = 0 ; nline = 0
        HouleRF%index = 1       ! HouleRF%Index fixed equal to 1.
        
        ! Creating of the geomtry of the domain
        call create_geom(fdomaine,nface,nline,0,InputData,0) ! 0 is useless.
        fgeom_vect%nface_vect(1) = nface
        
        ! Creating of the geometry of the floaters
        if(is_body)then
            do jj = 1,NBodies
                nface_old = nface ! Saving the previous value of nface.

                call create_geom(fgeom_vect%geom(jj),nface,nline,1,InputData,jj)

                fgeom_vect%nface_vect(jj+1) = nface - nface_old ! Number of faces for this geometry.
                
                
                call update_geom(rep0,fgeom_vect%geom(jj),InputData%Position(:,2,jj),InputData%Position(:,1,jj))
                
                
                !if (lineaireBody) then
                !    call update_geom(rep0,fgeom_vect%geom(jj),InputData%Position(:,2,jj),[0._rp,0._rp,0._rp]) !à vérifier (1111)
                !else
                !    call update_geom(rep0,fgeom_vect%geom(jj),InputData%Position(:,2,jj),InputData%Position(:,1,jj))
                !end if

            end do
        end if

        ! Intersection curves
        if(is_body)then
            call compute_intersection(t0,fgeom_vect,tab2,n_tab2,5,ierror,InputData,InputData%dx2(1),n_tab)
        else
            ierror = 0 ! Initialization of the error flag
            n_tab2 = 0 ! No body, no intersection curve.
        end if

        
        ! Flag to know if at least one floater pierces the free surface.
        is_immerged = is_body .and. n_tab2==0
        
        ! Update fgeom_vect%Active.
        do jj = 1,NBodies 
            if(RemeshFS)then
                call isGeom_aboveFS(fgeom_vect%geom(jj),t0,AboveFS)
            else
                AboveFS = .false.
            end if
                        
            if(not(AboveFS))then
                fgeom_vect%Active(jj) = .true. ! Body immerged or piercing.
            else
                fgeom_vect%Active(jj) = .false. ! Body above the free surface.
            end if
        end do
        
    else if(Mesh_type.eq.1)then ! Method of Lucas: floater immerged totally or Method of PYW.
        
        ierror = 0 ! Initialization of the error flag
        
        is_immerged = .true.
        
        ! Allocating of fgeom_vect
        call init_GeomVect(fgeom_vect)
        
        ! Update fgeom_vect%Active.
        do jj = 1,NBodies
            fgeom_vect%Active(jj) = .true.
        end do
        
    end if
    
    ! Flag to known if the free surface mesh needs to be remeshed or not.
    DeformFS    =   not(lineaireFS)     .or.    lineaireFS.and.not(lineaireBody).and. not(is_immerged)
    

    
end subroutine Generation_Geometry
    
    
subroutine create_geom(geom,nface,nline,iflag,InputData,NumBody)
  
    !f2py integer*1, dimension(1000)    :: geom     
    type(type_geom),intent(inout)       :: geom         ! Geometry.
    integer,intent(inout)               :: nface,nline  ! Number of faces and lines in the geometry.
    integer,intent(in)                  :: iflag        ! Flag.
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(in)    :: InputData    ! Input data.
    integer,intent(in)                  :: NumBody      ! Body number.
    

    !f2py integer*1, dimension(1000)    :: O            
    type(point)                         :: O            ! Origin.
    !f2py integer*1, dimension(1000)    :: VX,VY,VZ     
    type(vector)                        :: VX,VY,VZ     ! Unit axis of the frame.
    integer                             :: ierror       ! Error flag.
    !f2py integer*1, dimension(1000)    :: C
    type(point)                         :: C            !
    !f2py integer*1, dimension(1000)    :: rep
    type(repere3d)                      :: rep          ! Frame of the floater or the tank.
  
    
  
        
    ! This subroutine creates a geometry.
    
    ierror = 0
    
    O  = [0._rp,0._rp,0._rp]
    VX%coord  = [1._rp,0._rp,0._rp] ! Previously only VX
    VY%coord  = [0._rp,1._rp,0._rp] ! Previously only VY
    VZ%coord  = [0._rp,0._rp,1._rp] ! Previously only VZ
      
    select case(iflag)
      
    ! Geometry of the body.
    case(1)
    
    ! Center of the body.
    C = [0._rp,0._rp,0._rp]
    call assign_repere(1,C%coord,VX%coord,VY%coord,VZ%coord,0._rp,0._rp,rep)
        
    ! Bodies.
    if(InputData%igtype(NumBody)==1)then ! Cube.
        call cube_from_length(rep,InputData%Lgeom(1,NumBody),nface,geom)
        geom%cmd_pl(1:6)=1
    elseif(InputData%igtype(NumBody)==2)then ! Cylinder.
        call cylindre_from_dist(rep,InputData%Lgeom(1,NumBody),InputData%Lgeom(2,NumBody),nface,nline,1,iflag,geom)
        geom%cmd_cyl(1)=1
        geom%cmd_dis(1:2)=[1,1]
        if(abs(geom%disque(2)%repere%origine(3)+Ldom(3)).lt.Epsilon .and. cuve_ferme)then
            geom%cmd_dis(2)=0
        endif
    elseif(InputData%igtype(NumBody)==3)then ! Useless.
        call geom_wavestar0(rep,nface,nline,geom)
        geom%cmd_sph(1) = 1
        geom%cmd_cone(1:2) = 1
        geom%cmd_dis(1) = 1
    elseif(InputData%igtype(NumBody)==4)then ! Useless.
        call wavestar1(InputData%file_axisym(NumBody),rep,nface,nline,ierror,geom)
        geom%cmd_axi(1) = 1
    elseif(InputData%igtype(NumBody)==5)then ! Axisym.
        call wavestar2(InputData%file_axisym(NumBody),rep,nface,nline,ierror,geom)
        geom%cmd_axi(:) = 1
    elseif(InputData%igtype(NumBody)==6)then ! Wigley hull.
        call WigleyGeom(rep,InputData%Lgeom(1:3,NumBody),nface,nline,ierror,geom)
        geom%cmd_wigley(:) = 1
    elseif(InputData%igtype(NumBody) == 7)then ! Sphere.
        call Sphere_from_radius(rep,InputData%Lgeom(1,NumBody),nface,nline,geom)
        geom%cmd_sph(:) = 1
    else
        ierror = 100
        goto 9999
    endif
    
    ! Geometry of the domain.
    case(0) 
    
    ! Center of the domain.
    C = [0._rp,0._rp,-Ldom(3)/2._rp]
    call assign_repere(1,C%coord,VX%coord,VY%coord,VZ%coord,0._rp,0._rp,rep)
    
    if(idtype==1)then
        ! Cube.
        call cube_from_length(rep,Ldom(2),nface,geom)
        geom%cmd_pl(1:6)=1
    elseif(idtype==2)then
        ! Cylinder.
        call cylindre_from_dist(rep,Ldom(3),Ldom(5),nface,nline,-1,iflag,geom)
        geom%cmd_cyl(1)=1
        geom%cmd_dis(1:2)=1
        if(bottom_sym)then
            geom%cmd_dis(2) = 0
        endif
    elseif(idtype==3)then
        ! Rectangle.
        call rectangle_from_length(rep,Ldom(1),Ldom(2),Ldom(3),nface,cuve_ferme,geom)
        geom%cmd_pl = 1
        if(bottom_sym)then
            geom%cmd_pl(2) = 0
        endif
    elseif(idtype==4)then
        ! ?
        C = [0._rp,0._rp,0._rp]
        call assign_repere(1,C%coord,VX%coord,VY%coord,VZ%coord,0._rp,0._rp,rep)
        call disque_from_dist(rep,Ldom(5),nface,-1,geom)
        geom%cmd_dis(1) = 1
    else
        ierror = 110
        goto 9999
    endif  
    
   end select           
  
    9999 continue
      if(ierror/=0)then
        write(*,99),ierror,idtype
      endif
    99 format('** error #',i3,' : type de geometrie ',i3,' non correct') 
    

end subroutine create_geom

subroutine isGeom_aboveFS(geom,t,AboveFS)
    
    !f2py integer*1, dimension(1000)    :: geom
    type(type_geom),intent(in)          :: geom             ! Geometry.
    real(rp),intent(in)                 :: t                ! Current time.
    logical,intent(out)                 :: AboveFS          ! All parts of the geometry are above the free surface (True) or not (False).
    
    integer                             :: j,k              ! Loop parameters.
    real(rp)                            :: Eta0             ! Free surface elevation.
    real(rp),dimension(3)               :: Vect_1,Vect_2    ! Vectors.
    integer                             :: nPts             ! Number of points below the free surface.
        
    ! This subroutine tests if a geometry is above the free surface or not.
    
    ! Initialization
    AboveFS = .true.
    
    ! Edges
    do j = 1,geom%narete
        ! P1
        call CEta0(geom%arete(j)%P1%coord, t, Eta0)
        if(geom%arete(j)%P1%coord(3) .le. Eta0)then
            AboveFS = .false.
            go to 9999
        end if
        
        ! P2
        call CEta0(geom%arete(j)%P2%coord, t, Eta0)
        if(geom%arete(j)%P2%coord(3) .le. Eta0)then
            AboveFS = .false.
            go to 9999
        end if
    end do
    
    ! Circles
    do j = 1,geom%ncercle
        ! Origin of the frame.
        call CEta0(geom%cercle(j)%repere%origine, t, Eta0)        
        if(geom%cercle(j)%repere%origine(3) .le. Eta0)then
            AboveFS = .false.
            go to 9999
        end if
    end do
    
    ! Cylinders
    do j = 1,geom%ncylindre
        print*,"isGeom_aboveFS: cylinders are not taken into account."
    end do
    
    ! Discs
    do j = 1,geom%ndisque
        print*,"isGeom_aboveFS: discs are not taken into account."
    end do
    
    ! Plans
    do j = 1,geom%nplan
        print*,"isGeom_aboveFS: plans are not taken into account."
    end do
    
    ! Cones
    do j = 1,geom%ncone
        print*,"isGeom_aboveFS: cones are not taken into account."
    end do
    
    ! Spheres
    do j = 1,geom%nsphere
        call Computation_eq_sph3D(0._rp,geom%sphere(j)%vmax,geom%sphere(j),Vect_1)
        call CEta0(Vect_1, t, Eta0)
        if(Vect_1(3) .le. Eta0)then
            AboveFS = .false.
            go to 9999
        end if
    end do
    
    ! Axisym
    nPts = 0 ! Not in the loop over j because an axisym body is made of two parts. If there is one point in each part, the body is not above the FS anymore.
    do j = 1,geom%naxisym
        do k = 1,geom%axisym(j)%npoint
            ! Local to cartesian coordinates of a slice of the axisym geometry.
            Vect_1(1) = geom%axisym(j)%P(k,1) ! x
            Vect_1(2) = 0._RP                 ! y
            Vect_1(3) = geom%axisym(j)%P(k,2) ! z
            
            ! Should be better but sometime false.
            ! call loc2cart(Vect_1,geom%axisym(j)%repere,Vect_2)
            Vect_2 = Vect_1 + geom%axisym(j)%repere%origine ! Not very good.
                                    
            call CEta0(Vect_2, t, Eta0)
            
            if(Vect_2(3) .le. Eta0)then
                nPts = nPts + 1
                ! If at least two points are below the free surface, AboveFS = .false.
                ! If only one point is below the free surface, that can leads to problems in the remeshing process.
                if( (not(Symmetry) .and. nPts.ge.2) .or. (Symmetry .and. nPts.ge.1))then ! Only one point necessary in case of Symmetry.
                    AboveFS = .false.
                    go to 9999
                end if
            end if
        end do
    end do
    
    ! Polyline
    do j = 1,geom%npolyline
        print*,"isGeom_aboveFS: polylines are not taken into account."
    end do
    
    9999 continue
    
end subroutine isGeom_aboveFS

end module GeomDef    
