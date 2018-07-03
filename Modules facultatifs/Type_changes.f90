module type_changes
use GeomStruct
    
implicit none

contains

subroutine repere3d2Tables(repere,index,origine,e1,e2,e3,phi,gamma,psi)

    type(repere3d),intent(in)           :: repere   ! Geometry of the floater and the domain
    integer,intent(out)                 :: index    ! Index of the frame.
    real(rp),dimension(3),intent(out)   :: origine  ! Origine of the frame.
    real(rp),dimension(3),intent(out)   :: e1       ! 1st unit vector of the frame.
    real(rp),dimension(3),intent(out)   :: e2       ! 2nd unit vector of the frame.
    real(rp),dimension(3),intent(out)   :: e3       ! 3rd unit vector of the frame.
    real(rp),intent(out)                :: phi      ! Rotation around e3
    real(rp),intent(out)                :: gamma    ! Rotation around e2
    real(rp),intent(out)                :: psi      ! Rotation around e1


    
    ! This subroutine translates the type repere3d into no derived type structures.
    
    ! Index
    index = repere%index
    
    ! Origine
    origine = repere%origine
    
    ! e1
    e1 = repere%e1
    
    ! e1
    e2 = repere%e2
    
    ! e3
    e3 = repere%e3
    
    ! Phi
    phi = repere%phi
    
    ! Gamma
    gamma = repere%gamma
    
    ! Psi
    psi = repere%psi
        
end subroutine repere3d2Tables
   
subroutine fgeom2Tables(fgeom,index)

    type(type_geom),intent(in)  :: fgeom    ! Geometry of the floater and the domain
    integer,intent(out)         :: index    ! Index of the geometry

    ! This subroutine translates the type type_geom into no derived type structures.
    
    ! Index
    index = fgeom%index

end subroutine fgeom2Tables
    
end module type_changes