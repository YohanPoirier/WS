!dec$ attributes dllexport::Add
subroutine Add(a,b,c)
    implicit none

    real(8), intent(in)  :: a,b

    real(8), intent(out) :: c 

    c = a+ b
    print*,"Bonjour"

end subroutine Add

!dec$ attributes dllexport::Execution_API
subroutine Execution_API()
    use Parameters_API
    use Exec
    implicit none

    ! Reading of the input files.
    call Execution

end subroutine Execution_API
    
!dec$ attributes dllexport::Generation_Geometry_API
subroutine Generation_Geometry_API()
    use Parameters_API
    use geomDef
    implicit none

    ! This subroutine generates the geometry of both the floater and the domain and the intersection line.
    call Generation_Geometry(fgeom,fdomaine,nface,tab2,n_tab2,rep0,ierror)

end subroutine Generation_Geometry_API   

!dec$ attributes dllexport::Generation_Mesh_API
subroutine Generation_Mesh_API()
    use Parameters_API
    use MeshGen
    implicit none

    ! This subroutine generates mesh of the domain, the free surface and the floater.
    call Generation_Mesh(Maillage,fdomaine,fgeom,nface,rep0,mesh,ierror,tab2,n_tab2)

end subroutine Generation_Mesh_API

!dec$ attributes dllexport::BoucleTemporelle_RK4_API
subroutine BoucleTemporelle_RK4_API()
    use Parameters_API
    use BoucleTemp
    implicit none
    
    ! This subroutine is the temporal loop. It updates the state vector at each time step.
    call BoucleTemporelle_RK4(nt, Maillage, fgeom, nface, rep0, mesh)
    
end subroutine BoucleTemporelle_RK4_API