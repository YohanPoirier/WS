module Parameters_API
    use Structuresdonnees
    use GeomStruct
    use MeshStruct
    
    ! This module defines all the (global) variables which are used in the API.
    
    type(TMaillage)                         :: Maillage         ! Final mesh
    integer                                 :: nface,nline      ! Number of faces and lines in both the floater and the domain geometries
    integer                                 :: ierror           ! Error flag
    type(type_geom)                         :: fgeom, fdomaine  ! Geometry of the floater and the domain
    type(chaine_point_pt),dimension(100)    :: tab2             ! Table of intersection points ! Why 100? (PYW)
    integer                                 :: n_tab2           ! Number of intersection points
    type(repere3d)                          :: rep0             ! Inertial frame
    type(MGrid)                             :: mesh             ! Transitional mesh
    
end module Parameters_API