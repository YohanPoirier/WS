module Test_Lissage
    use Constantes
    use Parameters
    use Incident_mod
    use Spline
    use Structuresdonnees
    use FonctionsCommunes
    implicit none
    
    contains
    
    subroutine TestLissage(Mesh)
    ! Parameters
    type(TMaillage) :: Mesh
    ! Local
    integer :: j, Ncomp, ierror
    real(rp) :: w1, w2, k1, k2
    real(rp), allocatable :: Energies(:,:), Omega(:)
    type(TEcoulement) :: Ecoulement
    type(THoule) :: HouleObj, HouleInit
    ! Begin
    call NewEcoulement(Ecoulement, Mesh%Nnoeud)
    call IniEcoulement(Ecoulement, Mesh%Nnoeud, 0._RP)
    ! Initialisation of the solution
    Ncomp = 50
    !!! A modifier pour virer les variables globales
    if (Mesh_Type .eq. 1 ) k1 = 4._RP*Pi/Mesh%FS%DimFS(1)
    if (Mesh_Type .eq. 2 .and. idtype.eq.1) k1 = 4._RP*Pi/Ldom(1)
    if (Mesh_Type .eq. 2 .and. idtype.eq.2) k1 = 4._RP*Pi/Ldom(5)
    if (Mesh_Type .eq. 1 ) k2 = 2._RP*pi/(3._RP*Ldom(5))
    if (Mesh_Type .eq. 2 ) k2 = 2._RP*pi/(3._RP*dx2)
    w1 = sqrt(g*k1*tanh(k1*profondeur))
    w2 = sqrt(g*k2*tanh(k2*profondeur))
    print*, 'Lissage : ', w1,' < w < ',w2
    print*, 2._RP*Pi/k2,' < lambda < ', 2*Pi/k1
    !w1 = 0.5_RP; w2 = 4.5_RP
    allocate(Energies(4,Ncomp), Omega(Ncomp))
    do j=1,ncomp
        HouleObj%Htype = 1;      HouleObj%nhoule = 1
        allocate(HouleObj%w(HouleObj%nhoule), HouleObj%dir(HouleObj%nhoule), HouleObj%Aphi(HouleObj%nhoule), HouleObj%konde(HouleObj%nhoule), HouleObj%lambda(HouleObj%nhoule))
        !do j=1,HouleObj%nhoule
        !    HouleObj%w(j) = j;     HouleObj%dir(j) = 90._RP/pi*0._RP;     HouleObj%Aphi(j) = 0.01_RP/j
        !end do
        HouleObj%w(1) = w1+(w2-w1)*(j-1)/(Ncomp-1);     HouleObj%dir(1) = 90._RP/pi*0._RP;     HouleObj%Aphi(1) = 0.01_RP
        !HouleObj%w(2) = 4.5_RP;     HouleObj%dir(2) = 90._RP/pi*1._RP;     HouleObj%Aphi(2) = 0.01_RP
    
        call find_konde(HouleObj%konde,HouleObj%w,profondeur,HouleObj%NHoule)
        HouleObj%lambda(1:NHoule) = 2._rp*Pi/HouleObj%konde(1:NHoule)
        Omega(j) = HouleObj%konde(1)
    
        call copyHoule(HouleObj, HouleInit)
        call Incident(0._RP, Mesh, Ecoulement)
        call CopyBackHoule(HouleInit)
        Ecoulement%Eta(1:Mesh%Nnoeud)%Perturbation = Ecoulement%Eta(1:Mesh%Nnoeud)%Incident
        !  Lissage
        call LissageGaussien(Mesh,Ecoulement,ierror)
    
        call CheckEnergy(Mesh,Ecoulement,Energies(:,j))
        deallocate(HouleObj%w, HouleObj%dir, HouleObj%Aphi, HouleObj%konde, HouleObj%lambda)
    end do
    ! Plot et verification
    if (Ncomp .eq. 1) call PlotFSLissage(Mesh, Ecoulement)
    if (Ncomp .ge. 1) call PlotEnergyLissage(Ncomp, Omega, Energies)
    deallocate(Energies, Omega)
    ! End
    end subroutine TestLissage
    
    subroutine PlotFSLissage(Mesh, Ecoulement)
    ! Parameters
    type(TMaillage) :: Mesh
    type(TEcoulement) :: Ecoulement
    ! Locals
    integer :: j, i1, i2, ios
    ! Begin
    i1 = Mesh%FS%IndFS(3)
    i2 = Mesh%FS%IndFS(4)
    ! Plot
    open(unit=3,file='TestLissage_FS_'//filename, iostat=ios)
    if (ios/=0) stop "Erreur à l'ouverture du fichier de lissage"
    write(3,fmt='(50a)') 'Title= "Lissage sur la SL"'
    write(3,fmt='(50a)') 'VARIABLES = "X","Y","Z","EtaLisse","Eta0","Nvoisin" '
    write(3,fmt='(a,i,a,i,a)') 'Zone N =', i1, ', E=', i2, ' , ET=TRIANGLE, F=FEPOINT' 
    do j=1,i1
        write(3,'(5E,I)') Mesh%Tnoeud(j)%Pnoeud, Ecoulement%Eta(j)%perturbation, Ecoulement%Eta(j)%incident, Mesh%Tnoeud(j)%Nvoisin(2)
    end do
    do j=1,i2
        write(3,'(3I)') Mesh%Tfacette(j)%TNoeud
    end do 
    close(unit=3)
    ! End
    end subroutine PlotFSLissage
    
    subroutine PlotEnergyLissage(Ncomp, k, Energies)
    ! Parameters
    integer, intent(in) :: Ncomp
    real(rp), dimension(:), intent(in) :: k
    real(rp), dimension(:,:), intent(in) :: Energies
    ! Locals
    integer :: j, ios
    real(rp) :: dx, Bec, Bep
    ! Begin
    open(unit=3,file='TestLissage_NRJ_'//filename, iostat=ios)
    if (ios/=0) stop "Erreur à l'ouverture du fichier de lissage"
    write(3,fmt='(50a)') 'Title= "Lissage sur la SL"'
    write(3,fmt='(50a)') 'VARIABLES = "k","B_Ep","Ep","Ep0","B_Ec","Ec","Ec0"'
    if (Mesh_Type .eq. 1) then
        dx = Ldom(5)
    elseif (Mesh_Type .eq. 2) then
        dx = dx2
    end if
    
    do j=1,Ncomp
        Bep = 0._RP
        if (abs(Energies(4,j)).ge.0) Bep = Energies(2,j)/Energies(4,j)
        Bec = 0._RP
        if (abs(Energies(3,j)).ge.0) Bec = Energies(1,j)/Energies(3,j)
        write(3,'(7E)') 2._RP*PI/k(j), Bep, Energies(2,j), Energies(4,j), Bec, Energies(1,j), Energies(3,j)
    end do
    close(unit=3)
    ! End
    end subroutine PlotEnergyLissage
    
    subroutine CheckEnergy(Mesh, Ecoulement, Energies)
    ! Parameters
    type(TMaillage), intent(in) :: Mesh
    type(TEcoulement), intent(in) :: Ecoulement
    real(rp), dimension(:), intent(out) :: Energies
    ! Locals
    integer :: i1, i2
    real(rp) :: Ec, Ep, Ec0, Ep0
    real(rp), allocatable :: DPhiDn(:), Phi(:), Eta(:)
    ! Begin
    i1 = Mesh%FS%IndFS(1)
    i2 = Mesh%FS%IndFS(3)
    ! Allocate
    allocate(Eta(i1:i2), Phi(i1:i2), DPhiDn(i1:i2))
    ! Assign
    Phi(i1:i2) = Ecoulement%Phi(i1:i2)%incident
    DPhiDn(i1:i2) = Ecoulement%DPhiDn(i1:i2)%incident
    ! Energy calculation for the filtered solution
    Eta(i1:i2) = Ecoulement%Eta(i1:i2)%perturbation
    call EnergySL(i1,i2,Mesh,Eta,Phi,DPhiDn,Ec,Ep)
    ! Energy calculation for the original solution
    Eta(i1:i2) = Ecoulement%Eta(i1:i2)%incident
    call EnergySL(i1,i2,Mesh,Eta,Phi,DPhiDn,Ec0,Ep0)
    ! Deallocate
    deallocate(Eta,Phi,DPhiDn)
    Energies(:) = [Ec, Ep, Ec0, Ep0]
    ! End
    end subroutine CheckEnergy
    
    subroutine EnergySL(i1,i2,Mesh,Eta,Phi,DPhiDn,Ec,Ep)
    ! Parameters
    integer, intent(in) :: i1, i2                           ! Bornes d'intégration
    type(TMaillage), intent(in) :: Mesh                     ! Maillage
    real(rp), dimension(:), intent(in) :: Eta, Phi, DPhiDn  ! Grandeur à intégrer
    real(rp), intent(out) :: Ec, Ep                         ! Sorties
    ! Locals
    integer :: j
    real(rp) :: Dn
    real(rp), dimension(3) :: M
    ! Begin
    Ec = 0._RP      ! Kinetic energy
    Ep = 0._RP      ! Potential energy
    do j=i1,i2
        M = Mesh%Tnoeud(j)%Pnoeud     
        Dn = Mesh%Tnoeud(j)%Aire
        Ec = Ec - Phi(j)*DPhiDn(j)*Dn
        Ep = Ep + Eta(j)**2*Dn
    end do
    ! End
    end subroutine EnergySL
    end module Test_Lissage