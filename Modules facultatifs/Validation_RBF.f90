module Validation_RBF
    use Constantes
    use Parameters
    use FonctionsCommunes
    use Preplot
    use GeomStruct  
    use GeomFonct
    use GenMaillage
    use BodyMotion_mod
    use rbf_interp
    
    implicit none
    
contains 
    
! ###########################################################
subroutine TestValidation_RBF(Maillage,phi,fgeom,fdoss)

! Parameters
    type(TMaillage),intent(inout) :: Maillage
    external phi
    type(type_geom),intent(inout) :: fgeom
    character(len=50),intent(in) :: fdoss
! local
    type(TMaillage) :: Maillage0
    real(rp), dimension(:,:), allocatable :: xd, xi, fd, fi
    integer, dimension(:), allocatable :: tab_corres, tab_data, tab_data2
    real(rp),dimension(:), allocatable :: tab_theta, tab_dist, rdist
    integer :: nd, ni, nc, NFS1, NFS2, NFS, nd2
    real(rp) :: r0, rd, a1, b1, r1, jr
    real(rp) :: dlxmax, dlymax, dtheta, theta, Va
    real(rp) :: ti
    real(rp),dimension(3) :: y, yp, yp1, vx, vy, vz, O, P, Q
    integer :: j, k, jout, jn, j2, jmax
    real(rp),dimension(:),allocatable :: w
    real(rp),dimension(:,:),allocatable :: dlx,dly
    logical :: bool0,bool
    character(len=50) :: filemaill
    type(repere3d) :: repere
    real(rp),dimension(Maillage%NBody,6) :: VBody
    integer, parameter :: icase = 2
    integer, parameter :: itest = 2
    real(rp) ,parameter :: periode = 1.
    !!integer , parameter :: nmax = 70
    integer, parameter :: iout1=35, iout2=38
    integer, parameter :: icheckmesh = 37
    logical, parameter :: second_order = .true.
    integer :: nmax
    integer :: nbody, nfacette
    integer, dimension(:), allocatable :: tab_data_body
    real(rp), dimension(:), allocatable :: rdist_body, tab_dist_body, taux0, area0
    real(rp), dimension(:,:), allocatable :: TMetric
    real(rp), dimension(:), allocatable :: fshape0
    integer :: NBD1, NBD2, NBD
    real(rp) :: calf, AlfBody(3)

    print*,'Validation deplacement maillage'

    bool0 = .true.

    NFS1 = Maillage%FS%IndFS(1)
    NFS2 = Maillage%FS%IndFS(3)
    NFS = NFS2-NFS1+1

    nmax = int(1./dt)

    allocate(xd(2,NFS))
    allocate(xi(2,NFS))
    allocate(fd(2,NFS))
    allocate(fi(2,NFS))
    allocate(tab_corres(NFS))
    allocate(tab_theta(NFS))
    allocate(tab_data(NFS))

    allocate(TMetric(Maillage%Nnoeud,2))

    call NewMaillage(Maillage0,Maillage%Nnoeud)
    call CopyMaillage(Maillage0,Maillage)

    !nfacette = Maillage%FS%IndFS(4)-Maillage%FS%IndFS(2)+1
    nfacette = Maillage%Nfacette
    allocate(taux0(nfacette))
    allocate(area0(nfacette))
    allocate(fshape0(nfacette))

    allocate(dlx(NFS,nmax),dly(NFS,nmax))

    fi(1:2,1:NFS) = 0._rp
    dlx(1:NFS,1:nmax) = 0._rp
    dly(1:NFS,1:nmax) = 0._rp

    dlxmax = 1.5*0.361_rp
    dlymax = 0._rp
    
    dtheta = 2._rp*PI/dble(nmax) ! utilis\E9 dans le cas de deplacement en translation
    !dtheta = 0.5_rp*PI/dble(nmax) ! dans le cas de deplacement en rotation
    
    a1 = 0.5_rp*PI
    b1 = 2._rp*PI/dble(nmax)
    r1 = 0.361003861003861_rp
    r0 = 0.001_rp

    O  = [0._rp , 0._rp , 0._rp]
    vx = [1._rp , 0._rp , 0._rp]
    vy = [0._rp , 1._rp , 0._rp]
    vz = [0._rp , 0._rp , 1._rp]
    repere = assign_repere(1, O, vx, vy, vz, 0._rp, 0._rp)
! --------------------------------------------------------------------------
!   Recherche des points de contr\F4les
! --------------------------------------------------------------------------
    nd=0 ; ni=0
    do j=NFS1,NFS2
        rd = sqrt(Maillage%Tnoeud(j)%Pnoeud(1)**2 + Maillage%Tnoeud(j)%Pnoeud(2)**2)
        if(Maillage%Tnoeud(j)%Ndouble /= 0)then
            nd=nd+1
            xd(1:2,nd) = Maillage%Tnoeud(j)%Pnoeud(1:2)
            tab_data(nd) = j
        elseif(rd.gt.Ldom(5)-0.1)then
            nd=nd+1
            xd(1:2,nd) = Maillage%Tnoeud(j)%Pnoeud(1:2)
            tab_data(nd) = j
        else
            ni=ni+1
            xi(1:2,ni) = Maillage%Tnoeud(j)%Pnoeud(1:2)
            tab_corres(ni) = j
        endif
    enddo

    do j=1,nd
        y = cart2sph([xd(1,j),xd(2,j),0._rp],repere)
        tab_theta(j) = y(2)
    enddo

!   ## A utiliser pour les test de conservation de la g\E9om\E9trie

    open(unit=iout1,file=trim(fdoss)//"intersection_quality.dat")
    open(unit=iout2,file=trim(fdoss)//"immerged_surf_quality.dat")
    call CheckGeom(Maillage,fgeom,t0,iout1,iout2)
    
    !filemaill=trim(fdoss)//'MaillageInit2.dat'
    !open(unit=3,file=filemaill)
    !write(3,*) 'Title = "Maillage Evolution"'
    !write(3,'(a300)') 'VARIABLES = "X","Y","Z","Nx","Ny","Nz","IX","IY","IZ","JX","JY","JZ","typeFrontiere","Npanneau","Damping","NVoisin","Ndouble","VX","VY","VZ"'

    call PrePlots

    call PlotMaill(filemaill, Maillage,t0)


!    print*,"** Test de conservation de la geometrie actif"
!    allocate(tab_dist(nd))
!    allocate(rdist(nd))
!    allocate(tab_data2(nd))
!    nd2 = 0
!    do j=1,nd
!        j2 = tab_data(j)
!        if(Maillage%Tnoeud(j2)%Ndouble /=0)then
!            nd2 = nd2+1
!            tab_dist(nd2) = sqrt(xd(1,j)**2 + xd(2,j)**2)
!            tab_data2(nd2) = j2
!        endif
!    enddo
!    open(unit=32,file="intersection_quality.dat")

!   Crontrole conservation de la surface immergee
!    NBD1 = Maillage%Body(1)%IndBody(1)
!    NBD2 = Maillage%Body(1)%IndBody(3)
!    NBD = NBD2-NBD1
!    allocate(tab_dist_body(NBD), rdist_body(NBD), tab_data_body(NBD))
!    nbody = 0
!    do j=NBD1,NBD2
!        bool = .false.
!        do jn =1,Maillage%Tnoeud(j)%Ndouble
!            j2 = Maillage%Tnoeud(j)%double(jn)
!            bool = bool .or. Maillage%Tnoeud(j2)%typeNoeud == 0
!        enddo
!        if(.not.bool)then
!            nbody = nbody+1
!            tab_dist_body(nbody) = norm2(Maillage%Tnoeud(j)%Pnoeud)
!            tab_data_body(nbody) = j
!        endif
!    enddo
!    open(unit=35,file="immerged_surf_quality.dat")
!   ##

    print*," number of data point          : nd = ",nd
    print*," number of interpolation point : ni = ",ni 
    print*," number of iteration           : nmax = ",nmax

    allocate(w(nd+3))

! ----------------------------------------------------------------------------
! Calcul qualite du maillage initial
! ----------------------------------------------------------------------------

    open(unit=icheckmesh,file=trim(fdoss)//"MeshQuality.dat")
    call CheckMeshInit(Maillage,taux0,area0,fshape0)

! -----------------------------------------------------------------------------
! Definition des deplacements
! -----------------------------------------------------------------------------
! -----------------------------------------
    if(itest.eq.1)then
! -----------------------------------------

    select case(icase)
    case(1)
!   Deplacement en rotation
    do j=1,nd
        k = tab_data(j)
        if(Maillage%Tnoeud(k)%Ndouble /= 0)then
            theta = tab_theta(j)
            do jout=1,nmax
                jr = dble(jout)
                !dlx(j,jout) = r1*(cos(theta + a1*sin((jr+1._rp)*b1)) - cos(theta + a1*sin(jr*b1)))
                !dly(j,jout) = r1*(sin(theta + a1*sin((jr+1._rp)*b1)) - sin(theta + a1*sin(jr*b1)))
                dlx(j,jout) = r1*dcos(theta + (jr+1.)*dtheta) - r1*dcos(theta + jr*dtheta)
                dly(j,jout) = r1*dsin(theta + (jr+1.)*dtheta) - r1*dsin(theta + jr*dtheta)
            enddo
        else
            dlx(j,1:nmax) = 0._rp
            dly(j,1:nmax) = 0._rp
        endif
    enddo

    case(2)
!   Deplacement en translation
    do j=1,nd
        k = tab_data(j)
        if(Maillage%Tnoeud(k)%Ndouble /= 0)then
            do jout=1,nmax
                dlx(j,jout) = dlxmax*(sin(dtheta*(jout+1)) - sin(dtheta*jout))
                dly(j,jout) = dlymax*(sin(2.*dtheta*(jout+1)) - sin(2.*dtheta*jout))
            enddo
        else
            dlx(j,1:nmax) = 0._rp
            dly(j,1:nmax) = 0._rp
        endif
    enddo

    case(3)
!   Deplacement en rotation 2
    do j=1,nd
        k = tab_data(j)
        if(Maillage%Tnoeud(k)%Ndouble /= 0)then
            theta = tab_theta(j)
            r1 = sqrt(Maillage%Tnoeud(k)%Pnoeud(1)**2 + Maillage%Tnoeud(k)%Pnoeud(2)**2)
            do jout=1,nmax
                jr = dble(jout)
                yp  = [r1,theta+jr*dtheta,0.5_rp*PI]
                P = sph2cart(yp,repere)
                yp1 = [r1,theta+(jr+1._rp)*dtheta,0.5_rp*PI]
                Q = sph2cart(yp1,repere)
                dlx(j,jout) = Q(1)-P(1)
                dly(j,jout) = Q(2)-P(2)
            enddo
        else
            dlx(j,1:nmax) = 0._rp
            dly(j,1:nmax) = 0._rp
        endif
    enddo

    end select

! 
!   Boucle deplacement maillage
    do jout=1,nmax

        call rbf_weight(2, nd, xd, r0, phi, dlx(1:nd,jout), w)
        call rbf_funct(2, nd, xd, r0, phi, w, ni, xi, fi(1,1:ni))

        call rbf_weight(2, nd, xd, r0, phi, dly(1:nd,jout), w)
        call rbf_funct(2, nd, xd, r0, phi, w, ni, xi, fi(2,1:ni))

        do j=1,ni
            k = tab_corres(j)
            Maillage%Tnoeud(k)%Pnoeud(1:2) = xi(1:2,j)+fi(1:2,j)
        enddo 
!        do j=1,nd
!            k = tab_data(j)
!            if(Maillage%Tnoeud(k)%Ndouble /= 0)then
!                theta = tab_theta(j)
!                r1 = sqrt(Maillage%Tnoeud(k)%Pnoeud(1)**2 + Maillage%Tnoeud(k)%Pnoeud(2)**2)
!                jr = dble(jout)
!                yp1 = [r1,theta+(jr+1._rp)*dtheta,0.5_rp*PI]
!                Q = sph2cart(yp1,repere)
!                Maillage%Tnoeud(k)%Pnoeud(1:2) = Q(1:2)
!            endif
!        enddo
        do j=1,nd
           k = tab_data(j)
           Maillage%Tnoeud(k)%Pnoeud(1:2) = xd(1:2,j)+[dlx(j,jout),dly(j,jout)]
           do jn = 1,Maillage%Tnoeud(k)%Ndouble
                j2 = Maillage%Tnoeud(k)%double(jn)
                Maillage%Tnoeud(j2)%Pnoeud(1:2) = Maillage%Tnoeud(k)%Pnoeud(1:2)
           enddo
        enddo

        filemaill=trim(fdoss)//'MaillageInit2.dat'
        call PlotMaill(filemaill, Maillage,dble(jout))

        do j=1,nd
            k = tab_data(j)
            xd(1:2,j) = Maillage%Tnoeud(k)%Pnoeud(1:2)
        enddo
        do j=1,ni
            k = tab_corres(j)
            xi(1:2,j) = Maillage%Tnoeud(k)%Pnoeud(1:2)
        enddo

    enddo    
! ----------------------------------------------------------
    elseif(itest.eq.2)then
! ----------------------------------------------------------
    do nc=1,Maillage%NBody
        Maillage%Body(nc)%ABody(1:6) = 0._rp
        Maillage%Body(nc)%GBody(1:3,1:2) = 0._rp
        Maillage%Body(nc)%MBody(1:6) = 0._rp 
    enddo
    AlfBody(1:3) = 0._rp

    !jmax = int(periode)*nmax
    do jout = 1,nt

!       Vitesse de deplacement des points appartenant au flotteur
        !Va = dlxmax*(sin(dtheta*(jout+1))-sin(dtheta*jout))/dt
        !Va = 5._rp*0.361_rp/(dt*dble(nt))
        Va = 0._rp
        ti = t0+dble(jout-1)*dt
        !VBody = [0._rp , 0._rp , 0._rp , 0._rp , 8.*PI/4.*cos((jout-1.)*dtheta) , 0._rp]
        !VBody = [0.2*Va , 0._rp , 0._rp , 0.5*PI/4.*cos(4.*PI*(jout-1.)/dble(nmax)) , 0.5*PI/4.*cos(2.*PI*(jout-1.)/dble(nmax)) , 0._rp]
        !VBody = [0.5*Va , 0._rp , 0._rp , PI/4.*cos(2.*dtheta*(jout-1.)) , PI/4.*cos(dtheta*(jout-1.)) , PI/2.*cos(dtheta*(jout-1.))]
        !VBody = [0.5*Va , 0._rp , 0._rp , PI/4.*cos(2.*dtheta*(jout-1.)) , PI/4.*cos(dtheta*(jout-1.)) , 0._rp]
        !VBody = [0._rp , 0._rp , 0._rp , 0._rp, 0._rp, PI*cos((jout-1.)*dtheta)]
        !VBody = [0._rp , 0._rp , 0._rp , 0.5*PI/4.*cos(2.*dtheta*(jout-1.)) , 0.5*PI/4.*cos(dtheta*(jout-1.)) , 0._rp]
        !VBody = [0.5*Va , 0._rp , 0.d0 , 0._rp, 0._rp, 0._rp]
        !VBody = [0._rp , 0._rp , Va , 0._rp, 0._rp, 0._rp]
        VBody(Int_Body,1:6) = [0._rp , 0._rp , 0._rp , 0._rp, 0._rp, 0._rp]
        !VBody = [Va, Va, 0._rp, 0._rp,0._rp,PI/3./(dt*dble(nt))]
        AlfBody(1:3) = AlfBody(1:3) + VBody(Int_Body,4:6)*dt
        if(second_order)then
            call BodyVel_rk(Maillage,ti,dt,VBody,fgeom)
        else
            call BodyVel(Maillage,ti,dt,VBody)
        endif

         do j=1,nd
            k = tab_data(j)
            if(Maillage%Tnoeud(k)%Ndouble /=0 .and. norm2(Maillage%Tnoeud(k)%Pnoeud).lt.Ldom(5)-LAbs)then
                dlx(j,1) = Maillage%Tnoeud(k)%Velocity(1)*dt
                dly(j,1) = Maillage%Tnoeud(k)%Velocity(2)*dt
            else
                dlx(j,1) = 0._rp
                dly(j,1) = 0._rp
            endif
        enddo

!       Calcul des deplacements des noeuds de la surface libre
        call rbf_weight(2, nd, xd, r0, phi, dlx(1:nd,1), w)
        call rbf_funct(2, nd, xd, r0, phi, w, ni, xi, fi(1,1:ni))

        call rbf_weight(2, nd, xd, r0, phi, dly(1:nd,1), w)
        call rbf_funct(2, nd, xd, r0, phi, w, ni, xi, fi(2,1:ni))

        ! ## Test avec deplacement seulement le long du flotteur
        !fi(1:2,:) = 0._rp

!       Mise a jour vitesse des noeuds en fonction des deplacement
        do j=1,ni
            k = tab_corres(j)
            Maillage%Tnoeud(k)%Velocity(1:2) = fi(1:2,j)/dt
            !Maillage%Tnoeud(k)%Pnoeud(1:2) = Maillage%Tnoeud(k)%Pnoeud(1:2) + fi(1:2,j)
        enddo 

        ! ##Mise en commentaire car repositionnement trait\E9 dans Remesh
        !do j=1,nd
        !    k = tab_data(j)
        !    Maillage%Tnoeud(k)%Pnoeud(1:2) = xd(1:2,j) + [dlx(j,1),dly(j,1)]
        !    do jn = 1,Maillage%Tnoeud(k)%Ndouble
        !        j2 = Maillage%Tnoeud(k)%double(jn)
        !        Maillage%Tnoeud(j2)%Pnoeud(1:2) = Maillage%Tnoeud(k)%Pnoeud(1:2)
        !    enddo
        !enddo


!       Repositionnement des noeuds du maillage
        call CopyMaillage(Maillage0,Maillage)
        call Remesh(Maillage,Maillage0,ti,dt,fgeom)

        if(mod(jout,nout)==0)then
            call PlotMaill(filemaill, Maillage,ti+dt)
        endif

!       Mise a jour des tableaux (on suppose que l'on garde le mm nb de noeud)
        do j=1,nd
            k = tab_data(j)
            xd(1:2,j) = Maillage%Tnoeud(k)%Pnoeud(1:2)
        enddo
        do j=1,ni
            k = tab_corres(j)
            xi(1:2,j) = Maillage%Tnoeud(k)%Pnoeud(1:2)
        enddo

!       ## A utiliser pour les test de conservation de la g\E9om\E9trie
        call CheckGeom(Maillage,fgeom,ti+dt,iout1,iout2)
!       ##

!       ## Controle qualite maillage
        call CheckMesh(Maillage,taux0,area0,fshape0,ti+dt,TMetric,icheckmesh)

        if(iwmetric)then
            call PlotMetric(ti+dt,Maillage,TMetric)
        endif


    enddo
! ----------------------------------------------------------
    endif !(Fin test itest)
! ----------------------------------------------------------
    if(allocated(xd)) deallocate(xd)
    if(allocated(xi)) deallocate(xi)
    if(allocated(fd)) deallocate(fd)
    if(allocated(fi)) deallocate(fi)
    if(allocated(tab_corres)) deallocate(tab_corres)
    if(allocated(tab_data)) deallocate(tab_data)
    if(allocated(tab_theta)) deallocate(tab_theta)
    if(allocated(dlx)) deallocate(dlx)
    if(allocated(dly)) deallocate(dly)
    if(allocated(w)) deallocate(w)
    close(3)

    if(allocated(rdist)) deallocate(rdist)
    if(allocated(tab_dist)) deallocate(tab_dist)
    if(allocated(tab_data2)) deallocate(tab_data2)
    if(allocated(rdist_body)) deallocate(rdist_body)
    if(allocated(tab_dist_body)) deallocate(tab_dist_body)
    if(allocated(tab_data_body)) deallocate(tab_data_body)
    if(allocated(taux0)) deallocate(taux0)
    if(allocated(area0)) deallocate(area0)
    if(allocated(fshape0)) deallocate(fshape0)
    if(allocated(TMetric)) deallocate(TMetric)
    close(iout1)
    close(iout2)
    close(icheckmesh)
    close(iometric)

    call DelMaillage(Maillage0)
    

end subroutine TestValidation_RBF


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BodyVel_rk
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BodyVel_rk(Mesh,t,dtl,VBody,fgeom)
implicit none
type(TMaillage),intent(inout) :: Mesh
real(rp),intent(in) :: t
real(rp),intent(in) :: dtl
real(rp),dimension(:,:),intent(in) :: VBody
type(type_geom),intent(in) :: fgeom
! local
integer :: j, jn, j2, k, n, NBD1, NBD2, NBD
real(rp) :: h, h2
real(rp),dimension(:,:,:),allocatable :: RK
real(rp),dimension(6,4) :: RKVBody
type(TMaillage) :: MeshTemp
real(rp),dimension(4),parameter :: coeff_rk = [1./6.,1./3.,1./3.,1./6.]
!character(len=50),parameter :: filemaill = 'MaillageRK4.dat' 
logical :: bplot
integer :: iout

!bplot = .true.
!iout = 36

h = dtl
h2 = 0.5_rp*dtl

NBD1 = Mesh%Body(1)%IndBody(1)
NBD2 = Mesh%Body(1)%IndBody(3)
NBD = NBD2-NBD1+1

allocate(RK(3,NBD,4))

n = Mesh%Nnoeud
call NewMaillage(MeshTemp,n)

! -- RK1
call CopyMaillage(MeshTemp,Mesh)
call BodyVel(MeshTemp,t,h2,VBody)
!call PlotMaill(filemaill, MeshTemp,t+0.2*h,bplot,iout)

do j=1,NBD
    k = j+NBD1-1
    RK(1:3,j,1) = MeshTemp%Tnoeud(k)%Velocity
enddo
RKVBody(1:6,1) = MeshTemp%Body(Int_body)%VBody

! -- RK2
call Remesh(MeshTemp,Mesh,t+h2,h2)
call BodyVel(MeshTemp,t+h2,h2,VBody)
!call PlotMaill(filemaill, MeshTemp,t+0.4*h,bplot,iout)

do j=1,NBD
    k = j+NBD1-1
    RK(1:3,j,2) = MeshTemp%Tnoeud(k)%Velocity
enddo
RKVBody(1:6,2) = MeshTemp%Body(Int_body)%VBody

!call CopyMaillage(MeshTemp,Mesh)
!do j=1,NBD
!    k = j+NBD1-1
!    MeshTemp%Tnoeud(k)%Velocity = RK(1:3,j,2)
!    do jn=1,MeshTemp%Tnoeud(k)%Ndouble
!        j2 = MeshTemp%Tnoeud(k)%double(jn)
!        MeshTemp%Tnoeud(j2)%Velocity = RK(1:3,j,2)
!    enddo
!enddo
!MeshTemp%Body(Int_Body)%VBody = RKVBody(1:6,2)

call Remesh(MeshTemp,Mesh,t+h2,h2)

! -- RK3
call BodyVel(MeshTemp,t+h2,h2,VBody)
!call PlotMaill(filemaill, MeshTemp,t+0.6*h,bplot,iout)

do j=1,NBD
    k = j+NBD1-1
    RK(1:3,j,3) = MeshTemp%Tnoeud(k)%Velocity
enddo
RKVBody(1:6,3) = MeshTemp%Body(Int_Body)%VBody

!call CopyMaillage(MeshTemp,Mesh)
!do j=1,NBD
!    k = j+NBD1-1
!    MeshTemp%Tnoeud(k)%Velocity = RK(1:3,j,3)
!    do jn=1,MeshTemp%Tnoeud(k)%Ndouble
!        j2 = MeshTemp%Tnoeud(k)%double(jn)
!        MeshTemp%Tnoeud(j2)%Velocity = RK(1:3,j,3)
!    enddo
!enddo
!MeshTemp%Body(Int_Body)%VBody = RKVBody(1:6,3)

call Remesh(MeshTemp,Mesh,t+h,h)

! -- RK4
call BodyVel(MeshTemp,t+h,h,VBody)
!call PlotMaill(filemaill, MeshTemp,t+0.8*h,bplot,iout)

do j=1,NBD
    k = j+NBD1-1
    RK(1:3,j,4) = MeshTemp%Tnoeud(k)%Velocity
enddo
RKVBody(1:6,4) = MeshTemp%Body(Int_Body)%VBody

! -- N+1
do j=1,NBD
    k = j+NBD1-1
    Mesh%Tnoeud(k)%Velocity = coeff_rk(1)*RK(1:3,j,1) + coeff_rk(2)*RK(1:3,j,2) + &
& coeff_rk(3)*RK(1:3,j,3) + coeff_rk(4)*RK(1:3,j,4)
    do jn=1,Mesh%Tnoeud(k)%Ndouble
        j2 = Mesh%Tnoeud(k)%double(jn)
         Mesh%Tnoeud(j2)%Velocity = Mesh%Tnoeud(k)%Velocity
    enddo
enddo
Mesh%Body(Int_Body)%VBody = coeff_rk(1)*RKVBody(1:6,1) + coeff_rk(2)*RKVBody(1:6,2) + &
& coeff_rk(3)*RKVBody(1:6,3) + coeff_rk(4)*RKVBody(1:6,4)

deallocate(RK)
call DelMaillage(MeshTemp)
!close(iout)


end subroutine BodyVel_rk

subroutine TestValidation_Deformation(Mesh,fgeom)
! Parameters
    type(TMaillage),intent(inout) :: Mesh
    type(type_geom),intent(inout) :: fgeom
!   local
    integer :: j, jt, Nnoeud
    real(rp) :: ti, dt2, VBody(Mesh%NBody,6)
    real(rp),dimension(:),allocatable :: taux0, area0, fshape0
    real(rp),dimension(:,:),allocatable :: TMetric
    real(rp),dimension(:,:,:),allocatable :: RKVel
    type(TMaillage) :: MeshTemp 
    type(TEcoulement) :: Ecoulement
!   parameter
    real(rp),parameter :: denom = 1._rp/6._rp
 
    VBody = 0._rp
    Nnoeud = Mesh%Nnoeud
    Mesh%Body(Int_Body)%VBody = 0._rp
    Mesh%Body(Int_Body)%ABody = 0._rp

    dt2 = 0.5_rp*dt

!   Allocation dynamique
    allocate(taux0(Mesh%Nfacette))
    taux0(1:Mesh%Nfacette) = 1._rp
    allocate(area0(Mesh%Nfacette))
    allocate(fshape0(Mesh%Nfacette))
    fshape0(1:Mesh%Nfacette) = 1._rp
    allocate(TMetric(Mesh%NNoeud,2))
    allocate(RKVel(3,Nnoeud,4))

!   Tableaux locaux
    call NewMaillage(MeshTemp,PointMax)
    call CopyMaillage(MeshTemp,Mesh)
    call NewEcoulement(Ecoulement,Mesh%Nnoeud)
    !call Initialisation(Ecoulement,Mesh,t0)

!   Mesure de la qualit\E9 maillage Initial

    call PrePlots
    call CheckMesh3D_2(Mesh,taux0,area0,fshape0,t0,TMetric)
    !open(unit=iometric,file='metric.dat')
    call PlotMetric(t0,Mesh,TMetric)

    !open(unit=ioint,file="intersection_quality.dat")
    !open(unit=ioqual,file="immerged_surf_quality.dat")
    call CheckGeom(Mesh,fgeom,t0,ioint,ioqual)    

!   Boucle temporelle    
    do jt=1,nt

        ti = t0+(jt-1)*dt

!       RK1
        call MeshVel(MeshTemp,Ecoulement,ti,dt2,VBody)
        do j=1,Nnoeud
            RKVel(1:3,j,1) = MeshTemp%Tnoeud(j)%Velocity(1:3)
        enddo
        call Remesh(MeshTemp,Mesh,ti+dt2,dt2,fgeom)

!       RK2
        call MeshVel(MeshTemp,Ecoulement,ti+dt2,dt2,VBody)
        do j=1,Nnoeud
            RKVel(1:3,j,2) = MeshTemp%Tnoeud(j)%Velocity(1:3)
        enddo       
        call Remesh(MeshTemp,Mesh,ti+dt2,dt2,fgeom)

!       RK3
        call MeshVel(MeshTemp,Ecoulement,ti+dt2,dt,VBody)
        do j=1,Nnoeud
            RKVel(1:3,j,3) = MeshTemp%Tnoeud(j)%Velocity(1:3)
        enddo       
        call Remesh(MeshTemp,Mesh,ti+dt,dt,fgeom)

!       RK4
        call MeshVel(MeshTemp,Ecoulement,ti+dt,dt,VBody)
        do j=1,Nnoeud
            RKVel(1:3,j,4) = MeshTemp%Tnoeud(j)%Velocity(1:3)
        enddo

!       Mise \E0 jour
        do j = 1,Nnoeud
            MeshTemp%Tnoeud(j)%Velocity(1:3) = (RKVel(1:3,j,1) + 2._rp*RKVel(1:3,j,2) + 2._rp*RKVel(1:3,j,3) + RKVel(1:3,j,4))*denom
        enddo
        
!       Remaillage
        call Remesh(MeshTemp,Mesh,ti+dt,dt,fgeom)    
        call CopyMaillage(Mesh,MeshTemp)

!       Controle Maillage
        if(icheck)then
            call CheckGeom(Mesh,fgeom,t0,ioint,ioqual)!
            call CheckMesh3D_2(Mesh,taux0,area0,fshape0,ti+dt,TMetric)
            call PlotMetric(ti+dt,Mesh,TMetric)
        endif 

    enddo

    deallocate(taux0)
    deallocate(area0)
    deallocate(fshape0)
    deallocate(TMetric)
    deallocate(RKVel)
    close(iometric)
    close(ioint)
    close(ioqual)

    call DelEcoulement(Ecoulement)
    call DelMaillage(MeshTemp)

    end subroutine TestValidation_Deformation 
  
subroutine CheckMeshInit(Maillage,taux0,area0,fshape0)
! ----------------------------------------------
! Calcul de la taille de reference taux0
! ----------------------------------------------
    implicit none
    type(TMaillage),intent(in) :: Maillage
    real(rp),dimension(:),intent(inout) :: taux0
    real(rp),dimension(:),intent(inout) :: area0
    real(rp),dimension(:),intent(inout) :: fshape0
!   local
    integer :: j, k, NFS1, NFS2, NFS
    real(rp),dimension(2,2) :: Jb
    real(rp),dimension(2,3) :: PNoeud
    real(rp),dimension(2) :: L1,L2,L3
    real(rp) :: det
    integer :: ierror

    ierror = 0
    
    NFS1 = Maillage%FS%IndFS(2)
    NFS2 = Maillage%FS%IndFS(4)
    NFS = NFS2-NFS1+1

    do j=NFS1,NFS2
        k = j-NFS1+1
        PNoeud = Maillage%Tfacette(j)%PNoeud(1:2,1:3)
!	Edge
	L1 = PNoeud(1:2,3) - PNoeud(1:2,2)
        L2 = PNoeud(1:2,1) - PNoeud(1:2,3)
        L3 = PNoeud(1:2,2) - PNoeud(1:2,1)
!	Matrice Jacobian        
        Jb(1:2,1) = PNoeud(1:2,2) - PNoeud(1:2,1)
        Jb(1:2,2) = PNoeud(1:2,3) - PNoeud(1:2,1)
        det = Jb(2,1)*Jb(1,2)-Jb(1,1)*Jb(2,2)
        if(abs(det).gt.Epsilon)then
            taux0(k) = 1._rp/det
        else
            ierror = 100
            goto 9999
        endif
        area0(k) = Maillage%Tfacette(j)%aire
!	Metric shape initial        
        fshape0(k) = dot_product(L3,L3)+dot_product(L2,L2)+dot_product(L2,L3)
        fshape0(k) = (2._rp*sqrt(3._rp)*area0(k))/fshape0(k)
    enddo

9999 continue
    if(ierror.ne.0)then
        print*,"** error #",ierror," in CheskMeshInit"
    endif

end subroutine CheckMeshInit

subroutine CheckMesh(Maillage,taux0,area0,fshape0,ti,TMetric,iout)
! --------------------------------------------------
! Calcul des variations de taille des elements
! --------------------------------------------------
    implicit none
    type(TMaillage),intent(in) :: Maillage
    real(rp),dimension(:),intent(in) :: taux0
    real(rp),dimension(:),intent(in) :: area0
    real(rp),dimension(:),intent(in) :: fshape0
    real(rp),intent(in) :: ti
    real(rp),dimension(:,:),intent(inout) :: TMetric
    integer,intent(in), optional :: iout
!   local
    integer :: j, k, NFS1,NFS2,NFS,NF1,NFT
    integer :: nf, jf, kf
    real(rp),dimension(2,3) :: PNoeud
    real(rp),dimension(2,2) :: Jb1,Jb2,Jb3
    real(rp) :: det1, det2, det3, taux, fsize1, fsize2
    real(rp),dimension(:),allocatable :: fsize,fskeew,fss,fshape,fshs
    real(rp) :: area
	real(rp),dimension(2) :: L1,L2,L3

    NF1 = 1
    NFT = Maillage%Nfacette

    NFS1 = Maillage%FS%IndFS(2)
    NFS2 = Maillage%FS%IndFS(4)
    NFS = NFS2-NFS1+1

    allocate(fsize(NFT))
    allocate(fskeew(NFT))
    allocate(fss(NFT))
	allocate(fshape(NFT))
	allocate(fshs(NFT))

    do j=NF1,NFT
        !k = j-NFS1+1
        k=j
        PNoeud = Maillage%Tfacette(j)%Pnoeud(1:2,1:3)
!	Edge
        L1 = PNoeud(1:2,3) - PNoeud(1:2,2)
        L2 = PNoeud(1:2,1) - PNoeud(1:2,3)
        L3 = PNoeud(1:2,2) - PNoeud(1:2,1)
!       Jacobian matrix & determinant
        Jb1(1:2,1) = PNoeud(1:2,2) - PNoeud(1:2,1)
        Jb1(1:2,2) = PNoeud(1:2,3) - PNoeud(1:2,1)
        det1 = Jb1(2,1)*Jb1(1,2)-Jb1(1,1)*Jb1(2,2)
        Jb2(1:2,1) = PNoeud(1:2,3) - PNoeud(1:2,2)
        Jb2(1:2,2) = PNoeud(1:2,1) - Pnoeud(1:2,2)
        det2 = Jb2(2,1)*Jb2(1,2)-Jb2(1,1)*Jb2(2,2)
        Jb3(1:2,1) = PNoeud(1:2,1)-PNoeud(1:2,3)
        Jb3(1:2,2) = PNoeud(1:2,2)-PNoeud(1:2,3)
        det3 = Jb3(2,1)*Jb3(1,2)-Jb3(1,1)*Jb3(2,2)
!       Metric of relative size
        fsize1 = det1*taux0(k)
        if(abs(fsize1).gt.Epsilon)then
            fsize2 = 1._rp/fsize1
        else
            fsize2 = 9.
        endif
        fsize(k) = min(fsize1,fsize2)
!       Metric skeew
        area = Maillage%Tfacette(k)%Aire
        !if(abs(area).gt.Epsilon)then
        !    fskeew(k) = min(det1,det2,det3)*area0(k)/area
        !else
        !    print*,'warning : degenerating element (index = ',j,')'
        !endif
!       Metric skeew-size
        !fss(k) = sqrt(fsize(k))*fskeew(k)
!	Metric shape
	    fshape(k) = dot_product(L3,L3)+dot_product(L2,L2)+dot_product(L2,L3)
        fshape(k) = (2._rp*sqrt(3._rp)*area)/fshape(k)
	!fshapeR(k) = fshape(k)/fshape0(k)
!	Metric shape-size
	    fshs(k) = fshape(k)*fsize(k)
    enddo

    print*,"t = ",ti
    print*,"Fsize  : Min = ",minval(fsize(1:NFT))," ; Moy = ",sum(fsize(1:NFT))/dble(NFT)
    print*,"Fshape : Min = ",minval(fshape(1:NFT))," ; Moy = ",sum(fshape(1:NFT))/dble(NFT)
    print*,"Fmixt  : Min = ",minval(fshs(1:NFT))," ; Moy = ",sum(fshs(1:NFT))/dble(NFT)

    TMetric(:,:) = 0._rp
    !do j = Maillage%FS%IndFS(1),Maillage%FS%IndFS(3)
    do j = 1,Maillage%Nnoeud
        nf = Maillage%Tnoeud(j)%Nfacette
        do k = 1,nf
            jf = Maillage%Tnoeud(j)%Tfacette(k,1)
            !kf = jf-NFS1+1
            kf = jf
            !TMetric(j,1:2) = TMetric(j,1:2) + [fsize(kf),fshape(kf)/fshape0(kf)]
            TMetric(j,1:2) = TMetric(j,1:2) + [fsize(kf),fshape(kf)]
        enddo
       TMetric(j,1:2) = TMetric(j,1:2)/dble(nf)
    enddo
                
    
    deallocate(fsize)
    deallocate(fskeew)
    deallocate(fss)
    deallocate(fshape)
    deallocate(fshs)
    
end subroutine CheckMesh

subroutine CheckMesh3D_2(Maillage,taux0,area0,fshape0,ti,TMetric,iout)
! --------------------------------------------------
! Calcul des variations de taille des elements
! --------------------------------------------------
    implicit none
    type(TMaillage),intent(in) :: Maillage
    real(rp),dimension(:),intent(in) :: taux0
    real(rp),dimension(:),intent(in) :: area0
    real(rp),dimension(:),intent(in) :: fshape0
    real(rp),intent(in) :: ti
    real(rp),dimension(:,:),intent(inout) :: TMetric
    integer,intent(in), optional :: iout
!   local
    integer :: j, k, NFS1,NFS2,NFS,NF1,NFT
    integer :: nf, jf, kf
    real(rp),dimension(3,3) :: PNoeud
    real(rp),dimension(2,2) :: Jb1,Jb2,Jb3
    real(rp) :: det1, det2, det3, taux, fsize1, fsize2
    real(rp),dimension(:),allocatable :: fsize,fskeew,fss,fshape,fshs
    real(rp) :: area
	real(rp),dimension(3) :: L1,L2,L3

    NF1 = 1
    NFT = Maillage%Nfacette

    allocate(fsize(NFT))
    allocate(fskeew(NFT))
    allocate(fss(NFT))
	allocate(fshape(NFT))
	allocate(fshs(NFT))

    do j=NF1,NFT
        PNoeud = Maillage%Tfacette(j)%Pnoeud(1:3,1:3)
        area   = Maillage%Tfacette(j)%Aire
!	Edge
        L1 = PNoeud(1:3,3) - PNoeud(1:3,2)
        L2 = PNoeud(1:3,1) - PNoeud(1:3,3)
        L3 = PNoeud(1:3,2) - PNoeud(1:3,1)
!   Jacobian matrix & determinant

        !Jb1(1:2,1) = PNoeud(1:2,2) - PNoeud(1:2,1)
        !Jb1(1:2,2) = PNoeud(1:2,3) - PNoeud(1:2,1)
        !det1 = Jb1(2,1)*Jb1(1,2)-Jb1(1,1)*Jb1(2,2)

        !Jb2(1:2,1) = PNoeud(1:2,3) - PNoeud(1:2,2)
        !Jb2(1:2,2) = PNoeud(1:2,1) - Pnoeud(1:2,2)
        !det2 = Jb2(2,1)*Jb2(1,2)-Jb2(1,1)*Jb2(2,2)

        !Jb3(1:2,1) = PNoeud(1:2,1)-PNoeud(1:2,3)
        !Jb3(1:2,2) = PNoeud(1:2,2)-PNoeud(1:2,3)
        !det3 = Jb3(2,1)*Jb3(1,2)-Jb3(1,1)*Jb3(2,2)

!    Metric fsize
        !fsize1 = det1*taux0(j)
        !if(abs(fsize1).gt.Epsilon)then
        !    fsize2 = 1._rp/fsize1
        !else
        !    fsize2 = 9.
        !endif
        !fsize(k) = min(fsize1,fsize2)
!	Metric shape
	    fshape(j) = dot_product(L3,L3)+dot_product(L2,L2)+dot_product(L2,L3)
        fshape(j) = (2._rp*sqrt(3._rp)*area)/fshape(j)
	    fshape(j) = fshape(j)/fshape0(j)       ! fshape reduit (si cas deformation)
!	Metric shape-size
	    !fshs(j) = fshape(j)*fsize(j)
    enddo

    print*,"t = ",ti
    print*,"Fsize  : Min = ",minval(fsize(1:NFT))," ; Moy = ",sum(fsize(1:NFT))/dble(NFT)
    print*,"Fshape : Min = ",minval(fshape(1:NFT))," ; Moy = ",sum(fshape(1:NFT))/dble(NFT)
    print*,"Fmixt  : Min = ",minval(fshs(1:NFT))," ; Moy = ",sum(fshs(1:NFT))/dble(NFT)

    TMetric(:,:) = 0._rp
    !do j = Maillage%FS%IndFS(1),Maillage%FS%IndFS(3)
    do j = 1,Maillage%Nnoeud
        nf = Maillage%Tnoeud(j)%Nfacette
        do k = 1,nf
            jf = Maillage%Tnoeud(j)%Tfacette(k,1)
            !TMetric(j,1:2) = TMetric(j,1:2) + [fsize(kf),fshape(kf)/fshape0(kf)]
            TMetric(j,1:2) = TMetric(j,1:2) + [fsize(jf),fshape(jf)]
        enddo
        TMetric(j,1:2) = TMetric(j,1:2)/dble(nf)
    enddo
                
    
    deallocate(fsize)
    deallocate(fskeew)
    deallocate(fss)
    deallocate(fshape)
    deallocate(fshs)
    
end subroutine CheckMesh3D_2


subroutine phi1(n,r,r0,val)
    implicit none
    integer,intent(in) :: n
    real(rp),dimension(n),intent(in) :: r
    real(rp),intent(in) :: r0
    real(rp),dimension(n),intent(inout) :: val    
    val(1:n) = (1.-r(1:n)/Ldom(5))**2
end subroutine phi1

subroutine phi2(n,r,r0,val)
    implicit none
    integer,intent(in) :: n
    real(rp),dimension(n),intent(in) :: r
    real(rp),intent(in) :: r0
    real(rp),dimension(n),intent(inout) :: val    
    val(1:n) = (1.-r(1:n)/Ldom(5))**4*(4._rp*r(1:n)/Ldom(5)+1)
end subroutine phi2

subroutine phi6(n,r,r0,val)
    implicit none
    integer,intent(in) :: n
    real(rp),dimension(n),intent(in) :: r
    real(rp),intent(in) :: r0
    real(rp),dimension(n),intent(inout) :: val    
    real(rp),dimension(n) :: psi
    integer :: j
    psi(1:n) = r(1:n)/Ldom(5)
    do j=1,n
        if(abs(psi(j)).gt.Epsilon)then
            val(j) = 1._rp + 80./3.*psi(j)**2 - 40.*psi(j)**3 + 15.*psi(j)**4 - 8./3.*psi(j)**5 + 20.*psi(j)**2*log(psi(j))
        else
            val(j) = 1._rp
        endif
    enddo
end subroutine phi6

subroutine phi9(n,r,r0,val)
    implicit none
    integer,intent(in) :: n
    real(rp),dimension(n),intent(in) :: r
    real(rp),intent(in) :: r0
    real(rp),dimension(n),intent(inout) :: val    
    integer :: j
    do j=1,n
        if(r(j).gt.Epsilon)then
            val(j) = r(j)**2 * log(r(j))
        else
            val(j) = 0._rp
        endif
    enddo
end subroutine phi9

subroutine phi10(n,r,r0,val)
  implicit none
  integer,intent(in) :: n
  real(rp),dimension(n),intent(in) :: r
  real(rp),intent(in) :: r0
  real(rp),dimension(n),intent(inout) :: val
  val(1:n) = sqrt ( r(1:n)**2 + r0**2 )
end subroutine phi10

subroutine phi11(n,r,r0,val)
  implicit none
  integer,intent(in) :: n
  real(rp),dimension(n),intent(in) :: r
  real(rp),intent(in) :: r0
  real(rp),dimension(n),intent(inout) :: val
  val(1:n) = 1._rp / sqrt ( r(1:n)**2 + r0**2 )
end subroutine phi11

subroutine phi12(n,r,r0,val)
  implicit none
  integer,intent(in) :: n
  real(rp),dimension(n),intent(in) :: r
  real(rp),intent(in) :: r0
  real(rp),dimension(n),intent(inout) :: val
  val(1:n) = 1+r(1:n)**2
end subroutine phi12

subroutine phi13(n,r,r0,val)
  implicit none
  integer,intent(in) :: n
  real(rp),dimension(n),intent(in) :: r
  real(rp),intent(in) :: r0
  real(rp),dimension(n),intent(inout) :: val
  val(1:n) = 1._rp / (1._rp + r(1:n)**2)
end subroutine phi13

subroutine phi14(n,r,r0,val)
  implicit none
  integer,intent(in) :: n
  real(rp),dimension(n),intent(in) :: r
  real(rp),intent(in) :: r0
  real(rp),dimension(n),intent(inout) :: val
  val(1:n) = exp( - r(1:n)**2 )
end subroutine phi14
end module Validation_RBF