module CoeffInfluence_CL
use FonctionsCommunes

implicit none

contains

    ! Norme
    function norm_d(u)
        real(rp), dimension(3) :: u
        real(rp) :: norm_d
        integer :: i
    
        norm_d = 0
        do i = 1, 3
            norm_d = norm_d + u(i)*u(i)
        end do
        norm_d = sqrt(norm_d)
    
    end function


    ! Produit scalaire
    function dot_prod_d(u,v)
        real(rp), dimension(3) :: u,v
        real(rp) :: dot_prod_d
        integer :: i
    
        dot_prod_d = 0
        do i = 1, 3
            dot_prod_d = dot_prod_d + u(i)*v(i)
        end do
    end function
    
    !Produit vectoriel
    subroutine rot_d(res, u, v)
        real(rp), dimension(3) :: res, u, v
    
        res(1) = u(2)*v(3) - u(3)*v(2)
        res(2) = u(3)*v(1) - u(1)*v(3)
        res(3) = u(1)*v(2) - u(2)*v(1)
    
    end subroutine


    !Produit X.A : X(vector) et A(matrice 3x3)
    subroutine matmul_d(res,X,A)
        real(rp) :: X(3), A(3,3), res(3)
        integer :: i,j
    
        res = 0
        do i = 1, 3
            do j = 1, 3
                res(i) = res(i) + X(j)*A(j,i)
            end do
        end do    
    end subroutine



    function S_sigma_sing_d(t, L_P, M)

        real(rp) :: S_sigma_sing_d
        integer :: t
        real(rp) :: L_P(3,3), M(3)
        real(rp) :: P1(3), P2(3), P3(3)
        real(rp) :: i(3), j(3), h(3)
        real(rp) :: delta, a, b
        integer :: k
    
        do k=1,3
            P1(k) = L_P(t,k)
            P2(k) = L_P(mod(t,3)+1,k)
            P3(k) = L_P(mod(t+1,3)+1,k)
        
            i(k) = 0.5*(P2(k)-P1(k) + P3(k) - P1(k))
            j(k) = P3(k) - P2(k)
        end do
    
        call rot_d(h, i, j)
    
        delta = dot_prod_d(i,i)*dot_prod_d(j,j) - dot_prod_d(i,j)**2
        a = 0.5*dot_prod_d(j,j)/sqrt(delta)
        b = dot_prod_d(i,j)/sqrt(delta)
    
        S_sigma_sing_d = norm_d(h)*log((a+b+sqrt(1+(a+b)**2))/(-a+b+sqrt(1+(-a+b)**2)))/norm_d(j)
    
    
    end function
    
    
    subroutine S_sigma_mu_d(L_P, M, G, n, mu, sigma)
        real(rp) :: L_P(3,3), M(3), G(3), n(3)
        real(rp) :: mu, sigma
        real(rp) :: GM(3), PkPk1(3), PkM(3), Pk1M(3)
        real(rp) :: Rk, Rk1, dk, Z, Nk1, Nkt, Dk1, Dkt
        integer :: t, i, j
        real(rp) :: temp(3)
    
        mu = 0
        sigma = 0
    
        t = -1
        do i = 1, 3
            do j = 1, 3
                temp(j) = L_P(i,j) - M(j)
            end do
            if (norm_d(temp) < 1e-5) then
                t = i
            end if
        end do
    
    
        if (t==-1) then
    
            do i = 1, 3
                GM(i) = M(i) - G(i)
            end do
        
            Z = dot_product(GM,n)
        
            do i = 1, 3
                do j = 1, 3
                    PkM(j) = M(j) - L_P(i,j)
                    Pk1M(j) = M(j) - L_P(mod(i,3)+1,j)
                    PkPk1(j) = L_P(mod(i,3)+1,j) - L_P(i,j)
                end do
            
                Rk = norm_d(PkM)
                Rk1 = norm_d(Pk1M)
                dk = norm_d(PkPk1)
                Nk1 = Rk1 + Rk + dk
                Dk1 = Rk1 + Rk - dk
                call rot_d(temp, n, PkPk1)
                Nkt = 2*dot_prod_d(PkM,temp) ! Mettre un point ??
                Dkt = (Rk1+Rk)**2 - dk**2 + 2.*abs(Z)*(Rk1+Rk)
            
                if (abs(Z) > 1e-5) then
                    mu = mu + 2.*(Z/abs(Z))*atan(Nkt/Dkt)
                end if
            
                sigma = sigma + Nkt*log(Nk1/Dk1)/(2.*dk) - 2.*abs(Z)*atan(Nkt/Dkt)
            end do
        
        else
    
            sigma = S_sigma_sing_d(t, L_P, M)
        
        end if
    end subroutine
        

    subroutine coeff_asympt_d(L_P, M, G, n, CD, CS)
        real(rp) :: L_P(3,3), M(3), G(3), n(3)
        real(rp) :: CD(3), CS(3)
        real(rp) :: Cdd, Css
        real(rp) :: GM(3), M1M2(3), M1M3(3)
        real(rp) :: temp(3)
        real(rp) :: A
        real(rp) :: inv3 = 1._rp/3._rp
        integer :: j
    
        do j = 1, 3
            GM(j) = M(j) - G(j)
            M1M2(j) = L_P(2,j) - L_P(1,j)
            M1M3(j) = L_P(3,j) - L_P(1,j)
        end do
    
        call rot_d(temp, M1M2, M1M3)
        A = norm_d(temp)/2.
        
        Css = inv3*A/norm_d(GM)
        Cdd = inv3*dot_prod_d(GM,n)*A/norm_d(GM)**3
        do j = 1, 3
            CD(j) = Cdd
            CS(j) = Css
        end do
    
    end subroutine



    function sing_d(L_P, M)

        logical :: sing_d
        real(rp) :: L_P(3,3), M(3)
        real(rp) :: PM(3)
        integer :: i,j
    
        sing_d = .FALSE.
    
        do i = 1, 3
            do j = 1, 3
                PM(j) = M(j) - L_P(i,j)
            end do
        
            if (norm_d(PM) == 0) then!A changer ????
                sing_d = .TRUE.
            end if
        end do
    end function


    function asympt_d(G, M, R_max, N_seuil)
        logical :: asympt_d
        real(rp) :: G(3), M(3), R_max
        real(rp) :: GM(3)
        integer :: N_seuil
        integer :: j
    
        do j = 1, 3
            GM(j) = M(j) - G(j)
        end do
    
        asympt_d = (norm_d(GM) > N_seuil*R_max)
    
    end function
    
    
    subroutine I_sigma_mu_d(L_P, M, n, mu, sigma)
        real(rp) :: L_P(3,3), M(3), n(3)
        real(rp) :: mu(3), sigma(3)
        
        real(rp) :: a(3), b(3), AB(3), AM(3), BM(3)
        real(rp) :: res(3)
        real(rp) :: K2, a_m, b_m, AB_n, dp_AB_AM, coeff
        
        integer :: i, j
        
        do j = 1, 3
            mu(j) = 0
            sigma(j) = 0
        end do
        
        do i = 1, 3
            do j = 1, 3
                A(j) = L_P(i,j)
                B(j) = L_P(mod(i,3)+1, j)
                
                AB(j) = B(j) - A(j)
                AM(j) = M(j) - A(j)
                BM(j) = M(j) - B(j)
            end do
            
            dp_AB_AM = dot_prod_d(AB, AM)
            AB_n = norm_d(AB)
            K2 = dot_prod_d(AM,AM) - (dp_AB_AM/AB_n)**2
            
            if (K2 > 1e-5) then
                
                a_m = asinh(-dp_AB_AM/(AB_n*sqrt(K2)))
                b_m = asinh( (AB_n - dp_AB_AM/AB_n)/sqrt(K2))
                
                call rot_d(res, AM, AB)   
                    
                do j = 1 , 3             
                    mu(j) = mu(j) - res(j)*(b_m-a_m)/AB_n
                end do
                                
                call rot_d(res, n, AB)
                
                coeff = K2/(2.*AB_n)*(b_m-a_m+(sinh(2.*b_m) - sinh(2.*a_m))/2.)
                
                do j = 1, 3
                    sigma(j) = sigma(j) + coeff*res(j)
                end do
                
            else
                if (norm_d(AM) < 1e-5 .or. norm_d(BM) < 1e-5) then
                    call rot_d(res,n,AB)
                    
                    do j = 1, 3
                        sigma(j) = sigma(j) + 0.5*AB_n*res(j)
                    end do
                else
                    call rot_d(res, n, AB)
                    
                    coeff = abs(norm_d(AM) - 0.5*dp_AB_AM/norm_d(AM))
                
                    do j = 1, 3
                        sigma(j) = sigma(j) + coeff*res(j)
                    end do
                end if
            end if
        end do
    end subroutine
    
    
    subroutine coeff_inf_d(L_P, n, G, GS, R_max, M, CS, CD, i_n, i_f)
        real(rp) :: CS(3), CD(3)
        real(rp) :: L_P(3,3), n(3), G(3), GS(3,3), M(3)
        real(rp) :: R_max
        
        real(rp) :: GM(3)
        real(rp) :: Ssigma, Smu, Isigma(3), Imu(3)
        real(rp) :: temp1(3), temp2(3), temp3(3)
        real(rp) :: inv3
        
        integer :: N_seuil = 8
        
        integer :: j
        
        integer :: i_n, i_f ! Asupprimer
        
        inv3 = 1._RP/3._RP
        
     
        
        
        
        
                
        if (asympt_d(G, M, R_max, N_seuil)) then
            call coeff_asympt_d(L_P, M, G, n, CD, CS)
        else
        
            do j = 1, 3
                GM(j) = M(j) - G(j)
            end do
            
            call I_sigma_mu_d(L_P, M, n, Imu, Isigma)
            call S_sigma_mu_d(L_P, M, G, n, Smu, Ssigma)
            
            call matmul_d(temp1, GM,GS)
            call matmul_d(temp2, Imu, GS)
            call matmul_d(temp3, Isigma, GS)
        
            do j = 1, 3
                CD(j) = (inv3 + temp1(j))*Smu - temp2(j)
                CS(j) = (inv3 + temp1(j))*Ssigma - temp3(j)
            end do
        end if

        

        
        if (i_n == 2 .and. i_f < 4) then
            print *, i_n, i_f
            print *, "Asympt", asympt_d(G, M, R_max, N_seuil)
            print *, "Cdd", CD

        end if
        
      

    end subroutine
                
    
    subroutine mat_CI_kernel(L_X, L_T, L_n, L_GS, L_G, L_R_max, A_CD, A_CS, N_n, N_f, i_f_i, i_f_f, N_sym, prof)
        real(rp) :: L_X(N_n,3), L_n(N_f,3), L_GS(N_f,3,3), L_G(N_f,3), L_R_max(N_f)
        integer :: L_T(N_f,3)
        real(rp) :: A_CD(N_n, N_n), A_CS(N_n, N_n)
        integer :: N_n, N_f, N_sym, i_f_i, i_f_f
        real(rp) :: prof
        
        real(rp) :: CS(3), CD(3)
        real(rp) :: L_P(3,3), M(3), GS(3,3), n(3), G(3)
        integer :: i_n, i_f, j, k, i, i_sym
        
        
        do i_n = 1, N_n
            do k = 1, 2
                M(k) = L_X(i_n,k)
            end do
            
            do i_f = i_f_i, i_f_f
                do j = 1, 3
                    i = L_T(i_f, j)
                    
                    do k = 1, 3
                        L_P(j,k) = L_X(i,k)
                    end do
                end do
                
                do j = 1,3
                    do k = 1, 3
                        GS(j,k) = L_GS(i_f, j, k)
                    end do
                end do
                
                do k = 1, 3
                    n(k) = L_n(i_f, k)
                    G(k) = L_G(i_f, k)
                end do
                
                do i_sym = 0, N_sym
                    M(3) = (1-2*i_sym)*L_X(i_n,3) - 2.*i_sym*prof
                    

                    
                    call coeff_inf_d(L_P, n, G, GS, L_R_max(i_f), M, CS, CD, i_n, i_f)
                    
                    ! A modifier : ne faire l'addition qu'apres avoir fait la somme des deux symetriques
                    do k = 1, 3
                        i = L_T(i_f, k)
                        
                        A_CS(i_n, i) = A_CS(i_n, i) + CS(k)
                        A_CD(i_n, i) = A_CD(i_n, i) + CD(k)
                    end do
                end do
            end do
        end do
        
    end subroutine
    

    subroutine angle_solide_kernel(A_CD, N_n)
        real(rp) :: A_CD(N_n, N_n)
        integer :: N_n
        integer :: i_n, k
        
        
        do i_n = 1, N_n
            A_CD(i_n, i_n) = 0
            
            do k = 1, i_n-1
                A_CD(i_n, i_n) = A_CD(i_n, i_n) - A_CD(i_n, k)
            end do
            
            do k = i_n + 1, N_n
                A_CD(i_n, i_n) = A_CD(i_n, i_n) - A_CD(i_n, k)
            end do
        end do
    end subroutine
    
end module