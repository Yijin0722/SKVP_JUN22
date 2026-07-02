!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************
!
 SUBROUTINE basic_aux_mat_calcul
!
!*********************************************************************************************
! 
! This subroutine calculates matrices containing values of integrals
!
! BAM_1(i,j)    : Int[ (rhomin..rhomax)              u_i'(x)    * u_j'(x) dx ]
! BAM_3(i,j)    : Int[ (rhomin..rhomax)    u_i(x)  *    1       * u_j (x) dx ]
! BAM_4(i,j)    : Int[ (rhomin..rhomax)    u_i(x)  *  1/x**2    * u_j (x) dx ]

! BAM_x1(i,j)   : Int[ (rhomin..rhomax)    u_i(x)  * d2/dx2     * u_0j(x) dx ]
! BAM_x3(i,j)   : Int[ (rhomin..rhomax)    u_i(x)  *    1       * u_0j(x) dx ]
! BAM_x4(i,j)   : Int[ (rhomin..rhomax)    u_i(x)  *  1/x**2    * u_0j(x) dx ]
!
! BAM_xx1(i,j)  : Int[ (rhomin..rhomax)    u_0i(x) * d2/dx2     * u_0j(x) dx ]
! BAM_xx3(i,j)  : Int[ (rhomin..rhomax)    u_0i(x) *    1       * u_0j(x) dx ]
! BAM_xx4(i,j)  : Int[ (rhomin..rhomax)    u_0i(x) *  1/x**2    * u_0j(x) dx ]
        
! BAM_xx1b(i,j) : Int[ (rhomin..rhomax)    u_1i(x) * d2/dx2     * u_0j(x) dx ]
! BAM_xx3b(i,j) : Int[ (rhomin..rhomax)    u_1i(x) *    1       * u_0j(x) dx ]
! BAM_xx4b(i,j) : Int[ (rhomin..rhomax)    u_1i(x) *  1/x**2    * u_0j(x) dx ]
!
!=============================================================================================
!
        USE generateparameters
        USE AtomDiatomskvp
        USE omp_lib

        
!
        IMPLICIT NONE
!--------------------
        INTEGER :: i, j, i_x, j_x, k, n, file_size, n_prime
        
        REAL(8) :: tm11, tm12, xm, dx, xr
        !REAL(8) :: GetPotential, x1, x2, y1, y2 
        COMPLEX(8) :: sx1, sx3, sx4, sx1b, sx3b, sx4b
        COMPLEX(8) :: aux0_1, aux0_2, aux1_1, aux1_2, aux2_1, aux2_2, aux3_1, aux3_2
        COMPLEX(8) :: aux4_1, aux4_2, aux5_1, aux5_2, aux6_1, aux6_2, aux7_1, aux7_2
        LOGICAL :: file_exists, new_calc, same_param
!
        new_calc=.TRUE.
        INQUIRE(file='previous_calc.dat', exist=file_exists, size=file_size)
        IF(file_exists .eqv. .FALSE.)new_calc=.TRUE.
        IF((file_exists .eqv. .TRUE.) .AND. (file_size==0))new_calc=.TRUE.
        IF ((file_exists .eqv. .TRUE.) .AND. (file_size/=0)) THEN

        OPEN(unit=17, file='previous_calc.dat', form='unformatted', status='unknown')
        READ(17) n_xmb, k_xmb, jmxmb, n_openmb
        same_param = ((pbasst(1)%pb_nbr==n_xmb).and.(pbasst(1)%pb_pa1==k_xmb).and.(pbasst(2)%pb_nbr==jmxmb).and.(n_open==n_openmb))
        IF(same_param .eqv. .TRUE.)new_calc=.FALSE. 
        CLOSE(17) 
        ENDIF

        IF (new_calc .eqv. .TRUE.) THEN
        PRINT*, " "
        PRINT*, "running noc independent calc"
        OPEN(unit=17, file='previous_calc.dat', form='unformatted', status='unknown')
        WRITE(17) pbasst(1)%pb_nbr, pbasst(1)%pb_pa1, pbasst(2)%pb_nbr, n_open
        CLOSE(17)
!
        ALLOCATE ( BAM_1(1:pbasst(1)%pb_nbr,1:pbasst(1)%pb_nbr), BAM_3(1:pbasst(1)%pb_nbr,1:pbasst(1)%pb_nbr), &
                   BAM_4(1:pbasst(1)%pb_nbr,1:pbasst(1)%pb_nbr), &
                norm(1:pbasst(1)%pb_nbr))
!

               
! ===Calculation of normalization factors for B-spline basis functions===                
        norm = 0d0
        DO i_x = 1, pbasst(1)%pb_nbr
        DO j = i_x, i_x+pbasst(1)%pb_pa1-1
                xm = 5D-1*(knots_x(j+1)+knots_x(j))
                xr = 5D-1*(knots_x(j+1)-knots_x(j))
                sx1=complex(0D0,0D0)
!
!               Perform gaussian quadrature for x
!               ---------------------------------
                DO i = 1, ngqp_x/2
                        dx = xr * gq_root_x(i+ngqp_x/2)
                        sx1 = sx1 + gq_weight_x(i+ngqp_x/2) * &
                        (bsp_x(i_x,xm+dx)**2d0 + bsp_x(i_x,xm-dx)**2d0)
                ENDDO
                norm(i_x) = norm(i_x) + dble(sx1 * xr)
        ENDDO
        ENDDO
!=============================================================================
!
!Calculation of BAM_1, BAM_3, and BAM_4 matrices which are independent of the B-spline optimization
!BAM_1(i,j)    : Int[ (rhomin..rhomax)              u_i'(x)    * u_j'(x) dx ]
!BAM_3(i,j)    : Int[ (rhomin..rhomax)    u_i(x)  *    1       * u_j (x) dx ]
!BAM_4(i,j)    : Int[ (rhomin..rhomax)    u_i(x)  *  1/x**2    * u_j (x) dx ]
        BAM_1=complex(0D0,0D0)
        BAM_3=complex(0D0,0D0)
        BAM_4=complex(0D0,0D0)
!
        DO i_x = 1, pbasst(1)%pb_nbr
        DO j_x = 1, min(i_x+pbasst(1)%pb_pa1-1,pbasst(1)%pb_nbr) !n_x
        DO j = j_x, i_x+pbasst(1)%pb_pa1-1
                xm = 5D-1*(knots_x(j+1)+knots_x(j))
                xr = 5D-1*(knots_x(j+1)-knots_x(j))
                !PRINT*, 'xm', xm
                !PRINT*, 'xr', xr
                sx1=complex(0D0,0D0)
                sx3=complex(0D0,0D0)
                sx4=complex(0D0,0D0)
               
!
!               Perform gaussian quadrature for x
!               ---------------------------------
                DO i = 1, ngqp_x/2
                        dx = xr * gq_root_x(i+ngqp_x/2)
                        IF (abs(xm + dx) .le. 1e-12) THEN
                                aux4_1 = 0D0
                                aux4_2 = 0D0
                        ELSE
                                aux4_1  = (xm+dx)**(-2)
                                aux4_2  = (xm-dx)**(-2)
                        ENDIF
                        sx1 = sx1 + gq_weight_x(i+ngqp_x/2) * &
                             ( dbsp_x(i_x,xm+dx)*dbsp_x(j_x,xm+dx) + &
                               dbsp_x(i_x,xm-dx)*dbsp_x(j_x,xm-dx) )
                        sx3 = sx3 + gq_weight_x(i+ngqp_x/2) * &
                             ( bsp_x(i_x,xm+dx)*bsp_x(j_x,xm+dx) + &
                               bsp_x(i_x,xm-dx)*bsp_x(j_x,xm-dx) )
                        sx4 = sx4 + gq_weight_x(i+ngqp_x/2) * &
                             ( aux4_1*bsp_x(i_x,xm+dx)*bsp_x(j_x,xm+dx) + &
                               aux4_2*bsp_x(i_x,xm-dx)*bsp_x(j_x,xm-dx) )
                
                ENDDO

                BAM_1(i_x,j_x) = BAM_1(i_x,j_x) + sx1 * xr / dsqrt(norm(i_x)*norm(j_x))
                BAM_3(i_x,j_x) = BAM_3(i_x,j_x) + sx3 * xr / dsqrt(norm(i_x)*norm(j_x))
                BAM_4(i_x,j_x) = BAM_4(i_x,j_x) + sx4 * xr / dsqrt(norm(i_x)*norm(j_x))

               
        ENDDO
        ENDDO
        ENDDO
!
        ! Writing matrices to file

        
        OPEN(unit=10, file='matrices.bin', form='unformatted', access='stream', status='replace')
        WRITE(10) norm
        WRITE(10) BAM_1
        WRITE(10) BAM_3
        WRITE(10) BAM_4
        CLOSE(10)
      

        ELSE 
        ALLOCATE ( BAM_1(1:pbasst(1)%pb_nbr,1:pbasst(1)%pb_nbr), BAM_3(1:pbasst(1)%pb_nbr,1:pbasst(1)%pb_nbr), &
                   BAM_4(1:pbasst(1)%pb_nbr,1:pbasst(1)%pb_nbr), norm(1:pbasst(1)%pb_nbr))

        BAM_1=complex(0D0,0D0)
        BAM_3=complex(0D0,0D0)
        BAM_4=complex(0D0,0D0)
      
        ENDIF
!       
!
        OPEN(unit=10, file='matrices.bin', form='unformatted', access='stream', status='old')
        READ(10) norm
        READ(10) BAM_1
        READ(10) BAM_3
        READ(10) BAM_4
        CLOSE(10)


!
        CALL CPU_TIME(tm11)
!
        ! =====================================================================================================
        !! (j,k) -> n, Compacting this set of quantum numbers into one quantum number 
        !!ncf = (((2*Jtot + 1)*jmx)/2)+1
        !!PRINT*, "ncf", ncf
        !!PRINT*, "ENERGY: ", E
        !!PRINT*, "alpha: ", alpha
        !!PRINT*, "R0: ", r0
        !!PRINT*, "rhomin: ", rhomin
        !!PRINT*, "rhomax", rhomax
        !!PRINT*, "Mu: ", mu
        !!PRINT*, "B_rot", B_rot
        !!
        !!ALLOCATE(quant_mat(2,ncf))
        !!n = 1
        !!DO j = 0, jmx, 2
        !!        IF (j == 0) THEN
        !!                quant_mat(1,n) = 0
        !!                quant_mat(2,n) = 0
        !!                n = n + 1
        !!        ELSE
        !!        DO k = -Jtot,Jtot
        !!        quant_mat(1,n) = j
        !!        quant_mat(2,n) = k
        !!        n = n + 1
        !!        ENDDO
        !!        ENDIF
        !!END DO
        !!PRINT*," " 
        !!!PRINT*, "bsp(0, 2.1)", bsp_x(0, 2.0d0)
        !!PRINT*, "quant_mat"
        !!PRINT*, quant_mat(1, 1:10)
        !!PRINT*, quant_mat(2, 1:10)

!
!       Computation of basic auxiliary matrices
!       ---------------------------------------


        ALLOCATE ( BAM_x1(1:pbasst(1)%pb_nbr,1:n_open), BAM_x3(1:pbasst(1)%pb_nbr,1:n_open), BAM_x4(1:pbasst(1)%pb_nbr,1:n_open), &
                   norm_vec(1:pbasst(1)%pb_nbr), STAT = istatus )
        ALLOCATE ( BAM_xx1(1:n_open,1:n_open), BAM_xx3(1:n_open,1:n_open), BAM_xx4(1:n_open,1:n_open), STAT = istatus )
        ALLOCATE ( BAM_xx1b(1:n_open,1:n_open), BAM_xx3b(1:n_open,1:n_open), BAM_xx4b(1:n_open,1:n_open), STAT = istatus )
        
        !PRINT*, GetPotential(2d0, 0d0)
        !PRINT*, "V0"
        !PRINT*, (v0(i*1d0), i = 1, 100)

! 
!       norm=complex(0D0,0D0)

        BAM_x1 = (0D0,0D0)
        BAM_x3 = (0D0,0D0)
        BAM_x4 = (0D0,0D0)
        
        BAM_xx1 = (0D0,0D0)
        BAM_xx3 = (0D0,0D0)
        BAM_xx4 = (0D0,0D0)
       
        BAM_xx1b = (0D0,0D0)
        BAM_xx3b = (0D0,0D0)
        BAM_xx4b = (0D0,0D0)
        
        !PRINT*, "Potential: ", GetPotential(2d0,0d0)
        !PRINT*, "error:", DSQRT(715d0)/143d0 - wigner3j(4,6,2,0,0,0)

        ! PRINT*, "Norm"
        ! DO i = 1, n_x
        !         PRINT*, norm(i)
        ! ENDDO
        ! norm_vec = 0d0
        ! PRINT*, "Norm2"
        ! DO k = 1, n_x
        ! DO i = 1, ngqp_x
        !         norm_vec(k) = norm_vec(k) +  wx(i) * ((bsp_x(k, x(i)))**2d0)*(x(i)**2d0)
        ! ENDDO
        ! ENDDO
        ! PRINT*, "Norm_vec"
        ! DO i = 1, n_x
        !         PRINT*, norm_vec(i)
        ! ENDDO

        ! PRINT*, "d2DR2u0(1,1)", d2dR2u0(1,1d0)
!       
        IF (bsp_optimized) THEN

        PRINT*, " "
        PRINT*, "Running B-spline optimized program"
        PRINT*, " "


        !print*,'dimension BAM_x1:', size(BAM_x1, dim=1), size(BAM_x1, dim=2)
        !print*,'dimension BAM_x3:', size(BAM_x3, dim=1), size(BAM_x3, dim=2)
        !print*,'dimension BAM_x4:', size(BAM_x4, dim=1), size(BAM_x4, dim=2)

        
!Calculation of BAM_x1, BAM_x3, and BAM_x4 matrices for Close-Open coupling
!BAM_x1(i,j)   : Int[ (rhomin..rhomax)    u_i(x)  * d2/dx2     * u_0j(x) dx ]
!BAM_x3(i,j)   : Int[ (rhomin..rhomax)    u_i(x)  *    1       * u_0j(x) dx ]
!BAM_x4(i,j)   : Int[ (rhomin..rhomax)    u_i(x)  *  1/x**2    * u_0j(x) dx ]

        DO i_x = 1, pbasst(1)%pb_nbr
        DO j_x = 1, n_open
        !quant_j = quant_mat(1,j_x)
!
        DO j = pbasst(1)%pb_pa1, pbasst(1)%pb_nbr
                xm = 5D-1*(knots_x(j+1)+knots_x(j))
                xr = 5D-1*(knots_x(j+1)-knots_x(j))
!
                sx1=complex(0D0,0D0)
                sx3=complex(0D0,0D0)
                sx4=complex(0D0,0D0)
                

!
!               Perform gaussian quadrature for x
!               ---------------------------------
                DO i = 1, ngqp_x/2
                        dx = xr * gq_root_x(i+ngqp_x/2)
                        aux0_1  = 5d-1*(1d0+dtanh(alpha*(xm+dx-r0)))
                        aux0_2  = 5d-1*(1d0+dtanh(alpha*(xm-dx-r0)))
                        aux1_1  = 5d-1*(alpha/((dcosh(alpha*(xm+dx-r0)))**2))
                        aux1_2  = 5d-1*(alpha/((dcosh(alpha*(xm-dx-r0)))**2))
                        aux2_1  = -2d0*alpha*dtanh(alpha*(xm+dx-r0))*aux1_1
                        aux2_2  = -2d0*alpha*dtanh(alpha*(xm-dx-r0))*aux1_2
                        aux4_1  = zexp(-(0d0,1d0)*(kvec(j_x)*(xm+dx)))
                        aux4_2  = zexp(-(0d0,1d0)*(kvec(j_x)*(xm-dx)))
                        aux5_1  = aux4_1*(aux2_1 - 2d0*complex(0d0, 1d0)*kvec(j_x)*aux1_1 - &
                                  (kvec(j_x)**2)*aux0_1)
                        aux5_2  = aux4_2*(aux2_2 - 2d0*complex(0d0,1d0)*kvec(j_x)*aux1_2 - &
                                  (kvec(j_x)**2)*aux0_2)
                        !auxj = zexp((0d0,1d0)*pi*quant_j/2d0)
                        IF (abs(xm + dx) .le. 1e-12) THEN
                                aux6_1 = 0D0
                                aux6_2 = 0D0
                        ELSE
                                aux6_1  = (xm+dx)**(-2)
                                aux6_2  = (xm-dx)**(-2)
                        ENDIF
!
                        sx1 = sx1 + gq_weight_x(i+ngqp_x/2) * &
                              ( aux5_1*bsp_x(i_x,xm+dx) + aux5_2*bsp_x(i_x,xm-dx)) !*auxj
                        sx3 = sx3 + gq_weight_x(i+ngqp_x/2) * &
                              ( aux0_1*aux4_1*bsp_x(i_x,xm+dx) + aux0_2*aux4_2*bsp_x(i_x,xm-dx)) !*auxj  
                        sx4 = sx4 + gq_weight_x(i+ngqp_x/2) * ( aux0_1*aux6_1*aux4_1*bsp_x(i_x,xm+dx) + &
                              aux0_2*aux6_2*aux4_2*bsp_x(i_x,xm-dx) ) !*auxj
                      
                ENDDO
!
                
                BAM_x1(i_x,j_x) = BAM_x1(i_x,j_x) + sx1 * xr * dsqrt(mu_R) / dsqrt(kvec(j_x)*norm(i_x))
                BAM_x3(i_x,j_x) = BAM_x3(i_x,j_x) + sx3 * xr * dsqrt(mu_R) / dsqrt(kvec(j_x)*norm(i_x))
                BAM_x4(i_x,j_x) = BAM_x4(i_x,j_x) + sx4 * xr * dsqrt(mu_R) / dsqrt(kvec(j_x)*norm(i_x))
!
        ENDDO
!
        ENDDO
        ENDDO

!

!Calulation of BAM_xx1, BAM_xx3, and BAM_xx4 matrices for Open-Open coupling
!BAM_xx1(i,j)  : Int[ (rhomin..rhomax)    u_0i(x) * d2/dx2     * u_0j(x) dx ]
!BAM_xx3(i,j)  : Int[ (rhomin..rhomax)    u_0i(x) *    1       * u_0j(x) dx ]
!BAM_xx4(i,j)  : Int[ (rhomin..rhomax)    u_0i(x) *  1/x**2    * u_0j(x) dx ]    
        
!Calculation of BAM_xx1b, BAM_xx3b, and BAM_xx4b matrices for Open-Open coupling
!BAM_xx1b(i,j) : Int[ (rhomin..rhomax)    u_1i(x) * d2/dx2     * u_0        ]                       
!BAM_xx3b(i,j) : Int[ (rhomin..rhomax)    u_1i(x) *    1       * u_0j(x) dx ]
!BAM_xx4b(i,j) : Int[ (rhomin..rhomax)    u_1i(x) *  1/x**2    * u_0j(x) dx ]



        DO i_x = 1, n_open
        DO j_x = 1, n_open
        !quant_j = quant_mat(1, j_x)
        !quant_j_prime = quant_mat(1, i_x)
!
        DO j = pbasst(1)%pb_pa1, pbasst(1)%pb_nbr
                xm = 5D-1*(knots_x(j+1)+knots_x(j))
                xr = 5D-1*(knots_x(j+1)-knots_x(j))
!
                sx1=complex(0D0,0D0)
                sx3=complex(0D0,0D0)
                sx4=complex(0D0,0D0)
                
                sx1b=complex(0D0,0D0)
                sx3b=complex(0D0,0D0)
                sx4b=complex(0D0,0D0)
              
!
!               Perform gaussian quadrature for x
!               ---------------------------------
                DO i = 1, ngqp_x/2
                        dx = xr * gq_root_x(i+ngqp_x/2)
                        aux0_1  = 5d-1*(1d0+dtanh(alpha*(xm+dx-r0)))
                        aux0_2  = 5d-1*(1d0+dtanh(alpha*(xm-dx-r0)))
                        aux1_1  = 5d-1*(alpha/((dcosh(alpha*(xm+dx-r0)))**2))
                        aux1_2  = 5d-1*(alpha/((dcosh(alpha*(xm-dx-r0)))**2))
                        aux2_1  = -2d0*alpha*dtanh(alpha*(xm+dx-r0))*aux1_1
                        aux2_2  = -2d0*alpha*dtanh(alpha*(xm-dx-r0))*aux1_2
                        aux3_1  = zexp(-complex(0d0,1d0)*(kvec(i_x)+kvec(j_x))*(xm+dx))
                        aux3_2  = zexp(-complex(0d0,1d0)*(kvec(i_x)+kvec(j_x))*(xm-dx))
                        aux4_1  = zexp(-complex(0d0,1d0)*(kvec(j_x)-kvec(i_x))*(xm+dx))
                        aux4_2  = zexp(-complex(0d0,1d0)*(kvec(j_x)-kvec(i_x))*(xm-dx))
                        aux6_1  = aux0_1*(aux2_1 - 2d0*complex(0d0,1d0)*kvec(j_x)*aux1_1 - &
                                  (kvec(j_x)**2)*aux0_1)
                        aux6_2  = aux0_2*(aux2_2 - 2d0*complex(0d0,1d0)*kvec(j_x)*aux1_2 - &
                                  (kvec(j_x)**2)*aux0_2)
                        !auxj = zexp((0d0,1d0) * pi * (quant_j + quant_j_prime)/2d0)
                        !auxj_1 = zexp((0d0,1d0)*pi*(quant_j-quant_j_prime)/2d0)
                        IF (abs(xm + dx) .le. 1e-12) THEN
                                aux7_1 = 0D0
                                aux7_2 = 0D0
                        ELSE
                                aux7_1  = (xm+dx)**(-2)
                                aux7_2  = (xm-dx)**(-2)
                        ENDIF

!
                        sx1 = sx1 + gq_weight_x(i+ngqp_x/2) * ( aux3_1*aux6_1 + aux3_2*aux6_2 ) !*auxj
                        sx3 = sx3 + gq_weight_x(i+ngqp_x/2) * ( (aux0_1**2)*aux3_1 + (aux0_2**2)*aux3_2 ) !*auxj
                        sx4 = sx4 + gq_weight_x(i+ngqp_x/2) * ( (aux0_1**2)*aux3_1*aux7_1 + &
                                    (aux0_2**2)*aux3_2*aux7_2 ) !*auxj
                        

                        sx1b = sx1b + gq_weight_x(i+ngqp_x/2) * ( aux4_1*aux6_1 + aux4_2*aux6_2 ) !*auxj_1
                        sx3b = sx3b + gq_weight_x(i+ngqp_x/2) * ( (aux0_1**2)*aux4_1 + (aux0_2**2)*aux4_2 ) !*auxj_1  
                        sx4b = sx4b + gq_weight_x(i+ngqp_x/2) * ( (aux0_1**2)*aux4_1*aux7_1 + &
                                    (aux0_2**2)*aux4_2*aux7_2 ) !*auxj_1
                        
                ENDDO
!
                BAM_xx1(i_x,j_x) = BAM_xx1(i_x,j_x) + sx1 * xr * (mu_R) / dsqrt(kvec(i_x)*kvec(j_x))
                BAM_xx3(i_x,j_x) = BAM_xx3(i_x,j_x) + sx3 * xr * (mu_R) / dsqrt(kvec(i_x)*kvec(j_x))
                BAM_xx4(i_x,j_x) = BAM_xx4(i_x,j_x) + sx4 * xr * (mu_R) / dsqrt(kvec(i_x)*kvec(j_x))
                

                BAM_xx1b(i_x,j_x) = BAM_xx1b(i_x,j_x) + sx1b * xr * (mu_R) / dsqrt(kvec(i_x)*kvec(j_x))
                BAM_xx3b(i_x,j_x) = BAM_xx3b(i_x,j_x) + sx3b * xr * (mu_R) / dsqrt(kvec(i_x)*kvec(j_x))
                BAM_xx4b(i_x,j_x) = BAM_xx4b(i_x,j_x) + sx4b * xr * (mu_R) / dsqrt(kvec(i_x)*kvec(j_x))
                
        ENDDO
!
        ENDDO
        ENDDO

        ELSE
        DO i = 1, pbasst(1)%pb_nbr
                DO j = 1, pbasst(1)%pb_nbr
                DO k = 1, ngqp_x
                        BAM_1(i,j) = BAM_1(i,j) + &
                         wx(k) * dbsp_x(i,x(k)) * dbsp_x(j,x(k))/ DSQRT(norm(i)*norm(j))

                        BAM_3(i,j) = BAM_3(i,j) + &
                         wx(k) * bsp_x(i,x(k)) * bsp_x(j,x(k))/ DSQRT(norm(i)*norm(j))

                        BAM_4(i,j) = BAM_4(i,j) + &
                         wx(k) * bsp_x(i,x(k)) * (1d0/(x(k))**2d0) * bsp_x(j,x(k))/ DSQRT(norm(i)*norm(j))
                ENDDO
                ENDDO
        ENDDO


        DO i = 1, pbasst(1)%pb_nbr
                DO j = 1, n_open
                DO k = 1, ngqp_x
                        BAM_x1(i,j) = BAM_x1(i,j) + &
                         wx(k) * bsp_x(i,x(k)) * d2dR2u0(j,x(k))/ DSQRT(norm(i))

                        BAM_x3(i,j) = BAM_x3(i,j) + &
                         wx(k) * bsp_x(i,x(k)) * u0(j, x(k))/ DSQRT(norm(i))

                        BAM_x4(i,j) = BAM_x4(i,j) + &
                          wx(k) * bsp_x(i,x(k)) * (1d0/(x(k))**(2d0)) * u0(j,x(k))/ DSQRT(norm(i))
                ENDDO
                ENDDO
        ENDDO

        DO n = 1, n_open
                DO n_prime = 1, n_open
                        DO i = 1, ngqp_x
                        BAM_xx1(n,n_prime) = BAM_xx1(n,n_prime) + &
                          wx(i) * u0(n, x(i)) * d2dR2u0(n_prime, x(i))

                        BAM_xx3(n,n_prime) = BAM_xx3(n,n_prime) + &
                         wx(i) * u0(n, x(i)) * u0(n_prime, x(i))

                        BAM_xx4(n,n_prime) = BAM_xx4(n,n_prime) + &
                         wx(i) * u0(n, x(i)) * ((1d0/x(i))**(2d0)) * u0(n_prime, x(i))


                        BAM_xx1b(n,n_prime) = BAM_xx1b(n,n_prime) + &
                         wx(i) * CONJG(u0(n,x(i))) * d2DR2u0(n_prime, x(i))

                        BAM_xx3b(n,n_prime) = BAM_xx3b(n,n_prime) + &
                         wx(i) * CONJG(u0(n,x(i))) * u0(n_prime, x(i))

                        BAM_xx4b(n,n_prime) = BAM_xx4b(n,n_prime) + &
                        wx(i) * CONJG(u0(n,x(i))) * ((1d0/x(i))**(2d0)) * u0(n_prime, x(i))
                        ENDDO
                ENDDO
        ENDDO

!        PRINT*,'BAM_theta2'
!        !print*,'dimension BAM_theta2:', size(BAM_theta2, dim=1), size(BAM_theta2, dim=2)
!        do i=1, ncf
!        write(6,'(8E18.8)') (dble(BAM_theta2(i,j)), j=1, 5)
!        enddo

        !print*,'here 1'

!        PRINT*,'BAM_1'
!        do i=1, min(dim_x+2,5)
!        write(6,'(8E18.8)') (dble(BAM_1(i+1,j+1)), j=1, min(dim_x+2,5))
!        enddo
!        PRINT*,'BAM_3'
!        do i=1, min(dim_x+2,5)
!        write(6,'(8E18.8)') (dble(BAM_3(i+1,j+1)), j=1, min(dim_x+2,5))
!        enddo
!        PRINT*,'BAM_4'
!        do i=1, min(dim_x+2,5)
!        write(6,'(8E18.8)') (dble(BAM_4(i+1,j+1)), j=1, min(dim_x+2,5))
!        enddo

!        PRINT*,'BAM_r1'
!        do i=1, min(dim_x+2,5)
!        write(6,'(8E18.8)') (dble(BAM_r1(i+1,j+1)), j=1, min(dim_x+2,5))
!        enddo
!        PRINT*,'BAM_r2'
!        do i=1, min(dim_x+2,5)
!        write(6,'(8E18.8)') (dble(BAM_r2(i+1,j+1)), j=1, min(dim_x+2,5))
!        enddo
!        PRINT*,'BAM_r3'
!        do i=1, min(dim_x+2,5)
!        write(6,'(8E18.8)') (dble(BAM_r3(i+1,j+1)), j=1, min(dim_x+2,5))
!        enddo

!        PRINT*,'BAM_x1'
!        do i=1, min(N,5)
!        write(6,'(8E18.8)'),(dble(BAM_x1(i,j)), j=1, n_open)
!        enddo
!        PRINT*,'BAM_x3'
!        do i=1, min(N,5)
!        write(6,'(8E18.8)'),(dble(BAM_x3(i,j)), j=1, n_open)
!        enddo
!        PRINT*,'BAM_x4'
!        do i=1, min(N,5)
!        write(6,'(8E18.8)'),(dble(BAM_x4(i,j)), j=1, n_open)
!        enddo

!        PRINT*,'BAM_xx1', SHAPE(BAM_xx1)
!        do i=1, n_open
!        write(6,'(8E18.8)'),(dble(BAM_xx1(i,j)), j=1, n_open)
!        enddo
!         PRINT*,'BAM_xx3', SHAPE(BAM_xx3)
!        do i=1, n_open
!        write(6,'(8E18.8)'),(dble(BAM_xx3(i,j)), j=1, n_open)
!        enddo
!         PRINT*,'BAM_xx4', SHAPE(BAM_xx4)
!        do i=1, n_open
!        write(6,'(8E18.8)'),(dble(BAM_xx4(i,j)), j=1, n_open)
!        enddo

!        PRINT*,'BAM_xx1b', SHAPE(BAM_xx1b)
!        do i=1, n_open
!        write(6,'(8E18.8)'),(dble(BAM_xx1b(i,j)), j=1, n_open)
!        enddo
!         PRINT*,'BAM_xx3b', SHAPE(BAM_xx3b)
!        do i=1, n_open
!        write(6,'(8E18.8)'),(dble(BAM_xx3b(i,j)), j=1, n_open)
!        enddo
!         PRINT*,'BAM_xx4b', SHAPE(BAM_xx4b)
!        do i=1, n_open
!        write(6,'(8E18.8)'),(dble(BAM_xx4b(i,j)), j=1,n_open)
!        enddo


        ! PRINT*, "M_V"
        ! DO i = 1,5
        !         PRINT*, M_V(i, 1:5)
        ! ENDDO

        ! PRINT*, "M0_V"
        ! DO i = 1,3
        !         PRINT*, M0_V(i, :)
        ! ENDDO

        ! Do not deallocate kvec here.
        ! potential_mat_calcul still needs kvec through u0().
        ! DEALLOCATE ( kvec, STAT=istatus)

        ENDIF

        CALL CPU_TIME(tm12)
        PRINT*, " "
        write(6,'(1A,F8.3,1A)') 'time in basic_auxiliary_matrix calculations:',(tm12-tm11),'sec'
        PRINT*, " "
!
!
!=============================================================================================
!*********************************************************************************************
 END SUBROUTINE basic_aux_mat_calcul

 SUBROUTINE basic_aux_mat_calcul_sym_direct
!
!*********************************************************************************************
!
! Direct exchange-symmetric variant used by the direct-only program.
! It keeps the closed-channel B-spline matrices unchanged, but builds
! the open-channel auxiliary matrices in n_sym_open instead of n_open.
!
!=============================================================================================
!
        USE generateparameters
        USE AtomDiatomskvp
        USE omp_lib

        IMPLICIT NONE

        INTEGER :: i, j, i_x, j_x, file_size
        INTEGER :: unit_sum
        INTEGER(8) :: old_x_entries, sym_x_entries
        INTEGER(8) :: old_xx_entries, sym_xx_entries
        REAL(8) :: tm11, tm12, xm, dx, xr
        COMPLEX(8) :: sx1, sx3, sx4, sx1b, sx3b, sx4b
        COMPLEX(8) :: aux0_1, aux0_2, aux1_1, aux1_2, aux2_1, aux2_2
        COMPLEX(8) :: aux3_1, aux3_2, aux4_1, aux4_2
        COMPLEX(8) :: aux5_1, aux5_2, aux6_1, aux6_2, aux7_1, aux7_2
        LOGICAL :: file_exists, matrices_exist, new_calc, same_param

        CALL CPU_TIME(tm11)

        new_calc = .TRUE.
        INQUIRE(file='previous_calc_sym.dat', exist=file_exists, size=file_size)
        INQUIRE(file='matrices_sym.bin', exist=matrices_exist)

        IF (file_exists .eqv. .FALSE.) new_calc = .TRUE.
        IF ((file_exists .eqv. .TRUE.) .AND. (file_size == 0)) new_calc = .TRUE.

        IF ((file_exists .eqv. .TRUE.) .AND. (file_size /= 0) .AND. matrices_exist) THEN
                OPEN(unit=17, file='previous_calc_sym.dat', form='unformatted', status='old')
                READ(17) n_xmb, k_xmb, jmxmb, n_openmb
                same_param = ((pbasst(1)%pb_nbr == n_xmb) .AND. &
                              (pbasst(1)%pb_pa1 == k_xmb) .AND. &
                              (pbasst(2)%pb_nbr == jmxmb) .AND. &
                              (n_sym_open == n_openmb))
                IF (same_param .eqv. .TRUE.) new_calc = .FALSE.
                CLOSE(17)
        ENDIF

        IF (new_calc .eqv. .TRUE.) THEN
                PRINT*, " "
                PRINT*, "running direct symmetric noc independent calc"
                OPEN(unit=17, file='previous_calc_sym.dat', form='unformatted', status='replace')
                WRITE(17) pbasst(1)%pb_nbr, pbasst(1)%pb_pa1, pbasst(2)%pb_nbr, n_sym_open
                CLOSE(17)

                ALLOCATE(BAM_1(1:pbasst(1)%pb_nbr,1:pbasst(1)%pb_nbr), &
                         BAM_3(1:pbasst(1)%pb_nbr,1:pbasst(1)%pb_nbr), &
                         BAM_4(1:pbasst(1)%pb_nbr,1:pbasst(1)%pb_nbr), &
                         norm(1:pbasst(1)%pb_nbr), STAT=istatus)
                IF (istatus /= 0) STOP 'Allocation failed for closed aux matrices'

                norm = 0d0
                DO i_x = 1, pbasst(1)%pb_nbr
                DO j = i_x, i_x+pbasst(1)%pb_pa1-1
                        xm = 5D-1*(knots_x(j+1)+knots_x(j))
                        xr = 5D-1*(knots_x(j+1)-knots_x(j))
                        sx1 = (0D0,0D0)

                        DO i = 1, ngqp_x/2
                                dx = xr * gq_root_x(i+ngqp_x/2)
                                sx1 = sx1 + gq_weight_x(i+ngqp_x/2) * &
                                      (bsp_x(i_x,xm+dx)**2d0 + bsp_x(i_x,xm-dx)**2d0)
                        ENDDO
                        norm(i_x) = norm(i_x) + DBLE(sx1 * xr)
                ENDDO
                ENDDO

                BAM_1 = (0D0,0D0)
                BAM_3 = (0D0,0D0)
                BAM_4 = (0D0,0D0)

                DO i_x = 1, pbasst(1)%pb_nbr
                DO j_x = 1, MIN(i_x+pbasst(1)%pb_pa1-1,pbasst(1)%pb_nbr)
                DO j = j_x, i_x+pbasst(1)%pb_pa1-1
                        xm = 5D-1*(knots_x(j+1)+knots_x(j))
                        xr = 5D-1*(knots_x(j+1)-knots_x(j))
                        sx1 = (0D0,0D0)
                        sx3 = (0D0,0D0)
                        sx4 = (0D0,0D0)

                        DO i = 1, ngqp_x/2
                                dx = xr * gq_root_x(i+ngqp_x/2)
                                IF (ABS(xm + dx) .LE. 1d-12) THEN
                                        aux4_1 = 0D0
                                        aux4_2 = 0D0
                                ELSE
                                        aux4_1 = (xm+dx)**(-2)
                                        aux4_2 = (xm-dx)**(-2)
                                ENDIF

                                sx1 = sx1 + gq_weight_x(i+ngqp_x/2) * &
                                      (dbsp_x(i_x,xm+dx)*dbsp_x(j_x,xm+dx) + &
                                       dbsp_x(i_x,xm-dx)*dbsp_x(j_x,xm-dx))
                                sx3 = sx3 + gq_weight_x(i+ngqp_x/2) * &
                                      (bsp_x(i_x,xm+dx)*bsp_x(j_x,xm+dx) + &
                                       bsp_x(i_x,xm-dx)*bsp_x(j_x,xm-dx))
                                sx4 = sx4 + gq_weight_x(i+ngqp_x/2) * &
                                      (aux4_1*bsp_x(i_x,xm+dx)*bsp_x(j_x,xm+dx) + &
                                       aux4_2*bsp_x(i_x,xm-dx)*bsp_x(j_x,xm-dx))
                        ENDDO

                        BAM_1(i_x,j_x) = BAM_1(i_x,j_x) + sx1*xr/DSQRT(norm(i_x)*norm(j_x))
                        BAM_3(i_x,j_x) = BAM_3(i_x,j_x) + sx3*xr/DSQRT(norm(i_x)*norm(j_x))
                        BAM_4(i_x,j_x) = BAM_4(i_x,j_x) + sx4*xr/DSQRT(norm(i_x)*norm(j_x))
                ENDDO
                ENDDO
                ENDDO

                OPEN(unit=10, file='matrices_sym.bin', form='unformatted', &
                     access='stream', status='replace')
                WRITE(10) norm
                WRITE(10) BAM_1
                WRITE(10) BAM_3
                WRITE(10) BAM_4
                CLOSE(10)
        ELSE
                ALLOCATE(BAM_1(1:pbasst(1)%pb_nbr,1:pbasst(1)%pb_nbr), &
                         BAM_3(1:pbasst(1)%pb_nbr,1:pbasst(1)%pb_nbr), &
                         BAM_4(1:pbasst(1)%pb_nbr,1:pbasst(1)%pb_nbr), &
                         norm(1:pbasst(1)%pb_nbr), STAT=istatus)
                IF (istatus /= 0) STOP 'Allocation failed for cached closed aux matrices'
        ENDIF

        OPEN(unit=10, file='matrices_sym.bin', form='unformatted', access='stream', status='old')
        READ(10) norm
        READ(10) BAM_1
        READ(10) BAM_3
        READ(10) BAM_4
        CLOSE(10)

        IF (.NOT. ALLOCATED(sym_kvec)) THEN
                PRINT*, 'Error: sym_kvec is not allocated before basic_aux_mat_calcul_sym_direct.'
                STOP
        ENDIF

        IF (ALLOCATED(BAM_x1_sym)) THEN
                DEALLOCATE(BAM_x1_sym, BAM_x3_sym, BAM_x4_sym, STAT=istatus)
        ENDIF
        IF (ALLOCATED(BAM_xx1_sym)) THEN
                DEALLOCATE(BAM_xx1_sym, BAM_xx3_sym, BAM_xx4_sym, &
                           BAM_xx1b_sym, BAM_xx3b_sym, BAM_xx4b_sym, STAT=istatus)
        ENDIF

        ALLOCATE(BAM_x1_sym(1:pbasst(1)%pb_nbr,1:n_sym_open), &
                 BAM_x3_sym(1:pbasst(1)%pb_nbr,1:n_sym_open), &
                 BAM_x4_sym(1:pbasst(1)%pb_nbr,1:n_sym_open), STAT=istatus)
        IF (istatus /= 0) STOP 'Allocation failed for symmetric BAM_x matrices'

        ALLOCATE(BAM_xx1_sym(1:n_sym_open,1:n_sym_open), &
                 BAM_xx3_sym(1:n_sym_open,1:n_sym_open), &
                 BAM_xx4_sym(1:n_sym_open,1:n_sym_open), &
                 BAM_xx1b_sym(1:n_sym_open,1:n_sym_open), &
                 BAM_xx3b_sym(1:n_sym_open,1:n_sym_open), &
                 BAM_xx4b_sym(1:n_sym_open,1:n_sym_open), STAT=istatus)
        IF (istatus /= 0) STOP 'Allocation failed for symmetric BAM_xx matrices'

        BAM_x1_sym = (0D0,0D0)
        BAM_x3_sym = (0D0,0D0)
        BAM_x4_sym = (0D0,0D0)
        BAM_xx1_sym = (0D0,0D0)
        BAM_xx3_sym = (0D0,0D0)
        BAM_xx4_sym = (0D0,0D0)
        BAM_xx1b_sym = (0D0,0D0)
        BAM_xx3b_sym = (0D0,0D0)
        BAM_xx4b_sym = (0D0,0D0)

        old_x_entries = 3_8 * INT(pbasst(1)%pb_nbr,8) * INT(n_open,8)
        sym_x_entries = 3_8 * INT(pbasst(1)%pb_nbr,8) * INT(n_sym_open,8)
        old_xx_entries = 6_8 * INT(n_open,8) * INT(n_open,8)
        sym_xx_entries = 6_8 * INT(n_sym_open,8) * INT(n_sym_open,8)

        OPEN(newunit=unit_sum, file='basic_sym_summary.dat', &
             status='unknown', position='append', action='write')
        WRITE(unit_sum,'(A)') &
                '# Energy_eV Jtot n_open n_sym_open old_x sym_x ratio_x old_xx sym_xx ratio_xx'
        WRITE(unit_sum,'(ES20.10,1X,I5,2I8,2I14,ES20.10,2I14,ES20.10)') &
                E*27.211399d0, Jtot, n_open, n_sym_open, &
                old_x_entries, sym_x_entries, &
                DBLE(sym_x_entries)/DBLE(old_x_entries), &
                old_xx_entries, sym_xx_entries, &
                DBLE(sym_xx_entries)/DBLE(old_xx_entries)
        CLOSE(unit_sum)

        IF (.NOT. bsp_optimized) THEN
                PRINT*, 'Error: basic_aux_mat_calcul_sym_direct currently expects bsp_optimized=.TRUE.'
                STOP
        ENDIF

        PRINT*, " "
        PRINT*, "Running direct symmetric B-spline optimized program"
        PRINT*, " "

        DO i_x = 1, pbasst(1)%pb_nbr
        DO j_x = 1, n_sym_open
        DO j = pbasst(1)%pb_pa1, pbasst(1)%pb_nbr
                xm = 5D-1*(knots_x(j+1)+knots_x(j))
                xr = 5D-1*(knots_x(j+1)-knots_x(j))
                sx1 = (0D0,0D0)
                sx3 = (0D0,0D0)
                sx4 = (0D0,0D0)

                DO i = 1, ngqp_x/2
                        dx = xr * gq_root_x(i+ngqp_x/2)
                        aux0_1 = 5d-1*(1d0+DTANH(alpha*(xm+dx-r0)))
                        aux0_2 = 5d-1*(1d0+DTANH(alpha*(xm-dx-r0)))
                        aux1_1 = 5d-1*(alpha/((DCOSH(alpha*(xm+dx-r0)))**2))
                        aux1_2 = 5d-1*(alpha/((DCOSH(alpha*(xm-dx-r0)))**2))
                        aux2_1 = -2d0*alpha*DTANH(alpha*(xm+dx-r0))*aux1_1
                        aux2_2 = -2d0*alpha*DTANH(alpha*(xm-dx-r0))*aux1_2
                        aux4_1 = ZEXP(-(0d0,1d0)*(sym_kvec(j_x)*(xm+dx)))
                        aux4_2 = ZEXP(-(0d0,1d0)*(sym_kvec(j_x)*(xm-dx)))
                        aux5_1 = aux4_1*(aux2_1 - 2d0*(0d0,1d0)*sym_kvec(j_x)*aux1_1 - &
                                  (sym_kvec(j_x)**2)*aux0_1)
                        aux5_2 = aux4_2*(aux2_2 - 2d0*(0d0,1d0)*sym_kvec(j_x)*aux1_2 - &
                                  (sym_kvec(j_x)**2)*aux0_2)
                        IF (ABS(xm + dx) .LE. 1d-12) THEN
                                aux6_1 = 0D0
                                aux6_2 = 0D0
                        ELSE
                                aux6_1 = (xm+dx)**(-2)
                                aux6_2 = (xm-dx)**(-2)
                        ENDIF

                        sx1 = sx1 + gq_weight_x(i+ngqp_x/2) * &
                              (aux5_1*bsp_x(i_x,xm+dx) + aux5_2*bsp_x(i_x,xm-dx))
                        sx3 = sx3 + gq_weight_x(i+ngqp_x/2) * &
                              (aux0_1*aux4_1*bsp_x(i_x,xm+dx) + &
                               aux0_2*aux4_2*bsp_x(i_x,xm-dx))
                        sx4 = sx4 + gq_weight_x(i+ngqp_x/2) * &
                              (aux0_1*aux6_1*aux4_1*bsp_x(i_x,xm+dx) + &
                               aux0_2*aux6_2*aux4_2*bsp_x(i_x,xm-dx))
                ENDDO

                BAM_x1_sym(i_x,j_x) = BAM_x1_sym(i_x,j_x) + &
                        sx1*xr*DSQRT(mu_R)/DSQRT(sym_kvec(j_x)*norm(i_x))
                BAM_x3_sym(i_x,j_x) = BAM_x3_sym(i_x,j_x) + &
                        sx3*xr*DSQRT(mu_R)/DSQRT(sym_kvec(j_x)*norm(i_x))
                BAM_x4_sym(i_x,j_x) = BAM_x4_sym(i_x,j_x) + &
                        sx4*xr*DSQRT(mu_R)/DSQRT(sym_kvec(j_x)*norm(i_x))
        ENDDO
        ENDDO
        ENDDO

        DO i_x = 1, n_sym_open
        DO j_x = 1, n_sym_open
        DO j = pbasst(1)%pb_pa1, pbasst(1)%pb_nbr
                xm = 5D-1*(knots_x(j+1)+knots_x(j))
                xr = 5D-1*(knots_x(j+1)-knots_x(j))
                sx1 = (0D0,0D0)
                sx3 = (0D0,0D0)
                sx4 = (0D0,0D0)
                sx1b = (0D0,0D0)
                sx3b = (0D0,0D0)
                sx4b = (0D0,0D0)

                DO i = 1, ngqp_x/2
                        dx = xr * gq_root_x(i+ngqp_x/2)
                        aux0_1 = 5d-1*(1d0+DTANH(alpha*(xm+dx-r0)))
                        aux0_2 = 5d-1*(1d0+DTANH(alpha*(xm-dx-r0)))
                        aux1_1 = 5d-1*(alpha/((DCOSH(alpha*(xm+dx-r0)))**2))
                        aux1_2 = 5d-1*(alpha/((DCOSH(alpha*(xm-dx-r0)))**2))
                        aux2_1 = -2d0*alpha*DTANH(alpha*(xm+dx-r0))*aux1_1
                        aux2_2 = -2d0*alpha*DTANH(alpha*(xm-dx-r0))*aux1_2
                        aux3_1 = ZEXP(-(0d0,1d0)*(sym_kvec(i_x)+sym_kvec(j_x))*(xm+dx))
                        aux3_2 = ZEXP(-(0d0,1d0)*(sym_kvec(i_x)+sym_kvec(j_x))*(xm-dx))
                        aux4_1 = ZEXP(-(0d0,1d0)*(sym_kvec(j_x)-sym_kvec(i_x))*(xm+dx))
                        aux4_2 = ZEXP(-(0d0,1d0)*(sym_kvec(j_x)-sym_kvec(i_x))*(xm-dx))
                        aux6_1 = aux0_1*(aux2_1 - 2d0*(0d0,1d0)*sym_kvec(j_x)*aux1_1 - &
                                  (sym_kvec(j_x)**2)*aux0_1)
                        aux6_2 = aux0_2*(aux2_2 - 2d0*(0d0,1d0)*sym_kvec(j_x)*aux1_2 - &
                                  (sym_kvec(j_x)**2)*aux0_2)
                        IF (ABS(xm + dx) .LE. 1d-12) THEN
                                aux7_1 = 0D0
                                aux7_2 = 0D0
                        ELSE
                                aux7_1 = (xm+dx)**(-2)
                                aux7_2 = (xm-dx)**(-2)
                        ENDIF

                        sx1 = sx1 + gq_weight_x(i+ngqp_x/2) * (aux3_1*aux6_1 + aux3_2*aux6_2)
                        sx3 = sx3 + gq_weight_x(i+ngqp_x/2) * ((aux0_1**2)*aux3_1 + &
                              (aux0_2**2)*aux3_2)
                        sx4 = sx4 + gq_weight_x(i+ngqp_x/2) * ((aux0_1**2)*aux3_1*aux7_1 + &
                              (aux0_2**2)*aux3_2*aux7_2)
                        sx1b = sx1b + gq_weight_x(i+ngqp_x/2) * (aux4_1*aux6_1 + aux4_2*aux6_2)
                        sx3b = sx3b + gq_weight_x(i+ngqp_x/2) * ((aux0_1**2)*aux4_1 + &
                               (aux0_2**2)*aux4_2)
                        sx4b = sx4b + gq_weight_x(i+ngqp_x/2) * ((aux0_1**2)*aux4_1*aux7_1 + &
                               (aux0_2**2)*aux4_2*aux7_2)
                ENDDO

                BAM_xx1_sym(i_x,j_x) = BAM_xx1_sym(i_x,j_x) + &
                        sx1*xr*mu_R/DSQRT(sym_kvec(i_x)*sym_kvec(j_x))
                BAM_xx3_sym(i_x,j_x) = BAM_xx3_sym(i_x,j_x) + &
                        sx3*xr*mu_R/DSQRT(sym_kvec(i_x)*sym_kvec(j_x))
                BAM_xx4_sym(i_x,j_x) = BAM_xx4_sym(i_x,j_x) + &
                        sx4*xr*mu_R/DSQRT(sym_kvec(i_x)*sym_kvec(j_x))
                BAM_xx1b_sym(i_x,j_x) = BAM_xx1b_sym(i_x,j_x) + &
                        sx1b*xr*mu_R/DSQRT(sym_kvec(i_x)*sym_kvec(j_x))
                BAM_xx3b_sym(i_x,j_x) = BAM_xx3b_sym(i_x,j_x) + &
                        sx3b*xr*mu_R/DSQRT(sym_kvec(i_x)*sym_kvec(j_x))
                BAM_xx4b_sym(i_x,j_x) = BAM_xx4b_sym(i_x,j_x) + &
                        sx4b*xr*mu_R/DSQRT(sym_kvec(i_x)*sym_kvec(j_x))
        ENDDO
        ENDDO
        ENDDO

        CALL CPU_TIME(tm12)
        PRINT*, " "
        WRITE(6,'(1A,F8.3,1A)') &
                'time in direct symmetric basic_auxiliary_matrix calculations:', &
                (tm12-tm11),'sec'
        PRINT*, " "

 END SUBROUTINE basic_aux_mat_calcul_sym_direct
!*********************************************************************************************
!=============================================================================================
!
!
!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************
!
 SUBROUTINE solve_target_levels
!
!*********************************************************************************************


! This subroutine expresses the scattering wavefunction from a model case of SKVP
!
!=============================================================================================
!
        USE generateparameters
        USE AtomDiatomskvp
        USE omp_lib
!
        IMPLICIT NONE
!--------------------
        INTEGER :: i, quant_j1, quant_j2
        REAL(8) :: E_rot
!
!       Open Channels
!       -------------
        n_open = 0 
        DO i = 1, ncf
             quant_j1 = quant_mat(1, i)
             quant_j2 = quant_mat(3, i)

             E_rot = Brot*quant_j1*(quant_j1+1d0) + Brot*quant_j2*(quant_j2+1d0)

             IF (E_rot < E) THEN
                n_open = n_open + 1
             ENDIF
        ENDDO
!
        PRINT*, "  "
        PRINT*, "Number of open channels: ", n_open
        PRINT*, "  "
!
        IF (ALLOCATED(open_idx)) DEALLOCATE(open_idx, STAT=istatus)
       ALLOCATE(open_idx(n_open), STAT=istatus)

       n_open = 0

       DO i = 1, ncf
               quant_j1 = quant_mat(1,i)
               quant_j2 = quant_mat(3,i)
                E_rot = Brot*quant_j1*(quant_j1+1d0) + &
                                Brot*quant_j2*(quant_j2+1d0)


                        IF (E_rot < E) THEN
                                n_open = n_open + 1
                                open_idx(n_open) = i
                        ENDIF
        ENDDO


       ALLOCATE(kvec(1:n_open), STAT=istatus)


         DO i = 1, n_open
               quant_j1 = quant_mat(1,open_idx(i))
               quant_j2 = quant_mat(3,open_idx(i))


               E_rot = Brot*quant_j1*(quant_j1+1d0) + &
                       Brot*quant_j2*(quant_j2+1d0)


               kvec(i) = DSQRT(2d0*mu_R*(E-E_rot))
       ENDDO

        PRINT *, 'Open channel map: open index -> quant_mat index -> j1 k1 j2 k2'
       DO i = 1, n_open
               
        PRINT *, i, open_idx(i), &
                       quant_mat(1,open_idx(i)), quant_mat(2,open_idx(i)), &
                       quant_mat(3,open_idx(i)), quant_mat(4,open_idx(i)), &
                       ' kvec=', kvec(i)
       ENDDO


        !   
!
!=============================================================================================
!*********************************************************************************************
 END SUBROUTINE solve_target_levels
!*********************************************************************************************
!=============================================================================================
