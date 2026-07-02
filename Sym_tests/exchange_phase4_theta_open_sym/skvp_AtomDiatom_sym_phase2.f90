!*********************************************************************************************
!*********************************************************************************************
!
PROGRAM skvpAtomDiatom
!
!*********************************************************************************************
!
!
!=============================================================================================
!
!	Date		:	04.04.2022
!	Author		:	Steve Ndengue 
!	Last Change	:	03.29.2026
!
!	The present program is aimed at applying the S-Matrix Kohn Variational Principle 
!	       on the 2D scattering problem using the traditional method from Miller.
! 
!=============================================================================================		
! 
 USE omp_lib
 USE generateparameters
 USE AtomDiatomskvp
 USE Potential_Interface
!


 IMPLICIT NONE
!-------------
!
	 INTEGER :: i, j
	 INTEGER :: t
	 INTEGER :: j1, j2, j1low, j2low
	 INTEGER :: k1, k2, k1range, k2range
	 INTEGER                                  :: step, n_steps 
	INTEGER, PARAMETER                       :: proba_all_unit = 101
        INTEGER                                  :: N_phase2_sym
        COMPLEX(8), ALLOCATABLE                 :: Smat_sym_projected_phase2(:,:)
        COMPLEX(8), ALLOCATABLE                 :: Smat_sym_direct_phase2(:,:)

 REAL(8)                                  :: tm1, tm2, pas_x, yy=0.020d0
 REAL(8)                                  :: x1, x2              !, DLAMCH
 !COMPLEX(8)                               :: varc
 !REAL(8), ALLOCATABLE, DIMENSION(:)       :: var
 !COMPLEX(8), ALLOCATABLE, DIMENSION(:)    :: biatx
!
! Call for input paramaters reading routine
!------------------------------------------
!      
CALL read_input

!CALL set_potential_backend('BMKP')
!CALL set_bmkp_filename('/Users/yuan/Documents/skvp_diatomdiaton_firstdrfat/coefficients.dat')


CALL build_potential_index

CALL calculate_RMS_potential_expansion

                DO t = 1, n_pot
                WRITE(*,'(4I6)') t, pot_mat(1,t), pot_mat(2,t), pot_mat(3,t)
                ENDDO

        CALL CPU_TIME(tm1)
!
        OPEN (100,file ='proba.dat',status='unknown',action='write',position='append')
        OPEN (proba_all_unit,file='proba_all.dat',status='unknown', action='write',position='append')
!
!
! Call for input paramaters reading routine
!------------------------------------------
!       
        dim_x    = pbasst(1)%pb_nbr-2
        ngqp_x   = p*pbasst(1)%pb_pa1  !2*int(rhomax-rhomin)*p
! 
! Determination of the knots sequences
!-------------------------------------
        ! Allocation of matrices
        ALLOCATE ( knots_x(1:(pbasst(1)%pb_nbr+pbasst(1)%pb_pa1)), gq_root_x(1:ngqp_x), gq_weight_x(1:ngqp_x), STAT = istatus )
        ! Determination of the abcisses and weights for gaussian quadrature
        CALL gauleg(-1D0,1D0,gq_root_x,gq_weight_x,ngqp_x)
        WRITE (1, 590) ! Write gaussian quadrature points in global output file
        DO i = 1, ngqp_x
                WRITE (1, 595)  gq_root_x(i), gq_weight_x(i)
        ENDDO
        ! y = 1d0
        pas_x = 0d0
        knots_x = 0d0
        pas_x=(pbasst(1)%pb_max-pbasst(1)%pb_min)/(dexp(yy*dble(pbasst(1)%pb_nbr-pbasst(1)%pb_pa1+1)/ &
               dble(pbasst(1)%pb_nbr-pbasst(1)%pb_pa1+2))-dble(1))
        knots_x=(/(pbasst(1)%pb_min,i=1,pbasst(1)%pb_pa1), &
                  ((dexp(yy*dble(i-1)/dble(pbasst(1)%pb_nbr-pbasst(1)%pb_pa1+2))-dble(1))*pas_x+ &
                  pbasst(1)%pb_min,i=2,pbasst(1)%pb_nbr-pbasst(1)%pb_pa1+1), &
                  (pbasst(1)%pb_max,i=pbasst(1)%pb_nbr+1,pbasst(1)%pb_nbr+pbasst(1)%pb_pa1)/)
        ! pas_x=(pbasst(1)%pb_max-pbasst(1)%pb_min)/dble(pbasst(1)%pb_nbr-pbasst(1)%pb_pa1+1)
        ! knots_x=(/(pbasst(1)%pb_min,i=1,pbasst(1)%pb_pa1),(dble(i)*pas_x+pbasst(1)%pb_min,i=1,pbasst(1)%pb_nbr-pbasst(1)%pb_pa1), &
        !           (pbasst(1)%pb_max,i=pbasst(1)%pb_nbr+1,pbasst(1)%pb_nbr+pbasst(1)%pb_pa1)/)
        PRINT*,'knots_x'
        DO i = 1, SIZE(knots_x)
        IF (i < SIZE(knots_x)) THEN
            WRITE(*,'(F6.2,A)', ADVANCE='NO') knots_x(i), ','
        ELSE
            WRITE(*,'(F6.2)') knots_x(i)
        END IF
        END DO


	        ! quant_mat is generated inside the energy/J loops below.

!
! Calculate the number of steps (using Nearest Integer to be safe) in the energy range and loop
!----------------------------------------------------------------------------------------------
        n_steps = nint((en_final - en_start) / en_step)
        ! Energy loop -- loopnmax = 1001
        DO step = 0, n_steps
        E = en_start + (real(step,8)*en_step)
        !E = dble(loopn-1)*(0.0025d0/10d0)
        !!!E = 0.000d0 + dble(loopn-1)*(0.040d0/1000d0)
        !E = 0.0d0 + dble(loopnmax-loopn+1)*(0.0367493d0/10d0)

        IF (ALLOCATED(Xsec_jpair)) DEALLOCATE(Xsec_jpair, STAT=istatus)
        ALLOCATE(Xsec_jpair(0:pbasst(2)%pb_nbr, 0:pbasst(3)%pb_nbr), STAT=istatus)
        Xsec_jpair = 0d0
!
! Computation of the targets wavefunctions and loop over Jtot
!------------------------------------------------------------
        ! Jtot loop
        DO Jtot = Jmin, Jmax
        ! (j,k) -> n, Compacting this set of quantum numbers into one quantum number 
        !!!ncf = (((2*Jtot + 1)*pbasst(2)%pb_nbr)/2)+1
        ! Compute ncf

        ncf=0
        do j1=0, pbasst(2)%pb_nbr, 2 !=> Revisit later
                j1low=min(j1,pbasst(2)%pb_pa1)
                k1range=j1low
                !k1range=min(j1low,Jtot)
                do k1=-k1range, k1range
                        do j2=0, pbasst(3)%pb_nbr, 2 !=> Revisit later
                                j2low=min(j2,pbasst(3)%pb_pa1)
                                k2range=j2low
                                !k2range=min(j2low,Jtot)
                                do k2=-k2range, k2range
                                        ncf=ncf+1
                                enddo    
                        enddo
                enddo
        enddo
	        ! Compute quant_mat
	        ALLOCATE(quant_mat(4,ncf))
	        ncf = 0
        do j1=0, pbasst(2)%pb_nbr, 2 !=> Revisit later
                j1low=min(j1,pbasst(2)%pb_pa1)
                k1range=j1low
                !k1range=min(j1low,Jtot)
                do k1=-k1range, k1range
                        do j2=0, pbasst(3)%pb_nbr, 2 !=> Revisit later
                                j2low=min(j2,pbasst(3)%pb_pa1)
                                k2range=j2low
                                !k2range=min(j2low,Jtot)
	                                do k2=-k2range, k2range
	                                        ncf=ncf+1
	                                        quant_mat(1,ncf)=j1
	                                        quant_mat(2,ncf)=k1
	                                        quant_mat(3,ncf)=j2
	                                        quant_mat(4,ncf)=k2
	                                enddo
	                        enddo
	                enddo
	        enddo

	        CALL build_exchange_symmetric_basis


        ! Print parameters
        PRINT*, " "
        PRINT*, "ngqp_x:", ngqp_x
        PRINT*, "ncf", ncf
        PRINT*, "ENERGY: ", E
        PRINT*, "alpha: ", alpha
        PRINT*, "R0: ", r0
        PRINT*, "rhomin: ", pbasst(1)%pb_min
        PRINT*, "rhomax", pbasst(1)%pb_max
        PRINT*, "mu_R: ", mu_R
        PRINT*, "Brot", Brot
        PRINT*, "Jtot", Jtot
        PRINT*, "quant_mat"
!       PRINT*, quant_mat(1, :) !quant_mat(1, 1:10)
 !       PRINT*, quant_mat(2, :) !quant_mat(2, 1:10)
 !       PRINT*, quant_mat(3, :) !quant_mat(3, 1:10)
 !       PRINT*, quant_mat(4, :) !quant_mat(4, 1:10)
!
        !ALLOCATE(quant_mat(2,ncf))
        !nn = 1
        !DO j = 0, pbasst(2)%pb_nbr, 2
        !   IF (j == 0) THEN
        !      quant_mat(1,nn) = 0
        !      quant_mat(2,nn) = 0
        !      nn = nn + 1
        !   ELSE
        !      DO k = -Jtot,Jtot
        !         quant_mat(1,nn) = j
        !         quant_mat(2,nn) = k
        !         nn = nn + 1
        !      ENDDO
        !   ENDIF
        !ENDDO
        !PRINT*," " 
        !!PRINT*, "bsp(0, 2.1)", bsp_x(0, 2.0d0)
        !PRINT*, "quant_mat"
        !PRINT*, quant_mat(1, 1:10)
        !PRINT*, quant_mat(2, 1:10)
!
	        CALL solve_target_levels
	        CALL build_exchange_symmetric_open_basis
	        x1 = pbasst(1)%pb_min   
	        x2 = pbasst(1)%pb_max
        ALLOCATE(x(1:ngqp_x), wx(1:ngqp_x), STAT = istatus) 
        CALL gauleg(x1, x2, x, wx, ngqp_x)
!
	! Phase-2 exchange-symmetry validation: solve S matrices only.
	! -------------------------------------------------------------
	! First build the phase-1 projected reference from the old ordered
	! matrices. Then build the exchange-symmetric matrices directly and
	! compare the two symmetric-basis S matrices.
	        CALL basic_aux_mat_calcul

	        N_phase2_sym = dim_x * nsym

	        CALL potential_mat_calcul
	        CALL make_scatt_mat

	        CALL project_matrices_to_exchange_symmetric

	        ALLOCATE(Smat_sym_projected_phase2(1:n_sym_open,1:n_sym_open), &
	                 STAT=istatus)
	        IF (istatus /= 0) STOP 'Allocation failed for Smat_sym_projected_phase2'

	        CALL SolveSMatrixGeneric(N_phase2_sym, n_sym_open, mat_M_sym, &
	                                 mat_M0_sym, mat_M00_sym, mat_M10_sym, &
	                                 Smat_sym_projected_phase2, &
	                                 'projected exchange-symmetric reference')

	        CALL cleanup_phase2_matrix_storage(.TRUE.)

	        CALL potential_mat_calcul_sym_direct
	        CALL make_scatt_mat_sym_direct

	        ALLOCATE(Smat_sym_direct_phase2(1:n_sym_open,1:n_sym_open), &
	                 STAT=istatus)
	        IF (istatus /= 0) STOP 'Allocation failed for Smat_sym_direct_phase2'

	        CALL SolveSMatrixGeneric(N_phase2_sym, n_sym_open, mat_M_sym, &
	                                 mat_M0_sym, mat_M00_sym, mat_M10_sym, &
	                                 Smat_sym_direct_phase2, &
	                                 'direct exchange-symmetric')

	        CALL write_exchange_phase2_outputs(Smat_sym_projected_phase2, &
	                                           Smat_sym_direct_phase2)

	        DEALLOCATE(Smat_sym_projected_phase2, Smat_sym_direct_phase2, &
	                   STAT=istatus)

        !DEALLOCATE(knots_x, gq_root_x, gq_weight_x, knots_y, gq_root_y, gq_weight_y, STAT = istatus)

        DEALLOCATE(quant_mat, STAT = istatus)
        DEALLOCATE(BAM_1, BAM_3, BAM_4, norm, STAT = istatus)
        DEALLOCATE(BAM_x1, BAM_x3, BAM_x4, STAT = istatus)
        DEALLOCATE(BAM_xx1, BAM_xx3, BAM_xx4, BAM_xx1b, BAM_xx3b, BAM_xx4b, STAT = istatus)
!=============================================================================================
! OLD CODE (COMMENTED OUT): M_V has already been moved into mat_M, and M_K no
! longer exists. The old unconditional DEALLOCATE could therefore fail.
!        DEALLOCATE(M_V, M0_V, M00_V, M10_V, M_K, Smat, STAT = istatus)
!        DEALLOCATE(mat_M, mat_M0, mat_M00, mat_M10, STAT = istatus)
!=============================================================================================
! NEW OPTIMIZED CODE: deallocate each remaining matrix only when allocated.
        IF (ALLOCATED(M_V))   DEALLOCATE(M_V,   STAT = istatus)
        IF (ALLOCATED(M0_V))  DEALLOCATE(M0_V,  STAT = istatus)
        IF (ALLOCATED(M00_V)) DEALLOCATE(M00_V, STAT = istatus)
        IF (ALLOCATED(M10_V)) DEALLOCATE(M10_V, STAT = istatus)
        IF (ALLOCATED(Smat))  DEALLOCATE(Smat,  STAT = istatus)
        IF (ALLOCATED(mat_M))   DEALLOCATE(mat_M,   STAT = istatus)
        IF (ALLOCATED(mat_M0))  DEALLOCATE(mat_M0,  STAT = istatus)
        IF (ALLOCATED(mat_M00)) DEALLOCATE(mat_M00, STAT = istatus)
        IF (ALLOCATED(mat_M10)) DEALLOCATE(mat_M10, STAT = istatus)
        IF (ALLOCATED(mat_M_sym))   DEALLOCATE(mat_M_sym,   STAT = istatus)
        IF (ALLOCATED(mat_M0_sym))  DEALLOCATE(mat_M0_sym,  STAT = istatus)
        IF (ALLOCATED(mat_M00_sym)) DEALLOCATE(mat_M00_sym, STAT = istatus)
        IF (ALLOCATED(mat_M10_sym)) DEALLOCATE(mat_M10_sym, STAT = istatus)

        IF (ALLOCATED(kvec)) THEN
                DEALLOCATE(kvec, STAT = istatus)
        ENDIF

        IF (ALLOCATED(BAM_r)) THEN
        DEALLOCATE(BAM_r, STAT = istatus)
        ENDIF

        IF (ALLOCATED(BAM_r0)) THEN
                DEALLOCATE(BAM_r0, STAT = istatus)
        ENDIF

        IF (ALLOCATED(BAM_r00)) THEN
                DEALLOCATE(BAM_r00, STAT = istatus)
        ENDIF

        IF (ALLOCATED(BAM_r10)) THEN
                DEALLOCATE(BAM_r10, STAT = istatus)
        ENDIF

        IF (ALLOCATED(BAM_theta)) THEN
                DEALLOCATE(BAM_theta, STAT = istatus)
        ENDIF

! OLD CODE (COMMENTED OUT): duplicate deallocation of the same matrices.
!        DEALLOCATE(M_V, M0_V, M00_V, M10_V, M_K, Smat, STAT = istatus)
!        DEALLOCATE(mat_M, mat_M0, mat_M00, mat_M10, STAT = istatus)

        ENDDO ! Jtot

        PRINT*, 'GENERAL UNSYMMETRIZED CROSS SECTIONS at E = ', E*27.211399d0, ' eV'
        PRINT*, '   j1p   j2p        sigma'

        DO j1 = 0, pbasst(2)%pb_nbr, 2
        DO j2 = 0, pbasst(3)%pb_nbr, 2
                IF (Xsec_jpair(j1,j2) > 0d0) THEN
                        WRITE(6,'(2I6,1X,ES20.10)') j1, j2, Xsec_jpair(j1,j2)
                ENDIF
        ENDDO
        ENDDO

        IF (ALLOCATED(Xsec_jpair)) DEALLOCATE(Xsec_jpair, STAT=istatus)

        ENDDO ! E


        !print*,'legendre pol - l=2,0.5', plgndr(2,0,0.5d0), plgndr(2,-1,0.5d0), plgndr(2,1,0.5d0) 
        !print*,'legendre pol - l=2,-0.5', plgndr(2,0,-0.5d0), plgndr(2,-1,-0.5d0), plgndr(2,1,-0.5d0)
        !print*,'legendre pol - l=2,1', plgndr(2,0,1.0d0), plgndr(2,-1,1.0d0), plgndr(2,1,1.0d0)
        !print*,'legendre pol - l=2,-1', plgndr(2,0,-1.0d0), plgndr(2,-1,-1.0d0), plgndr(2,1,-1.0d0)
   

!
! Writing format
!---------------
! 500    FORMAT (8E15.7)
 590    FORMAT (10X,'INITIAL QUADRATURE ROOTS AND WEIGHTS FOR QUADRATURE:')
 595    FORMAT (8X, D12.5, D12.5)
! 600    FORMAT (6E15.7)                                                  
!
        CALL CPU_TIME(tm2)
        write(6,'(1A,F8.3,1A)') 'execution time:',(tm2-tm1)/6d1,'min'
        PRINT*," "
!
!*********************************************************************************************
!*********************************************************************************************
!
CLOSE(100)
CLOSE(proba_all_unit)
 END PROGRAM skvpAtomDiatom
!
!*********************************************************************************************
!=============================================================================================
!*********************************************************************************************
!
!
!
!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************
!
  SUBROUTINE CrossSection
!
!*********************************************************************************************
!
! This subroutine creates the total cross section using the probability result of the S-matrix
! as described in Ruthie (2013)
!
!=============================================================================================
!
        USE AtomDiatomskvp
        USE generateparameters
        USE omp_lib
!
        IMPLICIT NONE
!------------------------------
        INTEGER :: i, j
        REAL(8), ALLOCATABLE, DIMENSION(:)   :: Kj2, Coeff
        REAL(8), ALLOCATABLE, DIMENSION(:,:) :: Smat_prob, Dsec, Xsec
        REAL(8) :: E_rot

        
        !inquire(file='matrices.bin', exist=file_exists, size=file_size)
!
!        IF (file_size == 0) THEN
!                CALL noc_independent_calc
!        ELSE 
!        ALLOCATE ( BAM_1(1:n_x,1:n_x), BAM_3(1:n_x,1:n_x), BAM_4(1:n_x,1:n_x), &
!        BAM_r1(n_x, n_x), BAM_r2(n_x,n_x), BAM_r3(n_x,n_x), norm(n_x))
!
!        BAM_1=complex(0D0,0D0)
!        BAM_3=complex(0D0,0D0)
!        BAM_4=complex(0D0,0D0)
!        BAM_r1 = 0d0
!        BAM_r2 = 0d0
!        BAM_r3 = 0d0
!        norm = 0d0
!        ENDIF

!!! DO Jtot = Jmin, Jmax
!--------------------
! Creation of the basic auxiliary matrices
!-----------------------------------------
        CALL basic_aux_mat_calcul
!--------------------
! Creation of the potential matrices
!-----------------------------------------
        CALL potential_mat_calcul
!               
!               
! Creation of the main scattering matrices
!-----------------------------------------
        CALL make_scatt_mat

!=============================================================================================
! OLD CODE (COMMENTED OUT): after MOVE_ALLOC, M_V is already unallocated and
! the optimized implementation does not allocate M_K.
!        IF (ALLOCATED(M_V)) THEN
!                DEALLOCATE(M_V, STAT=istatus)
!        END IF
!
!        IF (ALLOCATED(M_K)) THEN
!                DEALLOCATE(M_K, STAT=istatus)
!        END IF
!=============================================================================================


!
! Sample test routines call
!--------------------------
        CALL PhaseShift
!
!-------------------------------
!       Allocation
!------------------------
        !ALLOCATE(Kj2(1:n_open), Coeff(1:n_open), STAT = istatus)
        !ALLOCATE(Smat_prob(1:n_open,1:n_open), Dsec(1:n_open,1:n_open), Xsec(1:n_open,1:n_open), STAT = istatus)
!----------------------------
!       Find Coefficient for Xsec
!-----------------------------     
    !DO j = 1, min(n_open,7)
        !DO i = 1, min(n_open,7)
            !Smat_prob(i,j) = (abs(Smat(i,j))**2.D0) !defining probability matrix
        !ENDDO
        !write(6,'(7E18.8)') (Smat_prob(i,j), i=1, min(n_open,7)) !!问一下史蒂夫
                
        !E_rot = Brot*quant_mat(1, open_idx(j))*(quant_mat(1, open_idx(j)) + 1d0) + &
               !Brot*quant_mat(3, open_idx(j))*(quant_mat(3, open_idx(j)) + 1d0)


       !Kj2(j) = mu_R * (E - E_rot) !why there is no 2
        ! make this for a vector E instead of single value
        !Coeff(j) = pi/((2*quant_mat(1,j)+1d0)*Kj2(j)) !coefficient for cross section
        !IF (Coeff(j) < 0) THEN
                !PRINT*, 'Warning: Negative Xsec Coefficient'
        !ENDIF
    !ENDDO

    !PRINT*, 'Kj2', Kj2
    !PRINT*, 'Coeff', Coeff
    !PRINT*, ' '

!       Angular momentum component of Xsec
!--------------------------------------------
    !Dsec = 0.0  
    !DO Jtotal = Jmin, Jmax !summing over Jtot
        !Dsec = Dsec + (2*Jtot+1)*Smat_prob
        !PRINT*, 'Dsec', Jtotal
        !DO i = 1, 3
            !DO j = 1, 3
                !PRINT*, Dsec(i,j)
            !ENDDO
        !ENDDO
    !ENDDO

!       Combining Results for Xsec result
!--------------------------------------------
    !DO j=1, min(n_open,7)
    !    Xsec(j,:) = Dsec(j,:) * Coeff(j) !Coeff number determined by which state we are coming from
    !ENDDO

    !PRINT*, 'Cross Section'
    !DO i = 1, min(n_open,7)
    !    write(6,'(7E18.8)') (Xsec(i,j), j=1, min(n_open,7))
    !ENDDO
    !PRINT*, ' '

!       Deallocating Matrices
! ---------------------------------
        !DEALLOCATE ( Kj2, Coeff, Smat_prob, Dsec, STAT = istatus )

 !!!ENDDO
!------------------------------------------------------------
! General unsymmetrized state-to-state cross section
!
! Accumulates:
!
! sigma(j1p,j2p <- j1,j2)
!   = pi/2*mu_R*E_coll * 1/[(2j1+1)(2j2+1)]
!     * sum_J (2J+1)
!     * sum_{k1,k2}
!     * sum_{k1p,k2p} |S^J(i -> f)|^2
!
! No exchange symmetrization, no primed sum, no delta factors.
!------------------------------------------------------------

        CALL AccumulateCrossSection

!=============================================================================================
!
 END SUBROUTINE CrossSection
!
!=============================================================================================
!
!=============================================================================================
!
 SUBROUTINE AccumulateCrossSection
!
!=============================================================================================
!
        USE AtomDiatomskvp
        USE generateparameters

        IMPLICIT NONE

        INTEGER :: i_open, f_open
        INTEGER :: iq, fq
        INTEGER :: j1i, k1i, j2i, k2i
        INTEGER :: j1f, k1f, j2f, k2f
        REAL(8) :: E_rot_i
        REAL(8) :: E_coll_i
        REAL(8) :: prefactor
        REAL(8) :: probability
        REAL(8) :: degeneracy_average
        REAL(8) :: J_weight

        IF (.NOT. ALLOCATED(Xsec_jpair)) THEN
                PRINT*, 'ERROR: Xsec_jpair is not allocated before CrossSection.'
                STOP
        ENDIF

        J_weight = DBLE(2*Jtot + 1)

        DO i_open = 1, n_open

                iq = open_idx(i_open)

                j1i = quant_mat(1,iq)
                k1i = quant_mat(2,iq)
                j2i = quant_mat(3,iq)
                k2i = quant_mat(4,iq)

                ! ONly Compute the cross section for first inital state 0,0,0,0
                IF (j1i /= 0 .OR. k1i /= 0 .OR. j2i /= 0 .OR. k2i /= 0) CYCLE

                E_rot_i = Brot*j1i*(j1i+1d0) + Brot*j2i*(j2i+1d0)

                E_coll_i = E - E_rot_i

                IF (E_coll_i <= 0d0) CYCLE

                degeneracy_average = 1d0 / DBLE((2*j1i + 1)*(2*j2i + 1))

                prefactor = pi / (2d0*mu_R*E_coll_i) * degeneracy_average

                DO f_open = 1, n_open

                        fq = open_idx(f_open)

                        j1f = quant_mat(1,fq)
                        k1f = quant_mat(2,fq)
                        j2f = quant_mat(3,fq)
                        k2f = quant_mat(4,fq)

                        probability = ABS(Smat(i_open,f_open))**2D0

                        Xsec_jpair(j1f,j2f) = Xsec_jpair(j1f,j2f) &
                                + prefactor * J_weight * probability

                ENDDO
        ENDDO

 END SUBROUTINE AccumulateCrossSection
!
!=============================================================================================
!
!
!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************
!
 SUBROUTINE make_scatt_mat
!
!*********************************************************************************************
!
! This subroutine creates the relevant scattering matrices M, M_0, M_0,0, M_1,0
!
!=============================================================================================
!
        USE AtomDiatomskvp
        USE generateparameters
        USE omp_lib
        USE Potential_Interface
!
        IMPLICIT NONE
!--------------------
        INTEGER i, j, k, i1, i2, j1, j2, N, quant_j, quant_j_prime, &
        quant_k, quant_k_prime
        INTEGER :: quant_j1, quant_k1, quant_j2, quant_k2
        INTEGER :: quant_j1_prime, quant_k1_prime, quant_j2_prime, quant_k2_prime
        REAL(8) :: E_rot
        REAL(8) :: channel_delta
        !REAL(8) modul, shift2delta
        !INTEGER, ALLOCATABLE, DIMENSION(:)      :: ipivx
!
        N = dim_x*ncf
!
!       Allocation
!       ----------
!=============================================================================================
! OLD CODE (COMMENTED OUT): this allocated mat_M and M_K in addition to M_V,
! so three full N x N real matrices could coexist.
!        ALLOCATE(mat_M(1:N,1:N), mat_M0(1:N,1:n_open), STAT=istatus)
!        ALLOCATE(mat_M00(1:n_open,1:n_open), &
!                 mat_M10(1:n_open,1:n_open), STAT=istatus)
!        ALLOCATE(M_K(N,N))
!        mat_M  = 0d0
!        mat_M0 = (0d0,0d0)
!        mat_M00 = (0d0,0d0)
!        mat_M10 = (0d0,0d0)
!        M_K = 0d0
!=============================================================================================
! NEW OPTIMIZED CODE: transfer ownership of the already-built potential
! matrix M_V to mat_M without copying its N x N data.
        IF (.NOT. ALLOCATED(M_V)) THEN
                PRINT*, 'Error: M_V is not allocated.'
                STOP
        ENDIF

        IF (ALLOCATED(mat_M)) DEALLOCATE(mat_M)

        CALL MOVE_ALLOC(M_V, mat_M)

        ALLOCATE(mat_M0(1:N,1:n_open), STAT=istatus)
        ALLOCATE(mat_M00(1:n_open,1:n_open), &
                mat_M10(1:n_open,1:n_open), STAT=istatus)

        mat_M0  = (0d0,0d0)
        mat_M00 = (0d0,0d0)
        mat_M10 = (0d0,0d0)
!
!       Formation of the various M matrices
!       -----------------------------------
        DO i1=1, dim_x   !1 refers to initial state
        DO i2=1, ncf     !2 refers to final state

        i = (i1-1)*ncf+i2
        quant_j1 = quant_mat(1,i2)
        quant_k1 = quant_mat(2,i2)
        quant_j2 = quant_mat(3,i2)
        quant_k2 = quant_mat(4,i2)

        DO j1=1, dim_x
        DO j2=1, ncf
        j = (j1-1)*ncf+j2
        quant_j1_prime = quant_mat(1,j2)
        quant_k1_prime = quant_mat(2,j2)
        quant_j2_prime = quant_mat(3,j2)
        quant_k2_prime = quant_mat(4,j2)
        

        channel_delta = delta(quant_j1,quant_j1_prime) * delta(quant_k1,quant_k1_prime) * &
                delta(quant_j2,quant_j2_prime) * delta(quant_k2,quant_k2_prime)

        E_rot = Brot*quant_j1*(quant_j1+1d0) + Brot*quant_j2*(quant_j2+1d0)

!=============================================================================================
! OLD CODE (COMMENTED OUT): form a separate full kinetic matrix M_K and then
! combine M_K, M_V, and the remaining terms into mat_M.
!        M_K(i,j) = (1d0/(2d0*mu_R)) * BAM_1(i1+1,j1+1) * channel_delta
!        mat_M(i,j) = M_K(i,j) + M_V(i,j) &
!                - E*BAM_3(i1+1,j1+1)*channel_delta &
!                + E_rot*BAM_3(i1+1,j1+1)*channel_delta &
!                + (1d0/(2d0*mu_R))*BAM_4(i1+1,j1+1) * &
!                  (quant_j1*(quant_j1+1d0)+quant_j2*(quant_j2+1d0))*channel_delta &
!                + (1d0/(2d0*mu_R))*BAM_4(i1+1,j1+1) * &
!                  Wdd(Jtot,quant_j1,quant_k1,quant_j2,quant_k2, &
!                      quant_j1_prime,quant_k1_prime,quant_j2_prime,quant_k2_prime)
!=============================================================================================
! NEW OPTIMIZED CODE: mat_M already contains the former M_V values. Add all
! kinetic, threshold, and rotational terms directly into that same storage.
        mat_M(i,j) = mat_M(i,j) &
                + (1d0/(2d0*mu_R)) * DBLE(BAM_1(i1+1,j1+1)) * channel_delta &
                - E * DBLE(BAM_3(i1+1,j1+1)) * channel_delta &
                + E_rot * DBLE(BAM_3(i1+1,j1+1)) * channel_delta &
                + (1d0/(2d0*mu_R)) * DBLE(BAM_4(i1+1,j1+1)) * &
                (quant_j1*(quant_j1+1d0) + quant_j2*(quant_j2+1d0)) * channel_delta &
                + (1d0/(2d0*mu_R)) * DBLE(BAM_4(i1+1,j1+1)) * &
                Wdd(Jtot, quant_j1, quant_k1, quant_j2, quant_k2, &
                quant_j1_prime, quant_k1_prime, &
                quant_j2_prime, quant_k2_prime) ! 1/R**2 term
        ENDDO
        ENDDO
!
        DO k=1, n_open

        quant_j1_prime = quant_mat(1,open_idx(k))
        quant_k1_prime = quant_mat(2,open_idx(k))
        quant_j2_prime = quant_mat(3,open_idx(k))
        quant_k2_prime = quant_mat(4,open_idx(k))

        channel_delta = delta(quant_j1,quant_j1_prime) * delta(quant_k1,quant_k1_prime) * &
                    delta(quant_j2,quant_j2_prime) * delta(quant_k2,quant_k2_prime)

        E_rot = Brot * quant_j1 * (quant_j1 + 1d0) + Brot * quant_j2 * (quant_j2 + 1d0)

        mat_M0(i,k) = (-1d0/(2d0*mu_R))*BAM_x1(i1+1,k)*channel_delta + M0_V(i,k) &
                      - E*BAM_x3(i1+1,k)*channel_delta & 
                      + E_rot*channel_delta*BAM_x3(i1+1,k) & 
                      + (1d0/(2d0*mu_R)) * BAM_x4(i1+1,k) * &
                        (quant_j1*(quant_j1 + 1d0)+quant_j2*(quant_j2 + 1d0))*channel_delta &
                      + (1d0/(2d0*mu_R)) * BAM_x4(i1+1,k) * &
                        Wdd(Jtot, quant_j1, quant_k1, quant_j2, quant_k2, &
                              quant_j1_prime, quant_k1_prime, quant_j2_prime, quant_k2_prime) ! 1/R**2 term
        ENDDO

        ENDDO
        ENDDO
!
        DO i = 1, n_open !n
        DO j = 1, n_open !n_prime
        quant_j1 = quant_mat(1,open_idx(i))
        quant_k1 = quant_mat(2,open_idx(i))
        quant_j2 = quant_mat(3,open_idx(i))
        quant_k2 = quant_mat(4,open_idx(i))

        quant_j1_prime = quant_mat(1,open_idx(j))
        quant_k1_prime = quant_mat(2,open_idx(j))
        quant_j2_prime = quant_mat(3,open_idx(j))
        quant_k2_prime = quant_mat(4,open_idx(j))

        channel_delta = delta(quant_j1,quant_j1_prime) * delta(quant_k1,quant_k1_prime) * &
                    delta(quant_j2,quant_j2_prime) * delta(quant_k2,quant_k2_prime)

        E_rot = Brot * quant_j1 * (quant_j1 + 1d0) + &
                Brot * quant_j2 * (quant_j2 + 1d0)

        mat_M00(i,j)= -(1d0/(2d0*mu_R))*BAM_xx1(i,j)*channel_delta+M00_V(i,j) &
                      - E*BAM_xx3(i,j)*channel_delta + &
                     E_rot*channel_delta*BAM_xx3(i,j) &
                      + (1d0/(2d0*mu_R)) * BAM_xx4(i,j) * &
                        (quant_j1*(quant_j1 + 1d0)+quant_j2*(quant_j2 + 1d0))*channel_delta &
                      + (1d0/(2d0*mu_R))*BAM_xx4(i,j)* &
                        Wdd(Jtot, quant_j1, quant_k1, quant_j2, quant_k2, &
                              quant_j1_prime, quant_k1_prime, quant_j2_prime, quant_k2_prime)

        mat_M10(i,j)= -(1d0/(2d0*mu_R))*BAM_xx1b(i,j)*channel_delta+M10_V(i,j) &
                         - E*BAM_xx3b(i,j)*channel_delta &
                       + E_rot*channel_delta*BAM_xx3b(i,j) &
                       + (1d0/(2d0*mu_R)) * BAM_xx4b(i,j) * &
                         (quant_j1*(quant_j1 + 1d0)+quant_j2*(quant_j2 + 1d0))*channel_delta &
                       + (1d0/(2d0*mu_R))*BAM_xx4b(i,j)* &
                         Wdd(Jtot, quant_j1, quant_k1, quant_j2, quant_k2, &
                              quant_j1_prime, quant_k1_prime, quant_j2_prime, quant_k2_prime)
        ENDDO
        ENDDO
!
!        PRINT*,'M_V'
!        do i=1, min(N,6)
!        write(6,'(8E18.8)') (dble(M_V(i,j)), j=1, min(N,6))
!        enddo
!        PRINT*,'M0_V'
!        do i=1, min(N,6)
!        write(6,'(8E18.8)') (dble(M0_V(i,j)), j=1, min(n_open,6))
!        enddo
!        PRINT*,'mat_M_V2'
!        do i=1, min(N,6)
!        write(6,'(8E18.8)'),(dble(M_V2(i,j)), j=1, min(N,6))
!        enddo

!        PRINT*,'M00_V'
!        do i=1, min(N,n_open)
!        write(6,'(8E18.8)') (dble(M00_V(i,:)))
!        enddo
!        PRINT*,'M10_V'
!        do i=1, min(N,n_open)
!        write(6,'(8E18.8)') (dble(M10_V(i,:)))
!        enddo
!         PRINT*,'mat_M00_V2'
!        do i=1, min(N,n_open)
!        write(6,'(8E18.8)'),(dble(M00_V2(i,:)))
!        enddo
!        PRINT*,'mat_M10_V2'
!        do i=1, min(N,n_open)
!        write(6,'(8E18.8)'),(dble(M10_V2(i,:)))
!       enddo


!       PRINT*,'mat_M'
!       do i=1, min(N,5)
!       write(6,'(8E18.8)') (dble(mat_M(i,j)), j=1, min(N,5))
!       enddo
!       PRINT*,'mat_M0'
!       do i= N-10, N
!       write(6,'(8E18.8)') (dble(mat_M0(i,j)), j=1, min(n_open,4))
!       enddo
!        PRINT*,'mat_M0_V'
!        do i=1, min(N,4)
!        write(6,'(8E18.8)'),(dble(M0_V(i,j)), j=1, min(n_open,4))
!        enddo
!       PRINT*,'mat_M_K'
!       do i=1, min(N,4)
!       write(6,'(8E18.8)') (dble(M_K(i,j)), j=1, min(n_open,4))
!       enddo
!       PRINT*,'mat_M00'
!       do i=1, min(n_open,5)
!       write(6,'(8E18.8)') (dble(mat_M00(i,:)))
!       enddo
!       PRINT*,'mat_M10'
!       do i=1, min(n_open,5)
!       write(6,'(8E18.8)') (dble(mat_M10(i,:)))
!       enddo


        ! OPEN(unit=10, file="M_Matrix.txt", status='unknown', position='append', action='write')
        ! WRITE(10,*) 'Matrix M:'
        ! DO i = 1, 5
        ! DO j = 1, 5
        !     WRITE(10, '(ES14.5)', advance='no') mat_M(i,j)
        ! END DO
        ! WRITE(10,*)
        ! END DO
        !  WRITE(10,*)

        !  WRITE(10,*) 'Matrix M_V:'
        !   DO i = 1, 15
        ! DO j = 1, 15
        !     WRITE(10,'(F15.3)', advance='no') M_V(i,j)
        ! END DO
        ! WRITE(10,*)
        ! END DO
        !  WRITE(10,*)

        !  WRITE(10,*) 'Matrix M_K:'
        !   DO i = 1, 15
        ! DO j = 1, 15
        !     WRITE(10,'(F15.3)', advance='no') M_K(i,j)
        ! END DO
        ! WRITE(10,*)
        ! END DO
        !  WRITE(10,*)
        ! CLOSE(10)        


!
! 600    FORMAT (6E15.7)
!
!PRINT*, '==== matrix NaN/debug check ===='
!PRINT*, 'max M_V     = ', MAXVAL(ABS(M_V))
!PRINT*, 'max M0_V    = ', MAXVAL(ABS(M0_V))
!PRINT*, 'max M00_V   = ', MAXVAL(ABS(M00_V))
!PRINT*, 'max M10_V   = ', MAXVAL(ABS(M10_V))
!PRINT*, 'max mat_M   = ', MAXVAL(ABS(mat_M))
!PRINT*, 'max mat_M0  = ', MAXVAL(ABS(mat_M0))
!PRINT*, 'max mat_M00 = ', MAXVAL(ABS(mat_M00))
!PRINT*, 'max mat_M10 = ', MAXVAL(ABS(mat_M10))

DO i = 1, n_open
        PRINT*, 'open channel ', i, &
                ' j1 k1 j2 k2 = ', quant_mat(1,open_idx(i)), &
                quant_mat(2,open_idx(i)), &
                quant_mat(3,open_idx(i)), &
                quant_mat(4,open_idx(i)), &
                ' kvec = ', kvec(i)
ENDDO

PRINT*, '==============================='
!=============================================================================================
!
	 END SUBROUTINE make_scatt_mat
	!
	!=============================================================================================
	!
	SUBROUTINE make_scatt_mat_sym_direct
	!
	!*********************************************************************************************
	!
	! Adds the kinetic, threshold, and rotational terms directly to the
	! exchange-symmetric matrices already initialized by
	! potential_mat_calcul_sym_direct.
	!
	!=============================================================================================
	!
	        USE AtomDiatomskvp
	        USE generateparameters
	        USE omp_lib
	        USE Potential_Interface
	!
	        IMPLICIT NONE
	!--------------------
	        INTEGER :: ir, jr
	        INTEGER :: s, t, os, ot
	        INTEGER :: a, b
	        INTEGER :: old_a, old_b
	        INTEGER :: old_open_a, old_open_b
	        INTEGER :: row_sym, col_sym
	        INTEGER :: Nsym_dim
	        INTEGER :: find_open_position
	        INTEGER :: quant_j1, quant_k1, quant_j2, quant_k2
	        INTEGER :: quant_j1_prime, quant_k1_prime
	        INTEGER :: quant_j2_prime, quant_k2_prime
	        REAL(8) :: coeff
	        REAL(8) :: E_rot
	        REAL(8) :: channel_delta

	        Nsym_dim = dim_x * nsym

	        IF (.NOT. ALLOCATED(mat_M_sym)) THEN
	                PRINT*, 'Error: mat_M_sym is not allocated before make_scatt_mat_sym_direct.'
	                STOP
	        ENDIF

	        IF (.NOT. ALLOCATED(mat_M0_sym) .OR. &
	            .NOT. ALLOCATED(mat_M00_sym) .OR. &
	            .NOT. ALLOCATED(mat_M10_sym)) THEN
	                PRINT*, 'Error: direct symmetric open matrices are not allocated.'
	                STOP
	        ENDIF

	!
	!       Closed-closed block.
	!       --------------------
	        DO ir = 1, dim_x
	        DO s = 1, nsym
	                row_sym = (ir-1)*nsym + s

	                DO jr = 1, dim_x
	                DO t = 1, nsym
	                        col_sym = (jr-1)*nsym + t

	                        DO a = 1, sym_ncomp(s)
	                                old_a = sym_old_idx(a,s)
	                                quant_j1 = quant_mat(1,old_a)
	                                quant_k1 = quant_mat(2,old_a)
	                                quant_j2 = quant_mat(3,old_a)
	                                quant_k2 = quant_mat(4,old_a)
	                                E_rot = Brot*quant_j1*(quant_j1+1d0) + &
	                                        Brot*quant_j2*(quant_j2+1d0)

	                                DO b = 1, sym_ncomp(t)
	                                        old_b = sym_old_idx(b,t)
	                                        quant_j1_prime = quant_mat(1,old_b)
	                                        quant_k1_prime = quant_mat(2,old_b)
	                                        quant_j2_prime = quant_mat(3,old_b)
	                                        quant_k2_prime = quant_mat(4,old_b)

	                                        coeff = sym_coeff(a,s) * sym_coeff(b,t)
	                                        channel_delta = &
	                                                delta(quant_j1,quant_j1_prime) * &
	                                                delta(quant_k1,quant_k1_prime) * &
	                                                delta(quant_j2,quant_j2_prime) * &
	                                                delta(quant_k2,quant_k2_prime)

	                                        mat_M_sym(row_sym,col_sym) = &
	                                                mat_M_sym(row_sym,col_sym) + coeff * ( &
	                                                (1d0/(2d0*mu_R)) * DBLE(BAM_1(ir+1,jr+1)) * channel_delta &
	                                                - E * DBLE(BAM_3(ir+1,jr+1)) * channel_delta &
	                                                + E_rot * DBLE(BAM_3(ir+1,jr+1)) * channel_delta &
	                                                + (1d0/(2d0*mu_R)) * DBLE(BAM_4(ir+1,jr+1)) * &
	                                                (quant_j1*(quant_j1+1d0) + &
	                                                 quant_j2*(quant_j2+1d0)) * channel_delta &
	                                                + (1d0/(2d0*mu_R)) * DBLE(BAM_4(ir+1,jr+1)) * &
	                                                Wdd(Jtot, quant_j1, quant_k1, quant_j2, quant_k2, &
	                                                quant_j1_prime, quant_k1_prime, &
	                                                quant_j2_prime, quant_k2_prime) )
	                                ENDDO
	                        ENDDO
	                ENDDO
	                ENDDO
	        ENDDO
	        ENDDO

	!
	!       Closed-open block.
	!       ------------------
	        DO ir = 1, dim_x
	        DO s = 1, nsym
	                row_sym = (ir-1)*nsym + s

	                DO os = 1, n_sym_open
	                        t = sym_open_idx(os)

	                        DO a = 1, sym_ncomp(s)
	                                old_a = sym_old_idx(a,s)
	                                quant_j1 = quant_mat(1,old_a)
	                                quant_k1 = quant_mat(2,old_a)
	                                quant_j2 = quant_mat(3,old_a)
	                                quant_k2 = quant_mat(4,old_a)
	                                E_rot = Brot*quant_j1*(quant_j1+1d0) + &
	                                        Brot*quant_j2*(quant_j2+1d0)

	                                DO b = 1, sym_ncomp(t)
	                                        old_b = sym_old_idx(b,t)
	                                        old_open_b = find_open_position(old_b)
	                                        IF (old_open_b <= 0) CYCLE

	                                        quant_j1_prime = quant_mat(1,old_b)
	                                        quant_k1_prime = quant_mat(2,old_b)
	                                        quant_j2_prime = quant_mat(3,old_b)
	                                        quant_k2_prime = quant_mat(4,old_b)

	                                        coeff = sym_coeff(a,s) * sym_coeff(b,t)
	                                        channel_delta = &
	                                                delta(quant_j1,quant_j1_prime) * &
	                                                delta(quant_k1,quant_k1_prime) * &
	                                                delta(quant_j2,quant_j2_prime) * &
	                                                delta(quant_k2,quant_k2_prime)

	                                        mat_M0_sym(row_sym,os) = mat_M0_sym(row_sym,os) + coeff * ( &
	                                                (-1d0/(2d0*mu_R))*BAM_x1(ir+1,old_open_b)*channel_delta &
	                                                - E*BAM_x3(ir+1,old_open_b)*channel_delta &
	                                                + E_rot*channel_delta*BAM_x3(ir+1,old_open_b) &
	                                                + (1d0/(2d0*mu_R)) * BAM_x4(ir+1,old_open_b) * &
	                                                (quant_j1*(quant_j1+1d0) + &
	                                                 quant_j2*(quant_j2+1d0))*channel_delta &
	                                                + (1d0/(2d0*mu_R)) * BAM_x4(ir+1,old_open_b) * &
	                                                Wdd(Jtot, quant_j1, quant_k1, quant_j2, quant_k2, &
	                                                quant_j1_prime, quant_k1_prime, &
	                                                quant_j2_prime, quant_k2_prime) )
	                                ENDDO
	                        ENDDO
	                ENDDO
	        ENDDO
	        ENDDO

	!
	!       Open-open blocks.
	!       -----------------
	        DO os = 1, n_sym_open
	                s = sym_open_idx(os)

	                DO ot = 1, n_sym_open
	                        t = sym_open_idx(ot)

	                        DO a = 1, sym_ncomp(s)
	                                old_a = sym_old_idx(a,s)
	                                old_open_a = find_open_position(old_a)
	                                IF (old_open_a <= 0) CYCLE

	                                quant_j1 = quant_mat(1,old_a)
	                                quant_k1 = quant_mat(2,old_a)
	                                quant_j2 = quant_mat(3,old_a)
	                                quant_k2 = quant_mat(4,old_a)
	                                E_rot = Brot*quant_j1*(quant_j1+1d0) + &
	                                        Brot*quant_j2*(quant_j2+1d0)

	                                DO b = 1, sym_ncomp(t)
	                                        old_b = sym_old_idx(b,t)
	                                        old_open_b = find_open_position(old_b)
	                                        IF (old_open_b <= 0) CYCLE

	                                        quant_j1_prime = quant_mat(1,old_b)
	                                        quant_k1_prime = quant_mat(2,old_b)
	                                        quant_j2_prime = quant_mat(3,old_b)
	                                        quant_k2_prime = quant_mat(4,old_b)

	                                        coeff = sym_coeff(a,s) * sym_coeff(b,t)
	                                        channel_delta = &
	                                                delta(quant_j1,quant_j1_prime) * &
	                                                delta(quant_k1,quant_k1_prime) * &
	                                                delta(quant_j2,quant_j2_prime) * &
	                                                delta(quant_k2,quant_k2_prime)

	                                        mat_M00_sym(os,ot) = mat_M00_sym(os,ot) + coeff * ( &
	                                                -(1d0/(2d0*mu_R))*BAM_xx1(old_open_a,old_open_b)*channel_delta &
	                                                - E*BAM_xx3(old_open_a,old_open_b)*channel_delta &
	                                                + E_rot*channel_delta*BAM_xx3(old_open_a,old_open_b) &
	                                                + (1d0/(2d0*mu_R)) * BAM_xx4(old_open_a,old_open_b) * &
	                                                (quant_j1*(quant_j1+1d0) + &
	                                                 quant_j2*(quant_j2+1d0))*channel_delta &
	                                                + (1d0/(2d0*mu_R))*BAM_xx4(old_open_a,old_open_b) * &
	                                                Wdd(Jtot, quant_j1, quant_k1, quant_j2, quant_k2, &
	                                                quant_j1_prime, quant_k1_prime, &
	                                                quant_j2_prime, quant_k2_prime) )

	                                        mat_M10_sym(os,ot) = mat_M10_sym(os,ot) + coeff * ( &
	                                                -(1d0/(2d0*mu_R))*BAM_xx1b(old_open_a,old_open_b)*channel_delta &
	                                                - E*BAM_xx3b(old_open_a,old_open_b)*channel_delta &
	                                                + E_rot*channel_delta*BAM_xx3b(old_open_a,old_open_b) &
	                                                + (1d0/(2d0*mu_R)) * BAM_xx4b(old_open_a,old_open_b) * &
	                                                (quant_j1*(quant_j1+1d0) + &
	                                                 quant_j2*(quant_j2+1d0))*channel_delta &
	                                                + (1d0/(2d0*mu_R))*BAM_xx4b(old_open_a,old_open_b) * &
	                                                Wdd(Jtot, quant_j1, quant_k1, quant_j2, quant_k2, &
	                                                quant_j1_prime, quant_k1_prime, &
	                                                quant_j2_prime, quant_k2_prime) )
	                                ENDDO
	                        ENDDO
	                ENDDO
	        ENDDO

	        PRINT*, 'Direct symmetric scattering matrices built.'
	        PRINT*, 'Direct N = ', Nsym_dim, ' direct n_open = ', n_sym_open
	        PRINT*, '==============================='

	END SUBROUTINE make_scatt_mat_sym_direct
	!
	!=============================================================================================
	!
	!
	!*********************************************************************************************
	!*********************************************************************************************
	!*********************************************************************************************!
 
 SUBROUTINE PhaseShift
!
!*********************************************************************************************
!
! This subroutine creates the relevant matrices and values and computes the S-matrix
!
!=============================================================================================
!
        USE AtomDiatomskvp
        USE generateparameters
        USE omp_lib
!
        IMPLICIT NONE
!--------------------
        INTEGER :: i, j, infox, lworkx, N
        INTEGER, ALLOCATABLE :: ipiv_M(:), ipiv_B(:)

        REAL(8), ALLOCATABLE :: sum_array(:), workx(:)
        REAL(8), ALLOCATABLE :: rhs_re(:,:), rhs_im(:,:), modul(:,:)
        REAL(8) :: work_query(1)

        COMPLEX(8), ALLOCATABLE :: solution_M0(:,:)
        COMPLEX(8), ALLOCATABLE :: Bsub(:,:), Csub(:,:), Smat_st(:,:)
        COMPLEX(8), ALLOCATABLE :: mat_B_st_inv(:,:), dum6vx(:,:)

!=============================================================================================
! OLD CODE (COMMENTED OUT): declarations used to store a complete inverse of
! mat_M plus several N x n_open temporary copies.
!        INTEGER, ALLOCATABLE :: ipivx(:)
!        REAL(8), ALLOCATABLE :: mat_M_inv(:,:), workx(:), sum_array(:), modul(:,:)
!        COMPLEX(8), ALLOCATABLE :: dum1vx(:,:), dum2vx(:,:), dum3vx(:,:)
!        COMPLEX(8), ALLOCATABLE :: dum4vx(:,:), dum5vx(:,:), dum6vx(:,:)
!=============================================================================================

        N = dim_x*ncf

!=============================================================================================
! OLD CODE (COMMENTED OUT): allocate and explicitly form inverse(mat_M).
!        lworkx = 50*(pbasst(1)%pb_nbr**2)
!        ALLOCATE(mat_M_inv(N,N), ipivx(N))
!        ALLOCATE(dum1vx(N,n_open), dum2vx(N,n_open))
!        ALLOCATE(dum3vx(n_open,N), dum4vx(n_open,N))
!        mat_M_inv = mat_M
!        CALL DSYTRF('L',N,mat_M_inv,N,ipivx,workx,N,infox)
!        CALL DSYTRI('L',N,mat_M_inv,N,ipivx,workx,infox)
!        DO i = 1, N
!                DO j = i+1, N
!                        mat_M_inv(i,j) = mat_M_inv(j,i)
!                ENDDO
!        ENDDO
!        dum1vx = MATMUL(mat_M_inv,mat_M0)
!        dum2vx = CONJG(mat_M0)
!        dum3vx = TRANSPOSE(mat_M0)
!        dum4vx = TRANSPOSE(dum2vx)
!        Bsub = MATMUL(dum3vx,dum1vx)
!        Csub = MATMUL(dum4vx,dum1vx)
!=============================================================================================
! NEW OPTIMIZED CODE: factor mat_M once and solve mat_M*X=mat_M0 directly.
! The real and imaginary parts are solved separately because mat_M is real.
        ALLOCATE(ipiv_M(N), rhs_re(N,n_open), rhs_im(N,n_open), &
                 solution_M0(N,n_open), STAT=istatus)
        IF (istatus /= 0) STOP 'Allocation failed for mat_M solve arrays'

        ALLOCATE(mat_B(n_open,n_open), mat_C(n_open,n_open), &
                 Bsub(n_open,n_open), Csub(n_open,n_open), &
                 mat_B_st_inv(n_open,n_open), dum6vx(n_open,n_open), &
                 Smat(n_open,n_open), Smat_st(n_open,n_open), &
                 modul(n_open,n_open), sum_array(n_open), STAT=istatus)
        IF (istatus /= 0) STOP 'Allocation failed for open-channel matrices'

        rhs_re = REAL(mat_M0,KIND=8)
        rhs_im = AIMAG(mat_M0)

        CALL DSYTRF('L',N,mat_M,N,ipiv_M,work_query,-1,infox)
        IF (infox /= 0) THEN
                PRINT*, 'DSYTRF workspace query failed: ', infox
                STOP
        ENDIF

        lworkx = MAX(1,INT(work_query(1)))
        ALLOCATE(workx(lworkx),STAT=istatus)
        IF (istatus /= 0) STOP 'Allocation failed for DSYTRF workspace'

        CALL DSYTRF('L',N,mat_M,N,ipiv_M,workx,lworkx,infox)
        IF (infox /= 0) THEN
                PRINT*, 'DSYTRF failed: ', infox
                STOP
        ENDIF

        CALL DSYTRS('L',N,n_open,mat_M,N,ipiv_M,rhs_re,N,infox)
        IF (infox /= 0) STOP 'DSYTRS failed for real part'

        CALL DSYTRS('L',N,n_open,mat_M,N,ipiv_M,rhs_im,N,infox)
        IF (infox /= 0) STOP 'DSYTRS failed for imaginary part'

        solution_M0 = CMPLX(rhs_re,rhs_im,KIND=8)

        Bsub = MATMUL(TRANSPOSE(mat_M0),solution_M0)
        Csub = MATMUL(TRANSPOSE(CONJG(mat_M0)),solution_M0)
        mat_B = mat_M00 - Bsub
        mat_C = mat_M10 - Csub

!=============================================================================================
! OLD CODE (COMMENTED OUT): explicitly invert CONJG(mat_B) with ZGETRI.
!        mat_B_st_inv = CONJG(mat_B)
!        CALL ZGETRF(n_open,n_open,mat_B_st_inv,n_open,ipivx,infox)
!        CALL ZGETRI(n_open,mat_B_st_inv,n_open,ipivx,workx,lworkx,infox)
!        dum6vx = MATMUL(mat_B_st_inv,mat_C)
!        Smat = (0d0,1d0) * &
!               (mat_B-MATMUL(TRANSPOSE(mat_C),dum6vx))
!=============================================================================================
! NEW OPTIMIZED CODE: solve CONJG(mat_B)*X=mat_C directly with ZGETRS.
        ALLOCATE(ipiv_B(n_open),STAT=istatus)
        IF (istatus /= 0) STOP 'Allocation failed for mat_B pivots'

        mat_B_st_inv = CONJG(mat_B)
        dum6vx = mat_C

        CALL ZGETRF(n_open,n_open,mat_B_st_inv,n_open,ipiv_B,infox)
        IF (infox /= 0) STOP 'ZGETRF failed'

        CALL ZGETRS('N',n_open,n_open,mat_B_st_inv,n_open,ipiv_B, &
                    dum6vx,n_open,infox)
        IF (infox /= 0) STOP 'ZGETRS failed'

        Smat = (0d0,1d0) * &
               (mat_B-MATMUL(TRANSPOSE(mat_C),dum6vx))

!=============================================================================================
! END NEW OPTIMIZED MATRIX SOLVE
!=============================================================================================

        Smat_st = TRANSPOSE(CONJG(Smat))

      


!       print Sdegger to check 
       ! DO i = 1, n_open
       !         dum5vx(i,i) = dum5vx(i,i) - (1d0,0d0)
       ! ENDDO

       ! PRINT*, 'max |Sdagger*S - I| = ', MAXVAL(ABS(dum5vx))


        modul=MATMUL(Smat_st,Smat)
        PRINT*,' '
!        PRINT*,'Module'
!        do i=1, min(n_open,7)
!        write(6,'(7E18.8)') ((modul(i,j)), j=1, min(n_open,7))
!        enddo
        PRINT*,' '
        PRINT*, 'Probabilities of first', min(7,n_open), 'channel'

        do i=1, min(n_open,7)
        write(6,'(7E18.8)') ((abs(Smat(i,j))**2D0), j=1, min(n_open,7))
        enddo

        PRINT*, ' '
!*!        shift2delta=dacos(-1d0)+datan2(dimag(Smat)/modul,dble(Smat)/modul)
!*!        PRINT*,'Phase Shift=',shift2delta/2d0,(shift2delta-dacos(-1d0))/2d0
!        PRINT*, 'Probabilities: j1,k1,j2,k2 -> j1p,k1p,j2p,k2p'
!
 !       DO i = 1, n_open
 !       DO j = 1, n_open
  !              write(6,'(E18.8)', advance='no') ABS(Smat(i,j))**2D0
   !             PRINT*, quant_mat(1,open_idx(i)), quant_mat(2,open_idx(i)), &
   !                     quant_mat(3,open_idx(i)), quant_mat(4,open_idx(i)), &
   !                     ' -> ', &
   !                     quant_mat(1,open_idx(j)), quant_mat(2,open_idx(j)), &
   !                     quant_mat(3,open_idx(j)), quant_mat(4,open_idx(j))
    !            write(*,*)
    !    ENDDO
    !    ENDDO

        sum_array = 0d0
        DO i = 1, n_open
        DO j = 1, n_open
             sum_array(i) = sum_array(i) + (abs(Smat(i,j))**2)
        ENDDO
        ENDDO
        PRINT*, "Checking probabilities sum to one: ", sum_array(:)
        PRINT*, ' '
        ! Adding data to output file
        ! IF (n_open >= 2) THEN
        !         open(unit=10, file="2_0.txt", status='unknown', position='append', action='write')
        !         write(10, '(ES24.16, 2X, ES24.16)') E, (abs(Smat(2,1))**2D0)
        !         close(10)
        ! ELSE 
        !         open(unit=10, file="2_0.txt", status='unknown', position='append', action='write')
        !         write(10, '(ES24.16, 2X, ES24.16)') E, 0d0
        !         close(10)
        ! ENDIF
        ! IF (n_open >= 3) THEN
        !         open(unit=10, file="4_0.txt", status='unknown', position='append', action='write')
        !         write(10, '(ES24.16, 2X, ES24.16)') E, (abs(Smat(3,1))**2D0)
        !         close(10)
        ! ELSE
        !         open(unit=10, file="4_0.txt", status='unknown', position='append', action='write')
        !         write(10, '(ES24.16, 2X, ES24.16)') E, 0d0
        !         close(10)
        ! ENDIF


        IF (n_open >= 1) THEN
                open(unit=10, file="0_0_all.txt", status='unknown', position='append', action='write')
                DO i = 1, n_open
                        write(10, '(4I5, 2X, ES24.16, 2X, ES24.16)') &
                                quant_mat(1,open_idx(i)), quant_mat(2,open_idx(i)), &
                                quant_mat(3,open_idx(i)), quant_mat(4,open_idx(i)), &
                                E, ABS(Smat(1,i))**2D0
                ENDDO
                close(10)
        ELSE
                open(unit=10, file="0_0_all.txt", status='unknown', position='append', action='write')
                write(10, '(4I5, 2X, ES24.16, 2X, ES24.16)') 0, 0, 0, 0, E, 0d0
                close(10)
        ENDIF


        IF (n_open >= 4) THEN
                open(unit=10, file="open4_all.txt", status='unknown', position='append', action='write')
                DO i = 1, n_open
                        write(10, '(4I5, 2X, ES24.16, 2X, ES24.16)') &
                                quant_mat(1,open_idx(i)), quant_mat(2,open_idx(i)), &
                                quant_mat(3,open_idx(i)), quant_mat(4,open_idx(i)), &
                                E, ABS(Smat(4,i))**2D0
                ENDDO
                close(10)
        ELSE
                open(unit=10, file="open4_all.txt", status='unknown', position='append', action='write')
                write(10, '(4I5, 2X, ES24.16, 2X, ES24.16)') 0, 0, 0, 0, E, 0d0
                close(10)
        ENDIF
!
!       Deallocating Matrices
!       ---------------------
!=============================================================================================
! OLD CODE (COMMENTED OUT): cleanup for explicit-inverse temporary arrays.
!        DEALLOCATE(workx,ipivx,STAT=istatus)
!        DEALLOCATE(mat_M_inv,STAT=istatus)
!        DEALLOCATE(dum1vx,dum2vx,Smat_st,STAT=istatus)
!        DEALLOCATE(dum3vx,dum4vx,dum5vx,dum6vx,STAT=istatus)
!        DEALLOCATE(mat_B,mat_C,mat_B_st_inv,modul,STAT=istatus)
!        DEALLOCATE(Bsub,Csub,sum_array,STAT=istatus)
!=============================================================================================
! NEW OPTIMIZED CODE: cleanup for direct-solve arrays.
        DEALLOCATE(workx,ipiv_M,ipiv_B,STAT=istatus)
        DEALLOCATE(rhs_re,rhs_im,solution_M0,STAT=istatus)
        DEALLOCATE(Smat_st,dum6vx,STAT=istatus)
        DEALLOCATE(mat_B,mat_C,mat_B_st_inv,modul,STAT=istatus)
        DEALLOCATE(Bsub,Csub,sum_array,STAT=istatus)
!   
!
! 600    FORMAT (6E15.7)
!
!=============================================================================================
!
	 END SUBROUTINE PhaseShift
	!
	!=============================================================================================
	!*********************************************************************************************
	!*********************************************************************************************

	SUBROUTINE cleanup_phase2_matrix_storage(include_sym)

	        USE AtomDiatomskvp

	        IMPLICIT NONE

	        LOGICAL, INTENT(IN) :: include_sym

	        IF (ALLOCATED(M_V))   DEALLOCATE(M_V,   STAT=istatus)
	        IF (ALLOCATED(M0_V))  DEALLOCATE(M0_V,  STAT=istatus)
	        IF (ALLOCATED(M00_V)) DEALLOCATE(M00_V, STAT=istatus)
	        IF (ALLOCATED(M10_V)) DEALLOCATE(M10_V, STAT=istatus)

	        IF (ALLOCATED(mat_M))   DEALLOCATE(mat_M,   STAT=istatus)
	        IF (ALLOCATED(mat_M0))  DEALLOCATE(mat_M0,  STAT=istatus)
	        IF (ALLOCATED(mat_M00)) DEALLOCATE(mat_M00, STAT=istatus)
	        IF (ALLOCATED(mat_M10)) DEALLOCATE(mat_M10, STAT=istatus)

	        IF (include_sym) THEN
	                IF (ALLOCATED(mat_M_sym))   DEALLOCATE(mat_M_sym,   STAT=istatus)
	                IF (ALLOCATED(mat_M0_sym))  DEALLOCATE(mat_M0_sym,  STAT=istatus)
	                IF (ALLOCATED(mat_M00_sym)) DEALLOCATE(mat_M00_sym, STAT=istatus)
	                IF (ALLOCATED(mat_M10_sym)) DEALLOCATE(mat_M10_sym, STAT=istatus)
	        ENDIF

	        IF (ALLOCATED(BAM_r))     DEALLOCATE(BAM_r,     STAT=istatus)
	        IF (ALLOCATED(BAM_r0))    DEALLOCATE(BAM_r0,    STAT=istatus)
	        IF (ALLOCATED(BAM_r00))   DEALLOCATE(BAM_r00,   STAT=istatus)
	        IF (ALLOCATED(BAM_r10))   DEALLOCATE(BAM_r10,   STAT=istatus)
	        IF (ALLOCATED(BAM_theta)) DEALLOCATE(BAM_theta, STAT=istatus)

	END SUBROUTINE cleanup_phase2_matrix_storage

	!=============================================================================================
	!*********************************************************************************************
	!*********************************************************************************************

	SUBROUTINE build_exchange_symmetric_basis

	        USE AtomDiatomskvp

	        IMPLICIT NONE

	        INTEGER :: old_idx, partner_idx, sym_idx
	        INTEGER :: find_exchange_partner
	        INTEGER, ALLOCATABLE :: seen(:)
	        REAL(8) :: inv_sqrt2

	        inv_sqrt2 = 1d0 / DSQRT(2d0)

	        IF (ALLOCATED(sym_quant_mat)) DEALLOCATE(sym_quant_mat, STAT=istatus)
	        IF (ALLOCATED(sym_old_idx)) DEALLOCATE(sym_old_idx, STAT=istatus)
	        IF (ALLOCATED(sym_coeff)) DEALLOCATE(sym_coeff, STAT=istatus)
	        IF (ALLOCATED(sym_ncomp)) DEALLOCATE(sym_ncomp, STAT=istatus)
	        IF (ALLOCATED(old_to_sym)) DEALLOCATE(old_to_sym, STAT=istatus)

	        ALLOCATE(seen(1:ncf), STAT=istatus)
	        IF (istatus /= 0) STOP 'Allocation failed for exchange seen array'
	        seen = 0
	        nsym = 0

	        DO old_idx = 1, ncf
	                IF (seen(old_idx) /= 0) CYCLE
	                partner_idx = find_exchange_partner(old_idx)
	                IF (partner_idx <= 0) THEN
	                        PRINT*, 'Exchange partner not found for old channel ', old_idx
	                        STOP
	                ENDIF
	                nsym = nsym + 1
	                seen(old_idx) = nsym
	                seen(partner_idx) = nsym
	        ENDDO

	        ALLOCATE(sym_quant_mat(4,nsym), sym_old_idx(2,nsym), &
	                 sym_coeff(2,nsym), sym_ncomp(nsym), old_to_sym(ncf), &
	                 STAT=istatus)
	        IF (istatus /= 0) STOP 'Allocation failed for exchange basis arrays'

	        sym_quant_mat = 0
	        sym_old_idx = 0
	        sym_coeff = 0d0
	        sym_ncomp = 0
	        old_to_sym = 0
	        seen = 0
	        sym_idx = 0

	        DO old_idx = 1, ncf
	                IF (seen(old_idx) /= 0) CYCLE
	                partner_idx = find_exchange_partner(old_idx)
	                sym_idx = sym_idx + 1

	                sym_quant_mat(:,sym_idx) = quant_mat(:,old_idx)
	                sym_old_idx(1,sym_idx) = old_idx
	                old_to_sym(old_idx) = sym_idx

	                IF (partner_idx == old_idx) THEN
	                        sym_ncomp(sym_idx) = 1
	                        sym_coeff(1,sym_idx) = 1d0
	                        seen(old_idx) = sym_idx
	                ELSE
	                        sym_ncomp(sym_idx) = 2
	                        sym_old_idx(2,sym_idx) = partner_idx
	                        sym_coeff(1,sym_idx) = inv_sqrt2
	                        sym_coeff(2,sym_idx) = inv_sqrt2
	                        old_to_sym(partner_idx) = sym_idx
	                        seen(old_idx) = sym_idx
	                        seen(partner_idx) = sym_idx
	                ENDIF
	        ENDDO

	        PRINT*, 'Exchange-symmetric basis: old ncf = ', ncf, ' nsym = ', nsym
	        CALL write_exchange_basis_map

	        DEALLOCATE(seen, STAT=istatus)

	END SUBROUTINE build_exchange_symmetric_basis

	INTEGER FUNCTION find_exchange_partner(old_idx)

	        USE AtomDiatomskvp

	        IMPLICIT NONE

	        INTEGER, INTENT(IN) :: old_idx
	        INTEGER :: trial
	        INTEGER :: j1, k1, j2, k2

	        j1 = quant_mat(1,old_idx)
	        k1 = quant_mat(2,old_idx)
	        j2 = quant_mat(3,old_idx)
	        k2 = quant_mat(4,old_idx)

	        find_exchange_partner = 0

	        DO trial = 1, ncf
	                IF (quant_mat(1,trial) == j2 .AND. &
	                    quant_mat(2,trial) == k2 .AND. &
	                    quant_mat(3,trial) == j1 .AND. &
	                    quant_mat(4,trial) == k1) THEN
	                        find_exchange_partner = trial
	                        RETURN
	                ENDIF
	        ENDDO

	END FUNCTION find_exchange_partner

	SUBROUTINE write_exchange_basis_map

	        USE AtomDiatomskvp

	        IMPLICIT NONE

	        INTEGER :: unit_map
	        INTEGER :: s, c, old_idx

	        OPEN(newunit=unit_map, file='phase2_sym_channel_map.dat', &
	             status='replace', action='write')

	        WRITE(unit_map,'(A)') '# sym_index ncomp representative_j1 k1 j2 k2 component old_index j1 k1 j2 k2 coeff'

	        DO s = 1, nsym
	                DO c = 1, sym_ncomp(s)
	                        old_idx = sym_old_idx(c,s)
	                        WRITE(unit_map,'(2I8,4I6,I8,I8,4I6,ES20.10)') &
	                                s, sym_ncomp(s), sym_quant_mat(1,s), &
	                                sym_quant_mat(2,s), sym_quant_mat(3,s), &
	                                sym_quant_mat(4,s), c, old_idx, &
	                                quant_mat(1,old_idx), quant_mat(2,old_idx), &
	                                quant_mat(3,old_idx), quant_mat(4,old_idx), &
	                                sym_coeff(c,s)
	                ENDDO
	        ENDDO

	        CLOSE(unit_map)

	END SUBROUTINE write_exchange_basis_map

	SUBROUTINE build_exchange_symmetric_open_basis

	        USE AtomDiatomskvp

	        IMPLICIT NONE

	        INTEGER :: old_open, sym_idx, sym_open
	        INTEGER :: find_open_position
	        INTEGER, ALLOCATABLE :: seen(:)
	        INTEGER :: unit_map

	        IF (ALLOCATED(sym_open_idx)) DEALLOCATE(sym_open_idx, STAT=istatus)
	        IF (ALLOCATED(old_open_pos_to_sym_open)) THEN
	                DEALLOCATE(old_open_pos_to_sym_open, STAT=istatus)
	        ENDIF

	        ALLOCATE(seen(1:nsym), STAT=istatus)
	        IF (istatus /= 0) STOP 'Allocation failed for sym open seen array'
	        seen = 0
	        n_sym_open = 0

	        DO old_open = 1, n_open
	                sym_idx = old_to_sym(open_idx(old_open))
	                IF (seen(sym_idx) == 0) THEN
	                        n_sym_open = n_sym_open + 1
	                        seen(sym_idx) = n_sym_open
	                ENDIF
	        ENDDO

	        ALLOCATE(sym_open_idx(1:n_sym_open), &
	                 old_open_pos_to_sym_open(1:n_open), STAT=istatus)
	        IF (istatus /= 0) STOP 'Allocation failed for sym open arrays'

	        sym_open_idx = 0
	        old_open_pos_to_sym_open = 0
	        seen = 0
	        sym_open = 0

	        DO old_open = 1, n_open
	                sym_idx = old_to_sym(open_idx(old_open))
	                IF (seen(sym_idx) == 0) THEN
	                        sym_open = sym_open + 1
	                        seen(sym_idx) = sym_open
	                        sym_open_idx(sym_open) = sym_idx
	                ENDIF
	                old_open_pos_to_sym_open(old_open) = seen(sym_idx)
	        ENDDO

	        PRINT*, 'Exchange-symmetric open channels: old n_open = ', &
	                n_open, ' n_sym_open = ', n_sym_open

	        OPEN(newunit=unit_map, file='phase2_sym_open_map.dat', &
	             status='replace', action='write')
	        WRITE(unit_map,'(A)') '# sym_open sym_index old_component old_open old_index j1 k1 j2 k2 coeff'

	        DO sym_open = 1, n_sym_open
	                sym_idx = sym_open_idx(sym_open)
	                DO old_open = 1, sym_ncomp(sym_idx)
	                        WRITE(unit_map,'(5I8,4I6,ES20.10)') sym_open, sym_idx, &
	                                old_open, &
	                                find_open_position(sym_old_idx(old_open,sym_idx)), &
	                                sym_old_idx(old_open,sym_idx), &
	                                quant_mat(1,sym_old_idx(old_open,sym_idx)), &
	                                quant_mat(2,sym_old_idx(old_open,sym_idx)), &
	                                quant_mat(3,sym_old_idx(old_open,sym_idx)), &
	                                quant_mat(4,sym_old_idx(old_open,sym_idx)), &
	                                sym_coeff(old_open,sym_idx)
	                ENDDO
	        ENDDO
	        CLOSE(unit_map)

	        DEALLOCATE(seen, STAT=istatus)

	END SUBROUTINE build_exchange_symmetric_open_basis

	INTEGER FUNCTION find_open_position(old_channel_idx)

	        USE AtomDiatomskvp

	        IMPLICIT NONE

	        INTEGER, INTENT(IN) :: old_channel_idx
	        INTEGER :: pos

	        find_open_position = 0

	        DO pos = 1, n_open
	                IF (open_idx(pos) == old_channel_idx) THEN
	                        find_open_position = pos
	                        RETURN
	                ENDIF
	        ENDDO

	END FUNCTION find_open_position

	SUBROUTINE project_matrices_to_exchange_symmetric

	        USE AtomDiatomskvp

	        IMPLICIT NONE

	        INTEGER :: ir, jr, s, t, os, ot
	        INTEGER :: a, b
	        INTEGER :: find_open_position
	        INTEGER :: old_a, old_b
	        INTEGER :: old_open_a, old_open_b
	        INTEGER :: row_sym, col_sym, row_old, col_old
	        INTEGER :: Nsym_dim

	        Nsym_dim = dim_x * nsym

	        IF (ALLOCATED(mat_M_sym)) DEALLOCATE(mat_M_sym, STAT=istatus)
	        IF (ALLOCATED(mat_M0_sym)) DEALLOCATE(mat_M0_sym, STAT=istatus)
	        IF (ALLOCATED(mat_M00_sym)) DEALLOCATE(mat_M00_sym, STAT=istatus)
	        IF (ALLOCATED(mat_M10_sym)) DEALLOCATE(mat_M10_sym, STAT=istatus)

	        ALLOCATE(mat_M_sym(1:Nsym_dim,1:Nsym_dim), &
	                 mat_M0_sym(1:Nsym_dim,1:n_sym_open), &
	                 mat_M00_sym(1:n_sym_open,1:n_sym_open), &
	                 mat_M10_sym(1:n_sym_open,1:n_sym_open), STAT=istatus)
	        IF (istatus /= 0) STOP 'Allocation failed for sym projected matrices'

	        mat_M_sym = 0d0
	        mat_M0_sym = (0d0,0d0)
	        mat_M00_sym = (0d0,0d0)
	        mat_M10_sym = (0d0,0d0)

	        DO ir = 1, dim_x
	        DO s = 1, nsym
	                row_sym = (ir-1)*nsym + s
	                DO jr = 1, dim_x
	                DO t = 1, nsym
	                        col_sym = (jr-1)*nsym + t
	                        DO a = 1, sym_ncomp(s)
	                                old_a = sym_old_idx(a,s)
	                                row_old = (ir-1)*ncf + old_a
	                                DO b = 1, sym_ncomp(t)
	                                        old_b = sym_old_idx(b,t)
	                                        col_old = (jr-1)*ncf + old_b
	                                        mat_M_sym(row_sym,col_sym) = &
	                                                mat_M_sym(row_sym,col_sym) + &
	                                                sym_coeff(a,s) * sym_coeff(b,t) * &
	                                                mat_M(row_old,col_old)
	                                ENDDO
	                        ENDDO
	                ENDDO
	                ENDDO
	        ENDDO
	        ENDDO

	        DO ir = 1, dim_x
	        DO s = 1, nsym
	                row_sym = (ir-1)*nsym + s
	                DO os = 1, n_sym_open
	                        t = sym_open_idx(os)
	                        DO a = 1, sym_ncomp(s)
	                                old_a = sym_old_idx(a,s)
	                                row_old = (ir-1)*ncf + old_a
	                                DO b = 1, sym_ncomp(t)
	                                        old_b = sym_old_idx(b,t)
	                                        old_open_b = find_open_position(old_b)
	                                        IF (old_open_b <= 0) CYCLE
	                                        mat_M0_sym(row_sym,os) = &
	                                                mat_M0_sym(row_sym,os) + &
	                                                sym_coeff(a,s) * sym_coeff(b,t) * &
	                                                mat_M0(row_old,old_open_b)
	                                ENDDO
	                        ENDDO
	                ENDDO
	        ENDDO
	        ENDDO

	        DO os = 1, n_sym_open
	                s = sym_open_idx(os)
	                DO ot = 1, n_sym_open
	                        t = sym_open_idx(ot)
	                        DO a = 1, sym_ncomp(s)
	                                old_a = sym_old_idx(a,s)
	                                old_open_a = find_open_position(old_a)
	                                IF (old_open_a <= 0) CYCLE
	                                DO b = 1, sym_ncomp(t)
	                                        old_b = sym_old_idx(b,t)
	                                        old_open_b = find_open_position(old_b)
	                                        IF (old_open_b <= 0) CYCLE
	                                        mat_M00_sym(os,ot) = mat_M00_sym(os,ot) + &
	                                                sym_coeff(a,s) * sym_coeff(b,t) * &
	                                                mat_M00(old_open_a,old_open_b)
	                                        mat_M10_sym(os,ot) = mat_M10_sym(os,ot) + &
	                                                sym_coeff(a,s) * sym_coeff(b,t) * &
	                                                mat_M10(old_open_a,old_open_b)
	                                ENDDO
	                        ENDDO
	                ENDDO
	        ENDDO

	        PRINT*, 'Projected ordered matrices to exchange-symmetric basis.'
	        PRINT*, 'Projected N = ', Nsym_dim, ' projected n_open = ', n_sym_open

	END SUBROUTINE project_matrices_to_exchange_symmetric

	SUBROUTINE SolveSMatrixGeneric(Nmat, NopenIn, matM_in, matM0_in, &
	                               matM00_in, matM10_in, Smat_out, label)

	        USE AtomDiatomskvp, ONLY: istatus

	        IMPLICIT NONE

	        INTEGER, INTENT(IN) :: Nmat, NopenIn
	        REAL(8), INTENT(IN) :: matM_in(Nmat,Nmat)
	        COMPLEX(8), INTENT(IN) :: matM0_in(Nmat,NopenIn)
	        COMPLEX(8), INTENT(IN) :: matM00_in(NopenIn,NopenIn)
	        COMPLEX(8), INTENT(IN) :: matM10_in(NopenIn,NopenIn)
	        COMPLEX(8), INTENT(OUT) :: Smat_out(NopenIn,NopenIn)
	        CHARACTER(LEN=*), INTENT(IN) :: label

	        INTEGER :: infox, lworkx
	        INTEGER, ALLOCATABLE :: ipiv_M(:), ipiv_B(:)
	        REAL(8), ALLOCATABLE :: matM_work(:,:), workx(:)
	        REAL(8), ALLOCATABLE :: rhs_re(:,:), rhs_im(:,:)
	        REAL(8) :: work_query(1)
	        COMPLEX(8), ALLOCATABLE :: solution_M0(:,:)
	        COMPLEX(8), ALLOCATABLE :: mat_B_local(:,:), mat_C_local(:,:)
	        COMPLEX(8), ALLOCATABLE :: Bsub(:,:), Csub(:,:)
	        COMPLEX(8), ALLOCATABLE :: mat_B_factor(:,:), solved_C(:,:)

	        ALLOCATE(matM_work(Nmat,Nmat), ipiv_M(Nmat), &
	                 rhs_re(Nmat,NopenIn), rhs_im(Nmat,NopenIn), &
	                 solution_M0(Nmat,NopenIn), STAT=istatus)
	        IF (istatus /= 0) STOP 'Allocation failed in SolveSMatrixGeneric'

	        ALLOCATE(mat_B_local(NopenIn,NopenIn), mat_C_local(NopenIn,NopenIn), &
	                 Bsub(NopenIn,NopenIn), Csub(NopenIn,NopenIn), &
	                 mat_B_factor(NopenIn,NopenIn), solved_C(NopenIn,NopenIn), &
	                 ipiv_B(NopenIn), STAT=istatus)
	        IF (istatus /= 0) STOP 'Open-matrix allocation failed in SolveSMatrixGeneric'

	        matM_work = matM_in
	        rhs_re = REAL(matM0_in,KIND=8)
	        rhs_im = AIMAG(matM0_in)

	        CALL DSYTRF('L',Nmat,matM_work,Nmat,ipiv_M,work_query,-1,infox)
	        IF (infox /= 0) THEN
	                PRINT*, TRIM(label), ' DSYTRF workspace query failed: ', infox
	                STOP
	        ENDIF

	        lworkx = MAX(1,INT(work_query(1)))
	        ALLOCATE(workx(lworkx), STAT=istatus)
	        IF (istatus /= 0) STOP 'Workspace allocation failed in SolveSMatrixGeneric'

	        CALL DSYTRF('L',Nmat,matM_work,Nmat,ipiv_M,workx,lworkx,infox)
	        IF (infox /= 0) THEN
	                PRINT*, TRIM(label), ' DSYTRF failed: ', infox
	                STOP
	        ENDIF

	        CALL DSYTRS('L',Nmat,NopenIn,matM_work,Nmat,ipiv_M,rhs_re,Nmat,infox)
	        IF (infox /= 0) THEN
	                PRINT*, TRIM(label), ' DSYTRS failed for real part: ', infox
	                STOP
	        ENDIF

	        CALL DSYTRS('L',Nmat,NopenIn,matM_work,Nmat,ipiv_M,rhs_im,Nmat,infox)
	        IF (infox /= 0) THEN
	                PRINT*, TRIM(label), ' DSYTRS failed for imaginary part: ', infox
	                STOP
	        ENDIF

	        solution_M0 = CMPLX(rhs_re,rhs_im,KIND=8)
	        Bsub = MATMUL(TRANSPOSE(matM0_in),solution_M0)
	        Csub = MATMUL(TRANSPOSE(CONJG(matM0_in)),solution_M0)
	        mat_B_local = matM00_in - Bsub
	        mat_C_local = matM10_in - Csub

	        mat_B_factor = CONJG(mat_B_local)
	        solved_C = mat_C_local

	        CALL ZGETRF(NopenIn,NopenIn,mat_B_factor,NopenIn,ipiv_B,infox)
	        IF (infox /= 0) THEN
	                PRINT*, TRIM(label), ' ZGETRF failed: ', infox
	                STOP
	        ENDIF

	        CALL ZGETRS('N',NopenIn,NopenIn,mat_B_factor,NopenIn,ipiv_B, &
	                    solved_C,NopenIn,infox)
	        IF (infox /= 0) THEN
	                PRINT*, TRIM(label), ' ZGETRS failed: ', infox
	                STOP
	        ENDIF

	        Smat_out = (0d0,1d0) * &
	                (mat_B_local - MATMUL(TRANSPOSE(mat_C_local),solved_C))

	        PRINT*, TRIM(label), ' S matrix solved. N = ', Nmat, &
	                ' n_open = ', NopenIn

	        DEALLOCATE(matM_work, ipiv_M, rhs_re, rhs_im, solution_M0, &
	                   mat_B_local, mat_C_local, Bsub, Csub, &
	                   mat_B_factor, solved_C, ipiv_B, workx, STAT=istatus)

	END SUBROUTINE SolveSMatrixGeneric

	SUBROUTINE write_exchange_phase1_outputs(Smat_ordered, Smat_sym_in)

	        USE AtomDiatomskvp
	        USE generateparameters, ONLY: Brot

	        IMPLICIT NONE

	        COMPLEX(8), INTENT(IN) :: Smat_ordered(n_open,n_open)
	        COMPLEX(8), INTENT(IN) :: Smat_sym_in(n_sym_open,n_sym_open)

	        INTEGER :: old_in_open, sym_in_open
	        INTEGER :: find_open_position
	        INTEGER :: old_idx, os, s, c, old_open
	        INTEGER :: unit_cmp, unit_sum
	        REAL(8) :: p_ordered_projected, p_sym, diff, max_diff
	        REAL(8) :: sum_ordered_projected, sum_sym
	        COMPLEX(8) :: amp_ordered_projected

	        old_in_open = 0
	        DO old_open = 1, n_open
	                old_idx = open_idx(old_open)
	                IF (quant_mat(1,old_idx) == 0 .AND. &
	                    quant_mat(2,old_idx) == 0 .AND. &
	                    quant_mat(3,old_idx) == 0 .AND. &
	                    quant_mat(4,old_idx) == 0) THEN
	                        old_in_open = old_open
	                        EXIT
	                ENDIF
	        ENDDO

	        IF (old_in_open <= 0) THEN
	                PRINT*, 'Could not find ordered initial channel (0,0,0,0).'
	                STOP
	        ENDIF

	        sym_in_open = old_open_pos_to_sym_open(old_in_open)

	        OPEN(newunit=unit_cmp, file='phase1_smatrix_compare.dat', &
	             status='unknown', position='append', action='write')

	        WRITE(unit_cmp,'(A)') '# Energy_eV Jtot sym_open sym_index rep_j1 rep_k1 rep_j2 rep_k2 P_ordered_projected P_sym abs_diff'

	        max_diff = 0d0
	        sum_ordered_projected = 0d0
	        sum_sym = 0d0

	        DO os = 1, n_sym_open
	                s = sym_open_idx(os)
	                amp_ordered_projected = (0d0,0d0)
	                DO c = 1, sym_ncomp(s)
	                        old_open = find_open_position(sym_old_idx(c,s))
	                        IF (old_open <= 0) CYCLE
	                        amp_ordered_projected = amp_ordered_projected + &
	                                sym_coeff(c,s) * Smat_ordered(old_in_open,old_open)
	                ENDDO

	                p_ordered_projected = ABS(amp_ordered_projected)**2
	                p_sym = ABS(Smat_sym_in(sym_in_open,os))**2
	                diff = ABS(p_ordered_projected - p_sym)
	                max_diff = MAX(max_diff,diff)
	                sum_ordered_projected = sum_ordered_projected + p_ordered_projected
	                sum_sym = sum_sym + p_sym

	                WRITE(unit_cmp,'(ES20.10,1X,I5,1X,I6,1X,I6,4I6,3ES20.10)') &
	                        E*27.211399d0, Jtot, os, s, &
	                        sym_quant_mat(1,s), sym_quant_mat(2,s), &
	                        sym_quant_mat(3,s), sym_quant_mat(4,s), &
	                        p_ordered_projected, p_sym, diff
	        ENDDO

	        CLOSE(unit_cmp)

	        OPEN(newunit=unit_sum, file='phase1_smatrix_summary.dat', &
	             status='unknown', position='append', action='write')

	        WRITE(unit_sum,'(A)') '# Energy_eV Jtot ncf nsym n_open n_sym_open old_in sym_in sum_ordered_projected sum_sym max_abs_diff'
	        WRITE(unit_sum,'(ES20.10,1X,I5,6I8,3ES20.10)') &
	                E*27.211399d0, Jtot, ncf, nsym, n_open, n_sym_open, &
	                old_in_open, sym_in_open, sum_ordered_projected, sum_sym, max_diff
	        CLOSE(unit_sum)

	        PRINT*, 'Phase-1 exchange-symmetry comparison max |dP| = ', max_diff
	        PRINT*, 'Projected ordered probability sum = ', sum_ordered_projected
	        PRINT*, 'Symmetric-basis probability sum   = ', sum_sym

	END SUBROUTINE write_exchange_phase1_outputs

	SUBROUTINE write_exchange_phase2_outputs(Smat_projected, Smat_direct)

	        USE AtomDiatomskvp

	        IMPLICIT NONE

	        COMPLEX(8), INTENT(IN) :: Smat_projected(n_sym_open,n_sym_open)
	        COMPLEX(8), INTENT(IN) :: Smat_direct(n_sym_open,n_sym_open)

	        INTEGER :: old_in_open, sym_in_open
	        INTEGER :: find_open_position
	        INTEGER :: old_idx, os, ot, s, old_open
	        INTEGER :: unit_cmp, unit_sum
	        REAL(8) :: p_projected, p_direct, diff, max_diff
	        REAL(8) :: sum_projected, sum_direct
	        REAL(8) :: max_sdiff

	        old_in_open = 0
	        DO old_open = 1, n_open
	                old_idx = open_idx(old_open)
	                IF (quant_mat(1,old_idx) == 0 .AND. &
	                    quant_mat(2,old_idx) == 0 .AND. &
	                    quant_mat(3,old_idx) == 0 .AND. &
	                    quant_mat(4,old_idx) == 0) THEN
	                        old_in_open = old_open
	                        EXIT
	                ENDIF
	        ENDDO

	        IF (old_in_open <= 0) THEN
	                PRINT*, 'Could not find ordered initial channel (0,0,0,0).'
	                STOP
	        ENDIF

	        sym_in_open = old_open_pos_to_sym_open(old_in_open)

	        max_sdiff = 0d0
	        DO os = 1, n_sym_open
	        DO ot = 1, n_sym_open
	                max_sdiff = MAX(max_sdiff, &
	                        ABS(Smat_projected(os,ot) - Smat_direct(os,ot)))
	        ENDDO
	        ENDDO

	        OPEN(newunit=unit_cmp, file='phase2_smatrix_compare.dat', &
	             status='unknown', position='append', action='write')

	        WRITE(unit_cmp,'(A)') '# Energy_eV Jtot sym_open sym_index rep_j1 rep_k1 rep_j2 rep_k2 P_projected P_direct abs_diff'

	        max_diff = 0d0
	        sum_projected = 0d0
	        sum_direct = 0d0

	        DO os = 1, n_sym_open
	                s = sym_open_idx(os)

	                p_projected = ABS(Smat_projected(sym_in_open,os))**2
	                p_direct = ABS(Smat_direct(sym_in_open,os))**2
	                diff = ABS(p_projected - p_direct)
	                max_diff = MAX(max_diff,diff)
	                sum_projected = sum_projected + p_projected
	                sum_direct = sum_direct + p_direct

	                WRITE(unit_cmp,'(ES20.10,1X,I5,1X,I6,1X,I6,4I6,3ES20.10)') &
	                        E*27.211399d0, Jtot, os, s, &
	                        sym_quant_mat(1,s), sym_quant_mat(2,s), &
	                        sym_quant_mat(3,s), sym_quant_mat(4,s), &
	                        p_projected, p_direct, diff
	        ENDDO

	        CLOSE(unit_cmp)

	        OPEN(newunit=unit_sum, file='phase2_smatrix_summary.dat', &
	             status='unknown', position='append', action='write')

	        WRITE(unit_sum,'(A)') '# Energy_eV Jtot ncf nsym n_open n_sym_open N_sym old_in sym_in sum_projected sum_direct max_abs_dP max_abs_dS'
	        WRITE(unit_sum,'(ES20.10,1X,I5,7I8,4ES20.10)') &
	                E*27.211399d0, Jtot, ncf, nsym, n_open, n_sym_open, &
	                dim_x*nsym, old_in_open, sym_in_open, &
	                sum_projected, sum_direct, max_diff, max_sdiff
	        CLOSE(unit_sum)

	        PRINT*, 'Phase-2 direct symmetric comparison max |dP| = ', max_diff
	        PRINT*, 'Phase-2 direct symmetric comparison max |dS| = ', max_sdiff
	        PRINT*, 'Projected symmetric probability sum = ', sum_projected
	        PRINT*, 'Direct symmetric probability sum    = ', sum_direct

	END SUBROUTINE write_exchange_phase2_outputs
