!*********************************************************************************************
!
! Test selected angular potential matrix elements for Potential_Interface_3_selective.f90.
!
! The program compares, as functions of R:
!
!   1. the selected expansion used by the current pot_mat window;
!   2. all expansion terms available up to the current lambda/m maxima;
!   3. a direct four-dimensional angular quadrature over complex spherical harmonics.
!
! Target matrix elements:
!
!   <2,0;2,0|V|2,-1;4,1>
!   <2,0;2,0|V|4,-1;2,1>
!   <2,0;2,0|V|4,-1;4,1>
!
! Output files:
!
!   matrix_element_test_selective.dat
!   matrix_element_terms_selective.dat
!
!*********************************************************************************************

PROGRAM test_selective_matrix_elements

        USE AtomDiatomskvp, ONLY: dp, pi, gauleg, gaunt_coeff, assoc_legendre
        USE generateparameters, ONLY: read_input, pbasst
        USE Potential_Interface

        IMPLICIT NONE

        INTEGER, PARAMETER :: n_elements = 3
        INTEGER, PARAMETER :: nR_test = 5
        INTEGER, PARAMETER :: nx_direct = 10
        INTEGER, PARAMETER :: nphi_direct = 12

        INTEGER :: elem
        INTEGER :: iR
        INTEGER :: ipot
        INTEGER :: lambda1, lambda2, m
        INTEGER :: summary_unit
        INTEGER :: terms_unit

        INTEGER :: left(4,n_elements)
        INTEGER :: right(4,n_elements)

        REAL(dp) :: R
        REAL(dp) :: Rmin, Rmax
        REAL(dp) :: selected_value
        REAL(dp) :: all_available_value
        REAL(dp) :: selected_error
        REAL(dp) :: all_available_error
        REAL(dp) :: selected_rel_error
        REAL(dp) :: all_available_rel_error
        REAL(dp) :: denom
        REAL(dp) :: angular_factor
        REAL(dp) :: contribution
        REAL(dp) :: A_all(0:lambda1_max_default, &
                          0:lambda2_max_default, &
                          0:m_max_default)
        COMPLEX(dp) :: direct_value

        CALL read_input
        CALL build_potential_index

        left(:,1)  = (/2, 0, 2, 0/)
        right(:,1) = (/2,-1, 4, 1/)

        left(:,2)  = (/2, 0, 2, 0/)
        right(:,2) = (/4,-1, 2, 1/)

        left(:,3)  = (/2, 0, 2, 0/)
        right(:,3) = (/4,-1, 4, 1/)

        Rmin = REAL(pbasst(1)%pb_min, dp) + 0.03_dp
        Rmax = REAL(pbasst(1)%pb_max, dp) - 0.03_dp

        IF (Rmax <= Rmin) THEN
                PRINT*, 'Error: invalid R test range.'
                PRINT*, 'Rmin = ', Rmin
                PRINT*, 'Rmax = ', Rmax
                STOP
        ENDIF

        OPEN(newunit=summary_unit, file='matrix_element_test_selective.dat', &
             status='replace', action='write')
        OPEN(newunit=terms_unit, file='matrix_element_terms_selective.dat', &
             status='replace', action='write')

        WRITE(summary_unit,'(A)') '# Matrix element test for Potential_Interface_3_selective.f90'
        WRITE(summary_unit,'(A,6(1X,I0))') '# lambda/m window:', &
                lambda1_min_default, lambda1_max_default, &
                lambda2_min_default, lambda2_max_default, &
                m_min_default, m_max_default
        WRITE(summary_unit,'(A,I0)') '# n_pot = ', n_pot
        WRITE(summary_unit,'(A,3(1X,I0))') '# quadrature nR nx_direct nphi_direct = ', &
                nR_test, nx_direct, nphi_direct
        WRITE(summary_unit,'(A)') '# columns: elem R left_j1 left_k1 left_j2 left_k2 right_j1 right_k1 right_j2 right_k2 selected_expansion all_terms_to_max direct_real direct_imag selected_minus_direct all_minus_direct rel_selected rel_all'

        WRITE(terms_unit,'(A)') '# Selected-term contributions only'
        WRITE(terms_unit,'(A)') '# columns: elem R lambda1 lambda2 m A_lambda1_lambda2_m angular_factor contribution'

        DO iR = 1, nR_test

                IF (nR_test == 1) THEN
                        R = 0.5_dp * (Rmin + Rmax)
                ELSE
                        R = Rmin + REAL(iR - 1, dp) * (Rmax - Rmin) / REAL(nR_test - 1, dp)
                ENDIF

                CALL calculate_A_from_potential_v(R, A_all)

                DO elem = 1, n_elements

                        selected_value = selected_expansion_element(A_all, left(:,elem), right(:,elem))
                        all_available_value = all_available_expansion_element(A_all, left(:,elem), right(:,elem))
                        direct_value = direct_matrix_element(R, left(:,elem), right(:,elem))

                        selected_error = selected_value - REAL(direct_value, dp)
                        all_available_error = all_available_value - REAL(direct_value, dp)

                        denom = MAX(ABS(REAL(direct_value, dp)), 1.0d-30)
                        selected_rel_error = selected_error / denom
                        all_available_rel_error = all_available_error / denom

                        WRITE(summary_unit,'(I4,1X,ES20.10,1X,8(I5,1X),8(ES20.10,1X))') &
                                elem, R, &
                                left(1,elem), left(2,elem), left(3,elem), left(4,elem), &
                                right(1,elem), right(2,elem), right(3,elem), right(4,elem), &
                                selected_value, all_available_value, &
                                REAL(direct_value, dp), AIMAG(direct_value), &
                                selected_error, all_available_error, &
                                selected_rel_error, all_available_rel_error

                        DO ipot = 1, n_pot
                                lambda1 = pot_mat(1,ipot)
                                lambda2 = pot_mat(2,ipot)
                                m       = pot_mat(3,ipot)

                                angular_factor = angular_element_from_gaunt(lambda1, lambda2, m, &
                                        left(:,elem), right(:,elem))
                                contribution = A_all(lambda1,lambda2,m) * angular_factor

                                WRITE(terms_unit,'(I4,1X,ES20.10,1X,3(I5,1X),3(ES20.10,1X))') &
                                        elem, R, lambda1, lambda2, m, &
                                        A_all(lambda1,lambda2,m), angular_factor, contribution
                        ENDDO

                ENDDO
        ENDDO

        CLOSE(summary_unit)
        CLOSE(terms_unit)

        PRINT*, 'Matrix-element test completed.'
        PRINT*, 'Summary file: matrix_element_test_selective.dat'
        PRINT*, 'Selected-term file: matrix_element_terms_selective.dat'

CONTAINS

        REAL(dp) FUNCTION selected_expansion_element(A_all_in, left_state, right_state)
                IMPLICIT NONE

                REAL(dp), INTENT(IN) :: A_all_in(0:lambda1_max_default, &
                                                 0:lambda2_max_default, &
                                                 0:m_max_default)
                INTEGER, INTENT(IN) :: left_state(4)
                INTEGER, INTENT(IN) :: right_state(4)

                INTEGER :: ipot_local
                INTEGER :: lambda1_local, lambda2_local, m_local
                REAL(dp) :: angular_factor_local

                selected_expansion_element = 0.0_dp

                DO ipot_local = 1, n_pot
                        lambda1_local = pot_mat(1,ipot_local)
                        lambda2_local = pot_mat(2,ipot_local)
                        m_local       = pot_mat(3,ipot_local)

                        angular_factor_local = angular_element_from_gaunt(lambda1_local, &
                                lambda2_local, m_local, left_state, right_state)

                        selected_expansion_element = selected_expansion_element + &
                                A_all_in(lambda1_local,lambda2_local,m_local) * angular_factor_local
                ENDDO

        END FUNCTION selected_expansion_element


        REAL(dp) FUNCTION all_available_expansion_element(A_all_in, left_state, right_state)
                IMPLICIT NONE

                REAL(dp), INTENT(IN) :: A_all_in(0:lambda1_max_default, &
                                                 0:lambda2_max_default, &
                                                 0:m_max_default)
                INTEGER, INTENT(IN) :: left_state(4)
                INTEGER, INTENT(IN) :: right_state(4)

                INTEGER :: lambda1_local, lambda2_local, m_local
                REAL(dp) :: angular_factor_local

                all_available_expansion_element = 0.0_dp

                DO lambda1_local = 0, lambda1_max_default
                DO lambda2_local = 0, lambda2_max_default
                DO m_local = 0, MIN(lambda1_local, lambda2_local, m_max_default)

                        angular_factor_local = angular_element_from_gaunt(lambda1_local, &
                                lambda2_local, m_local, left_state, right_state)

                        all_available_expansion_element = all_available_expansion_element + &
                                A_all_in(lambda1_local,lambda2_local,m_local) * angular_factor_local

                ENDDO
                ENDDO
                ENDDO

        END FUNCTION all_available_expansion_element


        REAL(dp) FUNCTION angular_element_from_gaunt(lambda1, lambda2, m, left_state, right_state)
                IMPLICIT NONE

                INTEGER, INTENT(IN) :: lambda1, lambda2, m
                INTEGER, INTENT(IN) :: left_state(4)
                INTEGER, INTENT(IN) :: right_state(4)

                INTEGER :: j1, k1, j2, k2
                INTEGER :: j1_prime, k1_prime, j2_prime, k2_prime
                REAL(dp) :: G1p, G2p, G1m, G2m
                REAL(dp) :: basis_factor

                j1 = left_state(1)
                k1 = left_state(2)
                j2 = left_state(3)
                k2 = left_state(4)

                j1_prime = right_state(1)
                k1_prime = right_state(2)
                j2_prime = right_state(3)
                k2_prime = right_state(4)

                G1p = gaunt_coeff(j1,k1,lambda1, m,j1_prime,k1_prime)
                G2p = gaunt_coeff(j2,k2,lambda2,-m,j2_prime,k2_prime)

                G1m = gaunt_coeff(j1,k1,lambda1,-m,j1_prime,k1_prime)
                G2m = gaunt_coeff(j2,k2,lambda2, m,j2_prime,k2_prime)

                IF (m == 0) THEN
                        basis_factor = SQRT(2.0_dp*pi)
                ELSE
                        basis_factor = (-1.0_dp)**m * 2.0_dp * SQRT(pi)
                ENDIF

                angular_element_from_gaunt = basis_factor * 0.5_dp * (G1p*G2p + G1m*G2m)

        END FUNCTION angular_element_from_gaunt


        COMPLEX(dp) FUNCTION direct_matrix_element(R, left_state, right_state)
                IMPLICIT NONE

                REAL(dp), INTENT(IN) :: R
                INTEGER, INTENT(IN) :: left_state(4)
                INTEGER, INTENT(IN) :: right_state(4)

                INTEGER :: ix1, ix2, iphi1, iphi2
                REAL(dp) :: xgrid(nx_direct), wx(nx_direct)
                REAL(dp) :: theta1, theta2
                REAL(dp) :: phi1, phi2, phi
                REAL(dp) :: dphi
                REAL(dp) :: weight
                REAL(dp) :: V
                COMPLEX(dp) :: bra1, bra2, ket1, ket2

                CALL gauleg(-1.0_dp, 1.0_dp, xgrid, wx, nx_direct)

                dphi = 2.0_dp*pi / REAL(nphi_direct, dp)
                direct_matrix_element = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

                DO ix1 = 1, nx_direct
                        theta1 = ACOS(MAX(-1.0_dp, MIN(1.0_dp, xgrid(ix1))))

                DO ix2 = 1, nx_direct
                        theta2 = ACOS(MAX(-1.0_dp, MIN(1.0_dp, xgrid(ix2))))

                DO iphi1 = 1, nphi_direct
                        phi1 = (REAL(iphi1, dp) - 0.5_dp) * dphi

                DO iphi2 = 1, nphi_direct
                        phi2 = (REAL(iphi2, dp) - 0.5_dp) * dphi
                        phi = wrap_phi(phi2 - phi1)

                        bra1 = CONJG(spherical_harmonic(left_state(1), left_state(2), &
                                xgrid(ix1), phi1))
                        bra2 = CONJG(spherical_harmonic(left_state(3), left_state(4), &
                                xgrid(ix2), phi2))

                        ket1 = spherical_harmonic(right_state(1), right_state(2), &
                                xgrid(ix1), phi1)
                        ket2 = spherical_harmonic(right_state(3), right_state(4), &
                                xgrid(ix2), phi2)

                        V = potential_v(R, theta1, theta2, phi)
                        weight = wx(ix1) * wx(ix2) * dphi * dphi

                        direct_matrix_element = direct_matrix_element + &
                                weight * bra1 * bra2 * CMPLX(V, 0.0_dp, KIND=dp) * ket1 * ket2

                ENDDO
                ENDDO
                ENDDO
                ENDDO

        END FUNCTION direct_matrix_element


        COMPLEX(dp) FUNCTION spherical_harmonic(j, k, x, phi)
                IMPLICIT NONE

                INTEGER, INTENT(IN) :: j, k
                REAL(dp), INTENT(IN) :: x
                REAL(dp), INTENT(IN) :: phi

                INTEGER :: abs_k
                REAL(dp) :: pjk
                REAL(dp) :: negative_k_phase
                COMPLEX(dp) :: azimuth

                IF (ABS(k) > j) THEN
                        spherical_harmonic = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
                        RETURN
                ENDIF

                abs_k = ABS(k)
                pjk = assoc_legendre(j, abs_k, x)

                IF (k < 0) THEN
                        negative_k_phase = (-1.0_dp)**abs_k
                ELSE
                        negative_k_phase = 1.0_dp
                ENDIF

                azimuth = CMPLX(COS(REAL(k,dp)*phi), SIN(REAL(k,dp)*phi), KIND=dp)

                spherical_harmonic = negative_k_phase * pjk * azimuth / SQRT(2.0_dp*pi)

        END FUNCTION spherical_harmonic


        REAL(dp) FUNCTION wrap_phi(phi)
                IMPLICIT NONE

                REAL(dp), INTENT(IN) :: phi

                wrap_phi = MODULO(phi, 2.0_dp*pi)

        END FUNCTION wrap_phi

END PROGRAM test_selective_matrix_elements
