MODULE TwoMoment_UtilitiesModule_OrderV

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Two, Three
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor, &
    EddingtonFactor

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputePrimitive_TwoMoment
  PUBLIC :: ComputeConserved_TwoMoment
  PUBLIC :: Flux_X1
  PUBLIC :: NumericalFlux_LLF

CONTAINS


  SUBROUTINE ComputePrimitive_TwoMoment &
    ( N, G_d_1, G_d_2, G_d_3, D, I_u_1, I_u_2, I_u_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, nIterations_Option )

    REAL(DP), INTENT(in)  :: N, G_d_1, G_d_2, G_d_3 ! --- Index Down
    REAL(DP), INTENT(out) :: D, I_u_1, I_u_2, I_u_3 ! --- Index Up
    REAL(DP), INTENT(in)  ::    V_u_1, V_u_2, V_u_3 ! --- Index Up
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    INTEGER, INTENT(out), OPTIONAL :: nIterations_Option

    ! --- Parameters ---

    INTEGER,  PARAMETER :: M = 2
    INTEGER,  PARAMETER :: LWORK = 2 * M
    INTEGER,  PARAMETER :: MaxIterations = 100
    REAL(DP), PARAMETER :: Rtol = 1.0d-08

    ! --- Local Variables ---

    LOGICAL  :: CONVERGED
    INTEGER  :: i, k, mk, INFO
    REAL(DP) :: I_d_1, I_d_2, I_d_3
    REAL(DP) :: k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33
    REAL(DP) :: UVEC(4), CVEC(4)
    REAL(DP) :: GVEC(4,M), GVECm(4)
    REAL(DP) :: FVEC(4,M), FVECm(4)
    REAL(DP) :: LMAT(4,4), DET, Alpha(M)
    REAL(DP) :: BVEC(4), AMAT(4,M), WORK(LWORK)

    CVEC = [ N, G_d_1, G_d_2, G_d_3 ]

    ! --- Initial Guess ---

    D     = N
    I_u_1 = Zero
    I_u_2 = Zero
    I_u_3 = Zero

    I_d_1 = Gm_dd_11 * I_u_1
    I_d_2 = Gm_dd_22 * I_u_2
    I_d_3 = Gm_dd_33 * I_u_3

    k = 0
    CONVERGED = .FALSE.
    DO WHILE( .NOT. CONVERGED .AND. k < MaxIterations )

      k  = k + 1
      mk = MIN( M, k )

      UVEC = [ D, I_d_1, I_d_2, I_d_3 ]

      CALL ComputeEddingtonTensorComponents_dd &
             ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
               k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33 )

      DET = One - (   V_u_1**2 * k_dd_11 &
                    + V_u_2**2 * k_dd_22 &
                    + V_u_3**2 * k_dd_33 &
                    + Two * (   V_u_1 * V_u_2 * k_dd_12 &
                              + V_u_1 * V_u_3 * k_dd_13 &
                              + V_u_2 * V_u_3 * k_dd_23 ) )

      LMAT(1,1) = One
      LMAT(2,1) = - ( V_u_1 * k_dd_11 + V_u_2 * k_dd_12 + V_u_3 * k_dd_13 )
      LMAT(3,1) = - ( V_u_1 * k_dd_12 + V_u_2 * k_dd_22 + V_u_3 * k_dd_23 )
      LMAT(4,1) = - ( V_u_1 * k_dd_13 + V_u_2 * k_dd_23 + V_u_3 * k_dd_33 )

      LMAT(1,2) = - V_u_1
      LMAT(2,2) = One
      LMAT(3,2) = Zero
      LMAT(4,2) = Zero

      LMAT(1,3) = - V_u_2
      LMAT(2,3) = Zero
      LMAT(3,3) = One
      LMAT(4,3) = Zero

      LMAT(1,4) = - V_u_3
      LMAT(2,4) = Zero
      LMAT(3,4) = Zero
      LMAT(4,4) = One

      LMAT = LMAT / DET

      CALL DGEMV( 'N', 4, 4, One, LMAT, 4, CVEC, 1, Zero, GVEC(:,mk), 1 )

      FVEC(:,mk) = GVEC(:,mk) - UVEC

      IF( mk == 1 )THEN

        ! --- Picard Iteration ---

        GVECm = GVEC(:,mk)

      ELSE

        ! --- Anderson Accelerated Fixed-Point ---

        BVEC = - FVEC(:,mk)

        AMAT(:,1:mk-1) &
          = FVEC(:,1:mk-1) - SPREAD( FVEC(:,mk), DIM = 2, NCOPIES = mk-1 )

        CALL DGELS( 'N', 4, mk-1, 1, AMAT(:,1:mk-1), 4, BVEC, 4, &
                    WORK, LWORK, INFO )

        Alpha(1:mk-1) = BVEC(1:mk-1)
        Alpha(mk)     = One - SUM( Alpha(1:mk-1) )

        GVECm = Zero
        DO i = 1, mk

          GVECm = GVECm + Alpha(i) * GVEC(:,i)

        END DO

      END IF

      FVECm = GVECm - UVEC

      IF( ALL( ABS( FVECm ) <= Rtol * ABS( CVEC ) ) )THEN

        CONVERGED = .TRUE.

      END IF

      UVEC = GVECm

      IF( mk == M .AND. .NOT. CONVERGED )THEN

        GVEC = CSHIFT( GVEC, SHIFT = + 1, DIM = 2 )
        FVEC = CSHIFT( FVEC, SHIFT = + 1, DIM = 2 )

      END IF

      D     = UVEC(1)
      I_d_1 = UVEC(2); I_u_1 = I_d_1 / Gm_dd_11
      I_d_2 = UVEC(3); I_u_2 = I_d_2 / Gm_dd_22
      I_d_3 = UVEC(4); I_u_3 = I_d_3 / Gm_dd_33

    END DO

    IF( PRESENT( nIterations_Option ) )THEN

      nIterations_Option = k

    END IF

    IF( k == MaxIterations )THEN

      PRINT*
      PRINT*, "ComputePrimitive_TwoMoment"
      PRINT*
      PRINT*, "  N     = ", N
      PRINT*, "  G_d_1 = ", G_d_1
      PRINT*, "  G_d_2 = ", G_d_2
      PRINT*, "  G_d_3 = ", G_d_3
      PRINT*
      PRINT*, "  V_u_1 = ", V_u_1
      PRINT*, "  V_u_2 = ", V_u_2
      PRINT*, "  V_u_3 = ", V_u_3
      PRINT*

      PRINT*, "  Converged with k = ", k

      PRINT*
      PRINT*, "  FVECm = ", FVECm
      PRINT*

      PRINT*
      PRINT*, "  D     = ", D
      PRINT*, "  I_u_1 = ", I_u_1
      PRINT*, "  I_u_2 = ", I_u_2
      PRINT*, "  I_u_3 = ", I_u_3
      PRINT*

    END IF

  END SUBROUTINE ComputePrimitive_TwoMoment


  SUBROUTINE ComputeConserved_TwoMoment &
    ( D, I_u_1, I_u_2, I_u_3, N, G_d_1, G_d_2, G_d_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP), INTENT(in)  :: D, I_u_1, I_u_2, I_u_3 ! --- Index Up
    REAL(DP), INTENT(out) :: N, G_d_1, G_d_2, G_d_3 ! --- Index Down
    REAL(DP), INTENT(in)  ::    V_u_1, V_u_2, V_u_3 ! --- Index Up
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    REAL(DP) :: k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33

    CALL ComputeEddingtonTensorComponents_dd &
           ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
             k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33 )

    ! --- Conserved Number Density ---

    N = D + Gm_dd_11 * V_u_1 * I_u_1 &
          + Gm_dd_22 * V_u_2 * I_u_2 &
          + Gm_dd_33 * V_u_3 * I_u_3

    ! --- Conserved Number Flux Density (1) ---

    G_d_1 = Gm_dd_11 * I_u_1 &
              + (   V_u_1 * k_dd_11 &
                  + V_u_2 * k_dd_12 &
                  + V_u_3 * k_dd_13 ) * D

    ! --- Conserved Number Flux Density (2) ---

    G_d_2 = Gm_dd_22 * I_u_2 &
              + (   V_u_1 * k_dd_12 &
                  + V_u_2 * k_dd_22 &
                  + V_u_3 * k_dd_23 ) * D

    ! --- Conserved Number Flux Density (3) ---

    G_d_3 = Gm_dd_33 * I_u_3 &
              + (   V_u_1 * k_dd_13 &
                  + V_u_2 * k_dd_23 &
                  + V_u_3 * k_dd_33 ) * D

  END SUBROUTINE ComputeConserved_TwoMoment


  FUNCTION Flux_X1 &
    ( D, I_u_1, I_u_2, I_u_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP)             :: Flux_X1(4)
    REAL(DP), INTENT(in) :: D, I_u_1, I_u_2, I_u_3
    REAL(DP), INTENT(in) ::    V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    REAL(DP) :: FF, EF, a, b
    REAL(DP) :: h_u_1, h_d_1, h_d_2, h_d_3
    REAL(DP) :: K_ud_11, K_ud_12, K_ud_13

    FF = FluxFactor( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    EF = EddingtonFactor( D, FF )

    a = Half * ( One - EF )
    b = Half * ( Three * EF - One )

    h_u_1 = I_u_1 / ( FF * D )

    h_d_1 = Gm_dd_11 * I_u_1 / ( FF * D )
    h_d_2 = Gm_dd_22 * I_u_2 / ( FF * D )
    h_d_3 = Gm_dd_33 * I_u_3 / ( FF * D )

    K_ud_11 = ( a + b * h_u_1 * h_d_1 ) * D
    K_ud_12 =     ( b * h_u_1 * h_d_2 ) * D
    K_ud_13 =     ( b * h_u_1 * h_d_3 ) * D

    Flux_X1(1) = I_u_1 + D * V_u_1

    Flux_X1(2) = K_ud_11 + V_u_1 * Gm_dd_11 * I_u_1

    Flux_X1(3) = K_ud_12 + V_u_1 * Gm_dd_22 * I_u_2

    Flux_X1(4) = K_ud_13 + V_u_1 * Gm_dd_33 * I_u_3

    RETURN
  END FUNCTION Flux_X1


  SUBROUTINE ComputeEddingtonTensorComponents_dd &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33 )

    REAL(DP), INTENT(in)  :: &
      D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(out) :: &
      k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33

    REAL(DP) :: FF, EF, a, b
    REAL(DP) :: h_d_1, h_d_2, h_d_3

    FF = FluxFactor( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    EF = EddingtonFactor( D, FF )

    h_d_1 = Gm_dd_11 * I_u_1 / ( FF * D )
    h_d_2 = Gm_dd_22 * I_u_2 / ( FF * D )
    h_d_3 = Gm_dd_33 * I_u_3 / ( FF * D )

    a = Half * ( One - EF )
    b = Half * ( Three * EF - One )

    ! --- Diagonal Eddington Tensor Components ---

    k_dd_11 = a * Gm_dd_11 + b * h_d_1 * h_d_1
    k_dd_22 = a * Gm_dd_22 + b * h_d_2 * h_d_2
    k_dd_33 = a * Gm_dd_33 + b * h_d_3 * h_d_3

    ! --- Off-Diagonal Eddington Tensor Components ---

    k_dd_12 = b * h_d_1 * h_d_2
    k_dd_13 = b * h_d_1 * h_d_3
    k_dd_23 = b * h_d_2 * h_d_3

  END SUBROUTINE ComputeEddingtonTensorComponents_dd


  FUNCTION NumericalFlux_LLF( u_L, u_R, Flux_L, Flux_R, alpha )

    REAL(DP)             :: NumericalFlux_LLF
    REAL(DP), INTENT(in) :: u_L, u_R, flux_L, flux_R, alpha

    NumericalFlux_LLF &
      = Half * ( flux_L + flux_R - alpha * ( u_R - u_L ) )

    RETURN
  END FUNCTION NumericalFlux_LLF


END MODULE TwoMoment_UtilitiesModule_OrderV