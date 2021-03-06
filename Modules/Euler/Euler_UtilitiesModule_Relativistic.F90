!> Perform computations related to the 3+1, CFA Euler equations.
!> Find the equations in Rezzolla & Zanotti, Relativistic Hydrodynamics, 2013,
!> Equation 7.234.
MODULE Euler_UtilitiesModule_Relativistic

  USE KindModule, ONLY: &
    DP,       &
    Zero,     &
    SqrtTiny, &
    Half,     &
    One,      &
    Two,      &
    Four,     &
    Fourth
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDimsX
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Alpha,    &
    iGF_Beta_1,   &
    iGF_Beta_2,   &
    iGF_Beta_3
  USE FluidFieldsModule, ONLY: &
    nCF,    &
    iCF_D,  &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E,  &
    iCF_Ne, &
    nPF,    &
    iPF_D,  &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E,  &
    iPF_Ne, &
    iAF_P,  &
    iAF_T,  &
    iAF_Ye, &
    iAF_S,  &
    iAF_E,  &
    iAF_Gm, &
    iAF_Cs
  USE EquationOfStateModule, ONLY: &
    ComputeSoundSpeedFromPrimitive, &
    ComputeAuxiliary_Fluid,         &
    ComputePressureFromSpecificInternalEnergy
  USE UnitsModule, ONLY: &
    AtomicMassUnit
  USE TimersModule_Euler,ONLY: &
    TimersStart_Euler, &
    TimersStop_Euler,  &
    Timer_Euler_ComputeTimeStep
  USE Euler_ErrorModule, ONLY: &
    DescribeError_Euler

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputePrimitive_Euler_Relativistic
  PUBLIC :: ComputeConserved_Euler_Relativistic
  PUBLIC :: ComputeFromConserved_Euler_Relativistic
  PUBLIC :: ComputeTimeStep_Euler_Relativistic
  PUBLIC :: Eigenvalues_Euler_Relativistic
  PUBLIC :: AlphaMiddle_Euler_Relativistic
  PUBLIC :: Flux_X1_Euler_Relativistic
  PUBLIC :: Flux_X2_Euler_Relativistic
  PUBLIC :: Flux_X3_Euler_Relativistic
  PUBLIC :: StressTensor_Diagonal_Euler_Relativistic
  PUBLIC :: NumericalFlux_LLF_Euler_Relativistic
  PUBLIC :: NumericalFlux_HLL_Euler_Relativistic
  PUBLIC :: NumericalFlux_X1_HLLC_Euler_Relativistic
  PUBLIC :: NumericalFlux_X2_HLLC_Euler_Relativistic
  PUBLIC :: NumericalFlux_X3_HLLC_Euler_Relativistic

  INTERFACE ComputePrimitive_Euler_Relativistic
    MODULE PROCEDURE ComputePrimitive_Scalar
    MODULE PROCEDURE ComputePrimitive_Vector
  END INTERFACE ComputePrimitive_Euler_Relativistic

  INTERFACE ComputeConserved_Euler_Relativistic
    MODULE PROCEDURE ComputeConserved_Scalar
    MODULE PROCEDURE ComputeConserved_Vector
  END INTERFACE ComputeConserved_Euler_Relativistic


CONTAINS


  !> Compute the primitive variables from the conserved variables,
  !> a la Galeazzi et al., (2013), Phys. Rev. D., 88, 064009, Appendix C
  !> @todo Modify for tabular EOS
  SUBROUTINE ComputePrimitive_Scalar &
    ( CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
      PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
      GF_Gm11, GF_Gm22, GF_Gm33 )

    REAL(DP), INTENT(in)  :: CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne
    REAL(DP), INTENT(in)  :: GF_Gm11, GF_Gm22, GF_Gm33
    REAL(DP), INTENT(out) :: PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne

    REAL(DP) :: S, q, r, k, z0
    REAL(DP) :: W, eps, p, h

    S = SQRT( CF_S1**2 / GF_Gm11 + CF_S2**2 / GF_Gm22 + CF_S3**2 / GF_Gm33 )

    ! --- Eq. C2 ---

    q = CF_E / CF_D
    r = S    / CF_D
    k = r    / ( One + q )

    CALL SolveZ_Bisection( CF_D, CF_Ne, q, r, k, z0 )

    ! --- Eq. C15 ---

    W     = SQRT( One + z0**2 )
    PF_D  = CF_D / W
    PF_Ne = CF_Ne / W

    ! --- Eq. C16 ---

    eps = W * q - z0 * r + z0**2 / ( One + W )

    CALL ComputePressureFromSpecificInternalEnergy &
           ( PF_D, eps, AtomicMassUnit * PF_Ne / PF_D, p )

    h = One + eps + p / PF_D

    PF_V1 = ( CF_S1 / GF_Gm11 ) / ( CF_D * W * h )
    PF_V2 = ( CF_S2 / GF_Gm22 ) / ( CF_D * W * h )
    PF_V3 = ( CF_S3 / GF_Gm33 ) / ( CF_D * W * h )
    PF_E  = CF_D * ( eps + p / PF_D ) / W - p

  END SUBROUTINE ComputePrimitive_Scalar


  SUBROUTINE ComputePrimitive_Vector &
    ( CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
      PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
      GF_Gm11, GF_Gm22, GF_Gm33 )

    REAL(DP), INTENT(in)  :: CF_D(:), CF_S1(:), CF_S2(:), &
                             CF_S3(:), CF_E(:), CF_Ne(:)
    REAL(DP), INTENT(in)  :: GF_Gm11(:), GF_Gm22(:), GF_Gm33(:)
    REAL(DP), INTENT(out) :: PF_D(:), PF_V1(:), PF_V2(:), &
                             PF_V3(:), PF_E(:), PF_Ne(:)

    INTEGER :: iNodeX

    DO iNodeX = 1, SIZE( CF_D )

      CALL ComputePrimitive_Scalar &
             ( CF_D   (iNodeX), &
               CF_S1  (iNodeX), &
               CF_S2  (iNodeX), &
               CF_S3  (iNodeX), &
               CF_E   (iNodeX), &
               CF_Ne  (iNodeX), &
               PF_D   (iNodeX), &
               PF_V1  (iNodeX), &
               PF_V2  (iNodeX), &
               PF_V3  (iNodeX), &
               PF_E   (iNodeX), &
               PF_Ne  (iNodeX), &
               GF_Gm11(iNodeX), &
               GF_Gm22(iNodeX), &
               GF_Gm33(iNodeX) )

    END DO

  END SUBROUTINE ComputePrimitive_Vector


  !> Compute conserved variables from primitive variables.
  SUBROUTINE ComputeConserved_Scalar &
    ( PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
      CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
      Gm11, Gm22, Gm33, &
      AF_P )

    REAL(DP), INTENT(in)  :: PF_D, PF_V1, PF_V2, PF_V3, &
                             PF_E, PF_Ne, AF_P, &
                             Gm11, Gm22, Gm33
    REAL(DP), INTENT(out) :: CF_D, CF_S1, CF_S2, CF_S3, &
                             CF_E, CF_Ne

    REAL(DP) :: VSq, W, h

    VSq = Gm11 * PF_V1**2 + Gm22 * PF_V2**2 + Gm33 * PF_V3**2

    W = One / SQRT( One - VSq )
    h = One + ( PF_E + AF_P ) / PF_D

    CF_D  = W * PF_D
    CF_S1 = h * W**2 * PF_D * Gm11 * PF_V1
    CF_S2 = h * W**2 * PF_D * Gm22 * PF_V2
    CF_S3 = h * W**2 * PF_D * Gm33 * PF_V3
    CF_E  = h * W**2 * PF_D - AF_P - W * PF_D
    CF_Ne = W * PF_Ne

  END SUBROUTINE ComputeConserved_Scalar


  SUBROUTINE ComputeConserved_Vector &
    ( PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
      CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
      Gm11, Gm22, Gm33, &
      AF_P )

    REAL(DP), INTENT(in)  :: PF_D(:), PF_V1(:), PF_V2(:), PF_V3(:), &
                             PF_E(:), PF_Ne(:), AF_P(:), &
                             Gm11(:), Gm22(:), Gm33(:)
    REAL(DP), INTENT(out) :: CF_D(:), CF_S1(:), CF_S2(:), CF_S3(:), &
                             CF_E(:), CF_Ne(:)

    INTEGER :: iNodeX

    DO iNodeX = 1, SIZE( PF_D )

      CALL ComputeConserved_Scalar &
             ( PF_D (iNodeX), &
               PF_V1(iNodeX), &
               PF_V2(iNodeX), &
               PF_V3(iNodeX), &
               PF_E (iNodeX), &
               PF_Ne(iNodeX), &
               CF_D (iNodeX), &
               CF_S1(iNodeX), &
               CF_S2(iNodeX), &
               CF_S3(iNodeX), &
               CF_E (iNodeX), &
               CF_Ne(iNodeX), &
               Gm11 (iNodeX), &
               Gm22 (iNodeX), &
               Gm33 (iNodeX), &
               AF_P (iNodeX) )

    END DO

  END SUBROUTINE ComputeConserved_Vector


  !> Compute primitive variables, pressure, and sound-speed from conserved
  !> variables for a data block.
  SUBROUTINE ComputeFromConserved_Euler_Relativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, P, A )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(out) :: &
      P(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      A(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iX1, iX2, iX3

    ! --- Update primitive variables, pressure, and sound speed ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      CALL ComputePrimitive_Euler_Relativistic &
             ( U(1:nDOFX,iX1,iX2,iX3,iCF_D),         &
               U(1:nDOFX,iX1,iX2,iX3,iCF_S1),        &
               U(1:nDOFX,iX1,iX2,iX3,iCF_S2),        &
               U(1:nDOFX,iX1,iX2,iX3,iCF_S3),        &
               U(1:nDOFX,iX1,iX2,iX3,iCF_E),         &
               U(1:nDOFX,iX1,iX2,iX3,iCF_Ne),        &
               P(1:nDOFX,iX1,iX2,iX3,iPF_D),         &
               P(1:nDOFX,iX1,iX2,iX3,iPF_V1),        &
               P(1:nDOFX,iX1,iX2,iX3,iPF_V2),        &
               P(1:nDOFX,iX1,iX2,iX3,iPF_V3),        &
               P(1:nDOFX,iX1,iX2,iX3,iPF_E),         &
               P(1:nDOFX,iX1,iX2,iX3,iPF_Ne),        &
               G(1:nDOFX,iX1,iX2,iX3,iGF_Gm_dd_11),  &
               G(1:nDOFX,iX1,iX2,iX3,iGF_Gm_dd_22),  &
               G(1:nDOFX,iX1,iX2,iX3,iGF_Gm_dd_33) )

      CALL ComputeAuxiliary_Fluid &
             ( P(1:nDOFX,iX1,iX2,iX3,iPF_D ), &
               P(1:nDOFX,iX1,iX2,iX3,iPF_E ), &
               P(1:nDOFX,iX1,iX2,iX3,iPF_Ne), &
               A(1:nDOFX,iX1,iX2,iX3,iAF_P ), &
               A(1:nDOFX,iX1,iX2,iX3,iAF_T ), &
               A(1:nDOFX,iX1,iX2,iX3,iAF_Ye), &
               A(1:nDOFX,iX1,iX2,iX3,iAF_S ), &
               A(1:nDOFX,iX1,iX2,iX3,iAF_E ), &
               A(1:nDOFX,iX1,iX2,iX3,iAF_Gm), &
               A(1:nDOFX,iX1,iX2,iX3,iAF_Cs) )

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeFromConserved_Euler_Relativistic


  !> Loop over all the elements in the spatial domain and compute the minimum
  !> required time-step for numerical stability.
  SUBROUTINE ComputeTimeStep_Euler_Relativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, CFL, TimeStep )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)  :: &
      CFL
    REAL(DP), INTENT(out) :: &
      TimeStep

    INTEGER  :: iX1, iX2, iX3, iNodeX
    REAL(DP) :: dX(3), dt(3)
    REAL(DP) :: P(nDOFX,nPF)
    REAL(DP) :: Cs(nDOFX)
    REAL(DP) :: EigVals_X1(nCF,nDOFX), alpha_X1, &
                EigVals_X2(nCF,nDOFX), alpha_X2, &
                EigVals_X3(nCF,nDOFX), alpha_X3

    CALL TimersStart_Euler( Timer_Euler_ComputeTimeStep )

    TimeStep = HUGE( One )
    dt       = HUGE( One )

    ! --- Maximum wave-speeds ---

    alpha_X1 = -HUGE( One )
    alpha_X2 = -HUGE( One )
    alpha_X3 = -HUGE( One )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      dX(1) = MeshX(1) % Width(iX1)
      dX(2) = MeshX(2) % Width(iX2)
      dX(3) = MeshX(3) % Width(iX3)

      CALL ComputePrimitive_Euler_Relativistic &
             ( U(1:nDOFX,iX1,iX2,iX3,iCF_D ), &
               U(1:nDOFX,iX1,iX2,iX3,iCF_S1), &
               U(1:nDOFX,iX1,iX2,iX3,iCF_S2), &
               U(1:nDOFX,iX1,iX2,iX3,iCF_S3), &
               U(1:nDOFX,iX1,iX2,iX3,iCF_E ), &
               U(1:nDOFX,iX1,iX2,iX3,iCF_Ne), &
               P(1:nDOFX,iPF_D ), &
               P(1:nDOFX,iPF_V1), &
               P(1:nDOFX,iPF_V2), &
               P(1:nDOFX,iPF_V3), &
               P(1:nDOFX,iPF_E ), &
               P(1:nDOFX,iPF_Ne), &
               G(1:nDOFX,iX1,iX2,iX3,iGF_Gm_dd_11), &
               G(1:nDOFX,iX1,iX2,iX3,iGF_Gm_dd_22), &
               G(1:nDOFX,iX1,iX2,iX3,iGF_Gm_dd_33) )

      CALL ComputeSoundSpeedFromPrimitive &
             ( P(1:nDOFX,iPF_D), P(1:nDOFX,iPF_E), P(1:nDOFX,iPF_Ne), &
               Cs(1:nDOFX) )

      DO iNodeX = 1, nDOFX

        EigVals_X1(1:nCF,iNodeX) &
          = Eigenvalues_Euler_Relativistic &
              ( P (iNodeX,iPF_V1), &
                Cs(iNodeX),        &
                G (iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                P (iNodeX,iPF_V1), &
                P (iNodeX,iPF_V2), &
                P (iNodeX,iPF_V3), &
                G (iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                G (iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                G (iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                G (iNodeX,iX1,iX2,iX3,iGF_Alpha),    &
                G (iNodeX,iX1,iX2,iX3,iGF_Beta_1) )

      END DO

      alpha_X1 = MAX( alpha_X1, MAXVAL( ABS( EigVals_X1(1:nCF,1:nDOFX) ) ) )

      dt(1) = dX(1) / alpha_X1

      IF( nDimsX .GT. 1 )THEN

        DO iNodeX = 1, nDOFX

          EigVals_X2(1:nCF,iNodeX) &
            = Eigenvalues_Euler_Relativistic &
                ( P (iNodeX,iPF_V2), &
                  Cs(iNodeX),        &
                  G (iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                  P (iNodeX,iPF_V1), &
                  P (iNodeX,iPF_V2), &
                  P (iNodeX,iPF_V3), &
                  G (iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                  G (iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                  G (iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                  G (iNodeX,iX1,iX2,iX3,iGF_Alpha),    &
                  G (iNodeX,iX1,iX2,iX3,iGF_Beta_2) )
        END DO

        alpha_X2 = MAX( alpha_X2, MAXVAL( ABS( EigVals_X2(1:nCF,1:nDOFX) ) ) )

        dt(2) = dX(2) / alpha_X2

      END IF

      IF( nDimsX .GT. 2 )THEN

        DO iNodeX = 1, nDOFX

          EigVals_X3(1:nCF,iNodeX) &
            = Eigenvalues_Euler_Relativistic &
                ( P (iNodeX,iPF_V3),   &
                  Cs(iNodeX),          &
                  G (iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                  P (iNodeX,iPF_V1),   &
                  P (iNodeX,iPF_V2),   &
                  P (iNodeX,iPF_V3),   &
                  G (iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                  G (iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                  G (iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                  G (iNodeX,iX1,iX2,iX3,iGF_Alpha),    &
                  G (iNodeX,iX1,iX2,iX3,iGF_Beta_3) )

        END DO

        alpha_X3 = MAX( alpha_X3, MAXVAL( ABS( EigVals_X3(1:nCF,1:nDOFX) ) ) )

        dt(3) = dX(3) / alpha_X3

      END IF

      TimeStep = MIN( TimeStep, MINVAL( dt(1:3) ) )

    END DO
    END DO
    END DO

    TimeStep = MAX( CFL * TimeStep, SqrtTiny )

    CALL TimersStop_Euler( Timer_Euler_ComputeTimeStep )

  END SUBROUTINE ComputeTimeStep_Euler_Relativistic


  !> Compute the eigenvalues of the flux-Jacobian.
  !> Find the expressions in Font et al., (1998), Eqs. (14) and (18).
  !> @param Vi The ith contravariant component of the three-velocity.
  !> @param Gmii The ith covariant component of the spatial three-metric.
  !> @param Shift The ith contravariant component of the shift-vector.
  PURE FUNCTION Eigenvalues_Euler_Relativistic &
    ( Vi, Cs, Gmii, V1, V2, V3, Gm11, Gm22, Gm33, Lapse, Shift )

    REAL(DP), INTENT(in) :: Vi, Cs, Gmii, V1, V2, V3, &
                            Gm11, Gm22, Gm33, Lapse, Shift

    REAL(DP) :: VSq, Eigenvalues_Euler_Relativistic(nCF)

    VSq = Gm11 * V1**2 + Gm22 * V2**2 + Gm33 * V3**2

    Eigenvalues_Euler_Relativistic(1) &
      = Lapse / ( One - VSq * Cs**2 ) * ( Vi * ( One - Cs**2 ) &
        - Cs * SQRT( ( One - VSq ) * ( ( One - VSq * Cs**2 ) / Gmii &
           - Vi**2 * ( One - Cs**2 ) ) ) ) - Shift

    Eigenvalues_Euler_Relativistic(2) &
      = Lapse * Vi - Shift

    Eigenvalues_Euler_Relativistic(3) &
      = Lapse / ( One - VSq * Cs**2 ) * ( Vi * ( One - Cs**2 ) &
        + Cs * SQRT( ( One - VSq ) * ( ( One - VSq * Cs**2 ) / Gmii &
           - Vi**2 * ( One - Cs**2 ) ) ) ) - Shift

    Eigenvalues_Euler_Relativistic(4) &
      = Lapse * Vi - Shift

    Eigenvalues_Euler_Relativistic(5) &
      = Lapse * Vi - Shift

    Eigenvalues_Euler_Relativistic(6) &
      = Lapse * Vi - Shift

    RETURN
  END FUNCTION Eigenvalues_Euler_Relativistic


  !> Estimate the contact wave-speed as suggested by
  !> Mignone & Bodo, (2005), MNRAS, 364, 126.
  !> @param Shift The ith contravariant component of the shift-vector.
  !> @param Gmii The ith covariant component of the spatial three-metric.
  !> @todo Optimize special cases of quadratic formula solutions.
  REAL(DP) FUNCTION AlphaMiddle_Euler_Relativistic &
    ( DL, SL, tauL, F_DL, F_SL, F_tauL, DR, SR, tauR, F_DR, F_SR, F_tauR, &
      Gmii, aP, aM, Lapse, Shift )

    REAL(DP), INTENT(in) :: DL, SL, tauL, F_DL, F_SL, F_tauL, &
                            DR, SR, tauR, F_DR, F_SR, F_tauR, &
                            Gmii, aP, aM, Lapse, Shift

    REAL(DP) :: EL, F_EL, ER, F_ER, a2, a1, a0
    REAL(DP) :: E_HLL, S_HLL, FE_HLL, FS_HLL

#if defined HYDRO_RIEMANN_SOLVER_HLL

    AlphaMiddle_Euler_Relativistic = 1.0e1_DP

    RETURN

#endif

    EL   = tauL + DL
    F_EL = F_tauL + F_DL
    ER   = tauR + DR
    F_ER = F_tauR + F_DR

    E_HLL  = aP * ER + aM * EL + Lapse * ( F_EL - F_ER )
    S_HLL  = aP * SR + aM * SL + Lapse * ( F_SL - F_SR )
    FE_HLL = Lapse * ( aP * F_EL + aM * F_ER ) - aM * aP * ( ER - EL )
    FS_HLL = Lapse * ( aP * F_SL + aM * F_SR ) - aM * aP * ( SR - SL )

    ! --- Coefficients in quadratic equation ---

    a2 = Gmii**2 * ( FE_HLL + Shift * E_HLL )
    a1 = -Gmii * ( Lapse * E_HLL + FS_HLL + Shift * S_HLL )
    a0 = Lapse * S_HLL

    ! --- Accounting for special cases of the solution to a
    !     quadratic equation when a2 = 0 ---

    IF     ( ( ABS( a2 ) .LT. SqrtTiny ) .AND. ( ABS( a1 ) .LT. SqrtTiny ) &
            .AND. ( ABS( a0 ) .LT. SqrtTiny ) )THEN

      CALL DescribeError_Euler( 09 )

    ELSE IF( ( ABS( a2 ) .LT. SqrtTiny ) .AND. ( ABS( a1 ) .LT. SqrtTiny ) )THEN

      CALL DescribeError_Euler( 09, 'a0 < 0' )

    ELSE IF( ( ABS( a2 ) .LT. SqrtTiny ) .AND. ( ABS( a0 ) .LT. SqrtTiny ) )THEN

      AlphaMiddle_Euler_Relativistic = Zero

    ELSE IF( ABS( a2 ) .LT. SqrtTiny )THEN

      AlphaMiddle_Euler_Relativistic = -a0 / a1

    ELSE IF( ABS( a0 ) .LT. SqrtTiny )THEN

       AlphaMiddle_Euler_Relativistic = Zero

    ELSE

      AlphaMiddle_Euler_Relativistic &
        = ( -a1 - SQRT( MAX( a1**2 - Four * a2 * a0, SqrtTiny ) ) ) &
            / ( Two * a2 )
    END IF

    RETURN
  END FUNCTION AlphaMiddle_Euler_Relativistic


  !> Compute the physical flux in the X1-direction.
  !> @param Vi The ith contravariant components of the three-velocity.
  !> @param Gmii The ith covariant components of the spatial three-metric.
  !> @param Shift The first contravariant component of the shift-vector.
  PURE FUNCTION Flux_X1_Euler_Relativistic &
    ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, Lapse, Shift )

    REAL(DP), INTENT(in) :: D, V1, V2, V3, E, Ne, P, &
                            Gm11, Gm22, Gm33, Lapse, Shift

    REAL(DP) :: VSq, W, h, Flux_X1_Euler_Relativistic(nCF)

    VSq = Gm11 * V1**2 + Gm22 * V2**2 + Gm33 * V3**2
    W   = One / SQRT( One - VSq )
    h   = One + ( E + P ) / D

    Flux_X1_Euler_Relativistic(iCF_D)  &
      = D * W * ( V1 - Shift / Lapse )

    Flux_X1_Euler_Relativistic(iCF_S1) &
      = D * h * W**2 * Gm11 * V1  * ( V1 - Shift / Lapse ) + P

    Flux_X1_Euler_Relativistic(iCF_S2) &
      = D * h * W**2 * Gm22 * V2  * ( V1 - Shift / Lapse )

    Flux_X1_Euler_Relativistic(iCF_S3) &
      = D * h * W**2 * Gm33 * V3  * ( V1 - Shift / Lapse )

    Flux_X1_Euler_Relativistic(iCF_E)  &
      = D * W * ( h * W - One ) * ( V1 - Shift / Lapse ) + Shift / Lapse * P

    Flux_X1_Euler_Relativistic(iCF_Ne) &
      = Ne * W * ( V1 - Shift / Lapse )

    RETURN
  END FUNCTION Flux_X1_Euler_Relativistic


  !> Compute the physical flux in the X2-direction.
  !> @param Vi The ith contravariant components of the three-velocity.
  !> @param Gmii The ith covariant components of the spatial three-metric.
  !> @param Shift The first contravariant component of the shift-vector.
  PURE FUNCTION Flux_X2_Euler_Relativistic &
    ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, Lapse, Shift )

    REAL(DP), INTENT(in) :: D, V1, V2, V3, E, Ne, P, &
                            Gm11, Gm22, Gm33, Lapse, Shift

    REAL(DP) :: VSq, W, h, Flux_X2_Euler_Relativistic(nCF)

    VSq = Gm11 * V1**2 + Gm22 * V2**2 + Gm33 * V3**2
    W   = One / SQRT( One - VSq )
    h   = One + ( E + P ) / D

    Flux_X2_Euler_Relativistic(iCF_D)  &
      = D * W * ( V2 - Shift / Lapse )

    Flux_X2_Euler_Relativistic(iCF_S1) &
      = D * h * W**2 * Gm11 * V1  * ( V2 - Shift / Lapse )

    Flux_X2_Euler_Relativistic(iCF_S2) &
      = D * h * W**2 * Gm22 * V2  * ( V2 - Shift / Lapse ) + P

    Flux_X2_Euler_Relativistic(iCF_S3) &
      = D * h * W**2 * Gm33 * V3  * ( V2 - Shift / Lapse )

    Flux_X2_Euler_Relativistic(iCF_E)  &
      = D * W * ( h * W - One ) * ( V2 - Shift / Lapse ) + Shift / Lapse * P

    Flux_X2_Euler_Relativistic(iCF_Ne) &
      = Ne * W * ( V2 - Shift / Lapse )

    RETURN
  END FUNCTION Flux_X2_Euler_Relativistic


  !> Compute the physical flux in the X3-direction.
  !> @param Vi The ith contravariant components of the three-velocity.
  !> @param Gmii The ith covariant components of the spatial three-metric.
  !> @param Shift The first contravariant component of the shift-vector.
  PURE FUNCTION Flux_X3_Euler_Relativistic &
    ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, Lapse, Shift )

    REAL(DP), INTENT(in) :: D, V1, V2, V3, E, Ne, P, &
                            Gm11, Gm22, Gm33, Lapse, Shift

    REAL(DP) :: VSq, W, h, Flux_X3_Euler_Relativistic(nCF)

    VSq = Gm11 * V1**2 + Gm22 * V2**2 + Gm33 * V3**2
    W   = One / SQRT( One - VSq )
    h   = One + ( E + P ) / D

    Flux_X3_Euler_Relativistic(iCF_D)  &
      = D * W * ( V3 - Shift / Lapse )

    Flux_X3_Euler_Relativistic(iCF_S1) &
      = D * h * W**2 * Gm11 * V1  * ( V3 - Shift / Lapse )

    Flux_X3_Euler_Relativistic(iCF_S2) &
      = D * h * W**2 * Gm22 * V2  * ( V3 - Shift / Lapse )

    Flux_X3_Euler_Relativistic(iCF_S3) &
      = D * h * W**2 * Gm33 * V3  * ( V3 - Shift / Lapse ) + P

    Flux_X3_Euler_Relativistic(iCF_E)  &
      = D * W * ( h * W - One ) * ( V3 - Shift / Lapse ) + Shift / Lapse * P

    Flux_X3_Euler_Relativistic(iCF_Ne) &
      = Ne * W * ( V3 - Shift / Lapse )

    RETURN
  END FUNCTION Flux_X3_Euler_Relativistic


  !> Compute the diagonal elements of the stress-tensor, needed for the
  !> source-terms in the hydro equations.
  !> @param Si The ith covariant components of the conserved momentum-density.
  !> @param Vi The ith contravavriant components of the three-velocity.
  PURE FUNCTION StressTensor_Diagonal_Euler_Relativistic &
    ( S1, S2, S3, V1, V2, V3, P )

    REAL(DP), INTENT(in) :: S1, S2, S3, V1, V2, V3, P

    REAL(DP) :: StressTensor_Diagonal_Euler_Relativistic(3)

    StressTensor_Diagonal_Euler_Relativistic(1) = S1 * V1 + P
    StressTensor_Diagonal_Euler_Relativistic(2) = S2 * V2 + P
    StressTensor_Diagonal_Euler_Relativistic(3) = S3 * V3 + P

    RETURN
  END FUNCTION StressTensor_Diagonal_Euler_Relativistic


  !> Compute the Local-Lax-Friedrichs numerical flux at a given element
  !> interface, in a given dimension.
  PURE FUNCTION NumericalFlux_LLF_Euler_Relativistic &
    ( uL, uR, fL, fR, aP, aM )

    ! --- Local Lax-Friedrichs Flux ---

    REAL(DP), INTENT(in) :: uL(nCF), uR(nCF), fL(nCF), fR(nCF), aP, aM

    REAL(DP) :: NumericalFlux_LLF_Euler_Relativistic(nCF)

    REAL(DP) :: alpha

    alpha = MAX( aM, aP )

    NumericalFlux_LLF_Euler_Relativistic &
      = Half * ( fL + fR - alpha * ( uR - uL ) )

    RETURN
  END FUNCTION NumericalFlux_LLF_Euler_Relativistic


  !> Compute the Harten-Lax-van-Leer numerical flux at a given element
  !> interface, in a given dimension.
  PURE FUNCTION NumericalFlux_HLL_Euler_Relativistic &
    ( uL, uR, fL, fR, aP, aM )

    REAL(DP), INTENT(in) :: uL(nCF), uR(nCF), fL(nCF), fR(nCF), aP, aM

    REAL(DP) :: NumericalFlux_HLL_Euler_Relativistic(nCF)

    NumericalFlux_HLL_Euler_Relativistic &
      = ( aP * fL + aM * fR - aP * aM * ( uR - uL ) ) / ( aP + aM )

    RETURN
  END FUNCTION NumericalFlux_HLL_Euler_Relativistic


  !> Compute the Harten-Lax-van-Leer-Contact numerical flux at a given element
  !> in the X1-direction.
  !> @param Shift The first contravariant component of the shift-vector.
  !> @param Gm11 The first covariant component of the spatial three-metric.
  FUNCTION NumericalFlux_X1_HLLC_Euler_Relativistic &
    ( uL, uR, fL, fR, aP, aM, aC, Gm11, vL, vR, pL, pR, Lapse, Shift )

    REAL(DP), INTENT(in) :: uL(nCF), uR(nCF), fL(nCF), fR(nCF), &
                            aP, aM, aC, Gm11, vL, vR, pL, pR, Lapse, Shift

    REAL(DP) :: p, D, S1, S2, S3, E, Ne, UE, FE, FS, VelocityRatio, &
                NumericalFlux_X1_HLLC_Euler_Relativistic(nCF)

    IF( aM .EQ. Zero )THEN

      NumericalFlux_X1_HLLC_Euler_Relativistic = fL

    ELSEIF( aP .EQ. Zero )THEN

      NumericalFlux_X1_HLLC_Euler_Relativistic = fR

    ELSE

      ! --- From Mignone & Bodo (2005)
      !     Note the sign change on aM which is due to it being
      !     read in as positive but the formulae assuming it is negative ---

      IF( aC .GE. Zero )THEN

        VelocityRatio = ( -aM - Lapse * vL + Shift ) &
                      / ( -aM - Lapse * aC + Shift )

        ! --- uL_star ---
        UE = uL(iCF_E) + uL(iCF_D)
        FE = uL(iCF_S1) / Gm11 - Shift / Lapse * UE
        FS = uL(iCF_S1) * ( vL - Shift / Lapse ) + pL

        p  = ( Gm11 * aC * ( -aM * UE - Lapse * FE ) &
               - ( -aM * uL(iCF_S1) - Lapse * FS ) ) &
             / ( Lapse - Gm11 * aC * ( -aM + Shift ) )

        D  = uL(iCF_D)  * VelocityRatio

        S1 = uL(iCF_S1) * VelocityRatio &
               + Lapse * ( p - pL ) / ( -aM - Lapse * aC + Shift )

        S2 = uL(iCF_S2) * VelocityRatio

        S3 = uL(iCF_S3) * VelocityRatio

        E  = UE         * VelocityRatio + Lapse * ( p * aC - pL * vL ) &
               / ( -aM - Lapse * aC + Shift )

        Ne = uL(iCF_Ne) * VelocityRatio

      ELSE

        VelocityRatio = ( aP - Lapse * vR + Shift ) &
                      / ( aP - Lapse * aC + Shift )

        ! --- uR_star ---
        UE = uR(iCF_E) + uR(iCF_D)
        FE = uR(iCF_S1) / Gm11 - Shift / Lapse * UE
        FS = uR(iCF_S1) * ( vR - Shift / Lapse ) + pR

        p  = ( Gm11 * aC * ( aP * UE - Lapse * FE ) &
               - ( aP * uR(iCF_S1) - Lapse * FS ) ) &
               / ( Lapse - Gm11 * aC * ( aP + Shift ) )

        D  = uR(iCF_D)  * VelocityRatio

        S1 = uR(iCF_S1) * VelocityRatio &
               + Lapse * ( p - pR ) / ( aP - Lapse * aC + Shift )

        S2 = uR(iCF_S2) * VelocityRatio

        S3 = uR(iCF_S3) * VelocityRatio

        E  = UE         * VelocityRatio + Lapse * ( p * aC - pR * vR ) &
               / ( aP - Lapse * aC + Shift )

        Ne = uR(iCF_Ne) * VelocityRatio

      END IF

      NumericalFlux_X1_HLLC_Euler_Relativistic(iCF_D)  &
        = D  * ( aC - Shift / Lapse )

      NumericalFlux_X1_HLLC_Euler_Relativistic(iCF_S1) &
        = S1 * ( aC - Shift / Lapse ) + p

      NumericalFlux_X1_HLLC_Euler_Relativistic(iCF_S2) &
        = S2 * ( aC - Shift / Lapse )

      NumericalFlux_X1_HLLC_Euler_Relativistic(iCF_S3) &
        = S3 * ( aC - Shift / Lapse )

      NumericalFlux_X1_HLLC_Euler_Relativistic(iCF_E)  &
        = S1 / Gm11 - D * aC - Shift / Lapse * ( E - D )

      NumericalFlux_X1_HLLC_Euler_Relativistic(iCF_Ne) &
        = Ne * ( aC - Shift / Lapse )

    END IF

    RETURN
  END FUNCTION NumericalFlux_X1_HLLC_Euler_Relativistic


  !> Compute the Harten-Lax-van-Leer-Contact numerical flux at a given element
  !> in the X2-direction.
  !> @param Shift The second contravariant component of the shift-vector.
  !> @param Gm22 The second covariant component of the spatial three-metric.
  FUNCTION NumericalFlux_X2_HLLC_Euler_Relativistic &
    ( uL, uR, fL, fR, aP, aM, aC, Gm22, vL, vR, pL, pR, Lapse, Shift )

    REAL(DP), INTENT(in) :: uL(nCF), uR(nCF), fL(nCF), fR(nCF), &
                            aP, aM, aC, Gm22, vL, vR, pL, pR, Lapse, Shift

    REAL(DP) :: p, D, S1, S2, S3, E, Ne, UE, FE, FS, VelocityRatio
    REAL(DP) :: NumericalFlux_X2_HLLC_Euler_Relativistic(nCF)

    IF( aM .EQ. Zero )THEN

      NumericalFlux_X2_HLLC_Euler_Relativistic = fL

    ELSEIF( aP .EQ. Zero )THEN

      NumericalFlux_X2_HLLC_Euler_Relativistic = fR

    ELSE

      ! --- From Mignone & Bodo (2005)
      !     Note the sign change on aM which is due to it being
      !     read in as positive but the formulae assuming it is negative ---

      IF( aC .GE. Zero )THEN

        VelocityRatio = ( -aM - Lapse * vL + Shift ) &
                      / ( -aM - Lapse * aC + Shift )

        ! --- uL_star ---
        UE = uL(iCF_E) + uL(iCF_D)
        FE = uL(iCF_S2) / Gm22 - Shift / Lapse * UE
        FS = uL(iCF_S2) * ( vL - Shift / Lapse ) + pL

        p  = ( Gm22 * aC * ( -aM * UE - Lapse * FE ) &
               - ( -aM * uL(iCF_S2) - Lapse * FS ) ) &
             / ( Lapse - Gm22 * aC * ( -aM + Shift ) )

        D  = uL(iCF_D)  * VelocityRatio

        S1 = uL(iCF_S1) * VelocityRatio

        S2 = uL(iCF_S2) * VelocityRatio &
               + Lapse * ( p - pL ) / ( -aM - Lapse * aC + Shift )

        S3 = uL(iCF_S3) * VelocityRatio

        E  = UE         * VelocityRatio + Lapse * ( p * aC - pL * vL ) &
               / ( -aM - Lapse * aC + Shift )

        Ne = uL(iCF_Ne) * VelocityRatio

      ELSE

        VelocityRatio = ( aP - Lapse * vR + Shift ) &
                      / ( aP - Lapse * aC + Shift )

        ! --- uR_star ---
        UE = uR(iCF_E) + uR(iCF_D)
        FE = uR(iCF_S2) / Gm22 - Shift / Lapse * UE
        FS = uR(iCF_S2) * ( vR - Shift / Lapse ) + pR

        p  = ( Gm22 * aC * ( aP * UE - Lapse * FE ) &
               - ( aP * uR(iCF_S2) - Lapse * FS ) ) &
               / ( Lapse - Gm22 * aC * ( aP + Shift ) )

        D  = uR(iCF_D)  * VelocityRatio

        S1 = uR(iCF_S1) * VelocityRatio

        S2 = uR(iCF_S2) * VelocityRatio &
               + Lapse * ( p - pR ) / ( aP - Lapse * aC + Shift )

        S3 = uR(iCF_S3) * VelocityRatio

        E  = UE         * VelocityRatio + Lapse * ( p * aC - pR * vR ) &
               / ( aP - Lapse * aC + Shift )

        Ne = uR(iCF_Ne) * VelocityRatio

      END IF

      NumericalFlux_X2_HLLC_Euler_Relativistic(iCF_D)  &
        = D  * ( aC - Shift / Lapse )

      NumericalFlux_X2_HLLC_Euler_Relativistic(iCF_S1) &
        = S1 * ( aC - Shift / Lapse )

      NumericalFlux_X2_HLLC_Euler_Relativistic(iCF_S2) &
        = S2 * ( aC - Shift / Lapse ) + p

      NumericalFlux_X2_HLLC_Euler_Relativistic(iCF_S3) &
        = S3 * ( aC - Shift / Lapse )

      NumericalFlux_X2_HLLC_Euler_Relativistic(iCF_E)  &
        = S2 / Gm22 - D * aC - Shift / Lapse * ( E - D )

      NumericalFlux_X2_HLLC_Euler_Relativistic(iCF_Ne) &
        = Ne * ( aC - Shift / Lapse )

    END IF

    RETURN
  END FUNCTION NumericalFlux_X2_HLLC_Euler_Relativistic


  !> Compute the Harten-Lax-van-Leer-Contact numerical flux at a given element
  !> in the X3-direction.
  !> @param Shift The third contravariant component of the shift-vector.
  !> @param Gm33 The third covariant component of the spatial three-metric.
  FUNCTION NumericalFlux_X3_HLLC_Euler_Relativistic &
    ( uL, uR, fL, fR, aP, aM, aC, Gm33, vL, vR, pL, pR, Lapse, Shift )

    ! --- Shift is the third contravariant component of the shift-vector
    !     Gm is the third covariant component of the spatial three-metric ---

    REAL(DP), INTENT(in) :: uL(nCF), uR(nCF), fL(nCF), fR(nCF), &
                            aP, aM, aC, Gm33, vL, vR, pL, pR, Lapse, Shift

    REAL(DP) :: p, D, S1, S2, S3, E, Ne, UE, FE, FS, VelocityRatio
    REAL(DP) :: NumericalFlux_X3_HLLC_Euler_Relativistic(nCF)

    IF( aM .EQ. Zero )THEN

      NumericalFlux_X3_HLLC_Euler_Relativistic = fL

    ELSEIF( aP .EQ. Zero )THEN

      NumericalFlux_X3_HLLC_Euler_Relativistic = fR

    ELSE

      ! --- From Mignone & Bodo (2005)
      !     Note the sign change on aM which is due to it being
      !     read in as positive but the formulae assuming it is negative ---

      IF( aC .GE. Zero )THEN

        VelocityRatio = ( -aM - Lapse * vL + Shift ) &
                      / ( -aM - Lapse * aC + Shift )

        ! --- uL_star ---
        UE = uL(iCF_E) + uL(iCF_D)
        FE = uL(iCF_S3) / Gm33 - Shift / Lapse * UE
        FS = uL(iCF_S3) * ( vL - Shift / Lapse ) + pL

        p  = ( Gm33 * aC * ( -aM * UE - Lapse * FE ) &
               - ( -aM * uL(iCF_S3) - Lapse * FS ) ) &
             / ( Lapse - Gm33 * aC * ( -aM + Shift ) )

        D  = uL(iCF_D)  * VelocityRatio

        S1 = uL(iCF_S1) * VelocityRatio

        S2 = uL(iCF_S2) * VelocityRatio

        S3 = uL(iCF_S3) * VelocityRatio &
               + Lapse * ( p - pL ) / ( -aM - Lapse * aC + Shift )

        E  = UE         * VelocityRatio + Lapse * ( p * aC - pL * vL ) &
               / ( -aM - Lapse * aC + Shift )

        Ne = uL(iCF_Ne) * VelocityRatio

      ELSE

        VelocityRatio = ( aP - Lapse * vR + Shift ) &
                      / ( aP - Lapse * aC + Shift )

        ! --- uR_star ---
        UE = uR(iCF_E) + uR(iCF_D)
        FE = uR(iCF_S3) / Gm33 - Shift / Lapse * UE
        FS = uR(iCF_S3) * ( vR - Shift / Lapse ) + pR

        p  = ( Gm33 * aC * ( aP * UE - Lapse * FE ) &
               - ( aP * uR(iCF_S3) - Lapse * FS ) ) &
               / ( Lapse - Gm33 * aC * ( aP + Shift ) )

        D  = uR(iCF_D)  * VelocityRatio

        S1 = uR(iCF_S1) * VelocityRatio

        S2 = uR(iCF_S2) * VelocityRatio

        S3 = uR(iCF_S3) * VelocityRatio &
               + Lapse * ( p - pR ) / ( aP - Lapse * aC + Shift )

        E  = UE         * VelocityRatio + Lapse * ( p * aC - pR * vR ) &
               / ( aP - Lapse * aC + Shift )

        Ne = uR(iCF_Ne) * VelocityRatio

      END IF

      NumericalFlux_X3_HLLC_Euler_Relativistic(iCF_D)  &
        = D  * ( aC - Shift / Lapse )

      NumericalFlux_X3_HLLC_Euler_Relativistic(iCF_S1) &
        = S1 * ( aC - Shift / Lapse )

      NumericalFlux_X3_HLLC_Euler_Relativistic(iCF_S2) &
        = S2 * ( aC - Shift / Lapse )

      NumericalFlux_X3_HLLC_Euler_Relativistic(iCF_S3) &
        = S3 * ( aC - Shift / Lapse ) + p

      NumericalFlux_X3_HLLC_Euler_Relativistic(iCF_E)  &
        = S3 / Gm33 - D * aC - Shift / Lapse * ( E - D )

      NumericalFlux_X3_HLLC_Euler_Relativistic(iCF_Ne) &
        = Ne * ( aC - Shift / Lapse )

    END IF

    RETURN
  END FUNCTION NumericalFlux_X3_HLLC_Euler_Relativistic


  ! --- Auxiliary utilities for ComputePrimitive ---


  REAL(DP) FUNCTION FunZ( z, D, Ne, q, r, k )

    REAL(DP), INTENT(in) :: z, D, Ne, q, r, k

    REAL(DP) :: Wt, rhot, epst, pt, at, Ye

    ! --- Eq. C15 ---

    Wt = SQRT( One + z**2 )
    rhot = D / Wt

    ! --- Eq. C16 ---

    epst = Wt * q - z * r + z**2 / ( One + Wt )

    Ye = Ne * AtomicMassUnit / D

    CALL ComputePressureFromSpecificInternalEnergy &
           ( rhot, epst, Ye, pt )

    ! --- Eq. C20 ---

    at = pt / ( rhot * ( One + epst ) )

    ! --- Eq. C24 ---

    FunZ = z - k / ( ( Wt - z * k ) * ( One + at ) )

    RETURN
  END FUNCTION FunZ


  SUBROUTINE SolveZ_Bisection( CF_D, CF_Ne, q, r, k, z0 )

    REAL(DP), INTENT(in)  :: CF_D, CF_Ne, q, r, k
    REAL(DP), INTENT(out) :: z0

    LOGICAL             :: CONVERGED
    INTEGER             :: ITERATION
    REAL(DP)            :: za, zb, zc, dz
    REAL(DP)            :: fa, fb, fc
    REAL(DP), PARAMETER :: dz_min = 1.0e-8_DP
    INTEGER,  PARAMETER :: MAX_IT = 4 - INT( LOG( dz_min) / LOG( Two ) )

    ! --- Eq. C23 ---

    za = Half * k / SQRT( One - Fourth * k**2 ) - SqrtTiny
    zb = k        / SQRT( One - k**2 )          + SqrtTiny

    ! --- Compute FunZ for upper and lower bounds ---

    fa = FunZ( za, CF_D, CF_Ne, q, r, k )
    fb = FunZ( zb, CF_D, CF_Ne, q, r, k )

    ! --- Check that sign of FunZ changes across bounds ---

    IF( .NOT. fa * fb .LT. 0 )THEN

      WRITE(*,'(6x,A)') &
        'SolveZ_Bisection:'
      WRITE(*,'(8x,A)') &
        'Error: No Root in Interval'
      WRITE(*,'(8x,A,ES24.16E3)') 'CF_D: ', CF_D
      WRITE(*,'(8x,A,ES24.16E3)') 'q:    ', q
      WRITE(*,'(8x,A,ES24.16E3)') 'r:    ', r
      WRITE(*,'(8x,A,ES24.16E3)') 'k:    ', k
      WRITE(*,'(8x,A,2ES15.6E3)') &
        'za, zb = ', za, zb
      WRITE(*,'(8x,A,2ES15.6E3)') &
        'fa, fb = ', fa, fb
      CALL DescribeError_Euler( 08 )

    END IF

    dz = zb - za

    ITERATION = 0
    CONVERGED = .FALSE.
    DO WHILE ( .NOT. CONVERGED .AND. ITERATION .LT. MAX_IT )

      ITERATION = ITERATION + 1

      ! --- Compute midpoint, zc ---

      dz = Half * dz

      ! --- Bisection ---

      zc = za + dz

!!$      ! --- Regula Falsi ---
!!$
!!$      zc = ( za * fb - zb * fa  ) / ( fb - fa )

!!$      ! --- Illinois ---
!!$
!!$      IF( fa * fc .GT. Zero )THEN
!!$
!!$        zc = ( fb * za - Half * fa * zb ) / ( fb - Half * fa )
!!$
!!$      ELSE
!!$
!!$        zc = ( Half * fb * za - fa * zb ) / ( Half * fb - fa )
!!$
!!$      ENDIF

      ! --- Compute f(zc) for midpoint zc ---

      fc = FunZ( zc, CF_D, CF_Ne, q, r, k )

      ! --- Change zc to za or zb, depending on sign of fc ---

      IF( fa * fc .LT. Zero )THEN

        zb = zc
        fb = fc

      ELSE IF( fa * fc .GT. Zero )THEN

        za = zc
        fa = fc

      ELSE

        CONVERGED = .TRUE.

      END IF

      IF( ABS( dz / za ) .LT. dz_min ) CONVERGED = .TRUE.

!!$      IF( ITERATION .GT. MAX_IT - 3 )THEN
!!$
!!$        WRITE(*,*) 'Iter   = ', ITERATION
!!$        WRITE(*,*) 'za, zb = ', za, zb
!!$        WRITE(*,*) 'dz     = ', dz
!!$        WRITE(*,*)
!!$
!!$      END IF

    END DO

    z0 = zc

  END SUBROUTINE SolveZ_Bisection


END MODULE Euler_UtilitiesModule_Relativistic
