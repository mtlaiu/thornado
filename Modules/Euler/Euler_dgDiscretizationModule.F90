#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_EULER
#endif

MODULE Euler_dgDiscretizationModule

  USE KindModule, ONLY: &
    DP, Zero, SqrtTiny, &
    Half, One
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nDimsX
  USE TimersModule_Euler,  ONLY:      &
    TimersStart_Euler,                &
    TimersStop_Euler,                 &
    Timer_Euler_dgDiscretization,     &
    Timer_Euler_Divergence,           &
    Timer_Euler_Geometry,             &
    Timer_Euler_Gravity,              &
    Timer_Euler_SurfaceTerm,          &
    Timer_Euler_NumericalFlux,        &
    Timer_Euler_VolumeTerm,           &
    Timer_Euler_Increment,            &
    Timer_Euler_ComputePrimitive,     &
    Timer_Euler_ComputeFromPrimitive, &
    Timer_Euler_CopyIn,               &
    Timer_Euler_Permute,              &
    Timer_Euler_Interpolate,          &
    Timer_Euler_CopyOut
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, WeightsX_X1, &
    nDOFX_X2, WeightsX_X2, &
    nDOFX_X3, WeightsX_X3, &
    WeightsX_q
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    dLXdX1_q, LX_X1_Dn, LX_X1_Up, &
    dLXdX2_q, LX_X2_Dn, LX_X2_Up, &
    dLXdX3_q, LX_X3_Dn, LX_X3_Up
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    nGF,          &
    iGF_h_1,      &
    iGF_h_2,      &
    iGF_h_3,      &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_SqrtGm,   &
    iGF_Alpha,    &
    iGF_Beta_1,   &
    iGF_Beta_2,   &
    iGF_Beta_3,   &
    iGF_Phi_N,    &
    CoordinateSystem
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX_FromScaleFactors
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive, &
    ComputeSoundSpeedFromPrimitive
  USE Euler_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Euler
  USE Euler_UtilitiesModule, ONLY: &
    ComputePrimitive_Euler,      &
    Eigenvalues_Euler,           &
    AlphaMiddle_Euler,           &
    Flux_X1_Euler,               &
    Flux_X2_Euler,               &
    Flux_X3_Euler,               &
    StressTensor_Diagonal_Euler, &
    Euler_NumericalFlux_X1,      &
    Euler_NumericalFlux_X2,      &
    Euler_NumericalFlux_X3

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: Euler_ComputeIncrement_DG_Explicit


CONTAINS


  SUBROUTINE Euler_ComputeIncrement_DG_Explicit &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU, SuppressBC_Option )

    INTEGER, INTENT(in)            :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)           :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout)        :: &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(out)          :: &
      dU(:,iX_B0(1):,iX_B0(2):,iX_B0(3):,:)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      SuppressBC_Option

    INTEGER :: iNodeX, iX1, iX2, iX3, iCF
    LOGICAL :: SuppressBC

    CALL TimersStart_Euler( Timer_Euler_dgDiscretization )

    ASSOCIATE( dX1 => MeshX(1) % Width, &
               dX2 => MeshX(2) % Width, &
               dX3 => MeshX(3) % Width )

    CALL TimersStart_Euler( Timer_Euler_CopyIn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: G, U, iX_B0, iX_E0, iX_B1, iX_E1, &
    !$OMP          dX1, dX2, dX3 ) &
    !$OMP MAP( alloc: dU )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( G, U, iX_B0, iX_E0, iX_B1, iX_E1, &
    !$ACC         dX1, dX2, dX3 ) &
    !$ACC CREATE( dU )
#endif

    CALL TimersStop_Euler( Timer_Euler_CopyIn )

    SuppressBC = .FALSE.
    IF( PRESENT( SuppressBC_Option ) ) &
      SuppressBC = SuppressBC_Option

    IF( .NOT. SuppressBC ) &
      CALL ApplyBoundaryConditions_Euler &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( dU, iX_B1, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iCF = 1, nCF
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX

      dU(iNodeX,iX1,iX2,iX3,iCF) = Zero

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStart_Euler( Timer_Euler_Divergence )

    CALL ComputeIncrement_Divergence_X1 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    CALL ComputeIncrement_Divergence_X2 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    CALL ComputeIncrement_Divergence_X3 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    CALL TimersStop_Euler( Timer_Euler_Divergence )

    ! --- Multiply Inverse Mass Matrix ---

    CALL TimersStart_Euler( Timer_Euler_Increment )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
    !$OMP PRIVATE( iNodeX )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( G, dU, dX1, dX2, dX3, WeightsX_q, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iCF = 1, nCF
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX

      dU(iNodeX,iX1,iX2,iX3,iCF) &
        = dU(iNodeX,iX1,iX2,iX3,iCF) &
            / ( WeightsX_q(iNodeX) * G(iNodeX,iX1,iX2,iX3,iGF_SqrtGm) &
                  * dX1(iX1) * dX2(iX2) * dX3(iX3) )

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_Increment )

    CALL TimersStart_Euler( Timer_Euler_Geometry )

    CALL ComputeIncrement_Geometry &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    CALL TimersStop_Euler( Timer_Euler_Geometry )

    CALL TimersStart_Euler( Timer_Euler_Gravity )

    CALL ComputeIncrement_Gravity &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    CALL TimersStop_Euler( Timer_Euler_Gravity )

#ifdef THORNADO_DEBUG_EULER
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( dU )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( dU )
#endif
    WRITE(*,'(A20,7I4)')     'MAXLOC(dU)', MAXLOC(dU)
    WRITE(*,'(A20,ES23.15)') 'MAXVAL(dU)', MAXVAL(dU)
#endif

    CALL TimersStart_Euler( Timer_Euler_CopyOut )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: dU, U ) &
    !$OMP MAP( release: G, iX_B0, iX_E0, iX_B1, iX_E1, &
    !$OMP               dX1, dX2, dX3 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT( dU, U ) &
    !$ACC DELETE( G, iX_B0, iX_E0, iX_B1, iX_E1, &
    !$ACC             dX1, dX2, dX3 )
#endif

    CALL TimersStop_Euler( Timer_Euler_CopyOut )

    END ASSOCIATE

    CALL TimersStop_Euler( Timer_Euler_dgDiscretization )

  END SUBROUTINE Euler_ComputeIncrement_DG_Explicit


  SUBROUTINE ComputeIncrement_Divergence_X1 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nGF), &
      U (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCF)
    REAL(DP), INTENT(inout) :: &
      dU(1:nDOFX,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),1:nCF)

    INTEGER  :: nK(3), nK_X1(3), nCF_K, nCF_F, nGF_F
    INTEGER  :: iNodeX, iNodeX_X1, iX1, iX2, iX3, iCF, iGF
    REAL(DP) :: Alpha, AlphaMns, AlphaPls, AlphaMdl
    REAL(DP) :: uPF_L(nPF), uPF_R(nPF), uPF_K(nPF)
    REAL(DP) :: Flux_X1_L(nCF), Flux_X1_R(nCF), Flux_X1_K(nCF)
    REAL(DP) :: EigVals_L(nCF), EigVals_R(nCF)
    REAL(DP) :: P_L , P_R, P_K
    REAL(DP) :: Cs_L, Cs_R

    REAL(DP) :: G_K          (nDOFX,   nGF,iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(1)-1:iX_E0(1)+1)
    REAL(DP) :: G_F          (nDOFX_X1,nGF,iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(1)  :iX_E0(1)+1)

    REAL(DP) :: uCF_K        (nDOFX,   nCF,iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(1)-1:iX_E0(1)+1)
    REAL(DP) :: uCF_L        (nDOFX_X1,nCF,iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(1)  :iX_E0(1)+1)
    REAL(DP) :: uCF_R        (nDOFX_X1,nCF,iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(1)  :iX_E0(1)+1)

    REAL(DP) :: dU_X1        (nDOFX,   nCF,iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(1)  :iX_E0(1)  )
    REAL(DP) :: Flux_X1_q    (nDOFX,   nCF,iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(1)  :iX_E0(1)  )
    REAL(DP) :: NumericalFlux(nDOFX_X1,nCF,iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(1)  :iX_E0(1)+1)

    IF( iX_E0(1) .EQ. iX_B0(1) ) RETURN

    nK    = iX_E0 - iX_B0 + 1      ! Number of Elements per Spatial Dimension
    nK_X1 = nK + [1,0,0]           ! Number of X1 Faces per Spatial Dimension
    nCF_K = nCF * PRODUCT( nK )    ! Number of Fluid Fields in Domain
    nCF_F = nCF * PRODUCT( nK_X1 ) ! Number of Fluid Fields on Interfaces
    nGF_F = nGF * PRODUCT( nK_X1 ) ! Number of Geometry Fields on Interfaces

    ASSOCIATE( dX2 => MeshX(2) % Width, dX3 => MeshX(3) % Width )

    CALL TimersStart_Euler( Timer_Euler_CopyIn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: dX2, dX3, iX_B0, iX_E0, iX_B1, iX_E1 ) &
    !$OMP MAP( alloc: G_K, G_F, uCF_K, uCF_L, uCF_R, &
    !$OMP             dU_X1, Flux_X1_q, NumericalFlux )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( dX2, dX3, iX_B0, iX_E0, iX_B1, iX_E1 ) &
    !$ACC CREATE( G_K, G_F, uCF_K, uCF_L, uCF_R, &
    !$ACC         dU_X1, Flux_X1_q, NumericalFlux )
#endif

    CALL TimersStop_Euler( Timer_Euler_CopyIn )

    ! --- Geometry Fields in Element Nodes ---

    CALL TimersStart_Euler( Timer_Euler_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( G_K, G, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX1 = iX_B0(1) - 1, iX_E0(1) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iGF = 1, nGF
    DO iNodeX = 1, nDOFX

      G_K(iNodeX,iGF,iX2,iX3,iX1) = G(iNodeX,iX1,iX2,iX3,iGF)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_Permute )

    ! --- Interpolate Geometry Fields on Shared Face ---

    CALL TimersStart_Euler( Timer_Euler_Interpolate )

    ! --- Face States (Average of Left and Right States) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nGF_F, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             G_K(1,1,iX_B0(2),iX_B0(3),iX_B0(1)-1), nDOFX, Zero, &
             G_F(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nGF_F, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             G_K(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX, Half, &
             G_F(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1 )

    CALL TimersStop_Euler( Timer_Euler_Interpolate )

    !---------------------
    ! --- Surface Term ---
    !---------------------

    CALL TimersStart_Euler( Timer_Euler_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( G_F, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX1 = iX_B0(1), iX_E0(1) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iGF = 1, nGF
    DO iNodeX_X1 = 1, nDOFX_X1

      G_F(iNodeX_X1,iGF,iX2,iX3,iX1) &
        = MAX( G_F(iNodeX_X1,iGF,iX2,iX3,iX1), SqrtTiny )

    END DO
    END DO
    END DO
    END DO
    END DO

    ! --- Compute metric quantities on interfaces from scale factors ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( G_F, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX1 = iX_B0(1), iX_E0(1) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iNodeX_X1 = 1, nDOFX_X1

      G_F         (iNodeX_X1,iGF_Gm_dd_11,iX2,iX3,iX1) &
        = MAX( G_F(iNodeX_X1,iGF_h_1     ,iX2,iX3,iX1)**2, SqrtTiny )

      G_F         (iNodeX_X1,iGF_Gm_dd_22,iX2,iX3,iX1) &
        = MAX( G_F(iNodeX_X1,iGF_h_2     ,iX2,iX3,iX1)**2, SqrtTiny )

      G_F         (iNodeX_X1,iGF_Gm_dd_33,iX2,iX3,iX1) &
        = MAX( G_F(iNodeX_X1,iGF_h_3     ,iX2,iX3,iX1)**2, SqrtTiny )

      G_F        (iNodeX_X1,iGF_SqrtGm,iX2,iX3,iX1) &
        = G_F    (iNodeX_X1,iGF_h_1   ,iX2,iX3,iX1) &
            * G_F(iNodeX_X1,iGF_h_2   ,iX2,iX3,iX1) &
            * G_F(iNodeX_X1,iGF_h_3   ,iX2,iX3,iX1)

    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( uCF_K, U, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX1 = iX_B0(1) - 1, iX_E0(1) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iCF = 1, nCF
    DO iNodeX = 1, nDOFX

      uCF_K(iNodeX,iCF,iX2,iX3,iX1) = U(iNodeX,iX1,iX2,iX3,iCF)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_Permute )

    ! --- Interpolate Fluid Fields ---

    CALL TimersStart_Euler( Timer_Euler_Interpolate )

    ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nCF_F, nDOFX, One, LX_X1_Up, nDOFX_X1, &
             uCF_K(1,1,iX_B0(2),iX_B0(3),iX_B0(1)-1), nDOFX, Zero, &
             uCF_L(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nCF_F, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
             uCF_K(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX, Zero, &
             uCF_R(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1 )

    CALL TimersStop_Euler( Timer_Euler_Interpolate )

    ! --- Numerical Flux ---

    CALL TimersStart_Euler( Timer_Euler_SurfaceTerm )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( Flux_X1_L, Flux_X1_R, P_L, P_R, Cs_L, Cs_R, &
    !$OMP          uPF_L, uPF_R, alpha )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( Flux_X1_L, Flux_X1_R, P_L, P_R, Cs_L, Cs_R, &
    !$ACC          uPF_L, uPF_R, alpha ) &
    !$ACC PRESENT( uCF_L, uCF_R, NumericalFlux, G_F, dX2, dX3, WeightsX_X1, &
    !$ACC          iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( Flux_X1_L, Flux_X1_R, P_L, P_R, Cs_L, Cs_R, &
    !$OMP          uPF_L, uPF_R, alpha )
#endif
    DO iX1 = iX_B0(1), iX_E0(1) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)

      DO iNodeX_X1 = 1, nDOFX_X1

        ! --- Left State Primitive, etc. ---

        CALL TimersStart_Euler( Timer_Euler_ComputePrimitive )

        CALL ComputePrimitive_Euler &
               ( uCF_L(iNodeX_X1,iCF_D ,iX2,iX3,iX1),     &
                 uCF_L(iNodeX_X1,iCF_S1,iX2,iX3,iX1),     &
                 uCF_L(iNodeX_X1,iCF_S2,iX2,iX3,iX1),     &
                 uCF_L(iNodeX_X1,iCF_S3,iX2,iX3,iX1),     &
                 uCF_L(iNodeX_X1,iCF_E ,iX2,iX3,iX1),     &
                 uCF_L(iNodeX_X1,iCF_Ne,iX2,iX3,iX1),     &
                 uPF_L(iPF_D ),                           &
                 uPF_L(iPF_V1),                           &
                 uPF_L(iPF_V2),                           &
                 uPF_L(iPF_V3),                           &
                 uPF_L(iPF_E ),                           &
                 uPF_L(iPF_Ne),                           &
                 G_F(iNodeX_X1,iGF_Gm_dd_11,iX2,iX3,iX1), &
                 G_F(iNodeX_X1,iGF_Gm_dd_22,iX2,iX3,iX1), &
                 G_F(iNodeX_X1,iGF_Gm_dd_33,iX2,iX3,iX1) )

        CALL TimersStop_Euler( Timer_Euler_ComputePrimitive )

        CALL TimersStart_Euler( Timer_Euler_ComputeFromPrimitive )

        CALL ComputePressureFromPrimitive &
               ( uPF_L(iPF_D), uPF_L(iPF_E), uPF_L(iPF_Ne), P_L  )

        CALL ComputeSoundSpeedFromPrimitive &
               ( uPF_L(iPF_D), uPF_L(iPF_E), uPF_L(iPF_Ne), Cs_L )

        CALL TimersStop_Euler( Timer_Euler_ComputeFromPrimitive )

        EigVals_L(1:nCF) = Eigenvalues_Euler &
                             ( uPF_L(iPF_V1),                             &
                               Cs_L,                                      &
                               G_F  (iNodeX_X1,iGF_Gm_dd_11,iX2,iX3,iX1), &
                               uPF_L(iPF_V1),                             &
                               uPF_L(iPF_V2),                             &
                               uPF_L(iPF_V3),                             &
                               G_F  (iNodeX_X1,iGF_Gm_dd_11,iX2,iX3,iX1), &
                               G_F  (iNodeX_X1,iGF_Gm_dd_22,iX2,iX3,iX1), &
                               G_F  (iNodeX_X1,iGF_Gm_dd_33,iX2,iX3,iX1), &
                               G_F  (iNodeX_X1,iGF_Alpha,   iX2,iX3,iX1), &
                               G_F  (iNodeX_X1,iGF_Beta_1,  iX2,iX3,iX1) )

        Flux_X1_L(1:nCF) &
          = Flux_X1_Euler &
              ( uPF_L(iPF_D ),                           &
                uPF_L(iPF_V1),                           &
                uPF_L(iPF_V2),                           &
                uPF_L(iPF_V3),                           &
                uPF_L(iPF_E ),                           &
                uPF_L(iPF_Ne),                           &
                P_L,                                     &
                G_F(iNodeX_X1,iGF_Gm_dd_11,iX2,iX3,iX1), &
                G_F(iNodeX_X1,iGF_Gm_dd_22,iX2,iX3,iX1), &
                G_F(iNodeX_X1,iGF_Gm_dd_33,iX2,iX3,iX1), &
                G_F(iNodeX_X1,iGF_Alpha,   iX2,iX3,iX1), &
                G_F(iNodeX_X1,iGF_Beta_1,  iX2,iX3,iX1) )

        ! --- Right State Primitive, etc. ---

        CALL TimersStart_Euler( Timer_Euler_ComputePrimitive )

        CALL ComputePrimitive_Euler &
               ( uCF_R(iNodeX_X1,iCF_D ,iX2,iX3,iX1),     &
                 uCF_R(iNodeX_X1,iCF_S1,iX2,iX3,iX1),     &
                 uCF_R(iNodeX_X1,iCF_S2,iX2,iX3,iX1),     &
                 uCF_R(iNodeX_X1,iCF_S3,iX2,iX3,iX1),     &
                 uCF_R(iNodeX_X1,iCF_E ,iX2,iX3,iX1),     &
                 uCF_R(iNodeX_X1,iCF_Ne,iX2,iX3,iX1),     &
                 uPF_R(iPF_D ),                           &
                 uPF_R(iPF_V1),                           &
                 uPF_R(iPF_V2),                           &
                 uPF_R(iPF_V3),                           &
                 uPF_R(iPF_E ),                           &
                 uPF_R(iPF_Ne),                           &
                 G_F(iNodeX_X1,iGF_Gm_dd_11,iX2,iX3,iX1), &
                 G_F(iNodeX_X1,iGF_Gm_dd_22,iX2,iX3,iX1), &
                 G_F(iNodeX_X1,iGF_Gm_dd_33,iX2,iX3,iX1) )

        CALL TimersStop_Euler( Timer_Euler_ComputePrimitive )

        CALL TimersStart_Euler( Timer_Euler_ComputeFromPrimitive )

        CALL ComputePressureFromPrimitive &
               ( uPF_R(iPF_D), uPF_R(iPF_E), uPF_R(iPF_Ne), P_R  )

        CALL ComputeSoundSpeedFromPrimitive &
               ( uPF_R(iPF_D), uPF_R(iPF_E), uPF_R(iPF_Ne), Cs_R )

        CALL TimersStop_Euler( Timer_Euler_ComputeFromPrimitive )

        EigVals_R(1:nCF) = Eigenvalues_Euler &
                             ( uPF_R(iPF_V1),                             &
                               Cs_R,                                      &
                               G_F  (iNodeX_X1,iGF_Gm_dd_11,iX2,iX3,iX1), &
                               uPF_R(iPF_V1),                             &
                               uPF_R(iPF_V2),                             &
                               uPF_R(iPF_V3),                             &
                               G_F  (iNodeX_X1,iGF_Gm_dd_11,iX2,iX3,iX1), &
                               G_F  (iNodeX_X1,iGF_Gm_dd_22,iX2,iX3,iX1), &
                               G_F  (iNodeX_X1,iGF_Gm_dd_33,iX2,iX3,iX1), &
                               G_F  (iNodeX_X1,iGF_Alpha,   iX2,iX3,iX1), &
                               G_F  (iNodeX_X1,iGF_Beta_1,  iX2,iX3,iX1) )

        Flux_X1_R(1:nCF) &
          = Flux_X1_Euler &
              ( uPF_R(iPF_D ),                           &
                uPF_R(iPF_V1),                           &
                uPF_R(iPF_V2),                           &
                uPF_R(iPF_V3),                           &
                uPF_R(iPF_E ),                           &
                uPF_R(iPF_Ne),                           &
                P_R,                                     &
                G_F(iNodeX_X1,iGF_Gm_dd_11,iX2,iX3,iX1), &
                G_F(iNodeX_X1,iGF_Gm_dd_22,iX2,iX3,iX1), &
                G_F(iNodeX_X1,iGF_Gm_dd_33,iX2,iX3,iX1), &
                G_F(iNodeX_X1,iGF_Alpha,   iX2,iX3,iX1), &
                G_F(iNodeX_X1,iGF_Beta_1,  iX2,iX3,iX1) )

        CALL TimersStart_Euler( Timer_Euler_NumericalFlux )

        AlphaMns &
          = MAX( Zero, &
                 MAXVAL( - EigVals_L(1:nCF) ), &
                 MAXVAL( - EigVals_R(1:nCF) ) )
        AlphaPls &
          = MAX( Zero, &
                 MAXVAL( + EigVals_L(1:nCF) ), &
                 MAXVAL( + EigVals_R(1:nCF) ) )

        AlphaMdl &
          = AlphaMiddle_Euler &
              ( uCF_L(iNodeX_X1,iCF_D     ,iX2,iX3,iX1), &
                uCF_L(iNodeX_X1,iCF_S1    ,iX2,iX3,iX1), &
                uCF_L(iNodeX_X1,iCF_E     ,iX2,iX3,iX1), &
                Flux_X1_L(iCF_D ),                       &
                Flux_X1_L(iCF_S1),                       &
                Flux_X1_L(iCF_E ),                       &
                uCF_R(iNodeX_X1,iCF_D     ,iX2,iX3,iX1), &
                uCF_R(iNodeX_X1,iCF_S1    ,iX2,iX3,iX1), &
                uCF_R(iNodeX_X1,iCF_E     ,iX2,iX3,iX1), &
                Flux_X1_R(iCF_D ),                       &
                Flux_X1_R(iCF_S1),                       &
                Flux_X1_R(iCF_E ),                       &
                G_F(iNodeX_X1,iGF_Gm_dd_11,iX2,iX3,iX1), &
                AlphaPls, AlphaMns,                      &
                G_F(iNodeX_X1,iGF_Alpha   ,iX2,iX3,iX1), &
                G_F(iNodeX_X1,iGF_Beta_1  ,iX2,iX3,iX1) )

        NumericalFlux(iNodeX_X1,1:nCF     ,iX2,iX3,iX1)  &
          = Euler_NumericalFlux_X1 &
              ( uCF_L(iNodeX_X1,1:nCF     ,iX2,iX3,iX1), &
                uCF_R(iNodeX_X1,1:nCF     ,iX2,iX3,iX1), &
                Flux_X1_L(1:nCF),                        &
                Flux_X1_R(1:nCF),                        &
                AlphaPls, AlphaMns, AlphaMdl,            &
                G_F(iNodeX_X1,iGF_Gm_dd_11,iX2,iX3,iX1), &
                uPF_L(iPF_V1),                           &
                uPF_R(iPF_V1),                           &
                P_L,                                     &
                P_R,                                     &
                G_F(iNodeX_X1,iGF_Alpha   ,iX2,iX3,iX1), &
                G_F(iNodeX_X1,iGF_Beta_1  ,iX2,iX3,iX1) )

        CALL TimersStop_Euler( Timer_Euler_NumericalFlux )

        DO iCF = 1, nCF

          NumericalFlux    (iNodeX_X1,iCF       ,iX2,iX3,iX1) &
            = NumericalFlux(iNodeX_X1,iCF       ,iX2,iX3,iX1) &
                * G_F      (iNodeX_X1,iGF_Alpha ,iX2,iX3,iX1) &
                * G_F      (iNodeX_X1,iGF_SqrtGm,iX2,iX3,iX1) &
                * dX2(iX2) * dX3(iX3) * WeightsX_X1(iNodeX_X1)

        END DO

      END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_SurfaceTerm )

    ! --- Surface Contribution ---

    CALL TimersStart_Euler( Timer_Euler_Interpolate )

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nCF_K, nDOFX_X1, + One, LX_X1_Dn, nDOFX_X1, &
             NumericalFlux(1,1,iX_B0(2),iX_B0(3),iX_B0(1)), nDOFX_X1, Zero, &
             dU_X1, nDOFX )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nCF_K, nDOFX_X1, - One, LX_X1_Up, nDOFX_X1, &
             NumericalFlux(1,1,iX_B0(2),iX_B0(3),iX_B0(1)+1), nDOFX_X1, One, &
             dU_X1, nDOFX )

    CALL TimersStop_Euler( Timer_Euler_Interpolate )

    !--------------------
    ! --- Volume Term ---
    !--------------------

    CALL TimersStart_Euler( Timer_Euler_VolumeTerm )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4), &
    !$OMP PRIVATE( uPF_K, Flux_X1_K )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4), &
    !$ACC PRIVATE( uPF_K, Flux_X1_K ) &
    !$ACC PRESENT( Flux_X1_q, uCF_K, G_K, dX2, dX3, WeightsX_q, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( uPF_K, Flux_X1_K, P_K )
#endif
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)

      DO iNodeX = 1, nDOFX

        CALL TimersStart_Euler( Timer_Euler_ComputePrimitive )

        CALL ComputePrimitive_Euler &
               ( uCF_K(iNodeX,iCF_D ,iX2,iX3,iX1),     &
                 uCF_K(iNodeX,iCF_S1,iX2,iX3,iX1),     &
                 uCF_K(iNodeX,iCF_S2,iX2,iX3,iX1),     &
                 uCF_K(iNodeX,iCF_S3,iX2,iX3,iX1),     &
                 uCF_K(iNodeX,iCF_E ,iX2,iX3,iX1),     &
                 uCF_K(iNodeX,iCF_Ne,iX2,iX3,iX1),     &
                 uPF_K(iPF_D ),                        &
                 uPF_K(iPF_V1),                        &
                 uPF_K(iPF_V2),                        &
                 uPF_K(iPF_V3),                        &
                 uPF_K(iPF_E ),                        &
                 uPF_K(iPF_Ne),                        &
                 G_K(iNodeX,iGF_Gm_dd_11,iX2,iX3,iX1), &
                 G_K(iNodeX,iGF_Gm_dd_22,iX2,iX3,iX1), &
                 G_K(iNodeX,iGF_Gm_dd_33,iX2,iX3,iX1) )

        CALL TimersStop_Euler( Timer_Euler_ComputePrimitive )

        CALL TimersStart_Euler( Timer_Euler_ComputeFromPrimitive )

        CALL ComputePressureFromPrimitive &
               ( uPF_K(iPF_D), uPF_K(iPF_E), uPF_K(iPF_Ne), P_K )

        CALL TimersStop_Euler( Timer_Euler_ComputeFromPrimitive )

        Flux_X1_K(1:nCF) &
          = Flux_X1_Euler &
            ( uPF_K(iPF_D ),                        &
              uPF_K(iPF_V1),                        &
              uPF_K(iPF_V2),                        &
              uPF_K(iPF_V3),                        &
              uPF_K(iPF_E ),                        &
              uPF_K(iPF_Ne),                        &
              P_K,                                  &
              G_K(iNodeX,iGF_Gm_dd_11,iX2,iX3,iX1), &
              G_K(iNodeX,iGF_Gm_dd_22,iX2,iX3,iX1), &
              G_K(iNodeX,iGF_Gm_dd_33,iX2,iX3,iX1), &
              G_K(iNodeX,iGF_Alpha,   iX2,iX3,iX1), &
              G_K(iNodeX,iGF_Beta_1,  iX2,iX3,iX1) )

        DO iCF = 1, nCF

          Flux_X1_q(iNodeX,iCF,iX2,iX3,iX1) = Flux_X1_K(iCF)

          Flux_X1_q    (iNodeX,iCF       ,iX2,iX3,iX1) &
            = Flux_X1_q(iNodeX,iCF       ,iX2,iX3,iX1) &
                * G_K  (iNodeX,iGF_Alpha ,iX2,iX3,iX1) &
                * G_K  (iNodeX,iGF_SqrtGm,iX2,iX3,iX1) &
                * dX2(iX2) * dX3(iX3) * WeightsX_q(iNodeX)

        END DO

      END DO

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_VolumeTerm )

    ! --- Contribution from Volume ---

    CALL TimersStart_Euler( Timer_Euler_Interpolate )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nCF_K, nDOFX, One, dLXdX1_q, nDOFX, &
             Flux_X1_q, nDOFX, One, dU_X1, nDOFX )

    CALL TimersStop_Euler( Timer_Euler_Interpolate )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( dU_X1, dU, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$ACC PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iCF = 1, nCF
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX

      dU    (iNodeX,iX1,iX2,iX3,iCF) &
        = dU(iNodeX,iX1,iX2,iX3,iCF) &
            + dU_X1(iNodeX,iCF,iX2,iX3,iX1)

    END DO
    END DO
    END DO
    END DO
    END DO

#ifdef THORNADO_DEBUG_EULER
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( dU_X1 )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( dU_X1 )
#endif
    WRITE(*,'(A20,7I4)')     'MAXLOC(dU_X1)', MAXLOC(dU_X1)
    WRITE(*,'(A20,ES23.15)') 'MAXVAL(dU_X1)', MAXVAL(dU_X1)
#endif

    CALL TimersStart_Euler( Timer_Euler_CopyOut )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: dX2, dX3, iX_B0, iX_E0, iX_B1, iX_E1, &
    !$OMP               G_K, G_F, uCF_K, uCF_L, uCF_R, &
    !$OMP               dU_X1, Flux_X1_q, NumericalFlux )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( dX2, dX3, iX_B0, iX_E0, iX_B1, iX_E1, &
    !$ACC         G_K, G_F, uCF_K, uCF_L, uCF_R, &
    !$ACC         dU_X1, Flux_X1_q, NumericalFlux )
#endif

    CALL TimersStop_Euler( Timer_Euler_CopyOut )

    END ASSOCIATE

  END SUBROUTINE ComputeIncrement_Divergence_X1


  SUBROUTINE ComputeIncrement_Divergence_X2 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nGF), &
      U (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCF)
    REAL(DP), INTENT(inout) :: &
      dU(1:nDOFX,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),1:nCF)

    INTEGER  :: nK(3), nK_X2(3), nCF_K, nCF_F, nGF_F
    INTEGER  :: iNodeX, iNodeX_X2, iX1, iX2, iX3, iCF, iGF
    REAL(DP) :: Alpha, AlphaMns, AlphaPls, AlphaMdl
    REAL(DP) :: uPF_L(nPF), uPF_R(nPF), uPF_K(nPF)
    REAL(DP) :: Flux_X2_L(nCF), Flux_X2_R(nCF), Flux_X2_K(nCF)
    REAL(DP) :: EigVals_L(nCF), EigVals_R(nCF)
    REAL(DP) :: P_L , P_R, P_K
    REAL(DP) :: Cs_L, Cs_R

    REAL(DP) :: G_K          (nDOFX,   nGF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(2)-1:iX_E0(2)+1)
    REAL(DP) :: G_F          (nDOFX_X2,nGF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(2)  :iX_E0(2)+1)

    REAL(DP) :: uCF_K        (nDOFX,   nCF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(2)-1:iX_E0(2)+1)
    REAL(DP) :: uCF_L        (nDOFX_X2,nCF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(2)  :iX_E0(2)+1)
    REAL(DP) :: uCF_R        (nDOFX_X2,nCF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(2)  :iX_E0(2)+1)

    REAL(DP) :: dU_X2        (nDOFX,   nCF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(2)  :iX_E0(2)  )
    REAL(DP) :: Flux_X2_q    (nDOFX,   nCF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(2)  :iX_E0(2)  )
    REAL(DP) :: NumericalFlux(nDOFX_X2,nCF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(2)  :iX_E0(2)+1)

    IF( iX_E0(2) .EQ. iX_B0(2) ) RETURN

    nK    = iX_E0 - iX_B0 + 1      ! Number of Elements per Spatial Dimension
    nK_X2 = nK + [0,1,0]           ! Number of X2 Faces per Spatial Dimension
    nCF_K = nCF * PRODUCT( nK )    ! Number of Fluid Fields in Domain
    nCF_F = nCF * PRODUCT( nK_X2 ) ! Number of Fluid Fields on Interfaces
    nGF_F = nGF * PRODUCT( nK_X2 ) ! Number of Geometry Fields on Interfaces

    ASSOCIATE( dX1 => MeshX(1) % Width, dX3 => MeshX(3) % Width )

    CALL TimersStart_Euler( Timer_Euler_CopyIn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: dX1, dX3, iX_B0, iX_E0, iX_B1, iX_E1 ) &
    !$OMP MAP( alloc: G_K, G_F, uCF_K, uCF_L, uCF_R, &
    !$OMP             dU_X2, Flux_X2_q, NumericalFlux )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( dX1, dX3, iX_B0, iX_E0, iX_B1, iX_E1 ) &
    !$ACC CREATE( G_K, G_F, uCF_K, uCF_L, uCF_R, &
    !$ACC         dU_X2, Flux_X2_q, NumericalFlux )
#endif

    CALL TimersStop_Euler( Timer_Euler_CopyIn )

    ! --- Geometry Fields in Element Nodes ---

    CALL TimersStart_Euler( Timer_Euler_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( G_K, G, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX2 = iX_B0(2) - 1, iX_E0(2) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iGF = 1, nGF
    DO iNodeX = 1, nDOFX

      G_K(iNodeX,iGF,iX1,iX3,iX2) = G(iNodeX,iX1,iX2,iX3,iGF)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_Permute )

    ! --- Interpolate Geometry Fields on Shared Face ---

    CALL TimersStart_Euler( Timer_Euler_Interpolate )

    ! --- Face States (Average of Left and Right States) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nGF_F, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
             G_K(1,1,iX_B0(1),iX_B0(3),iX_B0(2)-1), nDOFX, Zero, &
             G_F(1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX_X2 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nGF_F, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
             G_K(1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX, Half, &
             G_F(1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX_X2 )

    CALL TimersStop_Euler( Timer_Euler_Interpolate )

    !---------------------
    ! --- Surface Term ---
    !---------------------

    CALL TimersStart_Euler( Timer_Euler_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( G_F, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX2 = iX_B0(2), iX_E0(2) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iGF = 1, nGF
    DO iNodeX_X2 = 1, nDOFX_X2

      G_F         (iNodeX_X2,iGF,iX1,iX3,iX2) &
        = MAX( G_F(iNodeX_X2,iGF,iX1,iX3,iX2), SqrtTiny )

    END DO
    END DO
    END DO
    END DO
    END DO

    ! --- Compute metric quantities on interfaces from scale factors ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( G_F, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX2 = iX_B0(2), iX_E0(2) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX_X2 = 1, nDOFX_X2

      G_F         (iNodeX_X2,iGF_Gm_dd_11,iX1,iX3,iX2) &
        = MAX( G_F(iNodeX_X2,iGF_h_1     ,iX1,iX3,iX2)**2, SqrtTiny )

      G_F         (iNodeX_X2,iGF_Gm_dd_22,iX1,iX3,iX2) &
        = MAX( G_F(iNodeX_X2,iGF_h_2     ,iX1,iX3,iX2)**2, SqrtTiny )

      G_F         (iNodeX_X2,iGF_Gm_dd_33,iX1,iX3,iX2) &
        = MAX( G_F(iNodeX_X2,iGF_h_3     ,iX1,iX3,iX2)**2, SqrtTiny )

      G_F        (iNodeX_X2,iGF_SqrtGm,iX1,iX3,iX2) &
        = G_F    (iNodeX_X2,iGF_h_1   ,iX1,iX3,iX2) &
            * G_F(iNodeX_X2,iGF_h_2   ,iX1,iX3,iX2) &
            * G_F(iNodeX_X2,iGF_h_3   ,iX1,iX3,iX2)

    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( uCF_K, U, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX2 = iX_B0(2) - 1, iX_E0(2) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCF = 1, nCF
    DO iNodeX = 1, nDOFX

      uCF_K(iNodeX,iCF,iX1,iX3,iX2) = U(iNodeX,iX1,iX2,iX3,iCF)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_Permute )

    ! --- Interpolate Fluid Fields ---

    CALL TimersStart_Euler( Timer_Euler_Interpolate )

    ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nCF_F, nDOFX, One, LX_X2_Up, nDOFX_X2, &
             uCF_K(1,1,iX_B0(1),iX_B0(3),iX_B0(2)-1), nDOFX, Zero, &
             uCF_L(1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX_X2 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nCF_F, nDOFX, One, LX_X2_Dn, nDOFX_X2, &
             uCF_K(1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX, Zero, &
             uCF_R(1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX_X2 )

    CALL TimersStop_Euler( Timer_Euler_Interpolate )

    ! --- Numerical Flux ---

    CALL TimersStart_Euler( Timer_Euler_SurfaceTerm )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( Flux_X2_L, Flux_X2_R, P_L, P_R, Cs_L, Cs_R, &
    !$OMP          uPF_L, uPF_R, alpha )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( Flux_X2_L, Flux_X2_R, P_L, P_R, Cs_L, Cs_R, &
    !$ACC          uPF_L, uPF_R, alpha ) &
    !$ACC PRESENT( uCF_L, uCF_R, NumericalFlux, G_F, dX1, dX3, WeightsX_X2, &
    !$ACC          iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( Flux_X2_L, Flux_X2_R, P_L, P_R, Cs_L, Cs_R, &
    !$OMP          uPF_L, uPF_R, alpha )
#endif
    DO iX2 = iX_B0(2), iX_E0(2) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX_X2 = 1, nDOFX_X2

        ! --- Left State Primitive, etc. ---

        CALL TimersStart_Euler( Timer_Euler_ComputePrimitive )

        CALL ComputePrimitive_Euler &
               ( uCF_L(iNodeX_X2,iCF_D ,iX1,iX3,iX2),     &
                 uCF_L(iNodeX_X2,iCF_S1,iX1,iX3,iX2),     &
                 uCF_L(iNodeX_X2,iCF_S2,iX1,iX3,iX2),     &
                 uCF_L(iNodeX_X2,iCF_S3,iX1,iX3,iX2),     &
                 uCF_L(iNodeX_X2,iCF_E ,iX1,iX3,iX2),     &
                 uCF_L(iNodeX_X2,iCF_Ne,iX1,iX3,iX2),     &
                 uPF_L(iPF_D ),                           &
                 uPF_L(iPF_V1),                           &
                 uPF_L(iPF_V2),                           &
                 uPF_L(iPF_V3),                           &
                 uPF_L(iPF_E ),                           &
                 uPF_L(iPF_Ne),                           &
                 G_F(iNodeX_X2,iGF_Gm_dd_11,iX1,iX3,iX2), &
                 G_F(iNodeX_X2,iGF_Gm_dd_22,iX1,iX3,iX2), &
                 G_F(iNodeX_X2,iGF_Gm_dd_33,iX1,iX3,iX2) )

        CALL TimersStop_Euler( Timer_Euler_ComputePrimitive )

        CALL TimersStart_Euler( Timer_Euler_ComputeFromPrimitive )

        CALL ComputePressureFromPrimitive &
               ( uPF_L(iPF_D), uPF_L(iPF_E), uPF_L(iPF_Ne), P_L  )

        CALL ComputeSoundSpeedFromPrimitive &
               ( uPF_L(iPF_D), uPF_L(iPF_E), uPF_L(iPF_Ne), Cs_L )

        CALL TimersStop_Euler( Timer_Euler_ComputeFromPrimitive )

        EigVals_L(1:nCF) = Eigenvalues_Euler &
                             ( uPF_L(iPF_V2),                             &
                               Cs_L,                                      &
                               G_F  (iNodeX_X2,iGF_Gm_dd_22,iX1,iX3,iX2), &
                               uPF_L(iPF_V1),                             &
                               uPF_L(iPF_V2),                             &
                               uPF_L(iPF_V3),                             &
                               G_F  (iNodeX_X2,iGF_Gm_dd_11,iX1,iX3,iX2), &
                               G_F  (iNodeX_X2,iGF_Gm_dd_22,iX1,iX3,iX2), &
                               G_F  (iNodeX_X2,iGF_Gm_dd_33,iX1,iX3,iX2), &
                               G_F  (iNodeX_X2,iGF_Alpha,   iX1,iX3,iX2), &
                               G_F  (iNodeX_X2,iGF_Beta_2,  iX1,iX3,iX2) )

        Flux_X2_L(1:nCF) &
          = Flux_X2_Euler &
              ( uPF_L(iPF_D ),                           &
                uPF_L(iPF_V1),                           &
                uPF_L(iPF_V2),                           &
                uPF_L(iPF_V3),                           &
                uPF_L(iPF_E ),                           &
                uPF_L(iPF_Ne),                           &
                P_L,                                     &
                G_F(iNodeX_X2,iGF_Gm_dd_11,iX1,iX3,iX2), &
                G_F(iNodeX_X2,iGF_Gm_dd_22,iX1,iX3,iX2), &
                G_F(iNodeX_X2,iGF_Gm_dd_33,iX1,iX3,iX2), &
                G_F(iNodeX_X2,iGF_Alpha,   iX1,iX3,iX2), &
                G_F(iNodeX_X2,iGF_Beta_2,  iX1,iX3,iX2) )

        ! --- Right State Primitive, etc. ---

        CALL TimersStart_Euler( Timer_Euler_ComputePrimitive )

        CALL ComputePrimitive_Euler &
               ( uCF_R(iNodeX_X2,iCF_D ,iX1,iX3,iX2),     &
                 uCF_R(iNodeX_X2,iCF_S1,iX1,iX3,iX2),     &
                 uCF_R(iNodeX_X2,iCF_S2,iX1,iX3,iX2),     &
                 uCF_R(iNodeX_X2,iCF_S3,iX1,iX3,iX2),     &
                 uCF_R(iNodeX_X2,iCF_E ,iX1,iX3,iX2),     &
                 uCF_R(iNodeX_X2,iCF_Ne,iX1,iX3,iX2),     &
                 uPF_R(iPF_D ),                           &
                 uPF_R(iPF_V1),                           &
                 uPF_R(iPF_V2),                           &
                 uPF_R(iPF_V3),                           &
                 uPF_R(iPF_E ),                           &
                 uPF_R(iPF_Ne),                           &
                 G_F(iNodeX_X2,iGF_Gm_dd_11,iX1,iX3,iX2), &
                 G_F(iNodeX_X2,iGF_Gm_dd_22,iX1,iX3,iX2), &
                 G_F(iNodeX_X2,iGF_Gm_dd_33,iX1,iX3,iX2) )

        CALL TimersStop_Euler( Timer_Euler_ComputePrimitive )

        CALL TimersStart_Euler( Timer_Euler_ComputeFromPrimitive )

        CALL ComputePressureFromPrimitive &
               ( uPF_R(iPF_D), uPF_R(iPF_E), uPF_R(iPF_Ne), P_R  )

        CALL ComputeSoundSpeedFromPrimitive &
               ( uPF_R(iPF_D), uPF_R(iPF_E), uPF_R(iPF_Ne), Cs_R )

        CALL TimersStop_Euler( Timer_Euler_ComputeFromPrimitive )

        EigVals_R(1:nCF) = Eigenvalues_Euler &
                             ( uPF_R(iPF_V2),                             &
                               Cs_R,                                      &
                               G_F  (iNodeX_X2,iGF_Gm_dd_22,iX1,iX3,iX2), &
                               uPF_R(iPF_V1),                             &
                               uPF_R(iPF_V2),                             &
                               uPF_R(iPF_V3),                             &
                               G_F  (iNodeX_X2,iGF_Gm_dd_11,iX1,iX3,iX2), &
                               G_F  (iNodeX_X2,iGF_Gm_dd_22,iX1,iX3,iX2), &
                               G_F  (iNodeX_X2,iGF_Gm_dd_33,iX1,iX3,iX2), &
                               G_F  (iNodeX_X2,iGF_Alpha,   iX1,iX3,iX2), &
                               G_F  (iNodeX_X2,iGF_Beta_2,  iX1,iX3,iX2) )

        Flux_X2_R(1:nCF) &
          = Flux_X2_Euler &
              ( uPF_R(iPF_D ),                           &
                uPF_R(iPF_V1),                           &
                uPF_R(iPF_V2),                           &
                uPF_R(iPF_V3),                           &
                uPF_R(iPF_E ),                           &
                uPF_R(iPF_Ne),                           &
                P_R,                                     &
                G_F(iNodeX_X2,iGF_Gm_dd_11,iX1,iX3,iX2), &
                G_F(iNodeX_X2,iGF_Gm_dd_22,iX1,iX3,iX2), &
                G_F(iNodeX_X2,iGF_Gm_dd_33,iX1,iX3,iX2), &
                G_F(iNodeX_X2,iGF_Alpha,   iX1,iX3,iX2), &
                G_F(iNodeX_X2,iGF_Beta_2,  iX1,iX3,iX2) )

        CALL TimersStart_Euler( Timer_Euler_NumericalFlux )

        AlphaMns &
          = MAX( Zero, &
                 MAXVAL( - EigVals_L(1:nCF) ), &
                 MAXVAL( - EigVals_R(1:nCF) ) )
        AlphaPls &
          = MAX( Zero, &
                 MAXVAL( + EigVals_L(1:nCF) ), &
                 MAXVAL( + EigVals_R(1:nCF) ) )

        AlphaMdl &
          = AlphaMiddle_Euler &
              ( uCF_L(iNodeX_X2,iCF_D     ,iX1,iX3,iX2), &
                uCF_L(iNodeX_X2,iCF_S2    ,iX1,iX3,iX2), &
                uCF_L(iNodeX_X2,iCF_E     ,iX1,iX3,iX2), &
                Flux_X2_L(iCF_D ),                       &
                Flux_X2_L(iCF_S2),                       &
                Flux_X2_L(iCF_E ),                       &
                uCF_R(iNodeX_X2,iCF_D     ,iX1,iX3,iX2), &
                uCF_R(iNodeX_X2,iCF_S2    ,iX1,iX3,iX2), &
                uCF_R(iNodeX_X2,iCF_E     ,iX1,iX3,iX2), &
                Flux_X2_R(iCF_D ),                       &
                Flux_X2_R(iCF_S2),                       &
                Flux_X2_R(iCF_E ),                       &
                G_F(iNodeX_X2,iGF_Gm_dd_22,iX1,iX3,iX2), &
                AlphaPls, AlphaMns,                      &
                G_F(iNodeX_X2,iGF_Alpha   ,iX1,iX3,iX2), &
                G_F(iNodeX_X2,iGF_Beta_2  ,iX1,iX3,iX2) )

        NumericalFlux(iNodeX_X2,1:nCF,iX1,iX3,iX2) &
          = Euler_NumericalFlux_X2 &
              ( uCF_L(iNodeX_X2,1:nCF     ,iX1,iX3,iX2), &
                uCF_R(iNodeX_X2,1:nCF     ,iX1,iX3,iX2), &
                Flux_X2_L(1:nCF),                        &
                Flux_X2_R(1:nCF),                        &
                AlphaPls, AlphaMns, AlphaMdl,            &
                G_F(iNodeX_X2,iGF_Gm_dd_22,iX1,iX3,iX2), &
                uPF_L(iPF_V2),                           &
                uPF_R(iPF_V2),                           &
                P_L,                                     &
                P_R,                                     &
                G_F(iNodeX_X2,iGF_Alpha   ,iX1,iX3,iX2), &
                G_F(iNodeX_X2,iGF_Beta_2  ,iX1,iX3,iX2) )

        CALL TimersStop_Euler( Timer_Euler_NumericalFlux )

        DO iCF = 1, nCF

          NumericalFlux    (iNodeX_X2,iCF       ,iX1,iX3,iX2) &
            = NumericalFlux(iNodeX_X2,iCF       ,iX1,iX3,iX2) &
                * G_F      (iNodeX_X2,iGF_Alpha ,iX1,iX3,iX2) &
                * G_F      (iNodeX_X2,iGF_SqrtGm,iX1,iX3,iX2) &
                * dX1(iX1) * dX3(iX3) * WeightsX_X2(iNodeX_X2)


        END DO

      END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_SurfaceTerm )

    ! --- Surface Contribution ---

    CALL TimersStart_Euler( Timer_Euler_Interpolate )

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nCF_K, nDOFX_X2, + One, LX_X2_Dn, nDOFX_X2, &
             NumericalFlux(1,1,iX_B0(1),iX_B0(3),iX_B0(2)), nDOFX_X2, Zero, &
             dU_X2, nDOFX )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nCF_K, nDOFX_X2, - One, LX_X2_Up, nDOFX_X2, &
             NumericalFlux(1,1,iX_B0(1),iX_B0(3),iX_B0(2)+1), nDOFX_X2, One, &
             dU_X2, nDOFX )

    CALL TimersStop_Euler( Timer_Euler_Interpolate )

    !--------------------
    ! --- Volume Term ---
    !--------------------

    CALL TimersStart_Euler( Timer_Euler_VolumeTerm )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4), &
    !$OMP PRIVATE( uPF_K, Flux_X2_K )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4), &
    !$ACC PRIVATE( uPF_K, Flux_X2_K ) &
    !$ACC PRESENT( Flux_X2_q, uCF_K, G_K, dX1, dX3, WeightsX_q, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( uPF_K, Flux_X2_K, P_K )
#endif
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        CALL TimersStart_Euler( Timer_Euler_ComputePrimitive )

        CALL ComputePrimitive_Euler &
               ( uCF_K(iNodeX,iCF_D ,iX1,iX3,iX2),     &
                 uCF_K(iNodeX,iCF_S1,iX1,iX3,iX2),     &
                 uCF_K(iNodeX,iCF_S2,iX1,iX3,iX2),     &
                 uCF_K(iNodeX,iCF_S3,iX1,iX3,iX2),     &
                 uCF_K(iNodeX,iCF_E ,iX1,iX3,iX2),     &
                 uCF_K(iNodeX,iCF_Ne,iX1,iX3,iX2),     &
                 uPF_K(iPF_D ),                        &
                 uPF_K(iPF_V1),                        &
                 uPF_K(iPF_V2),                        &
                 uPF_K(iPF_V3),                        &
                 uPF_K(iPF_E ),                        &
                 uPF_K(iPF_Ne),                        &
                 G_K(iNodeX,iGF_Gm_dd_11,iX1,iX3,iX2), &
                 G_K(iNodeX,iGF_Gm_dd_22,iX1,iX3,iX2), &
                 G_K(iNodeX,iGF_Gm_dd_33,iX1,iX3,iX2) )

        CALL TimersStop_Euler( Timer_Euler_ComputePrimitive )

        CALL TimersStart_Euler( Timer_Euler_ComputeFromPrimitive )

        CALL ComputePressureFromPrimitive &
               ( uPF_K(iPF_D), uPF_K(iPF_E), uPF_K(iPF_Ne), P_K )

        CALL TimersStop_Euler( Timer_Euler_ComputeFromPrimitive )

        Flux_X2_K(1:nCF) &
          = Flux_X2_Euler &
            ( uPF_K(iPF_D ),                        &
              uPF_K(iPF_V1),                        &
              uPF_K(iPF_V2),                        &
              uPF_K(iPF_V3),                        &
              uPF_K(iPF_E ),                        &
              uPF_K(iPF_Ne),                        &
              P_K,                                  &
              G_K(iNodeX,iGF_Gm_dd_11,iX1,iX3,iX2), &
              G_K(iNodeX,iGF_Gm_dd_22,iX1,iX3,iX2), &
              G_K(iNodeX,iGF_Gm_dd_33,iX1,iX3,iX2), &
              G_K(iNodeX,iGF_Alpha,   iX1,iX3,iX2), &
              G_K(iNodeX,iGF_Beta_2,  iX1,iX3,iX2) )

        DO iCF = 1, nCF

          Flux_X2_q(iNodeX,iCF,iX1,iX3,iX2) = Flux_X2_K(iCF)

          Flux_X2_q    (iNodeX,iCF       ,iX1,iX3,iX2) &
            = Flux_X2_q(iNodeX,iCF       ,iX1,iX3,iX2) &
                * G_K  (iNodeX,iGF_Alpha ,iX1,iX3,iX2) &
                * G_K  (iNodeX,iGF_SqrtGm,iX1,iX3,iX2) &
                * dX1(iX1) * dX3(iX3) * WeightsX_q(iNodeX)

        END DO

      END DO

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_VolumeTerm )

    ! --- Contribution from Volume ---

    CALL TimersStart_Euler( Timer_Euler_Interpolate )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nCF_K, nDOFX, One, dLXdX2_q, nDOFX, &
             Flux_X2_q, nDOFX, One, dU_X2, nDOFX )

    CALL TimersStop_Euler( Timer_Euler_Interpolate )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( dU_X2, dU, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$ACC PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iCF = 1, nCF
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX

      dU    (iNodeX,iX1,iX2,iX3,iCF) &
        = dU(iNodeX,iX1,iX2,iX3,iCF) &
            + dU_X2(iNodeX,iCF,iX1,iX3,iX2)

    END DO
    END DO
    END DO
    END DO
    END DO

#ifdef THORNADO_DEBUG_EULER
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( dU_X2 )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( dU_X2 )
#endif
    WRITE(*,'(A20,7I4)')     'MAXLOC(dU_X2)', MAXLOC(dU_X2)
    WRITE(*,'(A20,ES23.15)') 'MAXVAL(dU_X2)', MAXVAL(dU_X2)
#endif

    CALL TimersStart_Euler( Timer_Euler_CopyOut )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: dX1, dX3, iX_B0, iX_E0, iX_B1, iX_E1, &
    !$OMP               G_K, G_F, uCF_K, uCF_L, uCF_R, &
    !$OMP               dU_X2, Flux_X2_q, NumericalFlux )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( dX1, dX3, iX_B0, iX_E0, iX_B1, iX_E1, &
    !$ACC         G_K, G_F, uCF_K, uCF_L, uCF_R, &
    !$ACC         dU_X2, Flux_X2_q, NumericalFlux )
#endif

    CALL TimersStop_Euler( Timer_Euler_CopyOut )

    END ASSOCIATE

  END SUBROUTINE ComputeIncrement_Divergence_X2


  SUBROUTINE ComputeIncrement_Divergence_X3 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nGF), &
      U (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCF)
    REAL(DP), INTENT(inout) :: &
      dU(1:nDOFX,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),1:nCF)

    INTEGER  :: nK(3), nK_X3(3), nCF_K, nCF_F, nGF_F
    INTEGER  :: iNodeX, iNodeX_X3, iX1, iX2, iX3, iCF, iGF
    REAL(DP) :: Alpha, AlphaMns, AlphaPls, AlphaMdl
    REAL(DP) :: uPF_L(nPF), uPF_R(nPF), uPF_K(nPF)
    REAL(DP) :: Flux_X3_L(nCF), Flux_X3_R(nCF), Flux_X3_K(nCF)
    REAL(DP) :: EigVals_L(nCF), EigVals_R(nCF)
    REAL(DP) :: P_L , P_R, P_K
    REAL(DP) :: Cs_L, Cs_R

    REAL(DP) :: G_K          (nDOFX,   nGF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)-1:iX_E0(3)+1)
    REAL(DP) :: G_F          (nDOFX_X3,nGF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3)+1)

    REAL(DP) :: uCF_K        (nDOFX,   nCF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)-1:iX_E0(3)+1)
    REAL(DP) :: uCF_L        (nDOFX_X3,nCF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3)+1)
    REAL(DP) :: uCF_R        (nDOFX_X3,nCF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3)+1)

    REAL(DP) :: dU_X3        (nDOFX,   nCF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3)  )
    REAL(DP) :: Flux_X3_q    (nDOFX,   nCF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3)  )
    REAL(DP) :: NumericalFlux(nDOFX_X3,nCF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3)+1)

    IF( iX_E0(3) .EQ. iX_B0(3) ) RETURN

    nK    = iX_E0 - iX_B0 + 1      ! Number of Elements per Spatial Dimension
    nK_X3 = nK + [0,0,1]           ! Number of X3 Faces per Spatial Dimension
    nCF_K = nCF * PRODUCT( nK )    ! Number of Fluid Fields in Domain
    nCF_F = nCF * PRODUCT( nK_X3 ) ! Number of Fluid Fields on Interfaces
    nGF_F = nGF * PRODUCT( nK_X3 ) ! Number of Geometry Fields on Interfaces

    ASSOCIATE( dX1 => MeshX(1) % Width, dX2 => MeshX(2) % Width )

    CALL TimersStart_Euler( Timer_Euler_CopyIn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: dX1, dX2, iX_B0, iX_E0, iX_B1, iX_E1 ) &
    !$OMP MAP( alloc: G_K, G_F, uCF_K, uCF_L, uCF_R, &
    !$OMP             dU_X3, Flux_X3_q, NumericalFlux )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( dX1, dX2, iX_B0, iX_E0, iX_B1, iX_E1 ) &
    !$ACC CREATE( G_K, G_F, uCF_K, uCF_L, uCF_R, &
    !$ACC         dU_X3, Flux_X3_q, NumericalFlux )
#endif

    CALL TimersStop_Euler( Timer_Euler_CopyIn )

    ! --- Geometry Fields in Element Nodes ---

    CALL TimersStart_Euler( Timer_Euler_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( G_K, G, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX3 = iX_B0(3) - 1, iX_E0(3) + 1
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iGF = 1, nGF
    DO iNodeX = 1, nDOFX

      G_K(iNodeX,iGF,iX1,iX2,iX3) = G(iNodeX,iX1,iX2,iX3,iGF)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_Permute )

    ! --- Interpolate Geometry Fields on Shared Face ---

    CALL TimersStart_Euler( Timer_Euler_Interpolate )

    ! --- Face States (Average of Left and Right States) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nGF_F, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
             G_K(1,1,iX_B0(1),iX_B0(2),iX_B0(3)-1), nDOFX, Zero, &
             G_F(1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX_X3 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nGF_F, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
             G_K(1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX, Half, &
             G_F(1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX_X3 )

    CALL TimersStop_Euler( Timer_Euler_Interpolate )

    !---------------------
    ! --- Surface Term ---
    !---------------------

    CALL TimersStart_Euler( Timer_Euler_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( G_F, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX3 = iX_B0(3), iX_E0(3) + 1
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iGF = 1, nGF
    DO iNodeX_X3 = 1, nDOFX_X3

      G_F         (iNodeX_X3,iGF,iX1,iX2,iX3) &
        = MAX( G_F(iNodeX_X3,iGF,iX1,iX2,iX3), SqrtTiny )

    END DO
    END DO
    END DO
    END DO
    END DO

    ! --- Compute metric quantities on interfaces from scale factors ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( G_F, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3) + 1
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX_X3 = 1, nDOFX_X3

      G_F         (iNodeX_X3,iGF_Gm_dd_11,iX1,iX2,iX3) &
        = MAX( G_F(iNodeX_X3,iGF_h_1     ,iX1,iX2,iX3)**2, SqrtTiny )

      G_F         (iNodeX_X3,iGF_Gm_dd_22,iX1,iX2,iX3) &
        = MAX( G_F(iNodeX_X3,iGF_h_2     ,iX1,iX2,iX3)**2, SqrtTiny )

      G_F         (iNodeX_X3,iGF_Gm_dd_33,iX1,iX2,iX3) &
        = MAX( G_F(iNodeX_X3,iGF_h_3     ,iX1,iX2,iX3)**2, SqrtTiny )

      G_F        (iNodeX_X3,iGF_SqrtGm,iX1,iX2,iX3) &
        = G_F    (iNodeX_X3,iGF_h_1   ,iX1,iX2,iX3) &
            * G_F(iNodeX_X3,iGF_h_2   ,iX1,iX2,iX3) &
            * G_F(iNodeX_X3,iGF_h_3   ,iX1,iX2,iX3)

    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( uCF_K, U, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX3 = iX_B0(3) - 1, iX_E0(3) + 1
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCF = 1, nCF
    DO iNodeX = 1, nDOFX

      uCF_K(iNodeX,iCF,iX1,iX2,iX3) = U(iNodeX,iX1,iX2,iX3,iCF)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_Permute )

    ! --- Interpolate Fluid Fields ---

    CALL TimersStart_Euler( Timer_Euler_Interpolate )

    ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nCF_F, nDOFX, One, LX_X3_Up, nDOFX_X3, &
             uCF_K(1,1,iX_B0(1),iX_B0(2),iX_B0(3)-1), nDOFX, Zero, &
             uCF_L(1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX_X3 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nCF_F, nDOFX, One, LX_X3_Dn, nDOFX_X3, &
             uCF_K(1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX, Zero, &
             uCF_R(1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX_X3 )

    CALL TimersStop_Euler( Timer_Euler_Interpolate )

    ! --- Numerical Flux ---

    CALL TimersStart_Euler( Timer_Euler_SurfaceTerm )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( Flux_X3_L, Flux_X3_R, P_L, P_R, Cs_L, Cs_R, &
    !$OMP          uPF_L, uPF_R, alpha )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( Flux_X3_L, Flux_X3_R, P_L, P_R, Cs_L, Cs_R, &
    !$ACC          uPF_L, uPF_R, alpha ) &
    !$ACC PRESENT( uCF_L, uCF_R, NumericalFlux, G_F, dX1, dX2, WeightsX_X3, &
    !$ACC          iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( Flux_X3_L, Flux_X3_R, P_L, P_R, Cs_L, Cs_R, &
    !$OMP          uPF_L, uPF_R, alpha )
#endif
    DO iX3 = iX_B0(3), iX_E0(3) + 1
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX_X3 = 1, nDOFX_X3

        ! --- Left State Primitive, etc. ---

        CALL TimersStart_Euler( Timer_Euler_ComputePrimitive )

        CALL ComputePrimitive_Euler &
               ( uCF_L(iNodeX_X3,iCF_D ,iX1,iX2,iX3),     &
                 uCF_L(iNodeX_X3,iCF_S1,iX1,iX2,iX3),     &
                 uCF_L(iNodeX_X3,iCF_S2,iX1,iX2,iX3),     &
                 uCF_L(iNodeX_X3,iCF_S3,iX1,iX2,iX3),     &
                 uCF_L(iNodeX_X3,iCF_E ,iX1,iX2,iX3),     &
                 uCF_L(iNodeX_X3,iCF_Ne,iX1,iX2,iX3),     &
                 uPF_L(iPF_D ),                           &
                 uPF_L(iPF_V1),                           &
                 uPF_L(iPF_V2),                           &
                 uPF_L(iPF_V3),                           &
                 uPF_L(iPF_E ),                           &
                 uPF_L(iPF_Ne),                           &
                 G_F(iNodeX_X3,iGF_Gm_dd_11,iX1,iX2,iX3), &
                 G_F(iNodeX_X3,iGF_Gm_dd_22,iX1,iX2,iX3), &
                 G_F(iNodeX_X3,iGF_Gm_dd_33,iX1,iX2,iX3) )

        CALL TimersStop_Euler( Timer_Euler_ComputePrimitive )

        CALL TimersStart_Euler( Timer_Euler_ComputeFromPrimitive )

        CALL ComputePressureFromPrimitive &
               ( uPF_L(iPF_D), uPF_L(iPF_E), uPF_L(iPF_Ne), P_L  )

        CALL ComputeSoundSpeedFromPrimitive &
               ( uPF_L(iPF_D), uPF_L(iPF_E), uPF_L(iPF_Ne), Cs_L )

        CALL TimersStop_Euler( Timer_Euler_ComputeFromPrimitive )

        EigVals_L(1:nCF) = Eigenvalues_Euler &
                             ( uPF_L(iPF_V3),                             &
                               Cs_L,                                      &
                               G_F  (iNodeX_X3,iGF_Gm_dd_33,iX1,iX2,iX3), &
                               uPF_L(iPF_V1),                             &
                               uPF_L(iPF_V2),                             &
                               uPF_L(iPF_V3),                             &
                               G_F  (iNodeX_X3,iGF_Gm_dd_11,iX1,iX2,iX3), &
                               G_F  (iNodeX_X3,iGF_Gm_dd_22,iX1,iX2,iX3), &
                               G_F  (iNodeX_X3,iGF_Gm_dd_33,iX1,iX2,iX3), &
                               G_F  (iNodeX_X3,iGF_Alpha,   iX1,iX2,iX3), &
                               G_F  (iNodeX_X3,iGF_Beta_3,  iX1,iX2,iX3) )

        Flux_X3_L(1:nCF) &
          = Flux_X3_Euler &
              ( uPF_L(iPF_D ),                           &
                uPF_L(iPF_V1),                           &
                uPF_L(iPF_V2),                           &
                uPF_L(iPF_V3),                           &
                uPF_L(iPF_E ),                           &
                uPF_L(iPF_Ne),                           &
                P_L,                                     &
                G_F(iNodeX_X3,iGF_Gm_dd_11,iX1,iX2,iX3), &
                G_F(iNodeX_X3,iGF_Gm_dd_22,iX1,iX2,iX3), &
                G_F(iNodeX_X3,iGF_Gm_dd_33,iX1,iX2,iX3), &
                G_F(iNodeX_X3,iGF_Alpha,   iX1,iX2,iX3), &
                G_F(iNodeX_X3,iGF_Beta_3,  iX1,iX2,iX3) )

        ! --- Right State Primitive, etc. ---

        CALL TimersStart_Euler( Timer_Euler_ComputePrimitive )

        CALL ComputePrimitive_Euler &
               ( uCF_R(iNodeX_X3,iCF_D ,iX1,iX2,iX3),     &
                 uCF_R(iNodeX_X3,iCF_S1,iX1,iX2,iX3),     &
                 uCF_R(iNodeX_X3,iCF_S2,iX1,iX2,iX3),     &
                 uCF_R(iNodeX_X3,iCF_S3,iX1,iX2,iX3),     &
                 uCF_R(iNodeX_X3,iCF_E ,iX1,iX2,iX3),     &
                 uCF_R(iNodeX_X3,iCF_Ne,iX1,iX2,iX3),     &
                 uPF_R(iPF_D ),                           &
                 uPF_R(iPF_V1),                           &
                 uPF_R(iPF_V2),                           &
                 uPF_R(iPF_V3),                           &
                 uPF_R(iPF_E ),                           &
                 uPF_R(iPF_Ne),                           &
                 G_F(iNodeX_X3,iGF_Gm_dd_11,iX1,iX2,iX3), &
                 G_F(iNodeX_X3,iGF_Gm_dd_22,iX1,iX2,iX3), &
                 G_F(iNodeX_X3,iGF_Gm_dd_33,iX1,iX2,iX3) )

        CALL TimersStop_Euler( Timer_Euler_ComputePrimitive )

        CALL TimersStart_Euler( Timer_Euler_ComputeFromPrimitive )

        CALL ComputePressureFromPrimitive &
               ( uPF_R(iPF_D), uPF_R(iPF_E), uPF_R(iPF_Ne), P_R  )

        CALL ComputeSoundSpeedFromPrimitive &
               ( uPF_R(iPF_D), uPF_R(iPF_E), uPF_R(iPF_Ne), Cs_R )

        CALL TimersStop_Euler( Timer_Euler_ComputeFromPrimitive )

        EigVals_R(1:nCF) = Eigenvalues_Euler &
                             ( uPF_R(iPF_V3),                             &
                               Cs_R,                                      &
                               G_F  (iNodeX_X3,iGF_Gm_dd_33,iX1,iX2,iX3), &
                               uPF_R(iPF_V1),                             &
                               uPF_R(iPF_V2),                             &
                               uPF_R(iPF_V3),                             &
                               G_F  (iNodeX_X3,iGF_Gm_dd_11,iX1,iX2,iX3), &
                               G_F  (iNodeX_X3,iGF_Gm_dd_22,iX1,iX2,iX3), &
                               G_F  (iNodeX_X3,iGF_Gm_dd_33,iX1,iX2,iX3), &
                               G_F  (iNodeX_X3,iGF_Alpha,   iX1,iX2,iX3), &
                               G_F  (iNodeX_X3,iGF_Beta_3,  iX1,iX2,iX3) )

        Flux_X3_R(1:nCF) &
          = Flux_X3_Euler &
              ( uPF_R(iPF_D ),                           &
                uPF_R(iPF_V1),                           &
                uPF_R(iPF_V2),                           &
                uPF_R(iPF_V3),                           &
                uPF_R(iPF_E ),                           &
                uPF_R(iPF_Ne),                           &
                P_R,                                     &
                G_F(iNodeX_X3,iGF_Gm_dd_11,iX1,iX2,iX3), &
                G_F(iNodeX_X3,iGF_Gm_dd_22,iX1,iX2,iX3), &
                G_F(iNodeX_X3,iGF_Gm_dd_33,iX1,iX2,iX3), &
                G_F(iNodeX_X3,iGF_Alpha,   iX1,iX2,iX3), &
                G_F(iNodeX_X3,iGF_Beta_3,  iX1,iX2,iX3) )

        CALL TimersStart_Euler( Timer_Euler_NumericalFlux )

        AlphaMns &
          = MAX( Zero, &
                 MAXVAL( - EigVals_L(1:nCF) ), &
                 MAXVAL( - EigVals_R(1:nCF) ) )
        AlphaPls &
          = MAX( Zero, &
                 MAXVAL( + EigVals_L(1:nCF) ), &
                 MAXVAL( + EigVals_R(1:nCF) ) )

        AlphaMdl &
          = AlphaMiddle_Euler &
              ( uCF_L(iNodeX_X3,iCF_D     ,iX1,iX2,iX3), &
                uCF_L(iNodeX_X3,iCF_S3    ,iX1,iX2,iX3), &
                uCF_L(iNodeX_X3,iCF_E     ,iX1,iX2,iX3), &
                Flux_X3_L(iCF_D ),                       &
                Flux_X3_L(iCF_S3),                       &
                Flux_X3_L(iCF_E ),                       &
                uCF_R(iNodeX_X3,iCF_D     ,iX1,iX2,iX3), &
                uCF_R(iNodeX_X3,iCF_S3    ,iX1,iX2,iX3), &
                uCF_R(iNodeX_X3,iCF_E     ,iX1,iX2,iX3), &
                Flux_X3_R(iCF_D ),                       &
                Flux_X3_R(iCF_S3),                       &
                Flux_X3_R(iCF_E ),                       &
                G_F(iNodeX_X3,iGF_Gm_dd_33,iX1,iX2,iX3), &
                AlphaPls, AlphaMns,                      &
                G_F(iNodeX_X3,iGF_Alpha   ,iX1,iX2,iX3), &
                G_F(iNodeX_X3,iGF_Beta_3  ,iX1,iX2,iX3) )

        NumericalFlux(iNodeX_X3,1:nCF,iX1,iX2,iX3) &
          = Euler_NumericalFlux_X3 &
              ( uCF_L(iNodeX_X3,1:nCF     ,iX1,iX2,iX3), &
                uCF_R(iNodeX_X3,1:nCF     ,iX1,iX2,iX3), &
                Flux_X3_L(1:nCF),                        &
                Flux_X3_R(1:nCF),                        &
                AlphaPls, AlphaMns, AlphaMdl,            &
                G_F(iNodeX_X3,iGF_Gm_dd_33,iX1,iX2,iX3), &
                uPF_L(iPF_V3),                           &
                uPF_R(iPF_V3),                           &
                P_L,                                     &
                P_R,                                     &
                G_F(iNodeX_X3,iGF_Alpha   ,iX1,iX2,iX3), &
                G_F(iNodeX_X3,iGF_Beta_3  ,iX1,iX2,iX3) )

        CALL TimersStop_Euler( Timer_Euler_NumericalFlux )

        DO iCF = 1, nCF

          NumericalFlux    (iNodeX_X3,iCF       ,iX1,iX2,iX3) &
            = NumericalFlux(iNodeX_X3,iCF       ,iX1,iX2,iX3) &
                * G_F      (iNodeX_X3,iGF_Alpha ,iX1,iX2,iX3) &
                * G_F      (iNodeX_X3,iGF_SqrtGm,iX1,iX2,iX3) &
                * dX1(iX1) * dX2(iX2) * WeightsX_X3(iNodeX_X3)


        END DO

      END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_SurfaceTerm )

    ! --- Surface Contribution ---

    CALL TimersStart_Euler( Timer_Euler_Interpolate )

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nCF_K, nDOFX_X3, + One, LX_X3_Dn, nDOFX_X3, &
             NumericalFlux(1,1,iX_B0(1),iX_B0(2),iX_B0(3)), nDOFX_X3, Zero, &
             dU_X3, nDOFX )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nCF_K, nDOFX_X3, - One, LX_X3_Up, nDOFX_X3, &
             NumericalFlux(1,1,iX_B0(1),iX_B0(2),iX_B0(3)+1), nDOFX_X3, One, &
             dU_X3, nDOFX )

    CALL TimersStop_Euler( Timer_Euler_Interpolate )

    !--------------------
    ! --- Volume Term ---
    !--------------------

    CALL TimersStart_Euler( Timer_Euler_VolumeTerm )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4), &
    !$OMP PRIVATE( uPF_K, Flux_X3_K )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4), &
    !$ACC PRIVATE( uPF_K, Flux_X3_K ) &
    !$ACC PRESENT( Flux_X3_q, uCF_K, G_K, dX1, dX2, WeightsX_q, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( uPF_K, Flux_X3_K, P_K )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        CALL TimersStart_Euler( Timer_Euler_ComputePrimitive )

        CALL ComputePrimitive_Euler &
               ( uCF_K(iNodeX,iCF_D ,iX1,iX2,iX3),     &
                 uCF_K(iNodeX,iCF_S1,iX1,iX2,iX3),     &
                 uCF_K(iNodeX,iCF_S2,iX1,iX2,iX3),     &
                 uCF_K(iNodeX,iCF_S3,iX1,iX2,iX3),     &
                 uCF_K(iNodeX,iCF_E ,iX1,iX2,iX3),     &
                 uCF_K(iNodeX,iCF_Ne,iX1,iX2,iX3),     &
                 uPF_K(iPF_D ),                        &
                 uPF_K(iPF_V1),                        &
                 uPF_K(iPF_V2),                        &
                 uPF_K(iPF_V3),                        &
                 uPF_K(iPF_E ),                        &
                 uPF_K(iPF_Ne),                        &
                 G_K(iNodeX,iGF_Gm_dd_11,iX1,iX2,iX3), &
                 G_K(iNodeX,iGF_Gm_dd_22,iX1,iX2,iX3), &
                 G_K(iNodeX,iGF_Gm_dd_33,iX1,iX2,iX3) )

        CALL TimersStop_Euler( Timer_Euler_ComputePrimitive )

        CALL TimersStart_Euler( Timer_Euler_ComputeFromPrimitive )

        CALL ComputePressureFromPrimitive &
               ( uPF_K(iPF_D), uPF_K(iPF_E), uPF_K(iPF_Ne), P_K )

        CALL TimersStop_Euler( Timer_Euler_ComputeFromPrimitive )

        Flux_X3_K(1:nCF) &
          = Flux_X3_Euler &
            ( uPF_K(iPF_D ),                        &
              uPF_K(iPF_V1),                        &
              uPF_K(iPF_V2),                        &
              uPF_K(iPF_V3),                        &
              uPF_K(iPF_E ),                        &
              uPF_K(iPF_Ne),                        &
              P_K,                                  &
              G_K(iNodeX,iGF_Gm_dd_11,iX1,iX2,iX3), &
              G_K(iNodeX,iGF_Gm_dd_22,iX1,iX2,iX3), &
              G_K(iNodeX,iGF_Gm_dd_33,iX1,iX2,iX3), &
              G_K(iNodeX,iGF_Alpha,   iX1,iX2,iX3), &
              G_K(iNodeX,iGF_Beta_3,  iX1,iX2,iX3) )

        DO iCF = 1, nCF

          Flux_X3_q(iNodeX,iCF,iX1,iX2,iX3) = Flux_X3_K(iCF)

          Flux_X3_q    (iNodeX,iCF       ,iX1,iX2,iX3) &
            = Flux_X3_q(iNodeX,iCF       ,iX1,iX2,iX3) &
                * G_K  (iNodeX,iGF_Alpha ,iX1,iX2,iX3) &
                * G_K  (iNodeX,iGF_SqrtGm,iX1,iX2,iX3) &
                * dX1(iX1) * dX2(iX2) * WeightsX_q(iNodeX)

        END DO

      END DO

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_VolumeTerm )

    ! --- Contribution from Volume ---

    CALL TimersStart_Euler( Timer_Euler_Interpolate )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nCF_K, nDOFX, One, dLXdX3_q, nDOFX, &
             Flux_X3_q, nDOFX, One, dU_X3, nDOFX )

    CALL TimersStop_Euler( Timer_Euler_Interpolate )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( dU_X3, dU, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$ACC PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iCF = 1, nCF
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX

      dU    (iNodeX,iX1,iX2,iX3,iCF) &
        = dU(iNodeX,iX1,iX2,iX3,iCF) &
            + dU_X3(iNodeX,iCF,iX1,iX2,iX3)

    END DO
    END DO
    END DO
    END DO
    END DO

#ifdef THORNADO_DEBUG_EULER
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( dU_X3 )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( dU_X3 )
#endif
    WRITE(*,'(A20,7I4)')     'MAXLOC(dU_X3)', MAXLOC(dU_X3)
    WRITE(*,'(A20,ES23.15)') 'MAXVAL(dU_X3)', MAXVAL(dU_X3)
#endif

    CALL TimersStart_Euler( Timer_Euler_CopyOut )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: dX1, dX2, iX_B0, iX_E0, iX_B1, iX_E1, &
    !$OMP               G_K, G_F, uCF_K, uCF_L, uCF_R, &
    !$OMP               dU_X3, Flux_X3_q, NumericalFlux )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( dX1, dX2, iX_B0, iX_E0, iX_B1, iX_E1, &
    !$ACC         G_K, G_F, uCF_K, uCF_L, uCF_R, &
    !$ACC         dU_X3, Flux_X3_q, NumericalFlux )
#endif

    CALL TimersStop_Euler( Timer_Euler_CopyOut )

    END ASSOCIATE

  END SUBROUTINE ComputeIncrement_Divergence_X3


  SUBROUTINE ComputeIncrement_Geometry &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      dU(:,iX_B0(1):,iX_B0(2):,iX_B0(3):,:)


#if defined HYDRO_NONRELATIVISTIC

    CALL ComputeIncrement_Geometry_NonRelativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

#elif defined HYDRO_RELATIVISTIC

    CALL ComputeIncrement_Geometry_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

#endif

  END SUBROUTINE ComputeIncrement_Geometry


  SUBROUTINE ComputeIncrement_Gravity &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

   INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      dU(:,iX_B0(1):,iX_B0(2):,iX_B0(3):,:)

    RETURN
  END SUBROUTINE ComputeIncrement_Gravity


  SUBROUTINE ComputeIncrement_Geometry_NonRelativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      dU(:,iX_B0(1):,iX_B0(2):,iX_B0(3):,:)

    INTEGER  :: iX1, iX2, iX3, iCF, iGF, iNodeX
    REAL(DP) :: dX1, dX2
    REAL(DP) :: P_K(nDOFX)
    REAL(DP) :: dh2dX1(nDOFX), dh3dX1(nDOFX), dh3dX2(nDOFX)
    REAL(DP) :: Stress(nDOFX,3)
    REAL(DP) :: uCF_K(nDOFX,nCF)
    REAL(DP) :: uPF_K(nDOFX,nPF)
    REAL(DP) :: G_K(nDOFX,nGF)
    REAL(DP) :: G_P_X1(nDOFX,nGF), G_N_X1(nDOFX,nGF)
    REAL(DP) :: G_P_X2(nDOFX,nGF), G_N_X2(nDOFX,nGF)
    REAL(DP) :: G_X1_Dn(nDOFX_X1,nGF), G_X1_Up(nDOFX_X1,nGF)
    REAL(DP) :: G_X2_Dn(nDOFX_X2,nGF), G_X2_Up(nDOFX_X2,nGF)

    IF( TRIM( CoordinateSystem ) == 'CARTESIAN' ) RETURN

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      dX2 = MeshX(2) % Width(iX2)
      dX1 = MeshX(1) % Width(iX1)

!      print*,"iX1, iX2, iX3 = ", iX1, iX2, iX3

      DO iCF = 1, nCF

        uCF_K(:,iCF) = U(:,iX1,iX2,iX3,iCF)

      END DO

      DO iGF = 1, nGF

        G_K   (:,iGF) = G(:,iX1,  iX2,iX3,iGF)
        G_P_X1(:,iGF) = G(:,iX1-1,iX2,iX3,iGF)
        G_N_X1(:,iGF) = G(:,iX1+1,iX2,iX3,iGF)

      END DO

      IF( nDimsX .GT. 1 )THEN
        DO iGF = 1, nGF
          G_P_X2(:,iGF) = G(:,iX1,iX2-1,iX3,iGF)
          G_N_X2(:,iGF) = G(:,iX1,iX2+1,iX3,iGF)
        END DO
      END IF

      CALL ComputePrimitive_Euler &
             ( uCF_K(:,iCF_D ), uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), &
               uCF_K(:,iCF_S3), uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
               uPF_K(:,iPF_D ), uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), &
               uPF_K(:,iPF_V3), uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
               G_K(:,iGF_Gm_dd_11), &
               G_K(:,iGF_Gm_dd_22), &
               G_K(:,iGF_Gm_dd_33) )

      CALL ComputePressureFromPrimitive &
             ( uPF_K(:,iPF_D ), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), P_K )

      DO iNodeX = 1, nDOFX

        Stress(iNodeX,1:3) &
          = StressTensor_Diagonal_Euler &
              ( uCF_K(iNodeX,iCF_S1), &
                uCF_K(iNodeX,iCF_S2), &
                uCF_K(iNodeX,iCF_S3), &
                uPF_K(iNodeX,iPF_V1), &
                uPF_K(iNodeX,iPF_V2), &
                uPF_K(iNodeX,iPF_V3), &
                P_K  (iNodeX) )

      END DO

      ! --- Scale Factor Derivatives wrt X1 ---

      ! --- Face States (Average of Left and Right States) ---

      DO iGF = iGF_h_2, iGF_h_3

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
                 G_P_X1(:,iGF), 1, Zero, G_X1_Dn(:,iGF), 1 )
        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
                 G_K   (:,iGF), 1, Half, G_X1_Dn(:,iGF), 1 )

        G_X1_Dn(1:nDOFX_X1,iGF) &
          = MAX( G_X1_Dn(1:nDOFX_X1,iGF), SqrtTiny )

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
                 G_K   (:,iGF), 1, Zero, G_X1_Up(:,iGF), 1 )
        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
                 G_N_X1(:,iGF), 1, Half, G_X1_Up(:,iGF), 1 )

        G_X1_Up(1:nDOFX_X1,iGF) &
          = MAX( G_X1_Up(1:nDOFX_X1,iGF), SqrtTiny )

      END DO

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                  WeightsX_X1(:) * G_X1_Up(:,iGF_h_2), 1, Zero, dh2dX1, 1 )
      CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                  WeightsX_X1(:) * G_X1_Dn(:,iGF_h_2), 1,  One, dh2dX1, 1 )
      CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                  WeightsX_q (:) * G_K    (:,iGF_h_2), 1,  One, dh2dX1, 1 )

      dh2dx1 = dh2dx1 / ( WeightsX_q(:) * dX1 )

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                  WeightsX_X1(:) * G_X1_Up(:,iGF_h_3), 1, Zero, dh3dX1, 1 )
      CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                  WeightsX_X1(:) * G_X1_Dn(:,iGF_h_3), 1,  One, dh3dX1, 1 )
      CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                  WeightsX_q (:) * G_K    (:,iGF_h_3), 1,  One, dh3dX1, 1 )

      dh3dx1 = dh3dx1 / ( WeightsX_q(:) * dX1 )

      dU(:,iX1,iX2,iX3,iCF_S1) &
        = dU(:,iX1,iX2,iX3,iCF_S1) &
            + ( Stress(:,2) * dh2dX1(:) ) / G_K(:,iGF_h_2)  &
            + ( Stress(:,3) * dh3dX1(:) ) / G_K(:,iGF_h_3)

      IF( nDimsX .GT. 1 )THEN

        ! --- Scale Factor Derivatives wrt X2 ---

        ! --- Face States (Average of Left and Right States) ---

        DO iGF = iGF_h_3, iGF_h_3

          CALL DGEMV &
                 ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
                   G_P_X2(:,iGF), 1, Zero, G_X2_Dn(:,iGF), 1 )
          CALL DGEMV &
                 ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
                   G_K   (:,iGF), 1, Half, G_X2_Dn(:,iGF), 1 )

          G_X2_Dn(1:nDOFX_X2,iGF) &
            = MAX( G_X2_Dn(1:nDOFX_X2,iGF), SqrtTiny )

          CALL DGEMV &
                 ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
                   G_K   (:,iGF), 1, Zero, G_X2_Up(:,iGF), 1 )
          CALL DGEMV &
                 ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
                   G_N_X2(:,iGF), 1, Half, G_X2_Up(:,iGF), 1 )

          G_X2_Up(1:nDOFX_X2,iGF) &
            = MAX( G_X2_Up(1:nDOFX_X2,iGF), SqrtTiny )

        END DO

        CALL DGEMV( 'T', nDOFX_X2, nDOFX, + One, LX_X2_Up, nDOFX_X2, &
                  WeightsX_X2(:) * G_X2_Up(:,iGF_h_3), 1, Zero, dh3dX2, 1 )
        CALL DGEMV( 'T', nDOFX_X2, nDOFX, - One, LX_X2_Dn, nDOFX_X2, &
                  WeightsX_X2(:) * G_X2_Dn(:,iGF_h_3), 1,  One, dh3dX2, 1 )
        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX2_q, nDOFX,    &
                  WeightsX_q (:) * G_K    (:,iGF_h_3), 1,  One, dh3dX2, 1 )

        dh3dx2 = dh3dx2 / ( WeightsX_q(:) * dX2 )

        dU(:,iX1,iX2,iX3,iCF_S2) &
          = dU(:,iX1,iX2,iX3,iCF_S2) &
              + ( Stress(:,3) * dh3dX2(:) ) / G_K(:,iGF_h_3)

      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeIncrement_Geometry_NonRelativistic


  SUBROUTINE ComputeIncrement_Geometry_Relativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      dU(:,iX_B0(1):,iX_B0(2):,iX_B0(3):,:)

    ! --- This subroutine currently assumes that the shift-vector
    !     is identically zero which, along with a stationary spacetime,
    !     implies that the extrinsic curvature is identically zero ---

    INTEGER  :: nX(3), nX_X1(3), nX_X2(3), nX_X3(3)
    INTEGER  :: nGF_K, nGF_X1, nGF_X2, nGF_X3
    INTEGER  :: iX1, iX2, iX3, iGF
    INTEGER  :: iNodeX, iNodeX_X1, iNodeX_X2, iNodeX_X3
    REAL(DP) :: dX1, dX2, dX3

    REAL(DP) :: CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne
    REAL(DP) :: PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, P_K, Stress(3)
    REAL(DP) :: h1, h2, h3, Gm11, Gm22, Gm33, Lapse, dU_S1, dU_S2, dU_S3, dU_E

    REAL(DP) :: G_K (nDOFX,nGF ,iX_B0(1):iX_E0(1), &
                                iX_B0(2):iX_E0(2), &
                                iX_B0(3):iX_E0(3))
    REAL(DP) :: G_X1(nDOFX     ,iX_B0(2):iX_E0(2), &
                                iX_B0(3):iX_E0(3), &
                                iX_B0(1)-1:iX_E0(1)+1,nGF)
    REAL(DP) :: G_X2(nDOFX     ,iX_B0(1):iX_E0(1), &
                                iX_B0(3):iX_E0(3), &
                                iX_B0(2)-1:iX_E0(2)+1,nGF)
    REAL(DP) :: G_X3(nDOFX     ,iX_B0(1):iX_E0(1), &
                                iX_B0(2):iX_E0(2), &
                                iX_B0(3)-1:iX_E0(3)+1,nGF)
    REAL(DP) :: h1_X1( nDOFX_X1,iX_B0(2):iX_E0(2), &
                                iX_B0(3):iX_E0(3), &
                                iX_B0(1):iX_E0(1)+1)
    REAL(DP) :: h2_X1( nDOFX_X1,iX_B0(2):iX_E0(2), &
                                iX_B0(3):iX_E0(3), &
                                iX_B0(1):iX_E0(1)+1)
    REAL(DP) :: h3_X1( nDOFX_X1,iX_B0(2):iX_E0(2), &
                                iX_B0(3):iX_E0(3), &
                                iX_B0(1):iX_E0(1)+1)
    REAL(DP) :: h1_X2( nDOFX_X2,iX_B0(1):iX_E0(1), &
                                iX_B0(3):iX_E0(3), &
                                iX_B0(2):iX_E0(2)+1)
    REAL(DP) :: h2_X2( nDOFX_X2,iX_B0(1):iX_E0(1), &
                                iX_B0(3):iX_E0(3), &
                                iX_B0(2):iX_E0(2)+1)
    REAL(DP) :: h3_X2( nDOFX_X2,iX_B0(1):iX_E0(1), &
                                iX_B0(3):iX_E0(3), &
                                iX_B0(2):iX_E0(2)+1)
    REAL(DP) :: h1_X3( nDOFX_X3,iX_B0(1):iX_E0(1), &
                                iX_B0(2):iX_E0(2), &
                                iX_B0(3):iX_E0(3)+1)
    REAL(DP) :: h2_X3( nDOFX_X3,iX_B0(1):iX_E0(1), &
                                iX_B0(2):iX_E0(2), &
                                iX_B0(3):iX_E0(3)+1)
    REAL(DP) :: h3_X3( nDOFX_X3,iX_B0(1):iX_E0(1), &
                                iX_B0(2):iX_E0(2), &
                                iX_B0(3):iX_E0(3)+1)
    REAL(DP) :: dh1dX1(nDOFX   ,iX_B0(2):iX_E0(2), &
                                iX_B0(3):iX_E0(3), &
                                iX_B0(1):iX_E0(1))
    REAL(DP) :: dh2dX1(nDOFX   ,iX_B0(2):iX_E0(2), &
                                iX_B0(3):iX_E0(3), &
                                iX_B0(1):iX_E0(1))
    REAL(DP) :: dh3dX1(nDOFX   ,iX_B0(2):iX_E0(2), &
                                iX_B0(3):iX_E0(3), &
                                iX_B0(1):iX_E0(1))
    REAL(DP) :: dh1dX2(nDOFX   ,iX_B0(1):iX_E0(1), &
                                iX_B0(3):iX_E0(3), &
                                iX_B0(2):iX_E0(2))
    REAL(DP) :: dh2dX2(nDOFX   ,iX_B0(1):iX_E0(1), &
                                iX_B0(3):iX_E0(3), &
                                iX_B0(2):iX_E0(2))
    REAL(DP) :: dh3dX2(nDOFX   ,iX_B0(1):iX_E0(1), &
                                iX_B0(3):iX_E0(3), &
                                iX_B0(2):iX_E0(2))
    REAL(DP) :: dh1dX3(nDOFX   ,iX_B0(1):iX_E0(1), &
                                iX_B0(2):iX_E0(2), &
                                iX_B0(3):iX_E0(3))
    REAL(DP) :: dh2dX3(nDOFX   ,iX_B0(1):iX_E0(1), &
                                iX_B0(2):iX_E0(2), &
                                iX_B0(3):iX_E0(3))
    REAL(DP) :: dh3dX3(nDOFX   ,iX_B0(1):iX_E0(1), &
                                iX_B0(2):iX_E0(2), &
                                iX_B0(3):iX_E0(3))
    REAL(DP) :: a_X1(  nDOFX_X1,iX_B0(2):iX_E0(2), &
                                iX_B0(3):iX_E0(3), &
                                iX_B0(1):iX_E0(1)+1)
    REAL(DP) :: a_X2(  nDOFX_X2,iX_B0(1):iX_E0(1), &
                                iX_B0(3):iX_E0(3), &
                                iX_B0(2):iX_E0(2)+1)
    REAL(DP) :: a_X3(  nDOFX_X3,iX_B0(1):iX_E0(1), &
                                iX_B0(2):iX_E0(2), &
                                iX_B0(3):iX_E0(3)+1)
    REAL(DP) :: dadX1( nDOFX   ,iX_B0(2):iX_E0(2), &
                                iX_B0(3):iX_E0(3), &
                                iX_B0(1):iX_E0(1))
    REAL(DP) :: dadX2( nDOFX   ,iX_B0(1):iX_E0(1), &
                                iX_B0(3):iX_E0(3), &
                                iX_B0(2):iX_E0(2))
    REAL(DP) :: dadX3( nDOFX   ,iX_B0(1):iX_E0(1), &
                                iX_B0(2):iX_E0(2), &
                                iX_B0(3):iX_E0(3))

    nX = iX_E0 - iX_B0 + 1
    nX_X1 = nX + [1,0,0]
    nX_X2 = nX + [0,1,0]
    nX_X3 = nX + [0,0,1]
    nGF_K  = PRODUCT( nX    )
    nGF_X1 = PRODUCT( nX_X1 )
    nGF_X2 = PRODUCT( nX_X2 )
    nGF_X3 = PRODUCT( nX_X3 )

    ASSOCIATE( dX1 => MeshX(1) % Width, &
               dX2 => MeshX(2) % Width, &
               dX3 => MeshX(3) % Width )

    CALL TimersStart_Euler( Timer_Euler_CopyIn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: dX1, dX2, dX3, iX_B0, iX_E0, iX_B1, iX_E1 ) &
    !$OMP MAP( alloc: h1_X1, h2_X1, h3_X1, dh1dX1, dh2dX1, dh3dX1, &
    !$OMP             h1_X2, h2_X2, h3_X2, dh1dX2, dh2dX2, dh3dX2, &
    !$OMP             h1_X3, h2_X3, h3_X3, dh1dX3, dh2dX3, dh3dX3, &
    !$OMP              a_X1,  a_X2,  a_X3, dadX1, dadX2, dadX3, &
    !$OMP             G_X1, G_X2, G_X3, G_K )
#elif defined(THORANDO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( dX1, dX2, dX3, iX_B0, iX_E0, iX_B1, iX_E1 ) &
    !$ACC CREATE( h1_X1, h2_X1, h3_X1, dh1dX1, dh2dX1, dh3dX1, &
    !$ACC         h1_X2, h2_X2, h3_X2, dh1dX2, dh2dX2, dh3dX2, &
    !$ACC         h1_X3, h2_X3, h3_X3, dh1dX3, dh2dX3, dh3dX3, &
    !$ACC          a_X1,  a_X2,  a_X3, dadX1, dadX2, dadX3, &
    !$ACC         G_X1, G_X2, G_X3, G_K )
#endif

    CALL TimersStop_Euler( Timer_Euler_CopyIn )

    !----------------
    ! --- X1 Face ---
    !----------------

    CALL TimersStart_Euler( Timer_Euler_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( G_X1, G, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iGF = 1, nGF
    DO iX1 = iX_B0(1) - 1, iX_E0(1) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iNodeX = 1, nDOFX

      G_X1(iNodeX,iX2,iX3,iX1,iGF) = G(iNodeX,iX1,iX2,iX3,iGF)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_Permute )

    ! --- Interpolate scale factors and lapse function to X1 faces ---

    CALL TimersStart_Euler( Timer_Euler_Interpolate )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nGF_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1,  &
             G_X1 (1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iGF_h_1  ), nDOFX, Zero, &
             h1_X1(1,iX_B0(2),iX_B0(3),iX_B0(1)), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nGF_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1,  &
             G_X1 (1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_1  ), nDOFX, Half, &
             h1_X1(1,iX_B0(2),iX_B0(3),iX_B0(1)), nDOFX_X1 )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nGF_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1,  &
             G_X1 (1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iGF_h_2  ), nDOFX, Zero, &
             h2_X1(1,iX_B0(2),iX_B0(3),iX_B0(1)), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nGF_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1,  &
             G_X1 (1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_2  ), nDOFX, Half, &
             h2_X1(1,iX_B0(2),iX_B0(3),iX_B0(1)), nDOFX_X1 )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nGF_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1,  &
             G_X1 (1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iGF_h_3  ), nDOFX, Zero, &
             h3_X1(1,iX_B0(2),iX_B0(3),iX_B0(1)), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nGF_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1,  &
             G_X1 (1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_3  ), nDOFX, Half, &
             h3_X1(1,iX_B0(2),iX_B0(3),iX_B0(1)), nDOFX_X1 )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nGF_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1,  &
             G_X1 (1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iGF_Alpha), nDOFX, Zero, &
             a_X1 (1,iX_B0(2),iX_B0(3),iX_B0(1)), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nGF_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1,  &
             G_X1 (1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_Alpha), nDOFX, Half, &
             a_X1 (1,iX_B0(2),iX_B0(3),iX_B0(1)), nDOFX_X1 )

    CALL TimersStop_Euler( Timer_Euler_Interpolate )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( h1_X1, h2_X1, h3_X1, a_X1, WeightsX_X1, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX1 = iX_B0(1), iX_E0(1) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iNodeX_X1 = 1, nDOFX_X1

      h1_X1(iNodeX_X1,iX2,iX3,iX1) &
        = WeightsX_X1(iNodeX_X1) * MAX( h1_X1(iNodeX_X1,iX2,iX3,iX1), SqrtTiny )

      h2_X1(iNodeX_X1,iX2,iX3,iX1) &
        = WeightsX_X1(iNodeX_X1) * MAX( h2_X1(iNodeX_X1,iX2,iX3,iX1), SqrtTiny )

      h3_X1(iNodeX_X1,iX2,iX3,iX1) &
        = WeightsX_X1(iNodeX_X1) * MAX( h3_X1(iNodeX_X1,iX2,iX3,iX1), SqrtTiny )

      a_X1 (iNodeX_X1,iX2,iX3,iX1) &
        = WeightsX_X1(iNodeX_X1) * MAX( a_X1 (iNodeX_X1,iX2,iX3,iX1), SqrtTiny )

    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( G_X1, WeightsX_q, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iGF = 1, nGF
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iNodeX = 1, nDOFX

      G_X1(iNodeX,iX2,iX3,iX1,iGF) &
        = WeightsX_q(iNodeX) * G_X1(iNodeX,iX2,iX3,iX1,iGF)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStart_Euler( Timer_Euler_Interpolate )

    ! --- Derivative of scale factors and lapse function wrt X1 ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX_X1, + One, LX_X1_Up, nDOFX_X1, &
             h1_X1 (1,iX_B0(2),iX_B0(3),iX_B0(1)+1), nDOFX_X1, Zero,      &
             dh1dX1(1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX_X1, - One, LX_X1_Dn, nDOFX_X1, &
             h1_X1 (1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1, One,       &
             dh1dX1(1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX,    - One, dLXdX1_q, nDOFX,    &
             G_X1  (1,iX_B0(2),iX_B0(3),iX_B0(1),iGF_h_1  ), nDOFX, One,  &
             dh1dX1(1,iX_B0(2),iX_B0(3),iX_B0(1)), nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX_X1, + One, LX_X1_Up, nDOFX_X1, &
             h2_X1 (1,iX_B0(2),iX_B0(3),iX_B0(1)+1), nDOFX_X1, Zero,      &
             dh2dX1(1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX_X1, - One, LX_X1_Dn, nDOFX_X1, &
             h2_X1 (1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1, One,       &
             dh2dX1(1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX,    - One, dLXdX1_q, nDOFX,    &
             G_X1  (1,iX_B0(2),iX_B0(3),iX_B0(1),iGF_h_2  ), nDOFX, One,  &
             dh2dX1(1,iX_B0(2),iX_B0(3),iX_B0(1)), nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX_X1, + One, LX_X1_Up, nDOFX_X1, &
             h3_X1 (1,iX_B0(2),iX_B0(3),iX_B0(1)+1), nDOFX_X1, Zero,      &
             dh3dX1(1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX_X1, - One, LX_X1_Dn, nDOFX_X1, &
             h3_X1 (1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1, One,       &
             dh3dX1(1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX,    - One, dLXdX1_q, nDOFX,    &
             G_X1  (1,iX_B0(2),iX_B0(3),iX_B0(1),iGF_h_3  ), nDOFX, One,  &
             dh3dX1(1,iX_B0(2),iX_B0(3),iX_B0(1)), nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX_X1, + One, LX_X1_Up, nDOFX_X1, &
             a_X1  (1,iX_B0(2),iX_B0(3),iX_B0(1)+1), nDOFX_X1, Zero,      &
             dadX1 (1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX_X1, - One, LX_X1_Dn, nDOFX_X1, &
             a_X1  (1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1, One,       &
             dadX1 (1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX,    - One, dLXdX1_q, nDOFX,    &
             G_X1  (1,iX_B0(2),iX_B0(3),iX_B0(1),iGF_Alpha), nDOFX, One,  &
             dadX1 (1,iX_B0(2),iX_B0(3),iX_B0(1)), nDOFX )

    CALL TimersStop_Euler( Timer_Euler_Interpolate )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( dh1dX1, dh2dX1, dh3dX1, dadX1, dX1, &
    !$ACC          WeightsX_q, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iNodeX = 1, nDOFX

      dh1dX1(iNodeX,iX2,iX3,iX1) &
        = dh1dX1(iNodeX,iX2,iX3,iX1) / ( WeightsX_q(iNodeX) * dX1(iX1) )

      dh2dX1(iNodeX,iX2,iX3,iX1) &
        = dh2dX1(iNodeX,iX2,iX3,iX1) / ( WeightsX_q(iNodeX) * dX1(iX1) )

      dh3dX1(iNodeX,iX2,iX3,iX1) &
        = dh3dX1(iNodeX,iX2,iX3,iX1) / ( WeightsX_q(iNodeX) * dX1(iX1) )

      dadX1 (iNodeX,iX2,iX3,iX1) &
        = dadX1 (iNodeX,iX2,iX3,iX1) / ( WeightsX_q(iNodeX) * dX1(iX1) )

    END DO
    END DO
    END DO
    END DO

    !----------------
    ! --- X2 Face ---
    !----------------

    CALL TimersStart_Euler( Timer_Euler_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( G_X2, G, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iGF = 1, nGF
    DO iX2 = iX_B0(2) - 1, iX_E0(2)  + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX

      G_X2(iNodeX,iX1,iX3,iX2,iGF) = G(iNodeX,iX1,iX2,iX3,iGF)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_Permute )

    ! --- Interpolate scale factors and lapse function to X2 faces ---

    CALL TimersStart_Euler( Timer_Euler_Interpolate )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nGF_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2,  &
             G_X2 (1,iX_B0(1),iX_B0(3),iX_B0(2)-1,iGF_h_1  ), nDOFX, Zero, &
             h1_X2(1,iX_B0(1),iX_B0(3),iX_B0(2)), nDOFX_X2 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nGF_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2,  &
             G_X2 (1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_1  ), nDOFX, Half, &
             h1_X2(1,iX_B0(1),iX_B0(3),iX_B0(2)), nDOFX_X2 )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nGF_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2,  &
             G_X2 (1,iX_B0(1),iX_B0(3),iX_B0(2)-1,iGF_h_2  ), nDOFX, Zero, &
             h2_X2(1,iX_B0(1),iX_B0(3),iX_B0(2)), nDOFX_X2 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nGF_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2,  &
             G_X2 (1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_2  ), nDOFX, Half, &
             h2_X2(1,iX_B0(1),iX_B0(3),iX_B0(2)), nDOFX_X2 )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nGF_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2,  &
             G_X2 (1,iX_B0(1),iX_B0(3),iX_B0(2)-1,iGF_h_3  ), nDOFX, Zero, &
             h3_X2(1,iX_B0(1),iX_B0(3),iX_B0(2)), nDOFX_X2 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nGF_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2,  &
             G_X2 (1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_3  ), nDOFX, Half, &
             h3_X2(1,iX_B0(1),iX_B0(3),iX_B0(2)), nDOFX_X2 )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nGF_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2,  &
             G_X2 (1,iX_B0(1),iX_B0(3),iX_B0(2)-1,iGF_Alpha), nDOFX, Zero, &
             a_X2 (1,iX_B0(1),iX_B0(3),iX_B0(2)), nDOFX_X2 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nGF_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2,  &
             G_X2 (1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_Alpha), nDOFX, Half, &
             a_X2 (1,iX_B0(1),iX_B0(3),iX_B0(2)), nDOFX_X2 )

    CALL TimersStop_Euler( Timer_Euler_Interpolate )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( h1_X2, h2_X2, h3_X2, a_X2, WeightsX_X2, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX2 = iX_B0(2), iX_E0(2) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX_X2 = 1, nDOFX_X2

      h1_X2(iNodeX_X2,iX1,iX3,iX2) &
        = WeightsX_X2(iNodeX_X2) * MAX( h1_X2(iNodeX_X2,iX1,iX3,iX2), SqrtTiny )

      h2_X2(iNodeX_X2,iX1,iX3,iX2) &
        = WeightsX_X2(iNodeX_X2) * MAX( h2_X2(iNodeX_X2,iX1,iX3,iX2), SqrtTiny )

      h3_X2(iNodeX_X2,iX1,iX3,iX2) &
        = WeightsX_X2(iNodeX_X2) * MAX( h3_X2(iNodeX_X2,iX1,iX3,iX2), SqrtTiny )

      a_X2 (iNodeX_X2,iX1,iX3,iX2) &
        = WeightsX_X2(iNodeX_X2) * MAX( a_X2 (iNodeX_X2,iX1,iX3,iX2), SqrtTiny )

    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( G_X2, WeightsX_q, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iGF = 1, nGF
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX

      G_X2(iNodeX,iX1,iX3,iX2,iGF) &
        = WeightsX_q(iNodeX) * G_X2(iNodeX,iX1,iX3,iX2,iGF)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStart_Euler( Timer_Euler_Interpolate )

    ! --- Derivatives of scale factors and lapse function wrt X2 ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX_X2, + One, LX_X2_Up, nDOFX_X2, &
             h1_X2 (1,iX_B0(1),iX_B0(3),iX_B0(2)+1), nDOFX_X2, Zero,      &
             dh1dX2(1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX_X2, - One, LX_X2_Dn, nDOFX_X2, &
             h1_X2 (1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX_X2, One,       &
             dh1dX2(1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX,    - One, dLXdX2_q, nDOFX,    &
             G_X2  (1,iX_B0(1),iX_B0(3),iX_B0(2),iGF_h_1  ), nDOFX, One,  &
             dh1dX2(1,iX_B0(1),iX_B0(3),iX_B0(2)), nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX_X2, + One, LX_X2_Up, nDOFX_X2, &
             h2_X2 (1,iX_B0(1),iX_B0(3),iX_B0(2)+1), nDOFX_X2, Zero,      &
             dh2dX2(1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX_X2, - One, LX_X2_Dn, nDOFX_X2, &
             h2_X2 (1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX_X2, One,       &
             dh2dX2(1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX,    - One, dLXdX2_q, nDOFX,    &
             G_X2  (1,iX_B0(1),iX_B0(3),iX_B0(2),iGF_h_2  ), nDOFX, One,  &
             dh2dX2(1,iX_B0(1),iX_B0(3),iX_B0(2)), nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX_X2, + One, LX_X2_Up, nDOFX_X2, &
             h3_X2 (1,iX_B0(1),iX_B0(3),iX_B0(2)+1), nDOFX_X2, Zero,      &
             dh3dX2(1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX_X2, - One, LX_X2_Dn, nDOFX_X2, &
             h3_X2 (1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX_X2, One,       &
             dh3dX2(1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX,    - One, dLXdX2_q, nDOFX,    &
             G_X2  (1,iX_B0(1),iX_B0(3),iX_B0(2),iGF_h_3  ), nDOFX, One,  &
             dh3dX2(1,iX_B0(1),iX_B0(3),iX_B0(2)), nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX_X2, + One, LX_X2_Up, nDOFX_X2, &
             a_X2  (1,iX_B0(1),iX_B0(3),iX_B0(2)+1), nDOFX_X2, Zero,      &
             dadX2 (1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX_X2, - One, LX_X2_Dn, nDOFX_X2, &
             a_X2  (1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX_X2, One,       &
             dadX2 (1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX,    - One, dLXdX2_q, nDOFX,    &
             G_X2  (1,iX_B0(1),iX_B0(3),iX_B0(2),iGF_Alpha), nDOFX, One,  &
             dadX2 (1,iX_B0(1),iX_B0(3),iX_B0(2)), nDOFX )

    CALL TimersStop_Euler( Timer_Euler_Interpolate )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( dh1dX2, dh2dX2, dh3dX2, dadX2, dX2, &
    !$ACC          WeightsX_q, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX

      dh1dX2(iNodeX,iX1,iX3,iX2) &
        = dh1dX2(iNodeX,iX1,iX3,iX2) / ( WeightsX_q(iNodeX) * dX2(iX2) )

      dh2dX2(iNodeX,iX1,iX3,iX2) &
        = dh2dX2(iNodeX,iX1,iX3,iX2) / ( WeightsX_q(iNodeX) * dX2(iX2) )

      dh3dX2(iNodeX,iX1,iX3,iX2) &
        = dh3dX2(iNodeX,iX1,iX3,iX2) / ( WeightsX_q(iNodeX) * dX2(iX2) )

      dadX2 (iNodeX,iX1,iX3,iX2) &
        = dadX2 (iNodeX,iX1,iX3,iX2) / ( WeightsX_q(iNodeX) * dX2(iX2) )

    END DO
    END DO
    END DO
    END DO

    !----------------
    ! --- X3 Face ---
    !----------------

    CALL TimersStart_Euler( Timer_Euler_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( G_X3, G, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iGF = 1, nGF
    DO iX3 = iX_B0(3) - 1, iX_E0(3) + 1
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX

      G_X3(iNodeX,iX1,iX2,iX3,iGF) = G(iNodeX,iX1,iX2,iX3,iGF)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_Permute )

    ! --- Interpolate scale factors and lapse function to X3 faces ---

    CALL TimersStart_Euler( Timer_Euler_Interpolate )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nGF_X3, nDOFX, One,  LX_X3_Up, nDOFX_X3,  &
             G_X3 (1,iX_B0(1),iX_B0(2),iX_B0(3)-1,iGF_h_1  ), nDOFX, Zero, &
             h1_X3(1,iX_B0(1),iX_B0(2),iX_B0(3)), nDOFX_X3 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nGF_X3, nDOFX, Half, LX_X3_Dn, nDOFX_X3,  &
             G_X3 (1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iGF_h_1  ), nDOFX, Half, &
             h1_X3(1,iX_B0(1),iX_B0(2),iX_B0(3)), nDOFX_X3 )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nGF_X3, nDOFX, One,  LX_X3_Up, nDOFX_X3,  &
             G_X3 (1,iX_B0(1),iX_B0(2),iX_B0(3)-1,iGF_h_2  ), nDOFX, Zero, &
             h2_X3(1,iX_B0(1),iX_B0(2),iX_B0(3)), nDOFX_X3 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nGF_X3, nDOFX, Half, LX_X3_Dn, nDOFX_X3,  &
             G_X3 (1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iGF_h_2  ), nDOFX, Half, &
             h2_X3(1,iX_B0(1),iX_B0(2),iX_B0(3)), nDOFX_X3 )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nGF_X3, nDOFX, One,  LX_X3_Up, nDOFX_X3,  &
             G_X3 (1,iX_B0(1),iX_B0(2),iX_B0(3)-1,iGF_h_3  ), nDOFX, Zero, &
             h3_X3(1,iX_B0(1),iX_B0(2),iX_B0(3)), nDOFX_X3 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nGF_X3, nDOFX, Half, LX_X3_Dn, nDOFX_X3,  &
             G_X3 (1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iGF_h_3  ), nDOFX, Half, &
             h3_X3(1,iX_B0(1),iX_B0(2),iX_B0(3)), nDOFX_X3 )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nGF_X3, nDOFX, One,  LX_X3_Up, nDOFX_X3,  &
             G_X3 (1,iX_B0(1),iX_B0(2),iX_B0(3)-1,iGF_Alpha), nDOFX, Zero, &
             a_X3 (1,iX_B0(1),iX_B0(2),iX_B0(3)), nDOFX_X3 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nGF_X3, nDOFX, Half, LX_X3_Dn, nDOFX_X3,  &
             G_X3 (1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iGF_Alpha), nDOFX, Half, &
             a_X3 (1,iX_B0(1),iX_B0(2),iX_B0(3)), nDOFX_X3 )

    CALL TimersStop_Euler( Timer_Euler_Interpolate )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( h1_X3, h2_X3, h3_X3, a_X3, WeightsX_X3, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3) + 1
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX_X3 = 1, nDOFX_X3

      h1_X3(iNodeX_X3,iX1,iX2,iX3) &
        = WeightsX_X3(iNodeX_X3) * MAX( h1_X3(iNodeX_X3,iX1,iX2,iX3), SqrtTiny )

      h2_X3(iNodeX_X3,iX1,iX2,iX3) &
        = WeightsX_X3(iNodeX_X3) * MAX( h2_X3(iNodeX_X3,iX1,iX2,iX3), SqrtTiny )

      h3_X3(iNodeX_X3,iX1,iX2,iX3) &
        = WeightsX_X3(iNodeX_X3) * MAX( h3_X3(iNodeX_X3,iX1,iX2,iX3), SqrtTiny )

      a_X3 (iNodeX_X3,iX1,iX2,iX3) &
        = WeightsX_X3(iNodeX_X3) * MAX( a_X3 (iNodeX_X3,iX1,iX2,iX3), SqrtTiny )

    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( G_X3, WeightsX_q, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iGF = 1, nGF
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX

      G_X3(iNodeX,iX1,iX2,iX3,iGF) &
        = WeightsX_q(iNodeX) * G_X3(iNodeX,iX1,iX2,iX3,iGF)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStart_Euler( Timer_Euler_Interpolate )

    ! --- Derivatives of scale factors and lapse function wrt X3 ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX_X3, + One, LX_X3_Up, nDOFX_X3, &
             h1_X3 (1,iX_B0(1),iX_B0(2),iX_B0(3)+1), nDOFX_X3, Zero,      &
             dh1dX3(1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX_X3, - One, LX_X3_Dn, nDOFX_X3, &
             h1_X3 (1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX_X3, One,       &
             dh1dX3(1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX,    - One, dLXdX3_q, nDOFX,    &
             G_X3  (1,iX_B0(1),iX_B0(2),iX_B0(3),iGF_h_1  ), nDOFX, One,  &
             dh1dX3(1,iX_B0(1),iX_B0(2),iX_B0(3)), nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX_X3, + One, LX_X3_Up, nDOFX_X3, &
             h2_X3 (1,iX_B0(1),iX_B0(2),iX_B0(3)+1), nDOFX_X3, Zero,      &
             dh2dX3(1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX_X3, - One, LX_X3_Dn, nDOFX_X3, &
             h2_X3 (1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX_X3, One,       &
             dh2dX3(1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX,    - One, dLXdX3_q, nDOFX,    &
             G_X3  (1,iX_B0(1),iX_B0(2),iX_B0(3),iGF_h_2  ), nDOFX, One,  &
             dh2dX3(1,iX_B0(1),iX_B0(2),iX_B0(3)), nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX_X3, + One, LX_X3_Up, nDOFX_X3, &
             h3_X3 (1,iX_B0(1),iX_B0(2),iX_B0(3)+1), nDOFX_X3, Zero,      &
             dh3dX3(1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX_X3, - One, LX_X3_Dn, nDOFX_X3, &
             h3_X3 (1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX_X3, One,       &
             dh3dX3(1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX,    - One, dLXdX3_q, nDOFX,    &
             G_X3  (1,iX_B0(1),iX_B0(2),iX_B0(3),iGF_h_3  ), nDOFX, One,  &
             dh3dX3(1,iX_B0(1),iX_B0(2),iX_B0(3)), nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX_X3, + One, LX_X3_Up, nDOFX_X3, &
             a_X3  (1,iX_B0(1),iX_B0(2),iX_B0(3)+1), nDOFX_X3, Zero,      &
             dadX3 (1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX_X3, - One, LX_X3_Dn, nDOFX_X3, &
             a_X3  (1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX_X3, One,       &
             dadX3 (1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX,    - One, dLXdX3_q, nDOFX,    &
             G_X3  (1,iX_B0(1),iX_B0(2),iX_B0(3),iGF_Alpha), nDOFX, One,  &
             dadX3 (1,iX_B0(1),iX_B0(2),iX_B0(3)), nDOFX )

    CALL TimersStop_Euler( Timer_Euler_Interpolate )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( dh1dX3, dh2dX3, dh3dX3, dadX3, dX3, &
    !$ACC          WeightsX_q, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX

      dh1dX3(iNodeX,iX1,iX2,iX3) &
        = dh1dX3(iNodeX,iX1,iX2,iX3) / ( WeightsX_q(iNodeX) * dX3(iX3) )

      dh2dX3(iNodeX,iX1,iX2,iX3) &
        = dh2dX3(iNodeX,iX1,iX2,iX3) / ( WeightsX_q(iNodeX) * dX3(iX3) )

      dh3dX3(iNodeX,iX1,iX2,iX3) &
        = dh3dX3(iNodeX,iX1,iX2,iX3) / ( WeightsX_q(iNodeX) * dX3(iX3) )

      dadX3 (iNodeX,iX1,iX2,iX3) &
        = dadX3 (iNodeX,iX1,iX2,iX3) / ( WeightsX_q(iNodeX) * dX3(iX3) )

    END DO
    END DO
    END DO
    END DO


    CALL TimersStart_Euler( Timer_Euler_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( G_K, G, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iGF = 1, nGF
    DO iNodeX = 1, nDOFX

      G_K(iNodeX,iGF,iX1,iX2,iX3) = G(iNodeX,iX1,iX2,iX3,iGF)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_Permute )

    ! --- Add to Increments ---

    CALL TimersStart_Euler( Timer_Euler_Increment )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( h1, h2, h3, Gm11, Gm22, Gm33, Lapse, &
    !$OMP          dU_S1, dU_S2, dU_S3, dU_E, &
    !$OMP          CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
    !$OMP          PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, P_K, Stress )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( h1, h2, h3, Gm11, Gm22, Gm33, Lapse, &
    !$ACC          dU_S1, dU_S2, dU_S3, dU_E, &
    !$ACC          CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
    !$ACC          PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, P_K, Stress ) &
    !$ACC PRESENT( dh1dX1, dh2dX1, dh3dX1, dadX1, &
    !$ACC          dh1dX2, dh2dX2, dh3dX2, dadX2, &
    !$ACC          dh1dX3, dh2dX3, dh3dX3, dadX3, G_K, U, dU, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
      DO iNodeX = 1, nDOFX

        h1    = G_K(iNodeX,iGF_h_1     ,iX1,iX2,iX3)
        h2    = G_K(iNodeX,iGF_h_2     ,iX1,iX2,iX3)
        h3    = G_K(iNodeX,iGF_h_3     ,iX1,iX2,iX3)
        Gm11  = G_K(iNodeX,iGF_Gm_dd_11,iX1,iX2,iX3)
        Gm22  = G_K(iNodeX,iGF_Gm_dd_22,iX1,iX2,iX3)
        Gm33  = G_K(iNodeX,iGF_Gm_dd_33,iX1,iX2,iX3)
        Lapse = G_K(iNodeX,iGF_Alpha   ,iX1,iX2,iX3)

        CF_D  = U(iNodeX,iX1,iX2,iX3,iCF_D )
        CF_S1 = U(iNodeX,iX1,iX2,iX3,iCF_S1)
        CF_S2 = U(iNodeX,iX1,iX2,iX3,iCF_S2)
        CF_S3 = U(iNodeX,iX1,iX2,iX3,iCF_S3)
        CF_E  = U(iNodeX,iX1,iX2,iX3,iCF_E )
        CF_Ne = U(iNodeX,iX1,iX2,iX3,iCF_Ne)

        CALL ComputePrimitive_Euler &
               ( CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
                 PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
                 Gm11, Gm22, Gm33 )

        CALL ComputePressureFromPrimitive &
               ( PF_D, PF_E, PF_Ne, P_K )

        Stress(1:3) &
          = StressTensor_Diagonal_Euler &
              ( CF_S1, CF_S2, CF_S3, &
                PF_V1, PF_V2, PF_V3, P_K )

        ! --- Compute S_1 Increment ---
        dU_S1 = Lapse * ( ( Stress(1) * dh1dX1(iNodeX,iX2,iX3,iX1) ) / h1   &
                        + ( Stress(2) * dh2dX1(iNodeX,iX2,iX3,iX1) ) / h2   &
                        + ( Stress(3) * dh3dX1(iNodeX,iX2,iX3,iX1) ) / h3 ) &
                  - ( CF_D + CF_E ) * dadx1(iNodeX,iX2,iX3,iX1)

        ! --- Compute S_2 Increment ---
        dU_S2 = Lapse * ( ( Stress(1) * dh1dX2(iNodeX,iX1,iX3,iX2) ) / h1   &
                        + ( Stress(2) * dh2dX2(iNodeX,iX1,iX3,iX2) ) / h2   &
                        + ( Stress(3) * dh3dX2(iNodeX,iX1,iX3,iX2) ) / h3 ) &
                  - ( CF_D + CF_E ) * dadx2(iNodeX,iX1,iX3,iX2)

        ! --- Compute S_3 Increment ---
        dU_S3 = Lapse * ( ( Stress(1) * dh1dX3(iNodeX,iX1,iX2,iX3) ) / h1   &
                        + ( Stress(2) * dh2dX3(iNodeX,iX1,iX2,iX3) ) / h2   &
                        + ( Stress(3) * dh3dX3(iNodeX,iX1,iX2,iX3) ) / h3 ) &
                  - ( CF_D + CF_E ) * dadx3(iNodeX,iX1,iX2,iX3)

        ! --- Compute Energy Increment (missing extrinsic curvature term) ---
        dU_E = - (  CF_S1 / Gm11 * dadx1(iNodeX,iX2,iX3,iX1) &
                  + CF_S2 / Gm22 * dadx2(iNodeX,iX1,iX3,iX2) &
                  + CF_S3 / Gm33 * dadx3(iNodeX,iX1,iX2,iX3) )

        dU(iNodeX,iX1,iX2,iX3,iCF_S1) = dU(iNodeX,iX1,iX2,iX3,iCF_S1) + dU_S1
        dU(iNodeX,iX1,iX2,iX3,iCF_S2) = dU(iNodeX,iX1,iX2,iX3,iCF_S2) + dU_S2
        dU(iNodeX,iX1,iX2,iX3,iCF_S3) = dU(iNodeX,iX1,iX2,iX3,iCF_S3) + dU_S3
        dU(iNodeX,iX1,iX2,iX3,iCF_E ) = dU(iNodeX,iX1,iX2,iX3,iCF_E ) + dU_E

      END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_Increment )

#ifdef THORNADO_DEBUG_EULER
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( dh2dX1, dh3dX1, dh3dX2 )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( dh2dX1, dh3dX1, dh3dX2 )
#endif
    WRITE(*,'(A,4I4,ES23.15)') &
      'MINLOC(dh2dX1), MINVAL(dh2dX1)', MINLOC(dh2dX1), MINVAL(dh2dX1)
    WRITE(*,'(A,4I4,ES23.15)') &
      'MAXLOC(dh2dX1), MAXVAL(dh2dX1)', MAXLOC(dh2dX1), MAXVAL(dh2dX1)
    WRITE(*,'(A,4I4,ES23.15)') &
      'MINLOC(dh3dX1), MINVAL(dh3dX1)', MINLOC(dh3dX1), MINVAL(dh3dX1)
    WRITE(*,'(A,4I4,ES23.15)') &
      'MAXLOC(dh3dX1), MAXVAL(dh3dX1)', MAXLOC(dh3dX1), MAXVAL(dh3dX1)
    WRITE(*,'(A,4I4,ES23.15)') &
      'MINLOC(dh3dX2), MINVAL(dh3dX2)', MINLOC(dh3dX2), MINVAL(dh3dX2)
    WRITE(*,'(A,4I4,ES23.15)') &
      'MAXLOC(dh3dX2), MAXVAL(dh3dX2)', MAXLOC(dh3dX2), MAXVAL(dh3dX2)
#endif

    CALL TimersStart_Euler( Timer_Euler_CopyOut )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: dX1, dX2, dX3, iX_B0, iX_E0, iX_B1, iX_E1, &
    !$OMP               h1_X1, h2_X1, h3_X1, dh1dX1, dh2dX1, dh3dX1, &
    !$OMP               h1_X2, h2_X2, h3_X2, dh1dX2, dh2dX2, dh3dX2, &
    !$OMP               h1_X3, h2_X3, h3_X3, dh1dX3, dh2dX3, dh3dX3, &
    !$OMP                a_X1,  a_X2,  a_X3, dadX1, dadX2, dadX3, &
    !$OMP               G_X1, G_X2, G_X3, G_K )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( dX1, dX2, dX3, iX_B0, iX_E0, iX_B1, iX_E1, &
    !$ACC               h1_X1, h2_X1, h3_X1, dh1dX1, dh2dX1, dh3dX1, &
    !$ACC               h1_X2, h2_X2, h3_X2, dh1dX2, dh2dX2, dh3dX2, &
    !$ACC               h1_X3, h2_X3, h3_X3, dh1dX3, dh2dX3, dh3dX3, &
    !$ACC                a_X1,  a_X2,  a_X3, dadX1, dadX2, dadX3, &
    !$ACC               G_X1, G_X2, G_X3, G_K )
#endif

    CALL TimersStop_Euler( Timer_Euler_CopyOut )


    END ASSOCIATE

  END SUBROUTINE ComputeIncrement_Geometry_Relativistic


END MODULE Euler_dgDiscretizationModule
