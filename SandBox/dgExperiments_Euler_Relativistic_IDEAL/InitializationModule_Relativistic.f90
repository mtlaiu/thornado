MODULE InitializationModule_Relativistic

  USE KindModule,                         ONLY: &
    DP,       &
    SqrtTiny, &
    Zero,     &
    Half,     &
    One,      &
    Two,      &
    Three,    &
    Four,     &
    Pi,       &
    TwoPi,    &
    FourPi
  USE ProgramHeaderModule,                ONLY: &
    ProgramName, &
    nNodesX,     &
    nDimsX,      &
    nDOFX,       &
    swX,         &
    iX_B0,       &
    iX_B1,       &
    iX_E0,       &
    iX_E1
  USE ReferenceElementModuleX,            ONLY: &
    NodeNumberTableX, &
    WeightsX_q
  USE MeshModule,                         ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule,               ONLY: &
    nGF,          &
    uGF,          &
    iGF_Phi_N,    &
    iGF_h_1,      &
    iGF_h_2,      &
    iGF_h_3,      &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Beta_1,   &
    iGF_Beta_2,   &
    iGF_Beta_3,   &
    iGF_SqrtGm,   &
    iGF_Alpha,    &
    iGF_Psi
  USE FluidFieldsModule,                  ONLY: &
    nPF,    &
    uPF,    &
    iPF_D,  &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E,  &
    iPF_Ne, &
    uCF,    &
    iCF_D,  &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E,  &
    iCF_Ne, &
    uAF,    &
    iAF_P
  USE EquationOfStateModule_IDEAL,        ONLY: &
    Gamma_IDEAL, &
    ComputePressureFromPrimitive_IDEAL
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_Euler_Relativistic
  USE PhysicalConstantsModule,            ONlY: &
    SpeedOfLightCGS
  USE UnitsModule,                        ONLY: &
    GravitationalConstant, &
    SpeedOfLight,          &
    Kilometer,             &
    SolarMass,             &
    Gram,                  &
    Centimeter,            &
    Erg,                   &
    Millisecond,           &
    Second,                &
    PlanckConstant,        &
    AtomicMassUnit
  USE UtilitiesModule,                    ONLY: &
    NodeNumberX, &
    Locate,      &
    Interpolate1D_Linear
  USE QuadratureModule,                   ONLY: &
    GetQuadrature
  USE PolynomialBasisModule_Lagrange,     ONLY: &
    LagrangeP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_Relativistic

  REAL(DP), PARAMETER :: &
    PolytropicConstant_TOV &
!!$      = 1.0_DP / 20.0_DP * ( 3.0_DP / Pi )**( 2.0_DP / 3.0_DP ) &
!!$          * ( PlanckConstant / ( Erg * Second ) )**2 &
!!$          / ( AtomicMassUnit / Gram )**( 8.0_DP / 3.0_DP ) &
!!$          * ( Erg / Centimeter**3 ) &
!!$          / ( Gram / Centimeter**3 )**( 5.0_DP/ 3.0_DP )
      = 1.455e5_DP * Erg * Centimeter**3 / Gram**2


CONTAINS


  SUBROUTINE InitializeFields_Relativistic &
               ( AdvectionProfile_Option, &
                 RiemannProblemName_Option, &
                 nDetCells_Option, Eblast_Option, &
                 MassPNS_Option, ShockRadius_Option, &
                 AccretionRate_Option, MachNumber_Option, &
                 ApplyPerturbation_Option, PerturbationOrder_Option, &
                 PerturbationAmplitude_Option, &
                 rPerturbationInner_Option, rPerturbationOuter_Option, &
                 CentralDensity_Option, CentralPressure_Option, &
                 CoreRadius_Option, CollapseTime_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: AdvectionProfile_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: RiemannProblemName_Option
    INTEGER,          INTENT(in), OPTIONAL :: nDetCells_Option
    REAL(DP),         INTENT(in), OPTIONAL :: Eblast_Option
    REAL(DP),         INTENT(in), OPTIONAL :: MassPNS_Option
    REAL(DP),         INTENT(in), OPTIONAL :: ShockRadius_Option
    REAL(DP),         INTENT(in), OPTIONAL :: AccretionRate_Option
    REAL(DP),         INTENT(in), OPTIONAL :: MachNumber_Option
    LOGICAL,          INTENT(in), OPTIONAL :: ApplyPerturbation_Option
    INTEGER,          INTENT(in), OPTIONAL :: PerturbationOrder_Option
    REAL(DP),         INTENT(in), OPTIONAL :: PerturbationAmplitude_Option
    REAL(DP),         INTENT(in), OPTIONAL :: rPerturbationInner_Option
    REAL(DP),         INTENT(in), OPTIONAL :: rPerturbationOuter_Option
    REAL(DP),         INTENT(in), OPTIONAL :: CentralDensity_Option
    REAL(DP),         INTENT(in), OPTIONAL :: CentralPressure_Option
    REAL(DP),         INTENT(in), OPTIONAL :: CoreRadius_Option
    REAL(DP),         INTENT(in), OPTIONAL :: CollapseTime_Option

    CHARACTER(LEN=64) :: AdvectionProfile = 'SineWave'
    CHARACTER(LEN=64) :: RiemannProblemName = 'Sod'

    ! --- Sedov-Taylor Blast Wave (Defaults) ---
    INTEGER  :: nDetCells = 1
    REAL(DP) :: Eblast    = 1.0d-3

    ! --- Standing Accretion Shock (Defaults) ---
    REAL(DP) :: MassPNS               = 1.4_DP   * SolarMass
    REAL(DP) :: ShockRadius           = 180.0_DP * Kilometer
    REAL(DP) :: AccretionRate         = 0.3_DP   * SolarMass / Second
    REAL(DP) :: MachNumber            = 10.0_DP
    LOGICAL  :: ApplyPerturbation     = .FALSE.
    INTEGER  :: PerturbationOrder     = 0
    REAL(DP) :: PerturbationAmplitude = 0.0_DP
    REAL(DP) :: rPerturbationInner    = 0.0_DP
    REAL(DP) :: rPerturbationOuter    = 0.0_DP

    ! --- Yahil Collapse ---
    REAL(DP) :: CentralDensity  = 7.0e9_DP  * Gram / Centimeter**3
    REAL(DP) :: CentralPressure = 6.0e27_DP * Erg  / Centimeter**3
    REAL(DP) :: CoreRadius      = 1.0e5_DP  * Kilometer
    REAL(DP) :: CollapseTime    = 1.50e2_DP * Millisecond

    IF( PRESENT( AdvectionProfile_Option ) ) &
      AdvectionProfile = TRIM( AdvectionProfile_Option )

    IF( PRESENT( RiemannProblemName_Option ) ) &
      RiemannProblemName = TRIM( RiemannProblemName_Option )

    IF( PRESENT( nDetCells_Option ) ) &
      nDetCells = nDetCells_Option
    IF( PRESENT( Eblast_Option ) ) &
      Eblast = Eblast_Option

    IF( PRESENT( MassPNS_Option ) ) &
      MassPNS = MassPNS_Option
    IF( PRESENT( ShockRadius_Option ) ) &
      ShockRadius = ShockRadius_Option
    IF( PRESENT( AccretionRate_Option ) ) &
      AccretionRate = AccretionRate_Option
    IF( PRESENT( MachNumber_Option ) ) &
      MachNumber = MachNumber_Option
    IF( PRESENT( ApplyPerturbation_Option ) ) &
      ApplyPerturbation = ApplyPerturbation_Option
    IF( PRESENT( PerturbationOrder_Option ) ) &
      PerturbationOrder = PerturbationOrder_Option
    IF( PRESENT( PerturbationAmplitude_Option ) ) &
      PerturbationAmplitude = PerturbationAmplitude_Option
    IF( PRESENT( rPerturbationInner_Option ) ) &
      rPerturbationInner = rPerturbationInner_Option
    IF( PRESENT( rPerturbationOuter_Option ) ) &
      rPerturbationOuter = rPerturbationOuter_Option

    IF( PRESENT( CentralDensity_Option ) ) &
      CentralDensity = CentralDensity_Option
    IF( PRESENT( CentralPressure_Option ) ) &
      CentralPressure = CentralPressure_Option
    IF( PRESENT( CoreRadius_Option ) ) &
      CoreRadius = CoreRadius_Option
    IF( PRESENT( CollapseTime_Option ) ) &
      CollapseTime = CollapseTime_Option

    WRITE(*,*)
    WRITE(*,'(A,A)') '    INFO: ', TRIM( ProgramName )

    SELECT CASE ( TRIM( ProgramName ) )

      CASE( 'SlopeLimiterTest' )

        CALL InitializeFields_SlopeLimiterTest

      CASE( 'Advection' )

        CALL InitializeFields_Advection &
               ( TRIM( AdvectionProfile ) )

      CASE( 'Advection2D' )

        CALL InitializeFields_Advection2D &
               ( TRIM( AdvectionProfile ) )

      CASE( 'RiemannProblem' )

        CALL InitializeFields_RiemannProblem &
               ( TRIM( RiemannProblemName ), &
                 nDetCells_Option = nDetCells, &
                 Eblast_Option    = Eblast )

      CASE( 'RiemannProblem2D' )

        CALL InitializeFields_RiemannProblem2D &
               ( TRIM( RiemannProblemName ) )

      CASE( 'RiemannProblemSpherical' )

        CALL InitializeFields_RiemannProblemSpherical &
               ( TRIM( RiemannProblemName ) )

      CASE( 'SedovTaylorBlastWave' )

        CALL InitializeFields_SedovTaylorBlastWave &
               ( nDetCells, Eblast )

      CASE( 'KelvinHelmholtzInstability' )

         CALL InitializeFields_KelvinHelmholtzInstability

      CASE( 'StandingAccretionShock' )

        CALL InitializeFields_StandingAccretionShock &
               ( MassPNS, ShockRadius, AccretionRate, MachNumber, &
                 ApplyPerturbation, PerturbationOrder, PerturbationAmplitude, &
                 rPerturbationInner, rPerturbationOuter )

      CASE( 'StaticTOV' )

         CALL InitializeFields_StaticTOV

      CASE( 'YahilCollapse' )

         CALL InitializeFields_YahilCollapse &
                ( CentralDensity, CentralPressure, CoreRadius, CollapseTime )

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A21,A)') 'Invalid ProgramName: ', ProgramName
        WRITE(*,'(A)')     'Stopping...'
        STOP

    END SELECT

  END SUBROUTINE InitializeFields_Relativistic


  SUBROUTINE InitializeFields_SlopeLimiterTest

    INTEGER       :: iX1, iX2, iX3, iPF
    INTEGER       :: iNodeX, iNodeX1, iNodeX2, iNodeX3
    REAL(DP)      :: X, X1, X2
    REAL(DP)      :: a0, a1, a2, a3, a4, a5
    REAL(DP)      :: X0, theta
    CHARACTER(32) :: Problem

    INTEGER  :: M, M1, M2, M3, nDOFQ
    INTEGER  :: qNodeX, qNodeX1, qNodeX2, qNodeX3
    REAL(DP) :: etaG(nNodesX(1)), xG(nNodesX(1)), wG(nNodesX(1))
    REAL(DP), ALLOCATABLE :: etaQ_X1(:), wQ_X1(:)
    REAL(DP), ALLOCATABLE :: etaQ_X2(:), wQ_X2(:)
    REAL(DP), ALLOCATABLE :: etaQ_X3(:), wQ_X3(:)
    REAL(DP), ALLOCATABLE :: etaQ(:),    wQ(:)
    REAL(DP), ALLOCATABLE :: u0(:,:)
    INTEGER,  ALLOCATABLE :: NodeNumberTableQ(:,:)

    REAL(DP), ALLOCATABLE :: InterpolationMatrix(:,:)

    M = 5
    nDOFQ = M**nDimsX

    M1 = M
    M2 = 1
    M3 = 1
    IF( nDimsX .GT. 1 ) M2 = M
    IF( nDimsX .GT. 2 ) M3 = M

    nDOFQ = M1 * M2 * M3

    ALLOCATE( etaQ_X1(M1), wQ_X1(M1) )
    ALLOCATE( etaQ_X2(M2), wQ_X2(M2) )
    ALLOCATE( etaQ_X3(M3), wQ_X3(M3) )
    ALLOCATE( etaQ(nDOFQ), wQ(nDOFQ) )
    ALLOCATE( u0(nDOFQ,nPF) )
    ALLOCATE( InterpolationMatrix(nDOFX,nDOFQ) )
    ALLOCATE( NodeNumberTableQ(3,nDOFQ) )

    CALL GetQuadrature( M1        , etaQ_X1, wQ_X1 )
    CALL GetQuadrature( M2        , etaQ_X2, wQ_X2 )
    CALL GetQuadrature( M3        , etaQ_X3, wQ_X3 )
    CALL GetQuadrature( nNodesX(1), etaG, wG )

    qNodeX = 0
    DO qNodeX3 = 1, M3
    DO qNodeX2 = 1, M2
    DO qNodeX1 = 1, M1

      qNodeX = qNodeX + 1

      NodeNumberTableQ(1:3,qNodeX) &
        = [ qNodeX1, qNodeX2, qNodeX3 ]

      etaQ(qNodeX) = etaQ_X1(qNodeX1) * etaQ_X2(qNodeX2) * etaQ_X3(qNodeX3)
      wQ  (qNodeX) = wQ_X1  (qNodeX1) * wQ_X2  (qNodeX2) * wQ_X3  (qNodeX3)

    END DO
    END DO
    END DO

    DO iNodeX = 1, nDOFX

      iNodeX1 = NodeNumberTableX(1,iNodeX)
      iNodeX2 = NodeNumberTableX(2,iNodeX)
      iNodeX3 = NodeNumberTableX(3,iNodeX)

      DO qNodeX = 1, nDOFQ

        qNodeX1 = NodeNumberTableQ(1,qNodeX)
        qNodeX2 = NodeNumberTableQ(2,qNodeX)
        qNodeX3 = NodeNumberTableQ(3,qNodeX)

        InterpolationMatrix(iNodeX,qNodeX) &
          = wQ(qNodeX) &
              * LagrangeP( etaQ_X1(qNodeX1), iNodeX1, etaG, nNodesX(1) ) &
              * LagrangeP( etaQ_X2(qNodeX2), iNodeX2, etaG, nNodesX(2) ) &
              * LagrangeP( etaQ_X3(qNodeX3), iNodeX3, etaG, nNodesX(3) )

      END DO

    END DO

    Problem = 'SmoothTANH'

    ! --- Coefficients for polynomial problem ---
    a0 = 1.0_DP
    a1 = 0.1_DP
    a2 = 4.0_DP
    a3 = -0.4_DP
    a4 = -1.0_DP

    ! --- Scale length for TANH problem ---
    X0 = 0.01_DP

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO qNodeX = 1, nDOFQ

        qNodeX1 = NodeNumberTableQ(1,qNodeX)
        qNodeX2 = NodeNumberTableQ(2,qNodeX)

        X1 = etaQ_X1(qNodeX1) * MeshX(1) % Width(iX1) + MeshX(1) % Center(iX1)
        X2 = etaQ_X2(qNodeX2) * MeshX(2) % Width(iX2) + MeshX(2) % Center(iX2)

        X = X1

        u0(qNodeX,iPF_V1) = Zero
        u0(qNodeX,iPF_V2) = Zero
        u0(qNodeX,iPF_V3) = Zero
        u0(qNodeX,iPF_E ) = One / ( Gamma_IDEAL - One )
        u0(qNodeX,iPF_Ne) = Zero

        SELECT CASE( Problem )

          CASE( 'Polynomial' )

            u0(qNodeX,iPF_D) = a0*X**0 + a1*X**1 + a2*X**2 + a3*X**3 + a4*X**4

          CASE( 'SmoothTANH' )

            u0(qNodeX,iPF_D) = 2.0d0 + TANH( X / X0 )

          CASE( 'SinCosTANH' )

            theta = Half * ( One - TANH( X / X0 ) )

            u0(qNodeX,iPF_D) = 1.0d0 + theta * COS( TwoPi * X ) &
                                 + ( One - theta ) * SIN( TwoPi * X )

          CASE( 'IsolatedContact' )

            IF( X .LT. Zero )THEN

               u0(qNodeX,iPF_D) = 0.1_DP

           ELSE

               u0(qNodeX,iPF_D) = 0.5_DP

           END IF

        END SELECT

      END DO

      DO iPF = 1, nPF

        uPF(:,iX1,iX2,iX3,iPF) &
          = MATMUL( InterpolationMatrix, u0(:,iPF) ) / WeightsX_q

      END DO

      CALL ComputePressureFromPrimitive_IDEAL &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_E), &
               uPF(:,iX1,iX2,iX3,iPF_Ne), uAF(:,iX1,iX2,iX3,iAF_P) )

      CALL ComputeConserved_Euler_Relativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uAF(:,iX1,iX2,iX3,iAF_P) )

    END DO
    END DO
    END DO

    DEALLOCATE( NodeNumberTableQ )
    DEALLOCATE( InterpolationMatrix )
    DEALLOCATE( u0 )
    DEALLOCATE( etaQ   , wQ )
    DEALLOCATE( etaQ_X3, wQ_X3 )
    DEALLOCATE( etaQ_X2, wQ_X2 )
    DEALLOCATE( etaQ_X1, wQ_X1 )

  END SUBROUTINE InitializeFields_SlopeLimiterTest


  SUBROUTINE InitializeFields_Advection( AdvectionProfile )

    CHARACTER(LEN=*), INTENT(in) :: AdvectionProfile

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1
    REAL(DP) :: X1

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Advection Profile: ', TRIM( AdvectionProfile )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        SELECT CASE( TRIM( AdvectionProfile ) )

          CASE( 'SineWave' )

            uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = One + 0.1_DP * SIN( TwoPi * X1 )
            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.1_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
            uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.0_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_E )  &
              = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

          CASE( 'TopHat' )

            IF( X1 .GT. 0.45 .AND. X1 .LT. 0.55 )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 2.0_DP

            ELSE

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP

            END IF

              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.1_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

          CASE DEFAULT

            WRITE(*,*)
            WRITE(*,'(A,A)') &
              'Invalid choice for AdvectionProfile: ', AdvectionProfile
            WRITE(*,'(A)') 'Valid choices:'
            WRITE(*,'(A)') '  SineWave'
            WRITE(*,'(A)') '  TopHat'
            WRITE(*,*)
            WRITE(*,'(A)') 'Stopping...'
            STOP

        END SELECT

      END DO

      CALL ComputeConserved_Euler_Relativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uAF(:,iX1,iX2,iX3,iAF_P) )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_Advection


  SUBROUTINE InitializeFields_Advection2D( AdvectionProfile )

    CHARACTER(LEN=*), INTENT(in) :: AdvectionProfile

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: X1, X2

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Advection Profile: ', TRIM( AdvectionProfile )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

        SELECT CASE( TRIM( AdvectionProfile ) )

          CASE( 'SineWaveX1' )

            uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = One + 0.1_DP * SIN( TwoPi * X1 )
            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.1_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
            uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.0_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_E )  &
              = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

          CASE( 'SineWaveX2' )

            uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = One + 0.1_DP * SIN( TwoPi * X2 )
            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.1_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
            uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.0_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_E )  &
              = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

          CASE DEFAULT

            WRITE(*,*)
            WRITE(*,'(A,A)') &
              'Invalid choice for AdvectionProfile: ', AdvectionProfile
            WRITE(*,'(A)') 'Valid choices:'
            WRITE(*,'(A)') '  SineWaveX1'
            WRITE(*,'(A)') '  SineWaveX2'
            WRITE(*,*)
            WRITE(*,'(A)') 'Stopping...'
            STOP

        END SELECT

      END DO

      CALL ComputeConserved_Euler_Relativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uAF(:,iX1,iX2,iX3,iAF_P) )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_Advection2D


  SUBROUTINE InitializeFields_RiemannProblem &
               ( RiemannProblemName, &
                 nDetCells_Option, Eblast_Option )

    CHARACTER(LEN=*), INTENT(in)           :: RiemannProblemName
    INTEGER,          INTENT(in), OPTIONAL :: nDetCells_Option
    REAL(DP),         INTENT(in), OPTIONAL :: Eblast_Option

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1
    REAL(DP) :: X1, XD, Vs

    INTEGER  :: nDetCells
    REAL(DP) :: Eblast

    REAL(DP) :: LeftState(nPF), RightState(nPF)

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Riemann Problem Name: ', TRIM( RiemannProblemName )
    WRITE(*,*)

    SELECT CASE( TRIM( RiemannProblemName ) )

      CASE( 'Sod' )

        XD = 0.5_DP

        LeftState(iPF_D ) = 1.0_DP
        LeftState(iPF_V1) = 0.0_DP
        LeftState(iPF_V2) = 0.0_DP
        LeftState(iPF_V3) = 0.0_DP
        LeftState(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

        RightState(iPF_D ) = 0.125_DP
        RightState(iPF_V1) = 0.0_DP
        RightState(iPF_V2) = 0.0_DP
        RightState(iPF_V3) = 0.0_DP
        RightState(iPF_E ) = 0.1_DP / ( Gamma_IDEAL - One )

      CASE( 'IsolatedShock' )

        Vs = 0.01_DP
        XD = Half

        RightState(iPF_D)  = 1.0_DP
        RightState(iPF_V1) = -0.9_DP
        RightState(iPF_V2) = 0.0_DP
        RightState(iPF_V3) = 0.0_DP
        RightState(iPF_E)  = 1.0_DP / ( Gamma_IDEAL - One )

        CALL ComputeLeftState &
               ( Vs,                 &
                 RightState(iPF_D ), &
                 RightState(iPF_V1), &
                 RightState(iPF_E ) * ( Gamma_IDEAL - One ), &
                 LeftState (iPF_D ), &
                 LeftState (iPF_V1), &
                 LeftState (iPF_E ) )

        LeftState(iPF_V2) = 0.0_DP
        LeftState(iPF_V3) = 0.0_DP

      CASE( 'IsolatedContact' )

        Vs = 0.01_DP
        XD = 0.5_DP

        LeftState(iPF_D ) = 5.9718209694880811e0_DP
        LeftState(iPF_V1) = Vs
        LeftState(iPF_V2) = 0.0_DP
        LeftState(iPF_V3) = 0.0_DP
        LeftState(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

        RightState(iPF_D ) = 1.0_DP
        RightState(iPF_V1) = Vs
        RightState(iPF_V2) = 0.0_DP
        RightState(iPF_V3) = 0.0_DP
        RightState(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

      CASE( 'MBProblem1' )

        XD = 0.5_DP

        LeftState(iPF_D ) = 1.0_DP
        LeftState(iPF_V1) = 0.9_DP
        LeftState(iPF_V2) = 0.0_DP
        LeftState(iPF_V3) = 0.0_DP
        LeftState(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

        RightState(iPF_D ) = 1.0_DP
        RightState(iPF_V1) = 0.0_DP
        RightState(iPF_V2) = 0.0_DP
        RightState(iPF_V3) = 0.0_DP
        RightState(iPF_E ) = 10.0_DP / ( Gamma_IDEAL - One )

      CASE( 'MBProblem4' )

        XD = 0.5_DP

        LeftState(iPF_D ) = 1.0_DP
        LeftState(iPF_V1) = 0.0_DP
        LeftState(iPF_V2) = 0.0_DP
        LeftState(iPF_V3) = 0.0_DP
        LeftState(iPF_E ) = 1.0e3_DP / ( Gamma_IDEAL - One )

        RightState(iPF_D ) = 1.0_DP
        RightState(iPF_V1) = 0.0_DP
        RightState(iPF_V2) = 0.0_DP
        RightState(iPF_V3) = 0.0_DP
        RightState(iPF_E ) = 1.0e-2_DP / ( Gamma_IDEAL - One )

      CASE( 'PerturbedShockTube' )

        XD = 0.5_DP

        LeftState(iPF_D ) = 5.0_DP
        LeftState(iPF_V1) = 0.0_DP
        LeftState(iPF_V2) = 0.0_DP
        LeftState(iPF_V3) = 0.0_DP
        LeftState(iPF_E ) = 50.0_DP / ( Gamma_IDEAL - One )

        RightState(iPF_D ) = 0.0_DP ! --- Dummy ---
        RightState(iPF_V1) = 0.0_DP
        RightState(iPF_V2) = 0.0_DP
        RightState(iPF_V3) = 0.0_DP
        RightState(iPF_E ) = 5.0_DP / ( Gamma_IDEAL - One )

      CASE( 'ShockReflection' )

        XD = 1.0_DP

        LeftState(iPF_D ) = 1.0_DP
        LeftState(iPF_V1) = 0.99999_DP
        LeftState(iPF_V2) = 0.0_DP
        LeftState(iPF_V3) = 0.0_DP
        LeftState(iPF_E ) = 0.01_DP / ( Gamma_IDEAL - One )

        ! --- All of these are dummies ---
        RightState(iPF_D ) = 0.0_DP
        RightState(iPF_V1) = 0.0_DP
        RightState(iPF_V2) = 0.0_DP
        RightState(iPF_V3) = 0.0_DP
        RightState(iPF_E ) = 0.0_DP

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A,A)') &
          'Invalid choice for RiemannProblemName: ', RiemannProblemName
        WRITE(*,'(A)') 'Valid choices:'
        WRITE(*,'(A)') &
          "  'Sod' - &
          Sod's shock tube"
        WRITE(*,'(A)') &
          "  'MBProblem1' - &
          Mignone & Bodo (2005) MNRAS, 364, 126, Problem 1"
        WRITE(*,'(A)') &
          "  'MBProblem4' - &
          Mignone & Bodo (2005) MNRAS, 364, 126, Problem 4"
        WRITE(*,'(A)') &
          "  'PerturbedShockTube' - &
          Del Zanna & Bucciantini (2002) AA, 390, 1177, &
          Sinusoidal density perturbation"
        WRITE(*,'(A)') &
          "  'ShockReflection' - &
          Del Zanna & Bucciantini (2002) AA, 390, 1177, &
          Planar shock reflection"
        WRITE(*,'(A)') 'Stopping...'
        STOP

    END SELECT

    IF( TRIM( RiemannProblemName ) .EQ. 'IsolatedShock' )THEN

      WRITE(*,'(6x,A,ES14.6E3)') 'Shock Velocity = ', Vs
      WRITE(*,*)

    END IF

    WRITE(*,'(6x,A,F8.6)') 'Gamma_IDEAL = ', Gamma_IDEAL
    WRITE(*,*)
    WRITE(*,'(6x,A,F8.6)') 'XD = ', XD
    WRITE(*,*)
    WRITE(*,'(6x,A)') 'Right State:'
    WRITE(*,*)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_D  = ', RightState(iPF_D )
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V1 = ', RightState(iPF_V1)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V2 = ', RightState(iPF_V2)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V3 = ', RightState(iPF_V3)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_E  = ', RightState(iPF_E )
    WRITE(*,*)
    WRITE(*,'(6x,A)') 'Left State:'
    WRITE(*,*)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_D  = ', LeftState(iPF_D )
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V1 = ', LeftState(iPF_V1)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V2 = ', LeftState(iPF_V2)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V3 = ', LeftState(iPF_V3)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_E  = ', LeftState(iPF_E )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        IF( X1 .LE. XD )THEN

          uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = LeftState(iPF_D )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = LeftState(iPF_V1)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = LeftState(iPF_V2)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = LeftState(iPF_V3)
          uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = LeftState(iPF_E )

        ELSE

          uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = RightState(iPF_D )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = RightState(iPF_V1)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = RightState(iPF_V2)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = RightState(iPF_V3)
          uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = RightState(iPF_E )

          IF( TRIM( RiemannProblemName ) .EQ. 'PerturbedShockTube' ) &
            uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
              = 2.0_DP + 0.3_DP * SIN( 50.0_DP * X1 )

        END IF

      END DO

      CALL ComputePressureFromPrimitive_IDEAL &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_E ), &
               uPF(:,iX1,iX2,iX3,iPF_Ne), uAF(:,iX1,iX2,iX3,iAF_P) )

      CALL ComputeConserved_Euler_Relativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uAF(:,iX1,iX2,iX3,iAF_P) )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_RiemannProblem


  SUBROUTINE InitializeFields_RiemannProblem2D( RiemannProblemName )

    CHARACTER(LEN=*), INTENT(in) :: RiemannProblemName

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: X1, X2, X1D, X2D, Vs

    REAL(DP) :: NE(nPF), NW(nPF), SW(nPF), SE(nPF)

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', '2D Riemann Problem Name: ', TRIM( RiemannProblemName )
    WRITE(*,*)

    SELECT CASE( TRIM( RiemannProblemName ) )

      CASE( 'DzB2002' )

        X1D = 0.5_DP
        X2D = 0.5_DP

        NE(iPF_D ) = 0.1_DP
        NE(iPF_V1) = 0.0_DP
        NE(iPF_V2) = 0.0_DP
        NE(iPF_V3) = 0.0_DP
        NE(iPF_E ) = 0.01_DP / ( Gamma_IDEAL - One )

        NW(iPF_D ) = 0.1_DP
        NW(iPF_V1) = 0.99_DP
        NW(iPF_V2) = 0.0_DP
        NW(iPF_V3) = 0.0_DP
        NW(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

        SW(iPF_D ) = 0.5_DP
        SW(iPF_V1) = 0.0_DP
        SW(iPF_V2) = 0.0_DP
        SW(iPF_V3) = 0.0_DP
        SW(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

        SE(iPF_D ) = 0.1_DP
        SE(iPF_V1) = 0.0_DP
        SE(iPF_V2) = 0.99_DP
        SE(iPF_V3) = 0.0_DP
        SE(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

      CASE( 'IsolatedShock' )

        Vs = 0.01_DP

        X1D = 0.5_DP
        X2D = 0.5_DP

        NE(iPF_D ) = 1.0_DP
        NE(iPF_V1) = -0.9_DP
        NE(iPF_V2) = 0.0_DP
        NE(iPF_V3) = 0.0_DP
        NE(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

        CALL ComputeLeftState &
               ( Vs, &
                 NE(iPF_D ), &
                 NE(iPF_V1), &
                 NE(iPF_E ) * ( Gamma_IDEAL - One ), &
                 NW(iPF_D ), &
                 NW(iPF_V1), &
                 NW(iPF_E ) )

        NW(iPF_V2) = 0.0_DP
        NW(iPF_V3) = 0.0_DP

        SE(iPF_D ) = 1.0_DP
        SE(iPF_V1) = -0.9_DP
        SE(iPF_V2) = 0.0_DP
        SE(iPF_V3) = 0.0_DP
        SE(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

        CALL ComputeLeftState &
               ( Vs, &
                 SE(iPF_D ), &
                 SE(iPF_V1), &
                 SE(iPF_E ) * ( Gamma_IDEAL - One ), &
                 SW(iPF_D ), &
                 SW(iPF_V1), &
                 SW(iPF_E ) )

        SW(iPF_V2) = 0.0_DP
        SW(iPF_V3) = 0.0_DP

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A,A)') &
          'Invalid choice for RiemannProblemName: ', TRIM( RiemannProblemName )
        WRITE(*,'(A)') 'Valid choices:'
        WRITE(*,'(A)') &
          "  'DzB2002' - &
          Blast wave from Del-Zanna & Bucciantini (2002)"
        WRITE(*,'(A)') &
          "  'IsolatedShock'"
        WRITE(*,'(A)') 'Stopping...'
        STOP

    END SELECT

    IF( TRIM( RiemannProblemName ) .EQ. 'IsolatedShock' )THEN

      WRITE(*,'(6x,A,ES14.6E3)') 'Shock Velocity = ', Vs
      WRITE(*,*)

    END IF

    WRITE(*,'(6x,A,F8.6)') 'Gamma_IDEAL = ', Gamma_IDEAL
    WRITE(*,*)
    WRITE(*,'(6x,A,F8.6)') 'X1D = ', X1D
    WRITE(*,'(6x,A,F8.6)') 'X2D = ', X2D
    WRITE(*,*)
    WRITE(*,'(6x,A)') 'NE:'
    WRITE(*,*)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_D  = ', NE(iPF_D )
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V1 = ', NE(iPF_V1)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V2 = ', NE(iPF_V2)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V3 = ', NE(iPF_V3)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_E  = ', NE(iPF_E )
    WRITE(*,*)
    WRITE(*,'(6x,A)') 'NW:'
    WRITE(*,*)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_D  = ', NW(iPF_D )
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V1 = ', NW(iPF_V1)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V2 = ', NW(iPF_V2)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V3 = ', NW(iPF_V3)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_E  = ', NW(iPF_E )
    WRITE(*,*)
    WRITE(*,'(6x,A)') 'SE:'
    WRITE(*,*)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_D  = ', SE(iPF_D )
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V1 = ', SE(iPF_V1)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V2 = ', SE(iPF_V2)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V3 = ', SE(iPF_V3)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_E  = ', SE(iPF_E )
    WRITE(*,*)
    WRITE(*,'(6x,A)') 'SW:'
    WRITE(*,*)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_D  = ', SW(iPF_D )
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V1 = ', SW(iPF_V1)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V2 = ', SW(iPF_V2)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V3 = ', SW(iPF_V3)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_E  = ', SW(iPF_E )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

        ! --- NE ---
        IF     ( X1 .GT. X1D .AND. X2 .GT. X2D )THEN

          uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = NE(iPF_D )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = NE(iPF_V1)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = NE(iPF_V2)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = NE(iPF_V3)
          uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = NE(iPF_E )

        ! --- NW ---
        ELSE IF( X1 .LE. X1D .AND. X2 .GT. X2D )THEN

          uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = NW(iPF_D )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = NW(iPF_V1)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = NW(iPF_V2)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = NW(iPF_V3)
          uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = NW(iPF_E )

        ! --- SW ---
        ELSE IF( X1 .LE. X1D .AND. X2 .LE. X2D )THEN

          uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = SW(iPF_D )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = SW(iPF_V1)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = SW(iPF_V2)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = SW(iPF_V3)
          uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = SW(iPF_E )

        ! --- SE ---
        ELSE

          uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = SE(iPF_D )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = SE(iPF_V1)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = SE(iPF_V2)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = SE(iPF_V3)
          uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = SE(iPF_E )

        END IF

      END DO

      CALL ComputePressureFromPrimitive_IDEAL &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_E ), &
               uPF(:,iX1,iX2,iX3,iPF_Ne), uAF(:,iX1,iX2,iX3,iAF_P) )

      CALL ComputeConserved_Euler_Relativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uAF(:,iX1,iX2,iX3,iAF_P) )

    END DO
    END DO
    END DO


  END SUBROUTINE InitializeFields_RiemannProblem2D


  SUBROUTINE InitializeFields_RiemannProblemSpherical( RiemannProblemName )

    CHARACTER(LEN=*), INTENT(in) :: RiemannProblemName

    INTEGER       :: iX1, iX2, iX3
    INTEGER       :: iNodeX, iNodeX1
    REAL(DP)      :: X1

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Spherical Riemann Problem Name: ', &
        TRIM( RiemannProblemName )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        SELECT CASE ( TRIM( RiemannProblemName ) )

          CASE( 'SphericalSod' )

            IF( X1 <= One )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

            ELSE

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 0.125_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 0.1_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

            END IF

         CASE DEFAULT

            WRITE(*,*)
            WRITE(*,*) &
              'Invalid choice for RiemannProblemName: ', &
              RiemannProblemName
            WRITE(*,*) 'Valid choices:'
            WRITE(*,*) &
              "'SphericalSod' - ", &
              "Spherical Sod's shock tube"
            STOP

          END SELECT

        END DO

      CALL ComputeConserved_Euler_Relativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uAF(:,iX1,iX2,iX3,iAF_P) )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_RiemannProblemSpherical


  SUBROUTINE InitializeFields_SedovTaylorBlastWave( nDetCells, Eblast )

    INTEGER,  INTENT(in) :: nDetCells
    REAL(DP), INTENT(in) :: Eblast

    INTEGER  :: iX1, iX2, iX3, iNodeX1, iNodeX
    REAL(DP) :: X1, X_D

    X_D = DBLE( nDetCells ) * MeshX(1) % Width(1)
    WRITE(*,*)
    WRITE(*,'(A,I4.4)')      '     nDetCells:              ', nDetCells
    WRITE(*,'(A,ES23.16E3)') '     Initial blast radius:   ', X_D
    WRITE(*,'(A,ES23.16E3)') '     Blast energy:           ', Eblast
    WRITE(*,'(A,ES23.16E3)') '     Initial blast pressure: ', &
                                     ( Gamma_IDEAL - One ) &
                                       * Eblast / ( FourPi / Three * X_D**3 )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        IF( X1 <= X_D)THEN

          uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
          uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
            = Eblast / ( FourPi / Three * X_D**3 )
          uAF(iNodeX,iX1,iX2,iX3,iAF_P)  &
            = ( Gamma_IDEAL - One ) * uPF(iNodeX,iX1,iX2,iX3,iPF_E)

        ELSE

          uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
          uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
            = 1.0d-5
          uAF(iNodeX,iX1,iX2,iX3,iAF_P)  &
            = ( Gamma_IDEAL - One ) * uPF(iNodeX,iX1,iX2,iX3,iPF_E)

        END IF

      END DO

      CALL ComputeConserved_Euler_Relativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uAF(:,iX1,iX2,iX3,iAF_P) )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_SedovTaylorBlastWave


  ! --- Relativistic 2D Kelvin-Helmholtz instability a la
  !     Beckwith & Stone (2011), ApjS, 193, 6 (typo in Eq. (63)) ---
  SUBROUTINE InitializeFields_KelvinHelmholtzInstability

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: X1, X2
    REAL(DP) :: rho0, rho1
    REAL(DP) :: Vshear, a, X2_Offset, sigma, A0

    rho0 = 0.505d0
    rho1 = 0.495d0

    Vshear    = 0.5d0
    a         = 0.01d0
    X2_Offset = 0.5d0
    sigma     = 0.1d0

    A0 = 0.1d0

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

        ! --- Top ---
        IF( X2 .GT. 0.0d0 )THEN
          uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
            = rho0 + rho1 * TANH( ( X2 - X2_Offset ) / a )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
            = Vshear      * TANH( ( X2 - X2_Offset ) / a )

          ! --- This is where the typo is. The following expression is
          !     taken from Radice & Rezzolla, 2012, AA, 547, A26, Eq. (48) ---
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
            = A0 * Vshear * SIN( 2.0d0 * Pi * X1 ) &
                * EXP( -( ( X2 - X2_Offset ) / sigma )**2 )

        ! --- Bottom ---
        ELSE
          uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
            = rho0 - rho1 * TANH( ( X2 + X2_Offset ) / a )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
            = -Vshear     * TANH( ( X2 + X2_Offset ) / a )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
            = -A0 * Vshear * SIN( 2.0d0 * Pi * X1 ) &
                * EXP( -( ( X2 + X2_Offset ) / sigma )**2 )

         END IF

        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0d0
        uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0d0
        uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
          = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

      END DO

      CALL ComputeConserved_Euler_Relativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uAF(:,iX1,iX2,iX3,iAF_P) )

    END DO
    END DO
    END DO


  END SUBROUTINE InitializeFields_KelvinHelmholtzInstability


  SUBROUTINE InitializeFields_StandingAccretionShock &
    ( MassPNS, ShockRadius, AccretionRate, MachNumber, &
      ApplyPerturbation, PerturbationOrder, PerturbationAmplitude, &
      rPerturbationInner, rPerturbationOuter )

    REAL(DP), INTENT(in) :: MassPNS, ShockRadius, AccretionRate, MachNumber
    LOGICAL,  INTENT(in) :: ApplyPerturbation
    INTEGER,  INTENT(in) :: PerturbationOrder
    REAL(DP), INTENT(in) :: PerturbationAmplitude
    REAL(DP), INTENT(in) :: rPerturbationInner
    REAL(DP), INTENT(in) :: rPerturbationOuter

    INTEGER  :: iX1, iX2, iX3, iNodeX1, iNodeX2, iNodeX3, iNodeX
    INTEGER  :: iX1_1, iX1_2, iNodeX1_1, iNodeX1_2
    REAL(DP) :: X1_1, X1_2, D_1, D_2, V_1, V_2, P_2
    REAL(DP) :: Alpha, Psi, V0, VSq, W
    REAL(DP) :: X1, X2, dX1, PolytropicConstant, MassConstant
    REAL(DP) :: D(1:nNodesX(1),iX_B1(1):iX_E1(1))
    REAL(DP) :: V(1:nNodesX(1),iX_B1(1):iX_E1(1))
    REAL(DP) :: P(1:nNodesX(1),iX_B1(1):iX_E1(1))
    LOGICAL  :: FirstPreShockElement = .FALSE.

    WRITE(*,*)
    WRITE(*,'(6x,A,ES9.2E3,A)') &
      'Shock radius:   ', ShockRadius / Kilometer, ' km'
    WRITE(*,'(6x,A,ES9.2E3,A)') &
      'PNS Mass:       ', MassPNS / SolarMass, ' Msun'
    WRITE(*,'(6x,A,ES9.2E3,A)') &
      'Accretion Rate: ', AccretionRate / ( SolarMass / Second ), &
      ' Msun/s'
    WRITE(*,'(6x,A,ES9.2E3)') &
      'Mach number:    ', MachNumber
    WRITE(*,*)
    WRITE(*,'(6x,A,L)') &
      'Apply Perturbation: ', ApplyPerturbation
    WRITE(*,'(6x,A,I1)') &
      'Perturbation order: ', PerturbationOrder
    WRITE(*,'(6x,A,ES9.2E3)') &
      'Perturbation amplitude: ', PerturbationAmplitude
    WRITE(*,'(6x,A,ES9.2E3,A)') &
      'Inner radius of perturbation: ', rPerturbationInner / Kilometer, ' km'
    WRITE(*,'(6x,A,ES9.2E3,A)') &
      'Outer radius of perturbation: ', rPerturbationOuter / Kilometer, ' km'

    !  --- Locate first element containing un-shocked fluid ---

    X1 = 0.0_DP

    DO iX1 = iX_B1(1), iX_E1(1)

      DO iNodeX1 = 1, nNodesX(1)

        dX1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 ) - X1
        X1  = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        IF( X1 .LE. ShockRadius ) CYCLE

        IF( X1 .GT. ShockRadius .AND. .NOT. FirstPreShockElement )THEN

          iX1_1     = iX1
          iNodeX1_1 = iNodeX1
          X1_1      = X1
          X1_2      = X1 - dX1

          IF( iNodeX1_1 .EQ. 1 )THEN

            iX1_2     = iX1_1 - 1
            iNodeX1_2 = nNodesX(1)

          ELSE

            iX1_2     = iX1_1
            iNodeX1_2 = iNodeX1_1 - 1

          END IF

          FirstPreShockElement = .TRUE.

        END IF

      END DO

    END DO

    ! --- Compute fields, pre-shock ---

    ! --- Compute polytropic constant for pre-shock flow ---

    iX1     = iX_E1(1)
    iNodeX1 = nNodesX(1)

    X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

    Alpha = LapseFunction  ( X1, MassPNS )
    Psi   = ConformalFactor( X1, MassPNS )

    V(iNodeX1,iX1) &
      = -Psi**(-2) * SpeedOfLight * SQRT( One - Alpha**2 )

    D(iNodeX1,iX1) &
      = Psi**(-6) * AccretionRate &
          / ( FourPi * X1**2 * ABS( V(iNodeX1,iX1) ) )

    VSq = Psi**4 * V(iNodeX1,iX1)**2

    P(iNodeX1,iX1) &
      = D(iNodeX1,iX1) * VSq &
          / ( MachNumber**2 * Gamma_IDEAL &
                - Gamma_IDEAL / ( Gamma_IDEAL - One ) &
                * VSq / SpeedOfLight**2 )

    PolytropicConstant &
      = P(iNodeX1,iX1) * D(iNodeX1,iX1)**( -Gamma_IDEAL )

    DO iX1 = iX_E1(1), iX1_2, -1

      DO iNodeX1 = nNodesX(1), 1, -1

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        IF( X1 .LT. X1_1 ) CYCLE

        Alpha = LapseFunction  ( X1, MassPNS )
        Psi   = ConformalFactor( X1, MassPNS )

        V(iNodeX1,iX1) &
          = -Psi**(-2) * SpeedOfLight * SQRT( One - Alpha**2 )

        D(iNodeX1,iX1) &
          = Psi**(-6) * AccretionRate &
              / ( FourPi * X1**2 * ABS( V(iNodeX1,iX1) ) )

        P(iNodeX1,iX1) &
          = PolytropicConstant * D(iNodeX1,iX1)**( Gamma_IDEAL )

!!$        VSq = Psi**4 * V(iNodeX1,iX1)**2
!!$
!!$        P(iNodeX1,iX1) &
!!$          = D(iNodeX1,iX1) * VSq &
!!$              / ( Gamma_IDEAL * MachNumber**2 ) &
!!$              / ( One - ( VSq / SpeedOfLight**2 ) &
!!$              / ( MachNumber**2 * ( Gamma_IDEAL - One ) ) )

      END DO

    END DO

    ! --- Apply jump conditions ---

    D_1 = D(iNodeX1_1,iX1_1)
    V_1 = V(iNodeX1_1,iX1_1)

    CALL ApplyJumpConditions &
           ( iX1_1, iNodeX1_1, X1_1, D_1, V_1, &
             iX1_2, iNodeX1_2, X1_2, &
             D_2, V_2, P_2, MassPNS, PolytropicConstant )

    WRITE(*,*)
    WRITE(*,'(6x,A)') 'Shock location:'
    WRITE(*,'(8x,A)') 'Pre-shock:'
    WRITE(*,'(10x,A,I4.4)')       'iX1     = ', iX1_1
    WRITE(*,'(10x,A,I2.2)')       'iNodeX1 = ', iNodeX1_1
    WRITE(*,'(10x,A,ES13.6E3,A)') 'X1      = ', X1_1 / Kilometer, ' km'
    WRITE(*,'(8x,A)') 'Post-shock:'
    WRITE(*,'(10x,A,I4.4)')       'iX1     = ', iX1_2
    WRITE(*,'(10x,A,I2.2)')       'iNodeX1 = ', iNodeX1_2
    WRITE(*,'(10x,A,ES13.6E3,A)') 'X1      = ', X1_2 / Kilometer, ' km'
    WRITE(*,*)
    WRITE(*,'(6x,A,ES13.6E3)') &
      'Compression Ratio LOG10(D_2/D_1) = ', LOG( D_2 / D_1 ) / LOG( 10.0_DP )
    WRITE(*,*)

    ! --- Compute fields, post-shock ---

    Alpha = LapseFunction  ( X1_1, MassPNS )
    Psi   = ConformalFactor( X1_1, MassPNS )
    W     = LorentzFactor( Psi, V_1 )

    MassConstant = Psi**6 * Alpha * X1_1**2 * D_1 * W * V_1

    V0 = V_2

    DO iX1 = iX1_2, iX_B1(1), -1

      DO iNodeX1 = nNodesX(1), 1, -1

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        IF( X1 .GT. X1_2 ) CYCLE

        Alpha = LapseFunction  ( X1, MassPNS )
        Psi   = ConformalFactor( X1, MassPNS )

        CALL NewtonRaphson_PostShockVelocity &
               ( Alpha, Psi, MassConstant, PolytropicConstant, &
                 MassPNS, AccretionRate, X1, V0  )

        V(iNodeX1,iX1) = V0

        W = LorentzFactor( Psi, V0 )

        D(iNodeX1,iX1) &
          = MassConstant / ( Psi**6 * Alpha * X1**2  * W * V0 )

        P(iNodeX1,iX1) &
          = PolytropicConstant * D(iNodeX1,iX1)**( Gamma_IDEAL )

      END DO

    END DO

    ! --- Map to 3D domain ---

    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
      DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
      DO iNodeX1 = 1, nNodesX(1)

        iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

        IF( ApplyPerturbation )THEN

          X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
          X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

          IF( X1 .GE. rPerturbationInner .AND. X1 .LE. rPerturbationOuter )THEN

            IF( PerturbationOrder .EQ. 0 ) &
              uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                = D(iNodeX1,iX1) * ( One + PerturbationAmplitude )

            IF( PerturbationOrder .EQ. 1 ) &
              uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                = D(iNodeX1,iX1) * ( One + PerturbationAmplitude * COS( X2 ) )

          ELSE

            uPF(iNodeX,iX1,iX2,iX3,iPF_D) = D(iNodeX1,iX1)

          END IF

        ELSE

          uPF(iNodeX,iX1,iX2,iX3,iPF_D) = D(iNodeX1,iX1)

        END IF

        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = V(iNodeX1,iX1)
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = Zero
        uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = P(iNodeX1,iX1)
        uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
          = uAF(iNodeX,iX1,iX2,iX3,iAF_P ) / ( Gamma_IDEAL - One )

        CALL ComputeConserved_Euler_Relativistic &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ),       &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V1),       &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V2),       &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V3),       &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ),       &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne),       &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_D ),       &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S1),       &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S2),       &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S3),       &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_E ),       &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_Ne),       &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_P ) )

      END DO
      END DO
      END DO
    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_StandingAccretionShock


  SUBROUTINE InitializeFields_StaticTOV

    REAL(DP), PARAMETER :: &
!!$      CentralDensity = 3.301e14_DP * ( Gram / Centimeter**3 ), &
      CentralDensity = 7.906e14_DP * ( Gram / Centimeter**3 ), &
      dX1            = 1.0e-4_DP * Kilometer, &
      TolF           = 1.0e-15_DP

    INTEGER  :: iX1, iX2, iX3, iNodeX, iNodeX1, iNodeX2, iNodeX3, &
                jNodeX, jNodeX1, iL, ITER, nX, iGF
    REAL(DP) :: X1, X2
    REAL(DP) :: Pressure , E1 , E2, Psi, Alpha, Phi
    REAL(DP) :: PressureN, E1N, E2N
    REAL(DP) :: CentralPressure, Psi0, Alpha0
    REAL(DP) :: GravitationalMass, Radius, dAlpha, dPsi, Alpha_A, Psi_A, dF

    REAL(DP), ALLOCATABLE :: PressureArr(:), DensityArr(:), &
                             AlphaArr(:), PsiArr(:), X1Arr(:)

    INTEGER, PARAMETER :: nMaxIter = 1000

    CentralPressure = PolytropicConstant_TOV * CentralDensity**( Gamma_IDEAL )
    Psi0            = 1.4_DP ! --- Initial guess ---
    Alpha0          = 0.8_DP ! --- Initial guess ---

    WRITE(*,*)
    WRITE(*,'(6x,A,ES10.3E3,A)' ) &
      'Polytropic Constant = ', PolytropicConstant_TOV &
                                  / ( Erg / Centimeter**3 &
                                  / ( Gram / Centimeter**3 )** &
                                    ( Gamma_IDEAL &
                                    ) ), &
      ' [ erg / cm^3 / ( g / cm^3 )^( Gamma ) ]'
    WRITE(*,'(6x,A,ES10.3E3,A)')  &
      'Central Density     = ', CentralDensity &
                                  / ( Gram / Centimeter**3 ), &
      ' [ g / cm^3 ]'
    WRITE(*,'(6x,A,ES10.3E3,A)')  &
      'Central Pressure    = ', CentralPressure &
                                  / ( Erg / Centimeter**3 ), &
      ' [ erg / cm^3 ]'
    WRITE(*,'(6x,A,ES10.3E3,A)')  &
      'dX                  = ', dX1 &
                                  / ( Kilometer ), &
      ' [ km ]'
    WRITE(*,'(6x,A,ES10.3E3)')    &
      'TolF                = ', TolF
    WRITE(*,*)

    ! --- Find geometry fields at center by iteratively integrating outward ---

    dF   = 1.1_DP * TolF
    ITER = 0
    DO WHILE( dF .GT. TolF .AND. ITER .LT. nMaxIter )

      ITER = ITER + 1

      IF( MOD( ITER, 100 ) .EQ. 0 ) PRINT*, 'Iteration ', ITER

      CALL IntegrateOutwards &
             ( dX1, CentralPressure, Psi0, Alpha0, &
               GravitationalMass, Radius, dAlpha, dPsi, Alpha_A, Psi_A, nX )

      ! --- Update guess for central values ---

      Alpha0 = Alpha0 + dAlpha
      Psi0   = Psi0   + dPsi

      dF = MAX( ABS( dAlpha / Alpha_A ), ABS( dPsi / Psi_A ) )

    END DO

    WRITE(*,'(6x,A,I4.4)') &
      'nIterations         = ', ITER
    WRITE(*,'(6x,A,ES13.6E3)'  )  &
      'dF                  = ', dF
    WRITE(*,'(6x,A,ES13.6E3)'  )  &
      'Alpha0              = ', Alpha0
    WRITE(*,'(6x,A,ES13.6E3)'  )  &
      'Psi0                = ', Psi0
    WRITE(*,'(6x,A,ES13.6E3,A)')  &
      'Radius              = ', Radius / Kilometer, ' km'
    WRITE(*,'(6x,A,ES13.6E3,A)')  &
      'Gravitational Mass  = ', GravitationalMass / SolarMass, ' Msun'

    ! --- Populate arrays ---

    ALLOCATE( PressureArr(nX), DensityArr(nX), AlphaArr(nX), &
              PsiArr(nX), X1Arr(nX) )

    ! --- Central values ---

    X1                = SqrtTiny * Kilometer
    Pressure          = CentralPressure
    E1                = Zero
    E2                = Zero
    Psi               = Psi0
    Phi               = Alpha0 * Psi0
    GravitationalMass = E1 / SpeedOfLight**2

    DO iX1 = 1, nX

      X1Arr      (iX1) = X1
      PressureArr(iX1) = Pressure
      DensityArr (iX1) = ( Pressure &
                             / PolytropicConstant_TOV )**( One / Gamma_IDEAL )

      ! --- Explicit steps ---

      PressureN = Pressure + dX1 * dpdr  ( Pressure, Phi, Psi, E1, E2, X1 )
      E1N       = E1       + dX1 * dE1dr ( Pressure, Psi, X1 )
      E2N       = E2       + dX1 * dE2dr ( Pressure, Phi, Psi, X1 )

      ! --- Implicit steps ---

      X1  = X1  + dX1
      Psi = Psi + dX1 * dPsidr( PressureN, E1N, X1 )
      Phi = Phi + dX1 * dPhidr( PressureN, E2N, X1 )

      PsiArr  (iX1) = Psi
      AlphaArr(iX1) = Phi / Psi

      Pressure = PressureN
      E1       = E1N
      E2       = E2N

    END DO

    ! --- Map to 3D domain ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E1(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

        iL = Locate( X1, X1Arr, nX )

        ! --- Geometry Fields ---

        uGF(iNodeX,iX1,iX2,iX3,iGF_Alpha) &
          = Interpolate1D_Linear( X1, X1Arr(iL), X1Arr(iL+1), &
                                  AlphaArr(iL), AlphaArr(iL+1) )

        uGF(iNodeX,iX1,iX2,iX3,iGF_Psi) &
          = Interpolate1D_Linear( X1, X1Arr(iL), X1Arr(iL+1), &
                                  PsiArr(iL), PsiArr(iL+1) )

        uGF(iNodeX,iX1,iX2,iX3,iGF_h_1) &
          = uGF(iNodeX,iX1,iX2,iX3,iGF_Psi)**2
        uGF(iNodeX,iX1,iX2,iX3,iGF_h_2) &
          = uGF(iNodeX,iX1,iX2,iX3,iGF_Psi)**2 * X1
        uGF(iNodeX,iX1,iX2,iX3,iGF_h_3) &
          = uGF(iNodeX,iX1,iX2,iX3,iGF_Psi)**2 * X1 * SIN( X2 )

        uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11) &
          = uGF(iNodeX,iX1,iX2,iX3,iGF_h_1)**2
        uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22) &
          = uGF(iNodeX,iX1,iX2,iX3,iGF_h_2)**2
        uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) &
          = uGF(iNodeX,iX1,iX2,iX3,iGF_h_3)**2

        uGF(iNodeX,iX1,iX2,iX3,iGF_SqrtGm) &
          = SQRT( uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11) &
                    * uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22) &
                    * uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

        uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_1) = Zero
        uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_2) = Zero
        uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_3) = Zero

        ! --- Fluid Fields ---

        uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
          = Interpolate1D_Linear( X1, X1Arr(iL), X1Arr(iL+1), &
                                  DensityArr(iL), DensityArr(iL+1) )

        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = Zero

        uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
          = Interpolate1D_Linear &
              ( X1, X1Arr(iL), X1Arr(iL+1), &
                PressureArr(iL), PressureArr(iL+1) ) / ( Gamma_IDEAL - One )

      END DO

      CALL ComputePressureFromPrimitive_IDEAL &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_E), &
               uPF(:,iX1,iX2,iX3,iPF_Ne), uAF(:,iX1,iX2,iX3,iAF_P) )

      CALL ComputeConserved_Euler_Relativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uAF(:,iX1,iX2,iX3,iAF_P) )

    END DO
    END DO
    END DO

    ! --- Apply reflecting boundary conditions to geometry fields (X1) ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = 1, swX(1)

        DO iNodeX3 = 1, nNodesX(3)
        DO iNodeX2 = 1, nNodesX(2)
        DO iNodeX1 = 1, nNodesX(1)

          jNodeX1 = ( nNodesX(1) - iNodeX1 ) + 1

          iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )
          jNodeX = NodeNumberX( jNodeX1, iNodeX2, iNodeX3 )

          DO iGF = 1, nGF

            uGF(iNodeX,iX_B0(1)-iX1,iX2,iX3,iGF) &
              = uGF(jNodeX,iX_B0(1),iX2,iX3,iGF)

          END DO

        END DO
        END DO
        END DO

    END DO
    END DO
    END DO

    DEALLOCATE( X1Arr, PsiArr, AlphaArr, DensityArr, PressureArr )

  END SUBROUTINE InitializeFields_StaticTOV


  SUBROUTINE InitializeFields_YahilCollapse &
    ( CentralDensity, CentralPressure, CoreRadius, CollapseTime )

    REAL(DP), INTENT(in) :: CentralDensity
    REAL(DP), INTENT(in) :: CentralPressure
    REAL(DP), INTENT(in) :: CoreRadius
    REAL(DP), INTENT(in) :: CollapseTime

    CHARACTER(LEN=64) :: FileName
    INTEGER           :: nLines

    LOGICAL               :: InitializeFromFile
    INTEGER               :: iX1, iX2, iX3, iNodeX, iNodeX1, iX_L
    REAL(DP)              :: K, C_X, C_D, C_V, C_M, R, XX
    REAL(DP), ALLOCATABLE :: X(:), D(:), V(:), M(:)

    InitializeFromFile = .TRUE.

    K = CentralPressure / CentralDensity**( Gamma_IDEAL )

    WRITE(*,*)
    WRITE(*,'(6x,A,L)') &
      'InitializeFromFile:  ', &
      InitializeFromFile
    WRITE(*,*)
    WRITE(*,'(6x,A,F5.3)') &
      'Adiabatic Gamma:     ', &
      Gamma_IDEAL
    WRITE(*,'(6x,A,ES10.3E3,A)') &
      'Central Density:     ', &
      CentralDensity / ( Gram / Centimeter**3  ), ' g/cm^3'
    WRITE(*,'(6x,A,ES10.3E3,A)') &
      'Central Pressure:    ', &
      CentralPressure / ( Erg / Centimeter**3  ), ' erg/cm^3'
    WRITE(*,'(6x,A,ES10.3E3,A)') &
      'Polytropic Constant: ', &
      K / ( Erg / Centimeter**3 &
              / ( Gram / Centimeter**3  )**( Gamma_IDEAL ) ), &
      ' ( erg/cm^3 ) / ( g/cm^3 )^( Gamma )'
    WRITE(*,'(6x,A,ES10.3E3,A)') &
      'Core Radius:         ', &
      CoreRadius / Kilometer, ' km'
    WRITE(*,'(6x,A,ES10.3E3,A)') &
      'Collapse Time:       ', &
      CollapseTime / Millisecond, ' ms'

    C_X = K**( -Half ) &
            * GravitationalConstant**( ( Gamma_IDEAL - One ) / Two ) &
            * CollapseTime**( Gamma_IDEAL - Two )
    C_D = GravitationalConstant**( -1 ) * CollapseTime**( -2 )
    C_V = K**( Half ) &
            * GravitationalConstant**( ( One - Gamma_IDEAL ) / Two ) &
            * CollapseTime**( One - Gamma_IDEAL )
    C_M = K**( Three / Two ) &
            * GravitationalConstant**( ( One - Three * Gamma_IDEAL ) / Two ) &
            * CollapseTime**( Four - Three * Gamma_IDEAL )

    IF( InitializeFromFile )THEN

      FileName = 'YahilHomologousCollapse_Gm_130.dat'

      ! --- https://stackoverflow.com/questions/30692424/
      !     how-to-read-number-of-lines-in-fortran-90-from-a-text-file
      nLines = 0
      OPEN(100,FILE=TRIM(FileName))
      READ(100,*)
      DO
        READ(100,*,END=10)
        nLines = nLines + 1
      END DO
      10 CLOSE(100)

      ALLOCATE( X(nLines) )
      ALLOCATE( D(nLines) )
      ALLOCATE( V(nLines) )
      ALLOCATE( M(nLines) )

      OPEN(100,FILE=TRIM(FileName))
      READ(100,*)

      DO iX1 = 1, nLines

        READ(100,*) X(iX1), D(iX1), V(iX1), M(iX1)

      END DO

      CLOSE(100)

      WRITE(*,'(6x,A,ES10.3E3,A)') &
        'Mass:                ', &
        M(nLines) * C_M / SolarMass, ' Msun'
      WRITE(*,*)

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E1(1)

       DO iNodeX = 1, nDOFX

         iNodeX1 = NodeNumberTableX(1,iNodeX)

         R = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
         XX = C_X * R

         iX_L = Locate( XX, X, nLines )

         uPF(iNodeX,iX1,iX2,iX3,iPF_D ) &
           = C_D * Interpolate1D_Linear( XX, X(iX_L), X(iX_L+1), &
                                         D(iX_L), D(iX_L+1) )

         uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
           = C_V * Interpolate1D_Linear( XX, X(iX_L), X(iX_L+1), &
                                         V(iX_L), V(iX_L+1) )

         uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = Zero

         uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = Zero

         uPF(iNodeX,iX1,iX2,iX3,iPF_E ) &
           = K * uPF(iNodeX,iX1,iX2,iX3,iPF_D)**( Gamma_IDEAL ) &
               / ( Gamma_IDEAL - One )

       END DO

        CALL ComputePressureFromPrimitive_IDEAL &
               ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_E), &
                 uPF(:,iX1,iX2,iX3,iPF_Ne), uAF(:,iX1,iX2,iX3,iAF_P) )

        CALL ComputeConserved_Euler_Relativistic &
               ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
                 uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
                 uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
                 uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
                 uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
                 uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
                 uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
                 uAF(:,iX1,iX2,iX3,iAF_P) )

      END DO
      END DO
      END DO

      DEALLOCATE( M )
      DEALLOCATE( V )
      DEALLOCATE( D )
      DEALLOCATE( X )

    ELSE
    END IF

  END SUBROUTINE InitializeFields_YahilCollapse

!!$  SUBROUTINE InitializeFields_YahilCollapse &
!!$    ( CentralDensity, CentralPressure, CoreRadius, CollapseTime )
!!$
!!$    REAL(DP), INTENT(in) :: CentralDensity
!!$    REAL(DP), INTENT(in) :: CentralPressure
!!$    REAL(DP), INTENT(in) :: CoreRadius
!!$    REAL(DP), INTENT(in) :: CollapseTime
!!$
!!$    INTEGER               :: N
!!$    REAL(DP)              :: dr, dX, K, C_X, C_D, C_V, C_M, dD
!!$    REAL(DP), ALLOCATABLE :: R(:), X(:), D(:), U(:), V(:), M(:)
!!$
!!$    INTEGER :: iX1
!!$real(dp)::constraint
!!$
!!$    K = CentralPressure / CentralDensity**( Gamma_IDEAL )
!!$
!!$    WRITE(*,*)
!!$    WRITE(*,'(6x,A,F5.3)') &
!!$      'Adiabatic Gamma:     ', &
!!$      Gamma_IDEAL
!!$    WRITE(*,'(6x,A,ES10.3E3,A)') &
!!$      'Central Density:     ', &
!!$      CentralDensity / ( Gram / Centimeter**3  ), ' g/cm^3'
!!$    WRITE(*,'(6x,A,ES10.3E3,A)') &
!!$      'Central Pressure:    ', &
!!$      CentralPressure / ( Erg / Centimeter**3  ), ' erg/cm^3'
!!$    WRITE(*,'(6x,A,ES10.3E3,A)') &
!!$      'Polytropic Constant: ', &
!!$      K / ( Erg / Centimeter**3 &
!!$              / ( Gram / Centimeter**3  )**( Gamma_IDEAL ) ), &
!!$      ' ( erg/cm^3 ) / ( g/cm^3 )^( Gamma )'
!!$    WRITE(*,'(6x,A,ES10.3E3,A)') &
!!$      'Core Radius:         ', &
!!$      CoreRadius / Kilometer, ' km'
!!$    WRITE(*,'(6x,A,ES10.3E3,A)') &
!!$      'Collapse Time:       ', &
!!$      CollapseTime / Millisecond, ' ms'
!!$    WRITE(*,*)
!!$
!!$    dr = 1.0e0_DP * Kilometer
!!$    N = CoreRadius / dr
!!$
!!$    C_X = K**( -Half ) &
!!$            * GravitationalConstant**( ( Gamma_IDEAL - One ) / Two ) &
!!$            * CollapseTime**( Gamma_IDEAL - Two )
!!$    C_D = GravitationalConstant**( -1 ) * CollapseTime**( -2 )
!!$    C_V = K**( Half ) &
!!$            * GravitationalConstant**( ( One - Gamma_IDEAL ) / Two ) &
!!$            * CollapseTime**( One - Gamma_IDEAL )
!!$    C_M = K**( Three / Two ) &
!!$            * GravitationalConstant**( ( One - Three * Gamma_IDEAL ) / Two ) &
!!$            * CollapseTime**( Four - Three * Gamma_IDEAL )
!!$
!!$    dX = C_X * dr
!!$
!!$    ALLOCATE( R(N) )
!!$    ALLOCATE( X(N) )
!!$    ALLOCATE( D(N) )
!!$    ALLOCATE( U(N) )
!!$    ALLOCATE( V(N) )
!!$    ALLOCATE( M(N) )
!!$
!!$    R(1) = SqrtTiny * Kilometer
!!$    X(1) = C_X * R(1)
!!$    D(1) = CentralDensity / C_D
!!$    U(1) = Zero!-( Gamma_IDEAL - Two ) * X(1)
!!$    M(1) = Zero
!!$    DO iX1 = 2, N
!!$
!!$      IF( MOD( iX1, 100 ) .EQ. 0 )THEN
!!$
!!$        WRITE(*,'(F4.2)') DBLE(iX1) / DBLE(N)
!!$
!!$      END IF
!!$
!!$      X(iX1) = X(iX1-1) + dX
!!$
!!$      dD     = dDdX( X(iX1-1), D(iX1-1), U(iX1-1), M(iX1-1) )
!!$
!!$      D(iX1) = D(iX1-1) + dX * dD
!!$
!!$      U(iX1) = U(iX1-1) + dX * dUdX( X(iX1-1), D(iX1-1), U(iX1-1), dD )
!!$
!!$      CALL IntegrateM( X(1:iX1), dX, D(1:iX1), M(iX1) )
!!$
!!$constraint = FourPi * X(iX1)**2 * D(iX1) * ABS( U(iX1) ) &
!!$               / ( ( Four - Three * Gamma_IDEAL ) * M(iX1) ) - One
!!$
!!$if(X(iX1).gt.1000.0d0)stop!!$if(abs(constraint).gt.1.0e-10_DP)then
!!$print*, 'X: ', X(iX1)
!!$print*, 'Radius: ', X(iX1) / C_X / Kilometer, ' km'
!!$print*, 'Mass:   ', C_M * M(iX1) / SolarMass, ' Msun'
!!$print*,'Constraint: ', Constraint
!!$print*
!!$!!$stop
!!$!!$endif
!!$
!!$    END DO
!!$
!!$    DEALLOCATE( M )
!!$    DEALLOCATE( V )
!!$    DEALLOCATE( U )
!!$    DEALLOCATE( D )
!!$    DEALLOCATE( X )
!!$    DEALLOCATE( R )
!!$
!!$stop
!!$
!!$  END SUBROUTINE InitializeFields_YahilCollapse


  ! --- Auxiliary functions for Yahil collapse problem ---

  REAL(DP) FUNCTION dDdX( X, D, U, M )

    REAL(DP), INTENT(in) :: X, D, U, M

    dDdX = ( -M / X**2 + Two * U**2 / X + ( Gamma_IDEAL - One ) * U &
             + ( Gamma_IDEAL - One ) * ( Two - Gamma_IDEAL ) * X ) &
             / ( Gamma_IDEAL * D**( Gamma_IDEAL - One ) - U**2 ) / D

  END FUNCTION dDdX


  REAL(DP) FUNCTION dUdX( X, D, U, dDdX )

    REAL(DP), INTENT(in) :: X, D, U, dDdX

    dUdX = Four - Three * Gamma_IDEAL - Two * U / X - U / D * dDdX

  END FUNCTION dUdX


  SUBROUTINE IntegrateM( X, dX, D, M )

    REAL(DP), INTENT(in)  :: X(:), dX, D(:)
    REAL(DP), INTENT(out) :: M

    INTEGER :: iX1, N

    N = SIZE( D )

    M = Half * ( D(1) * X(1)**2 + D(N) * X(N)**2 )

    IF( N .GT. 2 )THEN

      DO iX1 = 2, N-1

        M = M + D(iX1) * X(iX1)**2

      END DO

    END IF

    M = FourPi * dX * M

  END SUBROUTINE IntegrateM


  ! --- End of auxiliary functions for Yahil collapse problem ---


  ! --- Auxiliary utilities for TOV problem ---


  SUBROUTINE IntegrateOutwards &
    ( dX1, CentralPressure, Psi0, Alpha0, &
      GravitationalMass, Radius, dAlpha, dPsi, Alpha_A, Psi_A, nX )

    REAL(DP), INTENT(in)  :: dX1, CentralPressure, Psi0, Alpha0
    REAL(DP), INTENT(out) :: GravitationalMass, Radius, &
                             dAlpha, dPsi, Alpha_A, Psi_A
    INTEGER,  INTENT(out) :: nX

    REAL(DP) :: Pressure , E1 , E2 , Psi , Phi, X1
    REAL(DP) :: PressureN, E1N, E2N
    REAL(DP) :: Alpha

    ! --- Set inner boundary values ---

    X1                = SqrtTiny * Kilometer
    Pressure          = CentralPressure
    E1                = Zero
    E2                = Zero
    Psi               = Psi0
    Phi               = Alpha0 * Psi0
    GravitationalMass = E1 / SpeedOfLight**2

    nX = 1

    DO WHILE( Pressure .GT. 1.0e-8_DP * CentralPressure )

      ! --- Explicit steps ---

      PressureN = Pressure + dX1 * dpdr  ( Pressure, Phi, Psi, E1, E2, X1 )
      E1N       = E1       + dX1 * dE1dr ( Pressure, Psi, X1 )
      E2N       = E2       + dX1 * dE2dr ( Pressure, Phi, Psi, X1 )

      ! --- Implicit steps ---

      X1  = X1  + dX1
      Psi = Psi + dX1 * dPsidr( PressureN, E1N, X1 )
      Phi = Phi + dX1 * dPhidr( PressureN, E2N, X1 )

      Pressure = PressureN
      E1       = E1N
      E2       = E2N

      nX = nX + 1

    END DO

    GravitationalMass = E1 / SpeedOfLight**2
    Radius            = X1

    Alpha   &
      = Phi / Psi

    Alpha_A &
      =  ( One - GravitationalMass / ( Two * SpeedOfLight**2 * Radius ) ) &
       / ( One + GravitationalMass / ( Two * SpeedOfLight**2 * Radius ) )

    Psi     &
      = Psi

    Psi_A   &
      = One + GravitationalMass / ( Two * SpeedOfLight**2 * Radius )

    dAlpha = Alpha_A - Alpha
    dPsi   = Psi_A - Psi

  END SUBROUTINE IntegrateOutwards


  REAL(DP) FUNCTION dpdr( Pressure, Phi, Psi, E1, E2, X1  )

    REAL(DP), INTENT(in) :: Pressure, Phi, Psi, E1, E2, X1

    REAL(DP) :: Lapse

    Lapse = Phi / Psi

    dpdr = - Enthalpy( Pressure ) &
             * ( dPhidr( Pressure, E2, X1 ) &
                   - Lapse * dPsidr( Pressure, E1, X1 ) ) / Phi

    RETURN
  END FUNCTION dpdr


  REAL(DP) FUNCTION Enthalpy( Pressure )

    REAL(DP), INTENT(in) :: Pressure

    Enthalpy = ( Pressure / PolytropicConstant_TOV )**( One / Gamma_IDEAL ) &
                 * SpeedOfLight**2 + Pressure / ( Gamma_IDEAL - One ) + Pressure

    RETURN
  END FUNCTION Enthalpy


  REAL(DP) FUNCTION dPhidr( Pressure, E2, X1 )

    REAL(DP), INTENT(in) :: Pressure, E2, X1

    dPhidr = GravitationalConstant / ( Two * SpeedOfLight**4 ) * E2 / X1**2

    RETURN
  END FUNCTION dPhidr


  REAL(DP) FUNCTION dPsidr( Pressure, E1, X1 )

    REAL(DP), INTENT(in) :: Pressure, E1, X1

    dPsidr = -GravitationalConstant / ( Two * SpeedOfLight**4 ) * E1 / X1**2

    RETURN
  END FUNCTION dPsidr


  REAL(DP) FUNCTION dE1dr( Pressure, Psi, X1 )

    REAL(DP), INTENT(in) :: Pressure, Psi, X1

    dE1dr = FourPi * X1**2 * f_E1( Pressure ) * Psi**5

    RETURN
  END FUNCTION dE1dr


  REAL(DP) FUNCTION dE2dr( Pressure, Phi, Psi, X1 )

    REAL(DP), INTENT(in) :: Pressure, Phi, Psi, X1

    dE2dr = FourPi * X1**2 * f_E2( Pressure ) * Phi * Psi**4

    RETURN
  END FUNCTION dE2dr


  REAL(DP) FUNCTION f_E1( Pressure )

    REAL(DP), INTENT(in) :: Pressure

    f_E1 = Enthalpy( Pressure ) - Pressure

    RETURN
  END FUNCTION f_E1


  REAL(DP) FUNCTION f_E2( Pressure )

    REAL(DP), INTENT(in) :: Pressure

    f_E2 = f_E1( Pressure ) - 6.0_DP * Pressure

    RETURN
  END FUNCTION f_E2


  ! --- Auxiliary utilities for standing accretion shock problem ---


  SUBROUTINE ApplyJumpConditions &
    ( iX1_1, iNodeX1_1, X1_1, D_1, V_1, &
      iX1_2, iNodeX1_2, X1_2, &
      D_2, V_2, P_2, MassPNS, PolytropicConstant )

    INTEGER,  INTENT(in)  :: iX1_1, iNodeX1_1, iX1_2, iNodeX1_2
    REAL(DP), INTENT(in)  :: X1_1, X1_2, D_1, V_1, MassPNS
    REAL(DP), INTENT(out) :: D_2, V_2, P_2, PolytropicConstant

    REAL(DP) :: Alpha, Psi
    REAL(DP) :: C1, C2, C3, a0, a1, a2, a3, a4, X1
    REAL(DP) :: W

    REAL(DP), PARAMETER :: ShockTolerance = 0.1_DP
    LOGICAL             :: FoundShockVelocity = .FALSE.

    ! --- Constants from three jump conditions ---

    Alpha = LapseFunction  ( X1_1, MassPNS )
    Psi   = ConformalFactor( X1_1, MassPNS )

    C1 = D_1 * V_1 / Alpha

    C2 = D_1 * SpeedOfLight**2 / Alpha**2 * ( V_1 / SpeedOfLight )**2

    C3 = D_1 * SpeedOfLight**2 / Alpha**2 * V_1

    ! --- Five constants for post-shock fluid-velocity ---

    a4 = Psi**8 &
          * One / ( Gamma_IDEAL - One )**2 * C3**2 / SpeedOfLight**6
    a3 = -Two * Psi**8 &
          * Gamma_IDEAL / ( Gamma_IDEAL - One )**2 * C2 * C3 / SpeedOfLight**4
    a2 = Psi**4 &
          / SpeedOfLight**2 * ( Psi**4 * Gamma_IDEAL**2 &
          / ( Gamma_IDEAL - One )**2 * C2**2 + Two * One &
          / ( Gamma_IDEAL - One ) &
          * C3**2 / SpeedOfLight**2 + C1**2 * SpeedOfLight**2 )
    a1 = -Two * Psi**4 &
          * Gamma_IDEAL / ( Gamma_IDEAL - One ) * C2 * C3 / SpeedOfLight**2
    a0 = One / SpeedOfLight**2 * ( C3**2 - C1**2 * SpeedOfLight**4 )

    ! --- Newton-Raphson method for post-shock fluid-velocity ---

    V_2 = Two * V_1

    ! --- Ensure that shocked velocity is obtained ---

    FoundShockVelocity = .FALSE.
    DO WHILE( .NOT. FoundShockVelocity )

      V_2 = Half * V_2
      CALL NewtonRaphson_JumpConditions( a0, a1, a2, a3, a4, V_2 )
!!$      CALL Bisection_JumpConditions( a0, a1, a2, a3, a4, V_2 )

      IF( ABS( V_2 - V_1 ) / ABS( V_1 ) .GT. ShockTolerance ) &
        FoundShockVelocity = .TRUE.

    END DO

    ! --- Post-shock density, velocity, pressure, and polytropic constant ---

    Psi = ConformalFactor( X1_2, MassPNS )
    W   = LorentzFactor( Psi, V_2 )

    D_2 = ABS( C1 ) * SQRT( One / V_2**2 - Psi**4 / SpeedOfLight**2 )

    P_2 = ( Gamma_IDEAL - One ) / Gamma_IDEAL &
            * ( C3 - D_2 * SpeedOfLight**2 * W**2 * V_2 ) / ( W**2 * V_2 )

    PolytropicConstant = P_2 / D_2**( Gamma_IDEAL )

  END SUBROUTINE ApplyJumpConditions


  SUBROUTINE Bisection_JumpConditions( a0, a1, a2, a3, a4, V )

    REAL(DP), INTENT(in)    :: a0, a1, a2, a3, a4
    REAL(DP), INTENT(inout) :: V

    REAL(DP) :: F, dV, Vmin, Vmax, Fmin, Fmax
    LOGICAL  :: CONVERGED
    INTEGER  :: ITERATION

    INTEGER,  PARAMETER :: nMaxIter = 100
    REAL(DP), PARAMETER :: ToldV = 1.0d-15
    REAL(DP), PARAMETER :: EPS = 1.0e-15_DP

    IF( V .LT. Zero )THEN

      Vmin = V    + EPS
      Vmax = +One - EPS

    ELSE

      Vmin = -One + EPS
      Vmax = V    - EPS

    END IF

    Fmin = a4 * Vmin**4 + a3 * Vmin**3 + a2 * Vmin**2 + a1 * Vmin + a0
    Fmax = a4 * Vmax**4 + a3 * Vmax**3 + a2 * Vmax**2 + a1 * Vmax + a0

    IF( .NOT. Fmin * Fmax .LT. Zero )THEN

      WRITE(*,*) 'Root not bracketed. Stopping...'
      WRITE(*,*) 'Fmin = ', Fmin
      WRITE(*,*) 'Fmax = ', Fmax
      STOP

    END IF

    IF( Fmin .GT. Zero )THEN

      V = Vmax
      F = Fmax

      Vmax = Vmin
      Vmin = V

      Fmax = Fmin
      Fmin = F

    END IF

    ITERATION = 0
    DO WHILE( ITERATION .LT. nMaxIter )

      ITERATION = ITERATION + 1

      V = ( Vmin + Vmax ) / Two

      F = a4 * V**4 + a3 * V**3 + a2 * V**2 + a1 * V + a0

      IF( ABS( V - Vmin ) / MAX( ABS( Vmax ), ABS( Vmin ) ) .LT. ToldV ) EXIT

      IF( F .GT. Zero )THEN

        Vmax = V
        Fmax = F

     ELSE

        Vmin = V
        Fmin = F

     END IF

    END DO

!!$    WRITE(*,*) 'Converged at iteration ', ITERATION
!!$    WRITE(*,*) '|F|:  ' , ABS( F )
!!$    WRITE(*,*) 'dV/V: ', ABS( V - Vmax ) / ABS( Vmax )

  END SUBROUTINE Bisection_JumpConditions


  SUBROUTINE NewtonRaphson_JumpConditions( a0, a1, a2, a3, a4, V )

    REAL(DP), INTENT(in)    :: a0, a1, a2, a3, a4
    REAL(DP), INTENT(inout) :: V

    REAL(DP) :: f, df, dV
    LOGICAL  :: CONVERGED
    INTEGER  :: ITERATION

    INTEGER,  PARAMETER :: MAX_ITER = 10
    REAL(DP), PARAMETER :: TOLERANCE = 1.0d-15

    CONVERGED = .FALSE.
    ITERATION = 0
    DO WHILE( .NOT. CONVERGED .AND. ITERATION .LT. MAX_ITER )

      ITERATION = ITERATION + 1

      f  = a4 * V**4 + a3 * V**3 + a2 * V**2 + a1 * V + a0
      df = Four * a4 * V**3 + Three * a3 * V**2 + Two * a2 * V + a1

      dV = -f / df
      V = V + dV

      IF( ABS( dV / V ) .LT. TOLERANCE )THEN
        CONVERGED = .TRUE.
      END IF

    END DO

  END SUBROUTINE NewtonRaphson_JumpConditions


  SUBROUTINE NewtonRaphson_PostShockVelocity &
    ( Alpha, Psi, MassConstant, PolytropicConstant, &
      MassPNS, AccretionRate, X1, V )

    REAL(DP), INTENT(in)    :: Alpha, Psi, MassConstant, &
                               PolytropicConstant, MassPNS, AccretionRate, X1
    REAL(DP), INTENT(inout) :: V

    REAL(DP) :: f, df, dV, W
    INTEGER  :: ITERATION
    LOGICAL  :: CONVERGED

    INTEGER,  PARAMETER :: MAX_ITER = 20
    REAL(DP), PARAMETER :: TOLERANCE = 1.0d-15

    CONVERGED = .FALSE.
    ITERATION = 0
    DO WHILE( .NOT. CONVERGED .AND. ITERATION .LT. MAX_ITER )

      ITERATION = ITERATION + 1

      W = LorentzFactor( Psi, V )

      f  = Gamma_IDEAL / ( Gamma_IDEAL - One ) &
             * PolytropicConstant / SpeedOfLight**2 * ( MassConstant &
             / ( Psi**6 * Alpha * X1**2 * W * V ) )**( Gamma_IDEAL - One ) &
             - One / ( Alpha * W ) + One

      df = -Gamma_IDEAL * PolytropicConstant / SpeedOfLight**2 &
             * ( MassConstant &
                 / ( Psi**6 * Alpha * X1**2 * W * V ) )**( Gamma_IDEAL - One ) &
                 * ( Psi**4 * V / SpeedOfLight**2 * W**2 + One / V ) &
                 + W / Alpha * Psi**4 * V / SpeedOfLight**2

      dV = -f / df
      V = V + dV

      IF( ABS( dV / V ) .LT. TOLERANCE ) &
        CONVERGED = .TRUE.

    END DO

  END SUBROUTINE NewtonRaphson_PostShockVelocity


  REAL(DP) FUNCTION LapseFunction( R, M )

    REAL(DP), INTENT(in) :: R, M

    ! --- Schwarzschild Metric in Isotropic Coordinates ---

    LapseFunction = ABS( ( MAX( ABS( R ), SqrtTiny ) - Half * M ) &
                       / ( MAX( ABS( R ), SqrtTiny ) + Half * M ) )

    RETURN
  END FUNCTION LapseFunction


  REAL(DP) FUNCTION ConformalFactor( R, M )

    REAL(DP), INTENT(in) :: R, M

    ! --- Schwarzschild Metric in Isotropic Coordinates ---

    ConformalFactor = One + Half * M / MAX( ABS( R ), SqrtTiny )

    RETURN
  END FUNCTION ConformalFactor


  REAL(DP) FUNCTION LorentzFactor( Psi, V )

    REAL(DP), INTENT(in) :: Psi, V

    LorentzFactor = One / SQRT( One - Psi**4 * ( V / SpeedOfLight )**2 )

    RETURN
  END FUNCTION LorentzFactor


  ! --- Auxiliary functions/subroutines for computine left state ---


  SUBROUTINE ComputeLeftState( Vs, DR, VR, PR, DL, VL, PL )

    REAL(DP), INTENT(in)  :: Vs, DR, VR, PR
    REAL(DP), INTENT(out) ::     DL, VL, PL

    CALL ApplyJumpConditions_LeftState( Vs, DR, VR, PR, DL, VL, PL )

    ! --- Return energy-density instead of pressure ---
    PL = PL / ( Gamma_IDEAL - One )

  END SUBROUTINE ComputeLeftState


  SUBROUTINE ApplyJumpConditions_LeftState( Vs, DR, VR, PR, DL, VL, PL )

    REAL(DP), INTENT(in)  :: Vs, DR, VR, PR
    REAL(DP), INTENT(out) ::     DL, VL, PL

    REAL(DP), PARAMETER :: EPS = 1.0e-15_DP

    REAL(DP), PARAMETER :: ToldV = EPS
    REAL(DP), PARAMETER :: TolF  = EPS
    INTEGER,  PARAMETER :: nMaxIter = 1000

    INTEGER :: ITERATION
    REAL(DP) :: D, V, P, F
    REAL(DP) :: Vmin, Vmax, Fmin, Fmax, VV, FF

    IF( VR .LT. Zero )THEN

      Vmin = VR   + EPS
      Vmax = +One - EPS

    ELSE

      Vmin = -One + EPS
      Vmax = VR   - EPS

    END IF

    D = Density ( Vs, DR, VR, Vmin )
    P = Pressure( Vs, DR, VR, PR, D, Vmin )
    Fmin = PostShockVelocity( Vs, DR, VR, PR, D, Vmin, P )

    D = Density( Vs, DR, VR, Vmax )
    P = Pressure( Vs, DR, VR, PR, D, Vmax )
    Fmax = PostShockVelocity( Vs, DR, VR, PR, D, Vmax, P )

    IF( .NOT. Fmin * Fmax .LT. Zero )THEN

      WRITE(*,*) 'Root not bracketed. Stopping...'
      WRITE(*,*) 'Fmin = ', Fmin
      WRITE(*,*) 'Fmax = ', Fmax
      STOP

    END IF

    IF( Fmin .GT. Zero )THEN

      VV = Vmax
      FF = Fmax

      Vmax = Vmin
      Vmin = VV

      Fmax = Fmin
      Fmin = FF

    END IF

    ITERATION = 0
    DO WHILE( ITERATION .LT. nMaxIter )

      ITERATION = ITERATION + 1

      V = ( Vmin + Vmax ) / Two

      D = Density ( Vs, DR, VR, V )
      P = Pressure( Vs, DR, VR, PR, D, V )

      F = PostShockVelocity( Vs, DR, VR, PR, D, V, P )

      IF( ABS( V - Vmin ) / MAX( ABS( Vmax ), ABS( Vmin ) ) .LT. ToldV ) EXIT

      IF( F .GT. Zero )THEN

        Vmax = V
        Fmax = F

     ELSE

        Vmin = V
        Fmin = F

     END IF

    END DO

!!$    WRITE(*,*) 'Converged at iteration ', ITERATION
!!$    WRITE(*,*) '|F|:  ' , ABS( F )
!!$    WRITE(*,*) 'dV/V: ', ABS( V - Vmax ) / ABS( Vmax )

    VL = V
    DL = Density ( Vs, DR, VR, VL )
    PL = Pressure( Vs, DR, VR, PR, DL, VL )

  END SUBROUTINE ApplyJumpConditions_LeftState


  REAL(DP) FUNCTION Density( Vs, DR, VR, VL )

    REAL(DP), INTENT(in) :: Vs, DR, VR, VL

    REAL(DP) :: WR, WL

    WR = LorentzFactor( One, VR )
    WL = LorentzFactor( One, VL )

    Density = DR * ( WR * ( VR - Vs ) ) / ( WL * ( VL - Vs ) )

    RETURN
  END FUNCTION Density


  REAL(DP) FUNCTION Pressure( Vs, DR, VR, PR, DL, VL )

    REAL(DP), INTENT(in) :: Vs, DR, VR, PR, DL, VL

    REAL(DP) :: WR, WL, tau

    WR = LorentzFactor( One, VR )
    WL = LorentzFactor( One, VL )

    tau = Gamma_IDEAL / ( Gamma_IDEAL - One )

    Pressure = ( PR * ( One + tau * WR**2 * VR * ( VR - Vs ) ) &
                 - DL * WL**2 * VL**2 + DR * WR**2 * VR**2 &
                 + Vs * ( DL * WL**2 * VL - DR * WR**2 * VR ) ) &
               / ( One + tau * WL**2 * VL * ( VL - Vs ) )

    RETURN
  END FUNCTION Pressure


  REAL(DP) FUNCTION PostShockVelocity( Vs, DR, VR, PR, DL, VL, PL )

    REAL(DP), INTENT(in) :: Vs, DR, VR, PR, DL, VL, PL

    REAL(DP) :: WR, WL, tau

    WR = LorentzFactor( One, VR )
    WL = LorentzFactor( One, VL )

    tau = Gamma_IDEAL / ( Gamma_IDEAL - One )

    PostShockVelocity &
      = ( DL + tau * PL ) * WL**2 * ( VL - Vs ) &
          - ( DR + tau * PR ) * WR**2 * ( VR - Vs ) + Vs * ( PL - PR )

    RETURN
  END FUNCTION PostShockVelocity


END MODULE InitializationModule_Relativistic
