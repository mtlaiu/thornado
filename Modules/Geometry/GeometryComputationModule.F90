MODULE GeometryComputationModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    NodesLX_q
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    LX_L2G
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    CoordinateSystem, &
    iGF_Phi_N, &
    iGF_h_1,      iGF_h_2,      iGF_h_3,      &
    iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33, &
    iGF_Beta_1,   iGF_Beta_2,   iGF_Beta_3,   &
    iGF_SqrtGm, &
    iGF_Alpha, iGF_Psi, &
    nGF

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeGeometryX
  PUBLIC :: ComputeGeometryX_FromScaleFactors
  PUBLIC :: LapseFunction
  PUBLIC :: ConformalFactor


CONTAINS


  SUBROUTINE ComputeGeometryX( iX_B0, iX_E0, iX_B1, iX_E1, G, Mass_Option )

    INTEGER, INTENT(in)            :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout)        :: &
      G(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(in), OPTIONAL :: &
      Mass_Option

    REAL(DP) :: Mass

    IF( PRESENT( Mass_Option ) )THEN
      Mass = Mass_Option
    ELSE
      Mass = Zero
    END IF

    SELECT CASE ( TRIM( CoordinateSystem ) )

      CASE ( 'CARTESIAN' )

        CALL ComputeGeometryX_CARTESIAN &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mass )

      CASE ( 'SPHERICAL' )

        CALL ComputeGeometryX_SPHERICAL &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mass )

      CASE ( 'CYLINDRICAL' )

        CALL ComputeGeometryX_CYLINDRICAL &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mass )

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A5,A27,A)') &
          '', 'Invalid Coordinate System: ', TRIM( CoordinateSystem )
        STOP

    END SELECT

  END SUBROUTINE ComputeGeometryX


  SUBROUTINE ComputeGeometryX_CARTESIAN &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mass )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      G(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(in)    :: &
      Mass

    INTEGER :: iX1, iX2, iX3

    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)

      ! --- Initialize to flat spacetime ---
      G(:,iX1,iX2,iX3,iGF_Phi_N)    = Zero
      G(:,iX1,iX2,iX3,iGF_h_1)      = One
      G(:,iX1,iX2,iX3,iGF_h_2)      = One
      G(:,iX1,iX2,iX3,iGF_h_3)      = One
      G(:,iX1,iX2,iX3,iGF_Gm_dd_11) = One
      G(:,iX1,iX2,iX3,iGF_Gm_dd_22) = One
      G(:,iX1,iX2,iX3,iGF_Gm_dd_33) = One
      G(:,iX1,iX2,iX3,iGF_SqrtGm)   = One
      G(:,iX1,iX2,iX3,iGF_Alpha)    = One
      G(:,iX1,iX2,iX3,iGF_Beta_1)   = Zero
      G(:,iX1,iX2,iX3,iGF_Beta_2)   = Zero
      G(:,iX1,iX2,iX3,iGF_Beta_3)   = Zero
      G(:,iX1,iX2,iX3,iGF_Psi)      = One

      CALL ComputeGeometryX_FromScaleFactors( G(:,iX1,iX2,iX3,:) )

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeGeometryX_CARTESIAN


  SUBROUTINE ComputeGeometryX_SPHERICAL &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mass )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      G(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(in)    :: &
      Mass

    INTEGER  :: iX1, iX2, iX3, iNodeX
    REAL(DP) :: XC(3), dX(3), xL_q(3), xG_q(3)
    REAL(DP) :: G_L(nDOFX,nGF)

    XC = 0.0_DP
    dX = 0.0_DP

    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)

      XC(2) = MeshX(2) % Center(iX2)
      dX(2) = MeshX(2) % Width (iX2)

      XC(1) = MeshX(1) % Center(iX1)
      dX(1) = MeshX(1) % Width (iX1)

      ! --- Compute Geometry Fields in Lobatto Points ---

      DO iNodeX = 1, nDOFX

        ! --- Local Coordinates (Lobatto Points) ---

        xL_q = NodesLX_q(1:3,iNodeX)

        ! --- Global Coordinates ---

        xG_q = XC + dX * xL_q

        ! --- Compute Lapse Function and Conformal Factor ---

        G_L(iNodeX,iGF_Alpha) &
          = LapseFunction  ( xG_q(1), Mass )
        G_L(iNodeX,iGF_Psi) &
          = ConformalFactor( xG_q(1), Mass )

        ! --- Set Geometry in Lobatto Points ---

        G_L(iNodeX,iGF_h_1) &
          = G_L(iNodeX,iGF_Psi)**2
        G_L(iNodeX,iGF_h_2) &
          = G_L(iNodeX,iGF_Psi)**2 &
              * xG_q(1)
        G_L(iNodeX,iGF_h_3) &
          = G_L(iNodeX,iGF_Psi)**2 &
              * xG_q(1) * SIN( xG_q(2) )

      END DO

      ! --- Interpolate from Lobatto to Gaussian Points ---

      CALL DGEMV &
             ( 'N', nDOFX, nDOFX, One, LX_L2G, nDOFX, &
               G_L(:,iGF_h_1), 1, Zero, G(:,iX1,iX2,iX3,iGF_h_1), 1 )
      CALL DGEMV &
             ( 'N', nDOFX, nDOFX, One, LX_L2G, nDOFX, &
               G_L(:,iGF_h_2), 1, Zero, G(:,iX1,iX2,iX3,iGF_h_2), 1 )
      CALL DGEMV &
             ( 'N', nDOFX, nDOFX, One, LX_L2G, nDOFX, &
               G_L(:,iGF_h_3), 1, Zero, G(:,iX1,iX2,iX3,iGF_h_3), 1 )

      CALL ComputeGeometryX_FromScaleFactors( G(:,iX1,iX2,iX3,:) )

      CALL DGEMV &
             ( 'N', nDOFX, nDOFX, One, LX_L2G, nDOFX, &
               G_L(:,iGF_Alpha), 1, Zero, G(:,iX1,iX2,iX3,iGF_Alpha), 1 )

      CALL DGEMV &
             ( 'N', nDOFX, nDOFX, One, LX_L2G, nDOFX, &
               G_L(:,iGF_Psi),   1, Zero, G(:,iX1,iX2,iX3,iGF_Psi),   1 )

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeGeometryX_SPHERICAL


  SUBROUTINE ComputeGeometryX_CYLINDRICAL &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mass )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      G(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(in)    :: &
      Mass

    REAL(DP) :: XC(3), dX(3), xL_q(3), xG_q(3)
    REAL(DP) :: G_L(nDOFX,nGF)
    INTEGER  :: iX1, iX2, iX3, iNodeX

    XC = 0.0_DP
    dX = 0.0_DP

    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)

      XC(1) = MeshX(1) % Center(iX1)
      dX(1) = MeshX(1) % Width (iX1)

      DO iNodeX = 1, nDOFX

        ! --- Local Coordinates (Lobatto Points) ---

        xL_q = NodesLX_q(1:3,iNodeX)

        ! --- Global Coordinates ---

        xG_q = XC + dX * xL_q

        ! --- Compute Lapse Function and Conformal Factor ---

        G_L(iNodeX,iGF_Alpha) &
          = One
        G_L(iNodeX,iGF_Psi) &
          = One

        ! --- Set Geometry in Lobatto Points ---

        G_L(iNodeX,iGF_h_1) &
          = G_L(iNodeX,iGF_Psi)**2
        G_L(iNodeX,iGF_h_2) &
          = G_L(iNodeX,iGF_Psi)**2
        G_L(iNodeX,iGF_h_3) &
          = G_L(iNodeX,iGF_Psi)**2 * xG_q(1)

      END DO

      ! --- Interpolate from Lobatto to Gaussian Points ---

      CALL DGEMV &
             ( 'N', nDOFX, nDOFX, One, LX_L2G, nDOFX, &
               G_L(:,iGF_h_1), 1, Zero, G(:,iX1,iX2,iX3,iGF_h_1), 1 )
      CALL DGEMV &
             ( 'N', nDOFX, nDOFX, One, LX_L2G, nDOFX, &
               G_L(:,iGF_h_2), 1, Zero, G(:,iX1,iX2,iX3,iGF_h_2), 1 )
      CALL DGEMV &
             ( 'N', nDOFX, nDOFX, One, LX_L2G, nDOFX, &
               G_L(:,iGF_h_3), 1, Zero, G(:,iX1,iX2,iX3,iGF_h_3), 1 )

      CALL ComputeGeometryX_FromScaleFactors( G(:,iX1,iX2,iX3,:) )

      CALL DGEMV &
             ( 'N', nDOFX, nDOFX, One, LX_L2G, nDOFX, &
               G_L(:,iGF_Alpha), 1, Zero, G(:,iX1,iX2,iX3,iGF_Alpha), 1 )

      CALL DGEMV &
             ( 'N', nDOFX, nDOFX, One, LX_L2G, nDOFX, &
               G_L(:,iGF_Psi),   1, Zero, G(:,iX1,iX2,iX3,iGF_Psi),   1 )

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeGeometryX_CYLINDRICAL


  SUBROUTINE ComputeGeometryX_FromScaleFactors( G )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(inout) :: G(1:,1:)

    G(:,iGF_Gm_dd_11) = MAX( G(:,iGF_h_1)**2, SqrtTiny )
    G(:,iGF_Gm_dd_22) = MAX( G(:,iGF_h_2)**2, SqrtTiny )
    G(:,iGF_Gm_dd_33) = MAX( G(:,iGF_h_3)**2, SqrtTiny )

    G(:,iGF_SqrtGm) = G(:,iGF_h_1) * G(:,iGF_h_2) * G(:,iGF_h_3)

  END SUBROUTINE ComputeGeometryX_FromScaleFactors


  PURE REAL(DP) FUNCTION LapseFunction( R, M )

    REAL(DP), INTENT(in) :: R, M

    ! --- Schwarzschild Metric in Isotropic Coordinates ---

    LapseFunction = ABS( ( MAX( ABS( R ), SqrtTiny ) - Half * M ) &
                       / ( MAX( ABS( R ), SqrtTiny ) + Half * M ) )

    RETURN
  END FUNCTION LapseFunction


  PURE REAL(DP) FUNCTION ConformalFactor( R, M )

    REAL(DP), INTENT(in) :: R, M

    ! --- Schwarzschild Metric in Isotropic Coordinates ---

    ConformalFactor = One + Half * M / MAX( ABS( R ), SqrtTiny )

    RETURN
  END FUNCTION ConformalFactor


END MODULE GeometryComputationModule
