MODULE Euler_BoundaryConditionsModule

  USE KindModule, ONLY: &
    DP
  USE MeshModule, ONLY: &
    MeshX, NodeCoordinate
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX, WeightsX_q
  USE ProgramHeaderModule, ONLY: &
    bcX, swX, nDOFX, nNodesX
  USE UtilitiesModule, ONLY: &
    NodeNumberX
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, TimersStop_Euler, &
    Timer_Euler_BoundaryConditions

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ApplyBoundaryConditions_Euler

  INTEGER, PARAMETER, PUBLIC :: iEuler_ApplyBC_Both  = 0
  INTEGER, PARAMETER, PUBLIC :: iEuler_ApplyBC_Inner = 1
  INTEGER, PARAMETER, PUBLIC :: iEuler_ApplyBC_Outer = 2
  INTEGER, PARAMETER, PUBLIC :: iEuler_ApplyBC_None  = 3

#if defined(THORNADO_OMP_OL)
  !$OMP DECLARE TARGET &
  !$OMP ( iEuler_ApplyBC_Both, iEuler_ApplyBC_Inner, &
  !$OMP   iEuler_ApplyBC_Outer, iEuler_ApplyBC_None )
#elif defined(THORNADO_OACC)
  !$ACC DECLARE CREATE &
  !$ACC ( iEuler_ApplyBC_Both, iEuler_ApplyBC_Inner, &
  !$ACC   iEuler_ApplyBC_Outer, iEuler_ApplyBC_None )
#endif

CONTAINS


  LOGICAL FUNCTION ApplyInnerBC( iApplyBC )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER, INTENT(in) :: iApplyBC

    ApplyInnerBC = .FALSE.
    IF( iApplyBC .EQ. iEuler_ApplyBC_Inner .OR. &
        iApplyBC .EQ. iEuler_ApplyBC_Both ) &
    ApplyInnerBC = .TRUE.

  END FUNCTION ApplyInnerBC


  LOGICAL FUNCTION ApplyOuterBC( iApplyBC )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER, INTENT(in) :: iApplyBC

    ApplyOuterBC = .FALSE.
    IF( iApplyBC .EQ. iEuler_ApplyBC_Outer .OR. &
        iApplyBC .EQ. iEuler_ApplyBC_Both ) &
    ApplyOuterBC = .TRUE.

  END FUNCTION ApplyOuterBC


  SUBROUTINE ApplyBoundaryConditions_Euler &
    ( iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC_Option )

    INTEGER,  INTENT(in)           :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout)        :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    INTEGER,  INTENT(in), OPTIONAL :: &
      iApplyBC_Option(3)

    INTEGER :: iApplyBC(3)

    CALL TimersStart_Euler( Timer_Euler_BoundaryConditions )

    iApplyBC = iEuler_ApplyBC_Both

    IF( PRESENT( iApplyBC_Option ) ) &
      iApplyBC = iApplyBC_Option

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: U, iX_B0, iX_E0, iX_B1, iX_E1, iApplyBC )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( U, iX_B0, iX_E0, iX_B1, iX_E1, iApplyBC )
#endif

    CALL Euler_ApplyBC_X1 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, &
             U(1:nDOFX,iX_B1(1):iX_E1(1), &
                       iX_B1(2):iX_E1(2), &
                       iX_B1(3):iX_E1(3),1:nCF), iApplyBC(1) )

    CALL Euler_ApplyBC_X2 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, &
             U(1:nDOFX,iX_B1(1):iX_E1(1), &
                       iX_B1(2):iX_E1(2), &
                       iX_B1(3):iX_E1(3),1:nCF), iApplyBC(2) )

    CALL Euler_ApplyBC_X3 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, &
             U(1:nDOFX,iX_B1(1):iX_E1(1), &
                       iX_B1(2):iX_E1(2), &
                       iX_B1(3):iX_E1(3),1:nCF), iApplyBC(3) )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: U ) &
    !$OMP MAP( release: iX_B0, iX_E0, iX_B1, iX_E1, iApplyBC )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT( U ) &
    !$ACC DELETE( iX_B0, iX_E0, iX_B1, iX_E1, iApplyBC )
#endif

    CALL TimersStop_Euler( Timer_Euler_BoundaryConditions )

  END SUBROUTINE ApplyBoundaryConditions_Euler


  SUBROUTINE Euler_ApplyBC_X1( iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), &
      iApplyBC
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iCF, iX1, iX2, iX3
    INTEGER :: iNode, iNodeX, iNodeX_0
    INTEGER :: iNodeX1, iNodeX2, iNodeX3, jNodeX, jNodeX1
    REAL(DP) :: D_0, E_0, R_0, R_q

    SELECT CASE ( bcX(1) )

    CASE ( 0 ) ! No Boundary Condition

    CASE ( 1 ) ! Periodic

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( U, iX_B0, iX_E0, swX, iApplyBC )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
      DO iCF = 1, nCF
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = 1, swX(1)
      DO iNode = 1, nDOFX

        ! --- Inner Boundary ---
        IF( ApplyInnerBC( iApplyBC ) ) &
          U(iNode,iX_B0(1)-iX1,iX2,iX3,iCF) &
            = U(iNode,iX_E0(1)-(iX1-1),iX2,iX3,iCF)

        ! --- Outer Boundary ---
        IF( ApplyOuterBC( iApplyBC ) ) &
          U(iNode,iX_E0(1)+iX1,iX2,iX3,iCF) &
            = U(iNode,iX_B0(1)+(iX1-1),iX2,iX3,iCF)

      END DO
      END DO
      END DO
      END DO
      END DO

    CASE ( 2 ) ! Homogeneous

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( U, iX_B0, iX_E0, swX, iApplyBC )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
      DO iCF = 1, nCF
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = 1, swX(1)
      DO iNode = 1, nDOFX

        ! --- Inner Boundary ---
        IF( ApplyInnerBC( iApplyBC ) ) &
          U(iNode,iX_B0(1)-iX1,iX2,iX3,iCF) &
            = U(iNode,iX_B0(1),iX2,iX3,iCF)

        ! --- Outer Boundary ---
        IF( ApplyOuterBC( iApplyBC ) ) &
          U(iNode,iX_E0(1)+iX1,iX2,iX3,iCF) &
            = U(iNode,iX_E0(1),iX2,iX3,iCF)

      END DO
      END DO
      END DO
      END DO
      END DO

    CASE ( 3 ) ! Reflecting
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
      !$OMP PRIVATE( iNodeX, jNodeX, jNodeX1 )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
      !$ACC PRIVATE( iNodeX, jNodeX, jNodeX1 ) &
      !$ACC PRESENT( U, iX_B0, iX_E0, swX, nNodesX )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(4) &
      !$OMP PRIVATE( iNodeX, jNodeX, jNodeX1 )
#endif
      DO iCF = 1, nCF
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = 1, swX(1)

        IF( iCF .EQ. iCF_S1 )THEN

          DO iNodeX3 = 1, nNodesX(3)
          DO iNodeX2 = 1, nNodesX(2)
          DO iNodeX1 = 1, nNodesX(1)

            jNodeX1 = ( nNodesX(1) - iNodeX1 ) + 1

            iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )
            jNodeX = NodeNumberX( jNodeX1, iNodeX2, iNodeX3 )

            ! --- Inner Boundary ---
            IF( ApplyInnerBC( iApplyBC ) ) &
              U(iNodeX,iX_B0(1)-iX1,iX2,iX3,iCF) &
                = - U(jNodeX,iX_B0(1),iX2,iX3,iCF)

            ! --- Outer boundary ---
            IF( ApplyOuterBC( iApplyBC ) ) &
              U(iNodeX,iX_E0(1)+iX1,iX2,iX3,iCF) &
                = - U(jNodeX,iX_E0(1),iX2,iX3,iCF)

          END DO
          END DO
          END DO

        ELSE

          DO iNode = 1, nDOFX

            ! --- Inner Boundary ---
            IF( ApplyInnerBC( iApplyBC ) ) &
              U(iNode,iX_B0(1)-iX1,iX2,iX3,iCF) &
                = U(iNode,iX_B0(1),iX2,iX3,iCF)

            ! --- Outer Boundary ---
            IF( ApplyOuterBC( iApplyBC ) ) &
              U(iNode,iX_E0(1)+iX1,iX2,iX3,iCF) &
                = U(iNode,iX_E0(1),iX2,iX3,iCF)

          END DO

        END IF

      END DO
      END DO
      END DO
      END DO

    CASE ( 30 ) ! Reflecting (Inner), Zero (Outer)

      ! --- Inner Boundary ---
      IF( ApplyInnerBC( iApplyBC ) )THEN

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)

          DO iNodeX3 = 1, nNodesX(3)
          DO iNodeX2 = 1, nNodesX(2)
          DO iNodeX1 = 1, nNodesX(1)

            jNodeX1 = ( nNodesX(1) - iNodeX1 ) + 1

            iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )
            jNodeX = NodeNumberX( jNodeX1, iNodeX2, iNodeX3 )

            DO iCF = 1, nCF

              U(iNodeX,iX_B0(1)-iX1,iX2,iX3,iCF) &
                = U(jNodeX,iX_B0(1),iX2,iX3,iCF)

            END DO

            U(iNodeX,iX_B0(1)-iX1,iX2,iX3,iCF_S1) &
              = - U(jNodeX,iX_B0(1),iX2,iX3,iCF_S1)

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

    CASE ( 31 ) ! Reflecting (Inner), Homogeneous (Outer)

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = 1, swX(1)

        ! --- Inner Boundary ---
        IF( ApplyInnerBC( iApplyBC ) )THEN

          DO iNodeX3 = 1, nNodesX(3)
          DO iNodeX2 = 1, nNodesX(2)
          DO iNodeX1 = 1, nNodesX(1)

            jNodeX1 = ( nNodesX(1) - iNodeX1 ) + 1

            iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )
            jNodeX = NodeNumberX( jNodeX1, iNodeX2, iNodeX3 )

            DO iCF = 1, nCF

              U(iNodeX,iX_B0(1)-iX1,iX2,iX3,iCF) &
                = U(jNodeX,iX_B0(1),iX2,iX3,iCF)

            END DO

            U(iNodeX,iX_B0(1)-iX1,iX2,iX3,iCF_S1) &
              = - U(jNodeX,iX_B0(1),iX2,iX3,iCF_S1)

          END DO
          END DO
          END DO

        END IF

        ! --- Outer Boundary ---
        IF( ApplyOuterBC( iApplyBC ) )THEN
          DO iCF = 1, nCF
            U(:,iX_E0(1)+iX1,iX2,iX3,iCF) &
              = U(:,iX_E0(1),iX2,iX3,iCF)
          END DO
        END IF

      END DO
      END DO
      END DO

    CASE ( 11 ) ! Custom BCs for Accretion Problem

      ! --- Inner Boundary ---
      IF( ApplyInnerBC( iApplyBC ) )THEN

        R_0 = MeshX(1) % Center(1) &
                + MeshX(1) % Width(1) &
                    * MeshX(1) % Nodes(1)

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)

          DO iNodeX3 = 1, nNodesX(3)
          DO iNodeX2 = 1, nNodesX(2)
          DO iNodeX1 = 1, nNodesX(1)

            iNodeX   = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )
            iNodeX_0 = NodeNumberX( 1,       iNodeX2, iNodeX3 )

            D_0 = U(iNodeX_0,1,iX2,iX3,iCF_D)
            E_0 = U(iNodeX_0,1,iX2,iX3,iCF_E)

            R_q = NodeCoordinate( MeshX(1), iX_B0(1)-iX1, iNodeX1 )

            U(iNodeX,iX_B0(1)-iX1,iX2,iX3,iCF_D) &
              = D_0 * ( R_0 / R_q )**3

            U(iNodeX,iX_B0(1)-iX1,iX2,iX3,iCF_E) &
              = E_0 * ( R_0 / R_q )**4

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

    CASE DEFAULT

      WRITE(*,*)
      WRITE(*,'(A5,A45,I2.2)') &
        '', 'Invalid Boundary Condition for Fluid X1: ', bcX(1)
      STOP

    END SELECT

  END SUBROUTINE Euler_ApplyBC_X1


  SUBROUTINE Euler_ApplyBC_X2( iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), &
      iApplyBC
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iCF, iX1, iX2, iX3
    INTEGER :: iNode, iNodeX, iNodeX1, iNodeX2, iNodeX3, jNodeX, jNodeX2

    SELECT CASE ( bcX(2) )

    CASE ( 0 ) ! No Boundary Condition

    CASE ( 1 ) ! Periodic

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( U, iX_B0, iX_E0, swX, iApplyBC )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
      DO iCF = 1, nCF
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = 1, swX(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNode = 1, nDOFX

        ! --- Inner Boundary ---
        IF( ApplyInnerBC( iApplyBC ) ) &
          U(iNode,iX1,iX_B0(2)-iX2,iX3,iCF) &
            = U(iNode,iX1,iX_E0(2)-(iX2-1),iX3,iCF)

        ! --- Outer Boundary ---
        IF( ApplyOuterBC( iApplyBC ) ) &
          U(iNode,iX1,iX_E0(2)+iX2,iX3,iCF) &
            = U(iNode,iX1,iX_B0(2)+(iX2-1),iX3,iCF)

      END DO
      END DO
      END DO
      END DO
      END DO

    CASE ( 2 ) ! Homogeneous

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( U, iX_B0, iX_E0, swX, iApplyBC )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
      DO iCF = 1, nCF
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = 1, swX(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNode = 1, nDOFX

        ! --- Inner Boundary ---
        IF( ApplyInnerBC( iApplyBC ) ) &
          U(iNode,iX1,iX_B0(2)-iX2,iX3,iCF) &
            = U(iNode,iX1,iX_B0(2),iX3,iCF)

         ! --- Outer Boundary ---
         IF( ApplyOuterBC( iApplyBC ) ) &
          U(iNode,iX1,iX_E0(2)+iX2,iX3,iCF) &
            = U(iNode,iX1,iX_E0(2),iX3,iCF)

      END DO
      END DO
      END DO
      END DO
      END DO

    CASE ( 3 ) ! Reflecting

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
      !$OMP PRIVATE( iNodeX, jNodeX, jNodeX2 )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
      !$ACC PRIVATE( iNodeX, jNodeX, jNodeX2 ) &
      !$ACC PRESENT( U, iX_B0, iX_E0, swX, nNodesX )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(4) &
      !$OMP PRIVATE( iNodeX, jNodeX, jNodeX2 )
#endif
      DO iCF = 1, nCF
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = 1, swX(2)
      DO iX1 = iX_B0(1), iX_E0(1)

        IF( iCF .EQ. iCF_S2 )THEN

          DO iNodeX3 = 1, nNodesX(3)
          DO iNodeX2 = 1, nNodesX(2)
          DO iNodeX1 = 1, nNodesX(1)

            jNodeX2 = ( nNodesX(2) - iNodeX2 ) + 1

            iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )
            jNodeX = NodeNumberX( iNodeX1, jNodeX2, iNodeX3 )

            ! --- Inner Boundary ---
            IF( ApplyInnerBC( iApplyBC ) ) &
              U(iNodeX,iX1,iX_B0(2)-iX2,iX3,iCF) &
                = - U(jNodeX,iX1,iX_B0(2),iX3,iCF)

            ! --- Outer boundary ---
            IF( ApplyOuterBC( iApplyBC ) ) &
              U(iNodeX,iX1,iX_E0(2)+iX2,iX3,iCF) &
                = - U(jNodeX,iX1,iX_E0(2),iX3,iCF)

          END DO
          END DO
          END DO

        ELSE

          DO iNode = 1, nDOFX

            ! --- Inner Boundary ---
            IF( ApplyInnerBC( iApplyBC ) ) &
              U(iNode,iX1,iX_B0(2)-iX2,iX3,iCF) &
                = U(iNode,iX1,iX_B0(2),iX3,iCF)

            ! --- Outer Boundary ---
            IF( ApplyOuterBC( iApplyBC ) ) &
              U(iNode,iX1,iX_E0(2)+iX2,iX3,iCF) &
                = U(iNode,iX1,iX_E0(2),iX3,iCF)

          END DO

        END IF

      END DO
      END DO
      END DO
      END DO

    CASE ( 31 ) ! Reflecting (Inner), Homogeneous (Outer)

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          ! --- Inner Boundary ---
          IF( ApplyInnerBC( iApplyBC ) )THEN

            DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
            DO iNodeX1 = 1, nNodesX(1)

              jNodeX2 = ( nNodesX(2) - iNodeX2 ) + 1

              iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )
              jNodeX = NodeNumberX( iNodeX1, jNodeX2, iNodeX3 )

              DO iCF = 1, nCF

                U(iNodeX,iX1,iX_B0(2)-iX2,iX3,iCF) &
                  = U(jNodeX,iX1,iX_B0(2),iX3,iCF)

              END DO

              U(iNodeX,iX1,iX_B0(2)-iX2,iX3,iCF_S2) &
                = - U(jNodeX,iX1,iX_B0(2),iX3,iCF_S2)

            END DO
            END DO
            END DO

          END IF

          ! --- Outer Boundary ---
          IF( ApplyOuterBC( iApplyBC ) )THEN
            DO iCF = 1, nCF
              U(:,iX1,iX_E0(2)+iX2,iX3,iCF) &
                = U(:,iX1,iX_E0(2),iX3,iCF)
            END DO
          END IF

        END DO
        END DO
        END DO

    CASE DEFAULT

      WRITE(*,*)
      WRITE(*,'(A5,A45,I2.2)') &
        '', 'Invalid Boundary Condition for Fluid X2: ', bcX(2)
      STOP

    END SELECT

  END SUBROUTINE Euler_ApplyBC_X2


  SUBROUTINE Euler_ApplyBC_X3( iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), &
      iApplyBC
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iNode, iCF, iX1, iX2, iX3
    INTEGER :: iNodeX, jNodeX, jNodeX3, iNodeX1, iNodeX2, iNodeX3

    SELECT CASE ( bcX(3) )

    CASE ( 0 ) ! No Boundary Condition

    CASE ( 1 ) ! Periodic

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( U, iX_B0, iX_E0, swX, iApplyBC )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
      DO iCF = 1, nCF
      DO iX3 = 1, swX(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNode = 1, nDOFX

          ! --- Inner Boundary ---
          IF( ApplyInnerBC( iApplyBC ) ) &
            U(iNode,iX1,iX2,iX_B0(3)-iX3,iCF) &
              = U(iNode,iX1,iX2,iX_E0(3)-(iX3-1),iCF)

          ! --- Outer Boundary ---
          IF( ApplyOuterBC( iApplyBC ) ) &
            U(iNode,iX1,iX2,iX_E0(3)+iX3,iCF) &
              = U(iNode,iX1,iX2,iX_B0(3)+(iX3-1),iCF)

      END DO
      END DO
      END DO
      END DO
      END DO

    CASE ( 2 ) ! Homogeneous

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( U, iX_B0, iX_E0, swX, iApplyBC )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
      DO iCF = 1, nCF
      DO iX3 = 1, swX(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNode = 1, nDOFX

        ! --- Inner Boundary ---
        IF( ApplyInnerBC( iApplyBC ) ) &
          U(iNode,iX1,iX2,iX_B0(3)-iX3,iCF) &
            = U(iNode,iX1,iX2,iX_B0(3),iCF)

        ! --- Outer Boundary ---
        IF( ApplyOuterBC( iApplyBC ) ) &
          U(iNode,iX1,iX2,iX_E0(3)+iX3,iCF) &
            = U(iNode,iX1,iX2,iX_E0(3),iCF)

      END DO
      END DO
      END DO
      END DO
      END DO

    CASE ( 3 ) ! Reflecting

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
      !$OMP PRIVATE( iNodeX, jNodeX, jNodeX3 )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
      !$ACC PRIVATE( iNodeX, jNodeX, jNodeX3 ) &
      !$ACC PRESENT( U, iX_B0, iX_E0, swX, nNodesX )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(4) &
      !$OMP PRIVATE( iNodeX, jNodeX, jNodeX3 )
#endif
      DO iCF = 1, nCF
      DO iX3 = 1, swX(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)

        IF( iCF .EQ. iCF_S3 )THEN

          DO iNodeX3 = 1, nNodesX(3)
          DO iNodeX2 = 1, nNodesX(2)
          DO iNodeX1 = 1, nNodesX(1)

            jNodeX3 = ( nNodesX(3) - iNodeX3 ) + 1

            iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )
            jNodeX = NodeNumberX( iNodeX1, iNodeX2, jNodeX3 )

            ! --- Inner Boundary ---
            IF( ApplyInnerBC( iApplyBC ) ) &
              U(iNodeX,iX1,iX2,iX_B0(3)-iX3,iCF) &
                = - U(jNodeX,iX1,iX2,iX_B0(3),iCF)

            ! --- Outer boundary ---
            IF( ApplyOuterBC( iApplyBC ) ) &
              U(iNodeX,iX1,iX2,iX_E0(3)+iX3,iCF) &
                = - U(jNodeX,iX1,iX2,iX_E0(3),iCF)

          END DO
          END DO
          END DO

        ELSE

          DO iNode = 1, nDOFX

            ! --- Inner Boundary ---
            IF( ApplyInnerBC( iApplyBC ) ) &
              U(iNode,iX1,iX2,iX_B0(3)-iX3,iCF) &
                = U(iNode,iX1,iX2,iX_B0(3),iCF)

            ! --- Outer Boundary ---
            IF( ApplyOuterBC( iApplyBC ) ) &
              U(iNode,iX1,iX2,iX_E0(3)+iX3,iCF) &
                = U(iNode,iX1,iX2,iX_E0(3),iCF)

          END DO

        END IF

      END DO
      END DO
      END DO
      END DO

    CASE DEFAULT

      WRITE(*,*)
      WRITE(*,'(A5,A45,I2.2)') &
        '', 'Invalid Boundary Condition for Fluid X3: ', bcX(3)
      STOP

    END SELECT

  END SUBROUTINE Euler_ApplyBC_X3


END MODULE Euler_BoundaryConditionsModule
