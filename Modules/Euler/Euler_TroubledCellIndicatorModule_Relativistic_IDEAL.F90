MODULE Euler_TroubledCellIndicatorModule_Relativistic_IDEAL

  USE KindModule, ONLY: &
    DP, Zero, One
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nDimsX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX, &
    NodesX1, WeightsX1, &
    NodesX2, WeightsX2, &
    NodesX3, WeightsX3, &
    WeightsX_q
  USE LinearAlgebraModule, ONLY: &
    MatrixVectorMultiply
  USE PolynomialBasisModule_Lagrange, ONLY: &
    L_X1, L_X2, L_X3
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_SqrtGm
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_E, Shock
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, TimersStop_Euler, &
    Timer_Euler_TroubledCellIndicator, &
    Timer_Euler_CopyIn, Timer_Euler_CopyOut, &
    Timer_Euler_Permute

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeTroubledCellIndicator_Euler_Relativistic_IDEAL
  PUBLIC :: FinalizeTroubledCellIndicator_Euler_Relativistic_IDEAL
  PUBLIC :: DetectTroubledCells_Euler_Relativistic_IDEAL

  REAL(DP), ALLOCATABLE :: WeightsX_X1_P(:), WeightsX_X1_N(:)
  REAL(DP), ALLOCATABLE :: WeightsX_X2_P(:), WeightsX_X2_N(:)
  REAL(DP), ALLOCATABLE :: WeightsX_X3_P(:), WeightsX_X3_N(:)


CONTAINS


  SUBROUTINE InitializeTroubledCellIndicator_Euler_Relativistic_IDEAL

    INTEGER  :: iNode, iNodeX1, iNodeX2, iNodeX3
    INTEGER  :: jNode, jNodeX1, jNodeX2, jNodeX3
    REAL(DP) :: WeightX

    ALLOCATE( WeightsX_X1_P(nDOFX), WeightsX_X1_N(nDOFX) )
    ALLOCATE( WeightsX_X2_P(nDOFX), WeightsX_X2_N(nDOFX) )
    ALLOCATE( WeightsX_X3_P(nDOFX), WeightsX_X3_N(nDOFX) )

    ! --- Compute Weights for Extrapolating Neighbors into Target Cell ---

    DO jNode = 1, nDOFX

      jNodeX1 = NodeNumberTableX(1,jNode)
      jNodeX2 = NodeNumberTableX(2,jNode)
      jNodeX3 = NodeNumberTableX(3,jNode)

      WeightsX_X1_P(jNode) = Zero
      WeightsX_X1_N(jNode) = Zero
      WeightsX_X2_P(jNode) = Zero
      WeightsX_X2_N(jNode) = Zero
      WeightsX_X3_P(jNode) = Zero
      WeightsX_X3_N(jNode) = Zero

      DO iNode = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNode)
        iNodeX2 = NodeNumberTableX(2,iNode)
        iNodeX3 = NodeNumberTableX(3,iNode)

        WeightX = WeightsX1  (iNodeX1) &
                  * WeightsX2(iNodeX2) &
                  * WeightsX3(iNodeX3)

        WeightsX_X1_P(jNode) &
          = WeightsX_X1_P(jNode) &
              + WeightX &
                * ( L_X1  (jNodeX1) % P( NodesX1(iNodeX1) + One ) &
                    * L_X2(jNodeX2) % P( NodesX2(iNodeX2) ) &
                    * L_X3(jNodeX3) % P( NodesX3(iNodeX3) ) )

        WeightsX_X1_N(jNode) &
          = WeightsX_X1_N(jNode) &
              + WeightX &
                * ( L_X1  (jNodeX1) % P( NodesX1(iNodeX1) - One ) &
                    * L_X2(jNodeX2) % P( NodesX2(iNodeX2) ) &
                    * L_X3(jNodeX3) % P( NodesX3(iNodeX3) ) )

        WeightsX_X2_P(jNode) &
          = WeightsX_X2_P(jNode) &
              + WeightX &
                * ( L_X1  (jNodeX1) % P( NodesX1(iNodeX1) ) &
                    * L_X2(jNodeX2) % P( NodesX2(iNodeX2) + One ) &
                    * L_X3(jNodeX3) % P( NodesX3(iNodeX3) ) )

        WeightsX_X2_N(jNode) &
          = WeightsX_X2_N(jNode) &
              + WeightX &
                * ( L_X1  (jNodeX1) % P( NodesX1(iNodeX1) ) &
                    * L_X2(jNodeX2) % P( NodesX2(iNodeX2) - One ) &
                    * L_X3(jNodeX3) % P( NodesX3(iNodeX3) ) )

        WeightsX_X3_P(jNode) &
          = WeightsX_X3_P(jNode) &
              + WeightX &
                * ( L_X1  (jNodeX1) % P( NodesX1(iNodeX1) ) &
                    * L_X2(jNodeX2) % P( NodesX2(iNodeX2) ) &
                    * L_X3(jNodeX3) % P( NodesX3(iNodeX3) + One ) )

        WeightsX_X3_N(jNode) &
          = WeightsX_X3_N(jNode) &
              + WeightX &
                * ( L_X1  (jNodeX1) % P( NodesX1(iNodeX1) ) &
                    * L_X2(jNodeX2) % P( NodesX2(iNodeX2) ) &
                    * L_X3(jNodeX3) % P( NodesX3(iNodeX3) - One ) )

      END DO

    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: WeightsX_X1_P, WeightsX_X1_N, &
    !$OMP          WeightsX_X2_P, WeightsX_X2_N, &
    !$OMP          WeightsX_X3_P, WeightsX_X3_N )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(  WeightsX_X1_P, WeightsX_X1_N, &
    !$ACC          WeightsX_X2_P, WeightsX_X2_N, &
    !$ACC          WeightsX_X3_P, WeightsX_X3_N )
#endif

  END SUBROUTINE InitializeTroubledCellIndicator_Euler_Relativistic_IDEAL


  SUBROUTINE FinalizeTroubledCellIndicator_Euler_Relativistic_IDEAL

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: WeightsX_X1_P, WeightsX_X1_N, &
    !$OMP               WeightsX_X2_P, WeightsX_X2_N,  &
    !$OMP               WeightsX_X3_P, WeightsX_X3_N )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE(     WeightsX_X1_P, WeightsX_X1_N, &
    !$ACC             WeightsX_X2_P, WeightsX_X2_N,  &
    !$ACC             WeightsX_X3_P, WeightsX_X3_N )
#endif


    DEALLOCATE( WeightsX_X1_P, WeightsX_X1_N )
    DEALLOCATE( WeightsX_X2_P, WeightsX_X2_N )
    DEALLOCATE( WeightsX_X3_P, WeightsX_X3_N )

  END SUBROUTINE FinalizeTroubledCellIndicator_Euler_Relativistic_IDEAL


  SUBROUTINE DetectTroubledCells_Euler_Relativistic_IDEAL &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iNodeX, iX1, iX2, iX3, iCF
    INTEGER :: nX(3), nCF_K

    REAL(DP) :: U_X( 1:nDOFX,iX_B0(1)  :iX_E0(1), &
                             iX_B0(2)  :iX_E0(2), &
                             iX_B0(3)  :iX_E0(3),1:nCF)
    REAL(DP) :: U_X1(1:nDOFX,iX_B0(2)  :iX_E0(2), &
                             iX_B0(3)  :iX_E0(3),1:nCF, &
                             iX_B0(1)-1:iX_E0(1)+1)
    REAL(DP) :: U_X2(1:nDOFX,iX_B0(1)  :iX_E0(1), &
                             iX_B0(3)  :iX_E0(3),1:nCF, &
                             iX_B0(2)-1:iX_E0(2)+1)
    REAL(DP) :: U_X3(1:nDOFX,iX_B0(1)  :iX_E0(1), &
                             iX_B0(2)  :iX_E0(2),1:nCF, &
                             iX_B0(3)-1:iX_E0(3)+1)
    REAL(DP) :: U_K(         iX_B0(1)  :iX_E0(1), &
                             iX_B0(2)  :iX_E0(2), &
                             iX_B0(3)  :iX_E0(3),1:nCF)
    REAL(DP) :: U_K0_X1(     iX_B0(2)  :iX_E0(2), &
                             iX_B0(3)  :iX_E0(3),1:nCF, &
                             iX_B0(1)  :iX_E0(1),1:2)
    REAL(DP) :: U_K0_X2(     iX_B0(1)  :iX_E0(1), &
                             iX_B0(3)  :iX_E0(3),1:nCF, &
                             iX_B0(2)  :iX_E0(2),1:2)
    REAL(DP) :: U_K0_X3(     iX_B0(1)  :iX_E0(1), &
                             iX_B0(2)  :iX_E0(2),1:nCF, &
                             iX_B0(3)  :iX_E0(3),1:2)
    REAL(DP) :: Max_UK(      iX_B0(1)  :iX_E0(1), &
                             iX_B0(2)  :iX_E0(2), &
                             iX_B0(3)  :iX_E0(3),1:nCF)

    CALL TimersStart_Euler( Timer_Euler_TroubledCellIndicator )

    nX    = iX_E0 - iX_B0 + 1
    nCF_K = nCF * PRODUCT( nX )

    CALL TimersStart_Euler( Timer_Euler_CopyIn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: iX_B0, iX_E0, U ), &
    !$OMP MAP( alloc: U_X1   , U_X2   , U_X3, U_X, &
    !$OMP             U_K0_X1, U_K0_X2, U_K0_X3, U_K, Max_UK )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(  iX_B0, iX_E0, U ) &
    !$ACC CREATE(     U_X1   , U_X2   , U_X3   , U_X, &
    !$ACC             U_K0_X1, U_K0_X2, U_K0_X3, U_K, Max_UK )
#endif

    CALL TimersStop_Euler( Timer_Euler_CopyIn )

    CALL TimersStart_Euler( Timer_Euler_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5), &
    !$ACC PRESENT( U, U_X, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iCF = 1, nCF
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX

      U_X(iNodeX,iX1,iX2,iX3,iCF) = U(iNodeX,iX1,iX2,iX3,iCF)

    END DO
    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( U, U_X1, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX1 = iX_B0(1)-1, iX_E0(1)+1
    DO iCF = 1, nCF
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iNodeX = 1, nDOFX

      U_X1(iNodeX,iX2,iX3,iCF,iX1) = U(iNodeX,iX1,iX2,iX3,iCF)

    END DO
    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( U, U_X2, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX2 = iX_B0(2)-1, iX_E0(2)+1
    DO iCF = 1, nCF
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX

      U_X2(iNodeX,iX1,iX3,iCF,iX2) = U(iNodeX,iX1,iX2,iX3,iCF)

    END DO
    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( U, U_X3, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX3 = iX_B0(3)-1, iX_E0(3)+1
    DO iCF = 1, nCF
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX

      U_X3(iNodeX,iX1,iX2,iCF,iX3) = U(iNodeX,iX1,iX2,iX3,iCF)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_Permute )

    ! --- Compute cell-averages  ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nCF_K, One, U_X, nDOFX, WeightsX_q, 1, Zero, U_K, 1 )

    ! --- Compute cell-averages of neighbors (X1) ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nCF_K, One, &
             U_X1   (1,iX_B0(2),iX_B0(3),1,iX_B0(1)-1), &
             nDOFX, WeightsX_X1_P, 1, Zero, &
             U_K0_X1(  iX_B0(2),iX_B0(3),1,iX_B0(1)  ,1), 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nCF_K, One, &
             U_X1   (1,iX_B0(2),iX_B0(3),1,iX_B0(1)+1), &
             nDOFX, WeightsX_X1_N, 1, Zero, &
             U_K0_X1(  iX_B0(2),iX_B0(3),1,iX_B0(1)  ,2), 1 )

    ! --- Compute cell-averages of neighbors (X2) ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nCF_K, One, &
             U_X2   (1,iX_B0(1),iX_B0(3),1,iX_B0(2)-1), &
             nDOFX, WeightsX_X2_P, 1, Zero, &
             U_K0_X2(  iX_B0(1),iX_B0(3),1,iX_B0(2)  ,1), 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nCF_K, One, &
             U_X2   (1,iX_B0(1),iX_B0(3),1,iX_B0(2)+1), &
             nDOFX, WeightsX_X2_N, 1, Zero, &
             U_K0_X2(  iX_B0(1),iX_B0(3),1,iX_B0(2)  ,2), 1 )

    ! --- Compute cell-averages of neighbors (X3) ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nCF_K, One, &
             U_X3   (1,iX_B0(1),iX_B0(2),1,iX_B0(3)-1), &
             nDOFX, WeightsX_X3_P, 1, Zero, &
             U_K0_X3(  iX_B0(1),iX_B0(2),1,iX_B0(3)  ,1), 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nCF_K, One, &
             U_X3   (1,iX_B0(1),iX_B0(2),1,iX_B0(3)+1), &
             nDOFX, WeightsX_X3_N, 1, Zero, &
             U_K0_X3(  iX_B0(1),iX_B0(2),1,iX_B0(3)  ,2), 1 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4)
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iCF = 1, nCF
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      Max_UK(iX1,iX2,iX3,iCF) = MAX( ABS( U_K(    iX1,iX2,iX3,iCF) ), &
                                     ABS( U_K0_X1(iX2,iX3,iCF,iX1,1) ), &
                                     ABS( U_K0_X1(iX2,iX3,iCF,iX1,2) ), &
                                     ABS( U_K0_X2(iX1,iX3,iCF,iX2,1) ), &
                                     ABS( U_K0_X2(iX1,iX3,iCF,iX2,2) ), &
                                     ABS( U_K0_X3(iX1,iX2,iCF,iX3,1) ), &
                                     ABS( U_K0_X3(iX1,iX2,iCF,iX3,2) ) )


    END DO
    END DO
    END DO
    END DO

    ! --- Troubled-Cell Indicator from Fu & Shu (2017) ---
    ! --- JCP, 347, 305 - 327 ----------------------------

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( Shock, U_K, U_K0_X1, U_K0_X2, U_K0_X3, Max_UK, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(3)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      ! --- Use Conserved Density to Detect Troubled Cell ---

      Shock(iX1,iX2,iX3) &
        = ( ABS( U_K(iX1,iX2,iX3,iCF_D) - U_K0_X1(iX2,iX3,iCF_D,iX1,1) ) &
            + ABS( U_K(iX1,iX2,iX3,iCF_D) - U_K0_X1(iX2,iX3,iCF_D,iX1,2) ) &
            + ABS( U_K(iX1,iX2,iX3,iCF_D) - U_K0_X2(iX1,iX3,iCF_D,iX2,1) ) &
            + ABS( U_K(iX1,iX2,iX3,iCF_D) - U_K0_X2(iX1,iX3,iCF_D,iX2,2) ) &
            + ABS( U_K(iX1,iX2,iX3,iCF_D) - U_K0_X3(iX1,iX2,iCF_D,iX3,1) ) &
            + ABS( U_K(iX1,iX2,iX3,iCF_D) - U_K0_X3(iX1,iX2,iCF_D,iX3,2) ) ) &
          / Max_UK(iX1,iX2,iX3,iCF_D)

      ! --- Use Conserved Energy  to Detect Troubled Cell ---

      Shock(iX1,iX2,iX3) &
        = MAX( Shock(iX1,iX2,iX3), &
            ( ABS( U_K(iX1,iX2,iX3,iCF_E) - U_K0_X1(iX2,iX3,iCF_E,iX1,1) ) &
              + ABS( U_K(iX1,iX2,iX3,iCF_E) - U_K0_X1(iX2,iX3,iCF_E,iX1,2) ) &
              + ABS( U_K(iX1,iX2,iX3,iCF_E) - U_K0_X2(iX1,iX3,iCF_E,iX2,1) ) &
              + ABS( U_K(iX1,iX2,iX3,iCF_E) - U_K0_X2(iX1,iX3,iCF_E,iX2,2) ) &
              + ABS( U_K(iX1,iX2,iX3,iCF_E) - U_K0_X3(iX1,iX2,iCF_E,iX3,1) ) &
              + ABS( U_K(iX1,iX2,iX3,iCF_E) - U_K0_X3(iX1,iX2,iCF_E,iX3,2) ) ) &
            / Max_UK(iX1,iX2,iX3,iCF_E) )

    END DO
    END DO
    END DO

    CALL TimersStart_Euler( Timer_Euler_CopyOut )

!!$#if defined(THORNADO_OMP_OL)
!!$    !$OMP TARGET UPDATE FROM( Shock )
!!$#elif defined(THORNADO_OACC)
!!$    !$ACC UPDATE HOST( Shock )
!!$#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: U_X1   , U_X2   , U_X3, U_X, &
    !$OMP               U_K0_X1, U_K0_X2, U_K0_X3, U_K, Max_UK, &
    !$OMP               iX_B0, iX_E0, U )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE(       U_X1   , U_X2   , U_X3   , U_X, &
    !$ACC               U_K0_X1, U_K0_X2, U_K0_X3, U_K, Max_UK, &
    !$ACC               iX_B0, iX_E0, U )
#endif

    CALL TimersStop_Euler( Timer_Euler_CopyOut )

    CALL TimersStop_Euler( Timer_Euler_TroubledCellIndicator )

  END SUBROUTINE DetectTroubledCells_Euler_Relativistic_IDEAL


END MODULE Euler_TroubledCellIndicatorModule_Relativistic_IDEAL
