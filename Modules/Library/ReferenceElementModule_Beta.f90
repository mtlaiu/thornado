MODULE ReferenceElementModule_Beta

  USE KindModule, ONLY: &
    DP
  USE QuadratureModule, ONLY: &
    GetQuadrature
  USE ProgramHeaderModule, ONLY: &
    nNodesX, &
    nNodesE, &
    nNodesZ, &
    nDOF

  IMPLICIT NONE
  PRIVATE

  INTEGER,               PUBLIC :: nDOF_X1
  INTEGER,               PUBLIC :: nDOF_X2
  INTEGER,               PUBLIC :: nDOF_X3
  INTEGER,  ALLOCATABLE, PUBLIC :: NodeNumberTable(:,:)
  INTEGER,  ALLOCATABLE, PUBLIC :: NodeNumberTable_X1(:,:)
  INTEGER,  ALLOCATABLE, PUBLIC :: NodeNumberTable_X2(:,:)
  INTEGER,  ALLOCATABLE, PUBLIC :: NodeNumberTable_X3(:,:)
  INTEGER,  ALLOCATABLE, PUBLIC :: NodeNumberTable4D(:,:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: NodesE(:),  WeightsE(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: NodesX1(:), WeightsX1(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: NodesX2(:), WeightsX2(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: NodesX3(:), WeightsX3(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: Weights_q(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: Weights_X1(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: Weights_X2(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: Weights_X3(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: Nodes_q(:,:)

  PUBLIC :: InitializeReferenceElement
  PUBLIC :: FinalizeReferenceElement
  PUBLIC :: OuterProduct1D3D

CONTAINS


  SUBROUTINE InitializeReferenceElement

    INTEGER :: iNodeE, iNodeX1, iNodeX2, iNodeX3, iNode

    nDOF_X1 = nNodesX(2) * nNodesX(3) * nNodesE
    nDOF_X2 = nNodesX(1) * nNodesX(3) * nNodesE
    nDOF_X3 = nNodesX(1) * nNodesX(2) * nNodesE

    ALLOCATE( NodeNumberTable(4,nDOF) )
    ALLOCATE( NodeNumberTable4D(nNodesZ(1),nNodesZ(2),nNodesZ(3),nNodesZ(4)) )

    iNode = 0
    DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
        DO iNodeX1 = 1, nNodesX(1)
          DO iNodeE = 1, nNodesE

            iNode = iNode + 1

            NodeNumberTable(1:4,iNode) &
              = [ iNodeE, iNodeX1, iNodeX2, iNodeX3 ]

            NodeNumberTable4D(iNodeE,iNodeX1,iNodeX2,iNodeX3) &
              = iNode

          END DO
        END DO
      END DO
    END DO

    ALLOCATE( NodeNumberTable_X1(3,nDOF_X1) )

    iNode = 1
    DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
        DO iNodeE = 1, nNodesE

          NodeNumberTable_X1(1:3,iNode) &
            = [ iNodeE, iNodeX2, iNodeX3 ]

          iNode = iNode + 1

        END DO
      END DO
    END DO

    ALLOCATE( NodeNumberTable_X2(3,nDOF_X2) )

    iNode = 1
    DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX1 = 1, nNodesX(1)
        DO iNodeE = 1, nNodesE

          NodeNumberTable_X2(1:3,iNode) &
            = [ iNodeE, iNodeX1, iNodeX3 ]

          iNode = iNode + 1

        END DO
      END DO
    END DO

    ALLOCATE( NodeNumberTable_X3(3,nDOF_X3) )

    iNode = 1
    DO iNodeX2 = 1, nNodesX(2)
      DO iNodeX1 = 1, nNodesX(1)
        DO iNodeE = 1, nNodesE

          NodeNumberTable_X3(1:3,iNode) &
            = [ iNodeE, iNodeX1, iNodeX2 ]

          iNode = iNode + 1

        END DO
      END DO
    END DO

    ALLOCATE( NodesE (nNodesE),    WeightsE (nNodesE) )
    ALLOCATE( NodesX1(nNodesX(1)), WeightsX1(nNodesX(1)) )
    ALLOCATE( NodesX2(nNodesX(2)), WeightsX2(nNodesX(2)) )
    ALLOCATE( NodesX3(nNodesX(3)), WeightsX3(nNodesX(3)) )

    CALL GetQuadrature( nNodesE,    NodesE,  WeightsE )
    CALL GetQuadrature( nNodesX(1), NodesX1, WeightsX1 )
    CALL GetQuadrature( nNodesX(2), NodesX2, WeightsX2 )
    CALL GetQuadrature( nNodesX(3), NodesX3, WeightsX3 )

    ALLOCATE( Weights_q(nDOF) )
    ALLOCATE( Nodes_q(4,nDOF) )

    DO iNode = 1, nDOF

      iNodeE  = NodeNumberTable(1,iNode)
      iNodeX1 = NodeNumberTable(2,iNode)
      iNodeX2 = NodeNumberTable(3,iNode)
      iNodeX3 = NodeNumberTable(4,iNode)

      Weights_q(iNode) &
        = WeightsE(iNodeE) * WeightsX1(iNodeX1) &
            * WeightsX2(iNodeX2) * WeightsX3(iNodeX3)

      Nodes_q(1:4,iNode) &
        = [ NodesE (iNodeE),  NodesX1(iNodeX1), &
            NodesX2(iNodeX2), NodesX3(iNodeX3) ]

    END DO

    ALLOCATE( Weights_X1(nDOF_X1) )

    DO iNode = 1, nDOF_X1

      iNodeE  = NodeNumberTable_X1(1,iNode)
      iNodeX2 = NodeNumberTable_X1(2,iNode)
      iNodeX3 = NodeNumberTable_X1(3,iNode)

      Weights_X1(iNode) &
        = WeightsE(iNodeE) * WeightsX2(iNodeX2) * WeightsX3(iNodeX3)

    END DO

    ALLOCATE( Weights_X2(nDOF_X2) )

    DO iNode = 1, nDOF_X2

      iNodeE  = NodeNumberTable_X2(1,iNode)
      iNodeX1 = NodeNumberTable_X2(2,iNode)
      iNodeX3 = NodeNumberTable_X2(3,iNode)

      Weights_X2(iNode) &
        = WeightsE(iNodeE) * WeightsX1(iNodeX1) * WeightsX3(iNodeX3)

    END DO

    ALLOCATE( Weights_X3(nDOF_X3) )

    DO iNode = 1, nDOF_X3

      iNodeE  = NodeNumberTable_X3(1,iNode)
      iNodeX1 = NodeNumberTable_X3(2,iNode)
      iNodeX2 = NodeNumberTable_X3(3,iNode)

      Weights_X3(iNode) &
        = WeightsE(iNodeE) * WeightsX1(iNodeX1) * WeightsX2(iNodeX2)

    END DO

  END SUBROUTINE InitializeReferenceElement


  SUBROUTINE FinalizeReferenceElement

    DEALLOCATE( NodeNumberTable )
    DEALLOCATE( NodeNumberTable_X1 )
    DEALLOCATE( NodeNumberTable_X2 )
    DEALLOCATE( NodeNumberTable_X3 )
    DEALLOCATE( NodeNumberTable4D )
    DEALLOCATE( NodesE,  WeightsE )
    DEALLOCATE( NodesX1, WeightsX1 )
    DEALLOCATE( NodesX2, WeightsX2 )
    DEALLOCATE( NodesX3, WeightsX3 )
    DEALLOCATE( Nodes_q, Weights_q )
    DEALLOCATE( Weights_X1 )
    DEALLOCATE( Weights_X2 )
    DEALLOCATE( Weights_X3 )

  END SUBROUTINE FinalizeReferenceElement


  PURE FUNCTION OuterProduct1D3D( X1D, n1D, X3D, n3D )

    INTEGER,  INTENT(in) :: n1D, n3D
    REAL(DP), INTENT(in) :: X1D(n1D)
    REAL(DP), INTENT(in) :: X3D(n3D)
    REAL(DP)             :: OuterProduct1D3D(n1D*n3D)

    INTEGER :: i3D

    DO i3D = 1, n3D

      OuterProduct1D3D((i3D-1)*n1D+1:(i3D-1)*n1D+n1D) = X1D(:) * X3D(i3D)

    END DO

    RETURN
  END FUNCTION OuterProduct1D3D


END MODULE ReferenceElementModule_Beta