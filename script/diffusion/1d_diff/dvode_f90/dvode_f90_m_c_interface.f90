      SUBROUTINE DVODE_F90_SPARSE(NEQ,Y,T,TOUT,             &
          ITASK,ISTATE,RTOL,ATOL,F,JAC,REINIT)
! ..
! An interface to call DVODE_F90 with sparse matrix (MF = 126)
! ..
          USE DVODE_F90_M

          IMPLICIT NONE

          INTEGER, PARAMETER :: WP = KIND(1.0D0)

          REAL (WP), INTENT (INOUT) :: T, TOUT
          INTEGER, INTENT (INOUT) :: ISTATE
          INTEGER, INTENT (IN) :: ITASK, NEQ
          REAL (WP), INTENT (INOUT) :: Y(*)
          REAL (WP), INTENT (IN) :: RTOL, ATOL
          LOGICAL, INTENT (IN) :: REINIT

          EXTERNAL JAC
!       SUBROUTINE JAC (N, T, Y, IA, JA, NZ, PD)
!    -  If NZ = 0 on input:
!       Replace NZ by the number of nonzero elements in the Jacobian.
!       The diagonal of the Jacobian must be included.
!       Do NOT define the arrays IA, JA, PD at this time.
!       Once JAC has been called with NZ = 0 and you have defined the
!       value of NZ, future calls to JAC will use this value of NZ.
!    -  When a call is made with NZ unequal to 0, you must define the
!       sparsity structure arrays IA and JA, and the sparse Jacobian
!       PD.
!         - IA defines the number of nonzeros including the diagonal
!           in each column of the Jacobian. Define IA(1) = 1 and for
!           J = 1,..., N,
!           IA(J+1) = IA(J) + number of nonzeros in column J.
!           Diagonal elements must be include even if they are zero.
!           You should check to ensure that IA(N+1)-1 = NZ.
!         - JA defines the rows in which the nonzeros occur. For
!           I = 1,...,NZ, JA(I) is the row in which the Ith element
!           of the Jacobian occurs. JA must also include the diagonal
!           elements of the Jacobian.
!         - PD defines the numerical value of the Jacobian elements.
!           For I = 1,...,NZ, PD(I) is the numerical value of the
!           Ith element in the Jacobian. PD must also include the
!           diagonal elements of the Jacobian.

          INTERFACE
          SUBROUTINE F(NEQ,T,Y,YDOT)
              INTEGER, PARAMETER :: WP = KIND(1.0D0)
              INTEGER :: NEQ
              REAL(WP) :: T
              REAL(WP), DIMENSION(NEQ) :: Y, YDOT
              INTENT(IN)  :: NEQ, T, Y
              INTENT(OUT) :: YDOT
          END SUBROUTINE F
          END INTERFACE

          TYPE (VODE_OPTS), SAVE :: OPTS
          LOGICAL, SAVE :: FIRSTTIME = .TRUE.
          LOGICAL :: USER_JACOBIAN = .FALSE.

          IF (FIRSTTIME .OR. REINIT) THEN
              OPTS = SET_INTERMEDIATE_OPTS(SPARSE_J=.TRUE.,            &
                    ABSERR=ATOL, RELERR=RTOL,                          &
                    USER_SUPPLIED_JACOBIAN=USER_JACOBIAN               &
                    )
             FIRSTTIME = .FALSE.
          ENDIF

          IF (REINIT) ISTATE = 1
          CALL VODE_F90 (F,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTS,JAC)

      END SUBROUTINE DVODE_F90_SPARSE
