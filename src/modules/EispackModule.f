!> Numerical routines from the EISPACK library
!! See https://www.netlib.org/eispack/
!! \note REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
!!                Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
!!                SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
!!                1976.
      MODULE EISPACKMODULE

      PRIVATE
      PUBLIC :: FMS_CH, FMS_RS

      CONTAINS
!<
!!    Finds sqrt(A**2+B**2) without overflow or destructive underflow
!!    \author Modified to Real*8 by tjm
!!    @ingroup numerics
!>
      REAL*8 FUNCTION FMS_PYTHAG(A,B)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: A,B
      REAL*8 :: P,Q,R,S,T

      P = MAX(ABS(A),ABS(B))
      Q = MIN(ABS(A),ABS(B))
      IF (Q .EQ. 0.0D0) GO TO 20
   10 CONTINUE
         R = (Q/P)**2
         T = 4.0D0 + R
         IF (T .EQ. 4.0D0) GO TO 20
         S = R/T
         P = P + 2.0D0*P*S
         Q = Q*S
      GO TO 10
   20 FMS_PYTHAG = P
      RETURN
      END FUNCTION FMS_PYTHAG
!>
!!  Computes eigenvalues and eigenvectors of a symmetric tridiagonal matrix.
!!
!!  DATE WRITTEN   760101   (YYMMDD)
!!  REVISION DATE  861211   (YYMMDD)
!!  CATEGORY NO.  D4A5,D4C2A
!!  KEYWORDS  LIBRARY=SLATEC(EISPACK),TYPE=SINGLE PRECISION(TQL2-S),
!!            EIGENVALUES,EIGENVECTORS
!!  AUTHOR  SMITH, B. T., ET AL.
!!  DESCRIPTION
!!    This subroutine is a translation of the ALGOL procedure TQL2,
!!    NUM. MATH. 11, 293-306(1968) by Bowdler, Martin, Reinsch, and
!!    Wilkinson.
!!    HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
!!
!!    This subroutine finds the eigenvalues and eigenvectors
!!    of a SYMMETRIC TRIDIAGONAL matrix by the QL method.
!!    The eigenvectors of a FULL SYMMETRIC matrix can also
!!    be found if  TRED2  has been used to reduce this
!!    full matrix to tridiagonal form.
!!    On Input
!!    \param   NM must be set to the row dimension of two-dimensional
!!         array parameters as declared in the calling program
!!         dimension statement.
!!    \param   N is the order of the matrix.
!!    \param   D contains the diagonal elements of the input matrix.
!!    \param   E contains the subdiagonal elements of the input matrix
!!         in its last N-1 positions.  E(1) is arbitrary.
!!    \param   Z contains the transformation matrix produced in the
!!         reduction by  TRED2, if performed.  If the eigenvectors
!!         of the tridiagonal matrix are desired, Z must contain
!!         the identity matrix.
!!     On Output
!!     \param  D contains the eigenvalues in ascending order.  If an
!!         error exit is made, the eigenvalues are correct but
!!         unordered for indices 1,2,...,IERR-1.
!!     \param  E has been destroyed.
!!     \param  Z contains orthonormal eigenvectors of the symmetric
!!         tridiagonal (or full) matrix.  If an error exit is made,
!!         Z contains the eigenvectors associated with the stored
!!         eigenvalues.
!!   \param    IERR is set to
!!         Zero       for normal return,
!!         J          if the J-th eigenvalue has not been
!!                    determined after 30 iterations.
!!    Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
!!    Questions and comments should be directed to B. S. Garbow,
!!    APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!!    
!! \note REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
!!                Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
!!                SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
!!                1976.
!!    \author Modified to Real*8 by tjm
!<
      SUBROUTINE FMS_TQL2(NM,N,D,E,Z,IERR)
      implicit real*8 (a-h,o-z)
      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
      REAL*8 D(N),E(N),Z(NM,N)
      REAL*8 B,C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2

!     FIRST EXECUTABLE STATEMENT  TQL2
      IERR = 0
      IF (N .EQ. 1) GO TO 1001

      DO 100 I = 2, N
  100 E(I-1) = E(I)

      F = 0.0D0
      B = 0.0D0
      E(N) = 0.0D0

      DO 240 L = 1, N
         J = 0
         H = ABS(D(L)) + ABS(E(L))
         IF (B .LT. H) B = H
!     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
         DO 110 M = L, N
            IF (B + ABS(E(M)) .EQ. B) GO TO 120
!     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE

  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
!     .......... FORM SHIFT ..........
         L1 = L + 1
         L2 = L1 + 1
         G = D(L)
         P = (D(L1) - G) / (2.0D0 * E(L))
         R = FMS_PYTHAG(P,1.0D0)
         D(L) = E(L) / (P + SIGN(R,P))
         D(L1) = E(L) * (P + SIGN(R,P))
         DL1 = D(L1)
         H = G - D(L)
         IF (L2 .GT. N) GO TO 145

         DO 140 I = L2, N
  140    D(I) = D(I) - H

  145    F = F + H
!     .......... QL TRANSFORMATION ..........
         P = D(M)
         C = 1.0D0
         C2 = C
         EL1 = E(L1)
         S = 0.0D0
         MML = M - L
!     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            C3 = C2
            C2 = C
            S2 = S
            I = M - II
            G = C * E(I)
            H = C * P
            IF (ABS(P) .LT. ABS(E(I))) GO TO 150
            C = E(I) / P
            R = SQRT(C*C+1.0D0)
            E(I+1) = S * P * R
            S = C / R
            C = 1.0D0 / R
            GO TO 160
  150       C = P / E(I)
            R = SQRT(C*C+1.0D0)
            E(I+1) = S * E(I) * R
            S = 1.0D0 / R
            C = C * S
  160       P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
!     .......... FORM VECTOR ..........
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE

  200    CONTINUE

         P = -S * S2 * C3 * EL1 * E(L) / DL1
         E(L) = S * P
         D(L) = C * P
         IF (B + ABS(E(L)) .GT. B) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
!     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)

         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE

         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P

         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE

  300 CONTINUE

      GO TO 1001
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END SUBROUTINE FMS_TQL2
!>
!!    Computes eigenvalues of symmetric tridiagonal matrix
!!    using a rational variant of the QL method.
!!
!!  DATE WRITTEN   760101   (YYMMDD)
!!  REVISION DATE  861211   (YYMMDD)
!!  CATEGORY NO.  D4A5,D4C2A
!!  KEYWORDS  LIBRARY=SLATEC(EISPACK),TYPE=SINGLE PRECISION(TQLRAT-S),
!!           EIGENVALUES,EIGENVECTORS
!!  AUTHOR  SMITH, B. T., ET AL.
!!  DESCRIPTION
!!
!!   This subroutine is a translation of the ALGOL procedure TQLRAT,
!!   ALGORITHM 464, COMM. ACM 16, 689(1973) by Reinsch.
!!
!!   This subroutine finds the eigenvalues of a SYMMETRIC
!!   TRIDIAGONAL matrix by the rational QL method.
!!
!!   On Input
!!
!!   \param   N is the order of the matrix.
!!
!!   \param   D contains the diagonal elements of the input matrix.
!!
!!   \param   E2 contains the squares of the subdiagonal elements of the
!!        input matrix in its last N-1 positions.  E2(1) is arbitrary.
!!
!!    On Output
!!
!!   \param  D contains the eigenvalues in ascending order.  If an
!!        error exit is made, the eigenvalues are correct and
!!        ordered for indices 1,2,...IERR-1, but may not be
!!        the smallest eigenvalues.
!!
!!   \param E2 has been destroyed.
!!
!!   \param   IERR is set to
!!        Zero       for normal return,
!!        J          if the J-th eigenvalue has not been
!!                   determined after 30 iterations.
!!
!!   Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
!!   Questions and comments should be directed to B. S. Garbow,
!!   APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!!  
!! \note REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
!!               Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
!!               SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
!!               1976.
!! \author     Modified to real*8 by tjm
!<
      SUBROUTINE FMS_TQLRAT(N,D,E2,IERR)
      implicit real*8 (a-h,o-z)
      INTEGER I,J,L,M,N,II,L1,MML,IERR
      REAL*8 D(N),E2(N)
      REAL*8 B,C,F,G,H,P,R,S,MACHEP

      SAVE MACHEP
      DATA MACHEP/1.0D0/
!  FIRST EXECUTABLE STATEMENT  TQLRAT
      IF (MACHEP .NE. 1.0D0) GO TO 10
   05 MACHEP = 0.5D0*MACHEP
      IF (1.0D0 + MACHEP .GT. 1.0D0) GO TO 05
      MACHEP = 2.0D0*MACHEP

   10 IERR = 0
      IF (N .EQ. 1) GO TO 1001

      DO 100 I = 2, N
  100 E2(I-1) = E2(I)

      F = 0.0D0
      B = 0.0D0
      E2(N) = 0.0D0

      DO 290 L = 1, N
         J = 0
         H = MACHEP * (ABS(D(L)) + SQRT(E2(L)))
         IF (B .GT. H) GO TO 105
         B = H
         C = B * B
!   .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
            IF (E2(M) .LE. C) GO TO 120
!   .......... E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!              THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE

  120    IF (M .EQ. L) GO TO 210
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
!   .......... FORM SHIFT ..........
         L1 = L + 1
         S = SQRT(E2(L))
         G = D(L)
         P = (D(L1) - G) / (2.0D0 * S)
         R = FMS_PYTHAG(P,1.0D0)
         D(L) = S / (P + SIGN(R,P))
         H = G - D(L)

         DO 140 I = L1, N
  140    D(I) = D(I) - H

         F = F + H
!   .......... RATIONAL QL TRANSFORMATION ..........
         G = D(M)
         IF (G .EQ. 0.0D0) G = B
         H = G
         S = 0.0D0
         MML = M - L
!   .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            I = M - II
            P = G * H
            R = P + E2(I)
            E2(I+1) = S * R
            S = E2(I) / R
            D(I+1) = H + S * (H + D(I))
            G = D(I) - E2(I) / G
            IF (G .EQ. 0.0D0) G = B
            H = G * P / R
  200    CONTINUE

         E2(L) = S * G
         D(L) = H
!   .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ..........
         IF (H .EQ. 0.0D0) GO TO 210
         IF (ABS(E2(L)) .LE. ABS(C/H)) GO TO 210
         E2(L) = H * E2(L)
         IF (E2(L) .NE. 0.0D0) GO TO 130
  210    P = D(L) + F
!   .......... ORDER EIGENVALUES ..........
         IF (L .EQ. 1) GO TO 250
!   .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE

  250    I = 1
  270    D(I) = P
  290 CONTINUE

      GO TO 1001
!   .......... SET ERROR -- NO CONVERGENCE TO AN
!              EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END SUBROUTINE FMS_TQLRAT
!>
!!  Reduces a real symmetric matrix to symmetric tridiagonal
!!  matrix using orthogonal similarity transformations.
!!  
!!  DATE WRITTEN   760101   (YYMMDD)
!!  REVISION DATE  830518   (YYMMDD)
!!  CATEGORY NO.  D4C1B1
!!  KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
!!  AUTHOR  SMITH, B. T., ET AL.
!!  DESCRIPTION
!!
!!    This subroutine is a translation of the ALGOL procedure TRED1,
!!    NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
!!    HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!!
!!    This subroutine reduces a REAL SYMMETRIC matrix
!!    to a symmetric tridiagonal matrix using
!!    orthogonal similarity transformations.
!!
!!    On Input
!!
!!       NM must be set to the row dimension of two-dimensional
!!         array parameters as declared in the calling program
!!         dimension statement.
!!
!! \param      N is the order of the matrix.
!!
!! \param      A contains the real symmetric input matrix.  Only the
!!         lower triangle of the matrix need be supplied.
!!
!!    On Output
!!
!! \param      A contains information about the orthogonal trans-
!!         formations used in the reduction in its strict lower
!!         triangle.  The full upper triangle of A is unaltered.
!!
!! \param      D contains the diagonal elements of the tridiagonal matrix.
!!
!! \param      E contains the subdiagonal elements of the tridiagonal
!!         matrix in its last N-1 positions.  E(1) is set to zero.
!!
!! \param      E2 contains the squares of the corresponding elements of E.
!!         E2 may coincide with E if the squares are not needed.
!!
!!    Questions and comments should be directed to B. S. Garbow,
!!    APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!! 
!! \note REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
!!                Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
!!                SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
!!                1976.
!<
      SUBROUTINE FMS_TRED1(NM,N,A,D,E,E2)
   
      INTEGER I,J,K,L,N,II,NM,JP1
      REAL*8 A(NM,N),D(N),E(N),E2(N)
      REAL*8 F,G,H,SCALE

!  FIRST EXECUTABLE STATEMENT  TRED1
      DO 100 I = 1, N
  100 D(I) = A(I,I)
!    .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 300 II = 1, N
         I = N + 1 - II
         L = I - 1
         H = 0.0D0
         SCALE = 0.0D0
         IF (L .LT. 1) GO TO 130
!    .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(A(I,K))

         IF (SCALE .NE. 0.0D0) GO TO 140
  130    E(I) = 0.0D0
         E2(I) = 0.0D0
         GO TO 290

  140    DO 150 K = 1, L
            A(I,K) = A(I,K) / SCALE
            H = H + A(I,K) * A(I,K)
  150    CONTINUE

         E2(I) = SCALE * SCALE * H
         F = A(I,L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         A(I,L) = F - G
         IF (L .EQ. 1) GO TO 270
         F = 0.0D0

         DO 240 J = 1, L
            G = 0.0D0
!    .......... FORM ELEMENT OF A*U ..........
            DO 180 K = 1, J
  180       G = G + A(J,K) * A(I,K)

            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220

            DO 200 K = JP1, L
  200       G = G + A(K,J) * A(I,K)
!    .......... FORM ELEMENT OF P ..........
  220       E(J) = G / H
            F = F + E(J) * A(I,J)
  240    CONTINUE

         H = F / (H + H)
!    .......... FORM REDUCED A ..........
         DO 260 J = 1, L
            F = A(I,J)
            G = E(J) - H * F
            E(J) = G

            DO 260 K = 1, J
               A(J,K) = A(J,K) - F * E(K) - G * A(I,K)
  260    CONTINUE

  270    DO 280 K = 1, L
  280    A(I,K) = SCALE * A(I,K)

  290    H = D(I)
         D(I) = A(I,I)
         A(I,I) = H
  300 CONTINUE

      RETURN
      END SUBROUTINE FMS_TRED1
!>
!! Reduces a real symmetric matrix to symmetric tridiagonal
!!          matrix using and accumulating orthogonal transformation.
!!
!!   This subroutine is a translation of the ALGOL procedure TRED2,
!!   NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
!!   HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!!
!! DATE WRITTEN   760101   (YYMMDD)
!! REVISION DATE  830518   (YYMMDD)
!! CATEGORY NO.  D4C1B1
!! KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
!! AUTHOR  SMITH, B. T., ET AL.
!!
!!
!!   This subroutine FMS_reduces a REAL SYMMETRIC matrix to a
!!   symmetric tridiagonal matrix using and accumulating
!!   orthogonal similarity transformations.
!!
!!   On Input
!!
!! \param     NM must be set to the row dimension of two-dimensional
!!        array parameters as declared in the calling program
!!        dimension statement.
!!
!! \param     N is the order of the matrix.
!!
!! \param     A contains the real symmetric input matrix.  Only the
!!        lower triangle of the matrix need be supplied.
!!
!!   On Output
!!
!! \param     D contains the diagonal elements of the tridiagonal matrix.
!!
!! \param     E contains the subdiagonal elements of the tridiagonal
!!        matrix in its last N-1 positions.  E(1) is set to zero.
!!
!! \param     Z contains the orthogonal transformation matrix
!!        produced in the reduction.
!!
!! \param     A and Z may coincide.  If distinct, A is unaltered.
!!
!!   Questions and comments should be directed to B. S. Garbow,
!!   APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!!
!! \note REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
!!               Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
!!               SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
!!               1976.
!<
      SUBROUTINE FMS_TRED2(NM,N,A,D,E,Z)
   
      INTEGER I,J,K,L,N,II,NM,JP1
      REAL*8 A(NM,N),D(N),E(N),Z(NM,N)
      REAL*8 F,G,H,HH,SCALE

! FIRST EXECUTABLE STATEMENT  TRED2
      DO 100 I = 1, N

         DO 100 J = 1, I
            Z(I,J) = A(I,J)
  100 CONTINUE

      IF (N .EQ. 1) GO TO 320
!   .......... FOR I=N STEP -1 UNTIL 2 DO -- ..........
      DO 300 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0.0D0
         SCALE = 0.0D0
         IF (L .LT. 2) GO TO 130
!   .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(Z(I,K))

         IF (SCALE .NE. 0.0D0) GO TO 140
  130    E(I) = Z(I,L)
         GO TO 290

  140    DO 150 K = 1, L
            Z(I,K) = Z(I,K) / SCALE
            H = H + Z(I,K) * Z(I,K)
  150    CONTINUE

         F = Z(I,L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         Z(I,L) = F - G
         F = 0.0D0

         DO 240 J = 1, L
            Z(J,I) = Z(I,J) / H
            G = 0.0D0
!   .......... FORM ELEMENT OF A*U ..........
            DO 180 K = 1, J
  180       G = G + Z(J,K) * Z(I,K)

            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220

            DO 200 K = JP1, L
  200       G = G + Z(K,J) * Z(I,K)
!   .......... FORM ELEMENT OF P ..........
  220       E(J) = G / H
            F = F + E(J) * Z(I,J)
  240    CONTINUE

         HH = F / (H + H)
!   .......... FORM REDUCED A ..........
         DO 260 J = 1, L
            F = Z(I,J)
            G = E(J) - HH * F
            E(J) = G

            DO 260 K = 1, J
               Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)
  260    CONTINUE

  290    D(I) = H
  300 CONTINUE

  320 D(1) = 0.0D0
      E(1) = 0.0D0
!   .......... ACCUMULATION OF TRANSFORMATION MATRICES ..........
      DO 500 I = 1, N
         L = I - 1
         IF (D(I) .EQ. 0.0D0) GO TO 380

         DO 360 J = 1, L
            G = 0.0D0

            DO 340 K = 1, L
  340       G = G + Z(I,K) * Z(K,J)

            DO 360 K = 1, L
               Z(K,J) = Z(K,J) - G * Z(K,I)
  360    CONTINUE

  380    D(I) = Z(I,I)
         Z(I,I) = 1.0D0
         IF (L .LT. 1) GO TO 500

         DO 400 J = 1, L
            Z(I,J) = 0.0D0
            Z(J,I) = 0.0D0
  400    CONTINUE

  500 CONTINUE

      RETURN
      END SUBROUTINE FMS_TRED2
!>
!!  Forms eigenvectors of complex Hermitian matrix from
!!  eigenvectors of real symmetric tridiagonal matrix output
!!  from HTRIDI.
!!
!!  KEYWORDS  LIBRARY=SLATEC(EISPACK),TYPE=SINGLE PRECISION(HTRIBK-S),
!!            EIGENVALUES,EIGENVECTORS
!!  AUTHOR  SMITH, B. T., ET AL.
!!  DATE WRITTEN   760101   (YYMMDD)
!!  REVISION DATE  861211   (YYMMDD)
!!  CATEGORY NO.  D4C4
!!  DESCRIPTION
!!
!!    This subroutine is a translation of a complex analogue of
!!    the ALGOL procedure TRBAK1, NUM. MATH. 11, 181-195(1968)
!!    by Martin, Reinsch, and Wilkinson.
!!    HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!!
!!    This subroutine forms the eigenvectors of a COMPLEX HERMITIAN
!!    matrix by back transforming those of the corresponding
!!    real symmetric tridiagonal matrix determined by  HTRIDI.
!!
!!    On INPUT
!!
!! \param      NM must be set to the row dimension of two-dimensional
!!         array parameters as declared in the calling program
!!         dimension statement.
!!
!! \param      N is the order of the matrix.
!!
!! \param      AR and AI contain information about the unitary trans-
!!         formations used in the reduction by  HTRIDI  in their
!!         full lower triangles except for the diagonal of AR.
!!
!! \param      TAU contains further information about the transformations.
!!
!! \param      M is the number of eigenvectors to be back transformed.
!!
!! \param      ZR contains the eigenvectors to be back transformed
!!         in its first M columns.
!!
!!    On OUTPUT
!!
!! \param      ZR and ZI contain the real and imaginary parts,
!!         respectively, of the transformed eigenvectors
!!         in their first M columns.
!!
!!    Note that the last component of each returned vector
!!    is real and that vector Euclidean norms are preserved.
!!
!!    Questions and comments should be directed to B. S. Garbow,
!!    APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!!
!! \note REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
!!                Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
!!                SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
!!                1976.
!<
      SUBROUTINE FMS_HTRIBK(NM,N,AR,AI,TAU,M,ZR,ZI)
   
      implicit real*8 (a-h,o-z)
      INTEGER I,J,K,L,M,N,NM
      REAL*8 AR(NM,N),AI(NM,N),TAU(2,N),ZR(NM,M),ZI(NM,M)
      REAL*8 H,S,SI

!  FIRST EXECUTABLE STATEMENT  HTRIBK
      IF (M .EQ. 0) GO TO 200
!    .......... TRANSFORM THE EIGENVECTORS OF THE REAL SYMMETRIC
!               TRIDIAGONAL MATRIX TO THOSE OF THE HERMITIAN
!               TRIDIAGONAL MATRIX. ..........
      DO 50 K = 1, N

         DO 50 J = 1, M
            ZI(K,J) = -ZR(K,J) * TAU(2,K)
            ZR(K,J) = ZR(K,J) * TAU(1,K)
   50 CONTINUE

      IF (N .EQ. 1) GO TO 200
!    .......... RECOVER AND APPLY THE HOUSEHOLDER MATRICES ..........
      DO 140 I = 2, N
         L = I - 1
         H = AI(I,I)
         IF (H .EQ. 0.0D0) GO TO 140

         DO 130 J = 1, M
            S = 0.0D0
            SI = 0.0D0

            DO 110 K = 1, L
               S = S + AR(I,K) * ZR(K,J) - AI(I,K) * ZI(K,J)
               SI = SI + AR(I,K) * ZI(K,J) + AI(I,K) * ZR(K,J)
  110       CONTINUE
!    .......... DOUBLE DIVISIONS AVOID POSSIBLE UNDERFLOW ..........
            S = (S / H) / H
            SI = (SI / H) / H

            DO 120 K = 1, L
               ZR(K,J) = ZR(K,J) - S * AR(I,K) - SI * AI(I,K)
               ZI(K,J) = ZI(K,J) - SI * AR(I,K) + S * AI(I,K)
  120       CONTINUE

  130    CONTINUE

  140 CONTINUE

  200 RETURN
      END SUBROUTINE FMS_HTRIBK
!>
!! Reduces complex Hermitian matrix to real symmetric
!! tridiagonal matrix using unitary similarity
!! transformations.
!!
!!  DATE WRITTEN   760101   (YYMMDD)
!!  REVISION DATE  861211   (YYMMDD)
!!  CATEGORY NO.  D4C1B1
!!  KEYWORDS  LIBRARY=SLATEC(EISPACK),TYPE=SINGLE PRECISION(HTRIDI-S),
!!            EIGENVALUES,EIGENVECTORS
!!  AUTHOR  SMITH, B. T., ET AL.
!!  DESCRIPTION
!!
!!    This subroutine is a translation of a complex analogue of
!!    the ALGOL procedure TRED1, NUM. MATH. 11, 181-195(1968)
!!    by Martin, Reinsch, and Wilkinson.
!!    HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!!
!!    This subroutine reduces a COMPLEX HERMITIAN matrix
!!    to a real symmetric tridiagonal matrix using
!!    unitary similarity transformations.
!!
!!    On INPUT
!!
!!    \param   NM must be set to the row dimension of two-dimensional
!!         array parameters as declared in the calling program
!!         dimension statement.
!!
!!    \param   N is the order of the matrix.
!!
!!    \param   AR and AI contain the real and imaginary parts,
!!         respectively, of the complex hermitian input matrix.
!!         Only the lower triangle of the matrix need be supplied.
!!
!!    On OUTPUT
!!
!!    \param   AR and AI contain information about the unitary trans-
!!         formations used in the reduction in their full lower
!!         triangles.  Their strict upper triangles and the
!!         diagonal of AR are unaltered.
!!
!!    \param   D contains the diagonal elements of the tridiagonal matrix.
!!
!!    \param   E contains the subdiagonal elements of the tridiagonal
!!         matrix in its last N-1 positions.  E(1) is set to zero.
!!
!!    \param   E2 contains the squares of the corresponding elements of E.
!!         E2 may coincide with E if the squares are not needed.
!!
!!    \param   TAU contains further information about the transformations.
!!
!!    Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
!!
!!    Questions and comments should be directed to B. S. Garbow,
!!    APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!!
!! \note REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
!!                Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
!!                SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
!!                1976.
!<

      SUBROUTINE FMS_HTRIDI(NM,N,AR,AI,D,E,E2,TAU)
      implicit real*8 (a-h,o-z)
      INTEGER I,J,K,L,N,II,NM,JP1
      REAL*8 AR(NM,N),AI(NM,N),D(N),E(N),E2(N),TAU(2,N)
      REAL*8 F,G,H,FI,GI,HH,SI,SCALE

!     FIRST EXECUTABLE STATEMENT  HTRIDI
      TAU(1,N) = 1.0D0
      TAU(2,N) = 0.0D0

      DO 100 I = 1, N
  100 D(I) = AR(I,I)
!     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 300 II = 1, N
         I = N + 1 - II
         L = I - 1
         H = 0.0D0
         SCALE = 0.0D0
         IF (L .LT. 1) GO TO 130
!     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(AR(I,K)) + ABS(AI(I,K))

         IF (SCALE .NE. 0.0D0) GO TO 140
         TAU(1,L) = 1.0D0
         TAU(2,L) = 0.0D0
  130    E(I) = 0.0D0
         E2(I) = 0.0D0
         GO TO 290

  140    DO 150 K = 1, L
            AR(I,K) = AR(I,K) / SCALE
            AI(I,K) = AI(I,K) / SCALE
            H = H + AR(I,K) * AR(I,K) + AI(I,K) * AI(I,K)
  150    CONTINUE

         E2(I) = SCALE * SCALE * H
         G = SQRT(H)
         E(I) = SCALE * G
         F = FMS_PYTHAG(AR(I,L),AI(I,L))
!     .......... FORM NEXT DIAGONAL ELEMENT OF MATRIX T ..........
         IF (F .EQ. 0.0D0) GO TO 160
         TAU(1,L) = (AI(I,L) * TAU(2,I) - AR(I,L) * TAU(1,I)) / F
         SI = (AR(I,L) * TAU(2,I) + AI(I,L) * TAU(1,I)) / F
         H = H + F * G
         G = 1.0D0 + G / F
         AR(I,L) = G * AR(I,L)
         AI(I,L) = G * AI(I,L)
         IF (L .EQ. 1) GO TO 270
         GO TO 170
  160    TAU(1,L) = -TAU(1,I)
         SI = TAU(2,I)
         AR(I,L) = G
  170    F = 0.0D0

         DO 240 J = 1, L
            G = 0.0D0
            GI = 0.0D0
!     .......... FORM ELEMENT OF A*U ..........
            DO 180 K = 1, J
               G = G + AR(J,K) * AR(I,K) + AI(J,K) * AI(I,K)
               GI = GI - AR(J,K) * AI(I,K) + AI(J,K) * AR(I,K)
  180       CONTINUE

            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220

            DO 200 K = JP1, L
               G = G + AR(K,J) * AR(I,K) - AI(K,J) * AI(I,K)
               GI = GI - AR(K,J) * AI(I,K) - AI(K,J) * AR(I,K)
  200       CONTINUE
!     .......... FORM ELEMENT OF P ..........
  220       E(J) = G / H
            TAU(2,J) = GI / H
            F = F + E(J) * AR(I,J) - TAU(2,J) * AI(I,J)
  240    CONTINUE

         HH = F / (H + H)
!     .......... FORM REDUCED A ..........
         DO 260 J = 1, L
            F = AR(I,J)
            G = E(J) - HH * F
            E(J) = G
            FI = -AI(I,J)
            GI = TAU(2,J) - HH * FI
            TAU(2,J) = -GI

            DO 260 K = 1, J
               AR(J,K) = AR(J,K) - F * E(K) - G * AR(I,K)
     1                           + FI * TAU(2,K) + GI * AI(I,K)
               AI(J,K) = AI(J,K) - F * TAU(2,K) - G * AI(I,K)
     1                           - FI * E(K) - GI * AR(I,K)
  260    CONTINUE

  270    DO 280 K = 1, L
            AR(I,K) = SCALE * AR(I,K)
            AI(I,K) = SCALE * AI(I,K)
  280    CONTINUE

         TAU(2,L) = -SI
  290    HH = D(I)
         D(I) = AR(I,I)
         AR(I,I) = HH
         AI(I,I) = SCALE * SQRT(H)
  300 CONTINUE

      RETURN
      END SUBROUTINE FMS_HTRIDI
!>
!!    Computes eigenvalues and, optionally, eigenvectors of
!!    real symmetric matrix.
!!
!!  DATE WRITTEN   760101   (YYMMDD)
!!  REVISION DATE  830518   (YYMMDD)
!!  CATEGORY NO.  D4A1
!!  KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
!!  AUTHOR  SMITH, B. T., ET AL.
!!  DESCRIPTION
!!    This subroutine FMS_calls the recommended sequence of
!!    subroutines from the eigensystem subroutine FMS_package (EISPACK)
!!    to find the eigenvalues and eigenvectors (if desired)
!!    of a REAL SYMMETRIC matrix.
!!    On Input
!!    \param   NM  must be set to the row dimension of the two-dimensional
!!       array parameters as declared in the calling program
!!       dimension statement.
!!    \param   N  is the order of the matrix  A.
!!    \param   A  contains the real symmetric matrix.
!!    \param   MATZ  is an integer variable set equal to zero if
!!       only eigenvalues are desired.  Otherwise it is set to
!!       any non-zero integer for both eigenvalues and eigenvectors.
!!    On Output
!!    \param   W  contains the eigenvalues in ascending order.
!!    \param   Z  contains the eigenvectors if MATZ is not zero.
!!       IERR  is an integer output variable set equal to an
!!       error completion code described in section 2B of the
!!       documentation.  The normal completion code is zero.
!!    \param   FV1  and  FV2  are temporary storage arrays.
!!    Questions and comments should be directed to B. S. Garbow,
!!    APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!!    ------------------------------------------------------------------
!! \note REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
!!                Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
!!                SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
!!                1976.
!<
      SUBROUTINE FMS_RS(NM,N,A,W,MATZ,Z,FV1,FV2,IERR)
      INTEGER N,NM,IERR,MATZ
      REAL*8 A(NM,N),W(N),Z(NM,N),FV1(N),FV2(N)

!   FIRST EXECUTABLE STATEMENT  RS
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50

   10 IF (MATZ .NE. 0) GO TO 20
!     .......... FIND EIGENVALUES ONLY ..........
      CALL FMS_TRED1(NM,N,A,W,FV1,FV2)
      CALL FMS_TQLRAT(N,W,FV2,IERR)
      GO TO 50
!     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 CALL FMS_TRED2(NM,N,A,W,FV1,Z)
      CALL FMS_TQL2(NM,N,W,FV1,Z,IERR)
   50 RETURN
      END SUBROUTINE FMS_RS
!>
!!    This subroutine calls the recommended sequence of
!!    subroutines from the eigensystem subroutine FMS_package (EISPACK)
!!    to find the eigenvalues and eigenvectors (if desired)
!!    of a COMPLEX HERMITIAN matrix.
!!
!!  DATE WRITTEN   760101   (YYMMDD)
!!  REVISION DATE  861211   (YYMMDD)
!!  CATEGORY NO.  D4A3
!!  KEYWORDS  LIBRARY=SLATEC(EISPACK),TYPE=COMPLEX(RS-S CH-C),
!!            EIGENVALUES,EIGENVECTORS
!!  AUTHOR  SMITH, B. T., ET AL.
!!  PURPOSE  Computes the eigenvalues and, optionally, eigenvecto
!!           a complex Hermitian matrix.
!!  DESCRIPTION
!!
!!    On INPUT
!!
!! \param      NM  must be set to the row dimension of the two-dimensional
!!       array parameters as declared in the calling program
!!       dimension statement.
!!
!! \param      N  is the order of the matrix  A=(AR,AI).
!!
!! \param     AR  and  AI  contain the real and imaginary parts,
!!       respectively, of the complex hermitian matrix.
!!
!! \param      MATZ  is an integer variable set equal to zero if
!!       only eigenvalues are desired.  Otherwise it is set to
!!       any non-zero integer for both eigenvalues and eigenvectors.
!!
!!    On OUTPUT
!!
!! \param      W  contains the eigenvalues in ascending order.
!!
!! \param      ZR  and  ZI  contain the real and imaginary parts,
!!       respectively, of the eigenvectors if MATZ is not zero.
!!
!! \param      IERR  is an integer output variable set equal to an
!!       error completion code described in section 2B of the
!!       documentation.  The normal completion code is zero.
!!
!! \param      FV1, FV2, and  FM1  are temporary storage arrays.
!!
!!    Questions and comments should be directed to B. S. Garbow,
!!    APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!!    ------------------------------------------------------------------
!! \note REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
!!                Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
!!                SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
!!                1976.
!<
      SUBROUTINE FMS_CH(NM,N,AR,AI,W,MATZ,ZR,ZI,FV1,FV2,FM1,IERR)
   
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER I,J,N,NM,IERR,MATZ
      REAL*8 AR(NM,N),AI(NM,N),W(N),ZR(NM,N),ZI(NM,N)
      REAL*8 FV1(N),FV2(N),FM1(2,N)

!!  FIRST EXECUTABLE STATEMENT  CH
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50

   10 CALL FMS_HTRIDI(NM,N,AR,AI,W,FV1,FV2,FM1)
      IF (MATZ .NE. 0) GO TO 20
!!    .......... FIND EIGENVALUES ONLY ..........
      CALL FMS_TQLRAT(N,W,FV2,IERR)
      GO TO 50
!!    .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 DO 40 I = 1, N

         DO 30 J = 1, N
            ZR(J,I) = 0.0D0
   30    CONTINUE

         ZR(I,I) = 1.0D0
   40 CONTINUE

      CALL FMS_TQL2(NM,N,W,FV1,ZR,IERR)
      IF (IERR .NE. 0) GO TO 50
      CALL FMS_HTRIBK(NM,N,AR,AI,FM1,N,ZR,ZI)
   50 RETURN
      END SUBROUTINE FMS_CH

      END MODULE EISPACKMODULE
