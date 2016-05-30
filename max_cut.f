      PROGRAM MAC_CUT
      IMPLICIT NONE
*     PARAMETER
      INTEGER LDA, LDU, LDVT
      PARAMETER (LDA=10, LDU=10, LDVT=10)
*     To keep real dimension of matrix
      INTEGER M,N
      INTEGER K
      PARAMETER (K=1)
*     ARRAY
      DOUBLE PRECISION A(LDA,LDA), U(LDU,LDU), VT(LDVT,LDVT), S(LDA)
*     VARIABLE RESULTAT
      
      WRITE (*,*) 'Entrer la matrice triangulaire inférieure A'
      CALL READ_MPL (M,N, A, LDA)

      CALL CALC_SVD(M,N,A,LDA, U, LDU, VT, LDVT, S, LDA)
      
      CALL CALC_TILDE(K, U, M, M, LDU)
      CALL CALC_TILDE(K, VT, N, N, LDVT)
      CALL SOLVE_PB (M, U, LDU ,VT, LDVT,
     $               S, LDA, K)
      STOP
      END
************************
*      SUBROUTINE SOLVE_PB 
************************
      SUBROUTINE SOLVE_PB(MA,UTILDE,LDUTILDE,
     $       VTTILDE,LDVTTILDE, S, LDS, K)
      IMPLICIT NONE
      INTEGER MA
      INTEGER LDUTILDE, LDVTTILDE, LDS
      INTEGER K
      INTEGER I,J,Y,Z
      INTEGER NBELEM, NBELEMNEXT
      DOUBLE PRECISION UTILDE(LDUTILDE,LDUTILDE)
     $        , VTTILDE(LDVTTILDE,LDVTTILDE), S(LDS)
      DOUBLE PRECISION LW(2,2**MA, K,2)
      DOUBLE PRECISION LX(2,2**MA, MA)
      DOUBLE PRECISION W0(K,2), X0(MA)
      DOUBLE PRECISION W1(K,2), X1(MA)
*      Initialisation des tableaux
      DO I = 1, MA, 1
              LX(2,1,I) = 0
      END DO
      DO I = 1, K, 1
              DO J=1, 2, 1
                     LW(2,1,I,J) = 0
              END DO
      END DO
*     On initialise la liste LW avec un seul élement
      NBELEM = 1
*      Boucle principale
      DO I = 1, MA, 1
              DO J = 1, NBELEM,1
                     DO Y=1, K, 1
                            DO Z = 1, MA, 1
                                   LW(1,J,Y,Z) = LW(2,J,Y,Z)
                            END DO
                     END DO
              END DO
              DO J=1,NBELEM,1
                     DO Y=1, MA, 1
                            LX(1,J, Y) = LX(2,J,Y)
                     END DO
              END DO
              NBELEMNEXT = 0
              DO J = 1, NBELEM, 1
                    DO Y = 1, K, 1
                            DO Z = 1, 2, 1
                                   W0(Y,Z) = LW(1,J, Y, Z) 
                            END DO
                    END DO
                    DO Y = 1, MA, 1
                            X0(Y) = LX(1,J,Y)
                    END DO
                    X0(I) = 0
                    CALL W_UPDATE(K, W0, X0, MA, I, K, UTILDE
     $                      , LDUTILDE, MA*K**2d0, VTTILDE
     $                      , LDVTTILDE)
                     
                    DO Y = 1, K, 1
                            DO Z = 1, 2, 1
                                   W1(Y,Z) = LW(1,J, Y, Z) 
                            END DO
                    END DO
                    DO Y = 1, MA, 1
                            X1(Y) = LX(1,J,Y)
                    END DO
                    X1(I) = 1
                    CALL W_UPDATE(K, W1, X1, MA, I, K, UTILDE
     $                      , LDUTILDE, MA*K**2d0, VTTILDE
     $                      , LDVTTILDE)
                    CALL ADD (2, 2**MA, K,2, LW, K, 2, W0, MA,
     $                    X0,2, 2**MA,MA, LX,NBELEMNEXT) 
                    CALL ADD (2, 2**MA, K,2, LW, K, 2, W1, MA,
     $                    X1,2, 2**MA,MA, LX,NBELEMNEXT) 
              END DO
              NBELEM = NBELEMNEXT
      END DO
      CALL CALC_SOL(2, 2**MA,K,2,LW, 2, 2**MA,MA, LX,S,LDS
     $               ,NBELEMNEXT)
      END
************************
*      SUBROUTINE CALC_SOL
************************
      SUBROUTINE CALC_SOL(MLW,NLW,OLW,PLW,LW, MLX,NLX,OLX,LX,S,LDS,
     $                     NBELEM)
      IMPLICIT NONE
*     Parameters
      INTEGER MLW,NLW,OLW,PLW
      DOUBLE PRECISION LW(MLW,NLW,OLW,PLW)
      INTEGER MLX,NLX,OLX
      DOUBLE PRECISION LX(MLX,NLX,OLX)
      INTEGER LDS
      DOUBLE PRECISION S(LDS)
      INTEGER NBELEM
*     Local
      DOUBLE PRECISION MAP(NBELEM)
      INTEGER BEST, I
      DOUBLE PRECISION VAL
*     Code      
      CALL VALEUR(MAP,NBELEM, OLW,S,LDS, 
     $        MLW, NLW,OLW,PLW,LW)
      BEST = 1
      VAL = MAP(1)
      DO I = 2, NBELEM,1
              IF (MAP(I) > VAL) THEN
                     BEST = I
                     VAL = MAP(I)
              END IF
      END DO
      PRINT *, ''
      PRINT *, 'BEST := ', BEST
      PRINT *, 'Solution'
      PRINT *, ''
      DO I = 1,OLX,1
             PRINT *, NINT(LX(2,BEST,I))
      END DO
      PRINT *, ''
      PRINT *, 'Sommet de S1'
      DO I = 1,OLX,1
             IF (NINT(LX(2,BEST,I)) == 1) THEN
                     PRINT *,I
             END IF 
      END DO

      PRINT *, 'Sommet de S2'
      DO I = 1,OLX,1
             IF (NINT(LX(2,BEST,I)) == 0) THEN
                     PRINT *,I
             END IF
      END DO 
      END
************************
*      SUBROUTINE VALEUR
************************
      SUBROUTINE VALEUR(MAP,NBELEM, K,S,LDS, 
     $        MLW, NLW,OLW,PLW,LW)
      IMPLICIT NONE
*     Parameter
      INTEGER MLW,NLW,OLW,PLW
      DOUBLE PRECISION LW(MLW,NLW,OLW,PLW)
      INTEGER NBELEM
      DOUBLE PRECISION MAP(NBELEM)
      INTEGER LDS
      DOUBLE PRECISION S(LDS)
      INTEGER K
*     Local
      INTEGER I, J
*     Code
      DO I = 1,NBELEM, 1
          MAP(I) = 0
      END DO
      DO I =1,NBELEM,1
          DO J=1, K,1
              MAP(I) = MAP(I) + LW(2,I,J,1) * S(J) * LW(2,I,J,2)
          END DO
      END DO
      END 
************************
*      SUBROUTINE ADD
************************
      SUBROUTINE ADD (SLW, MLW, NLW, OLW, LW ,MW, NW, W,NX, X
     $           ,MLX,NLX,OLX,LX, NBELEMNEXT)
      IMPLICIT NONE
*     Parameter
      INTEGER MLW,NLW,MW,NW, OLW, SLW, MLX,NLX,OLX
      INTEGER NX
      DOUBLE PRECISION LW(SLW, MLW,NLW,OLW), W(MW,NW), X(NX)
      DOUBLE PRECISION LX(MLX,NLX,OLX)
      INTEGER NBELEMNEXT
*     Local
      INTEGER I, J, K
      LOGICAL PRE
*     Code
      PRE = .FALSE.
      DO I = 1, NBELEMNEXT, 1
          DO J = 1, NLW, 1
              DO K = 1, OLW, 1
                 IF (LW(2, I,J,K) .ne. W(J,K)) THEN
                     PRE = .TRUE.
                 END IF 
              END DO
          END DO
          IF (PRE .eqv. .FALSE.) THEN
              GO TO 1
          END IF
          PRE = .FALSE.
      END DO
      IF (.NOT. PRE) THEN
              NBELEMNEXT = NBELEMNEXT + 1
              DO J = 1, NLW, 1
                     DO K = 1, OLW, 1
                            LW(2,NBELEMNEXT,J,K) = W(J,K)
                     END DO
              END DO
              DO J = 1, OLX, 1
                     LX(2, NBELEMNEXT,J) = X(J)
              END DO
      END IF
 1    END
************************
*      SUBROUTINE W_UPDATE
************************
      SUBROUTINE W_UPDATE (LDW, W, X, SX, N, K,UTILDE,LDUTILDE, 
     $  FACT, VTTILDE,LDVTTILDE)
      IMPLICIT NONE
      INTEGER SX, LDW, N, K
      DOUBLE PRECISION FACT
      DOUBLE PRECISION X(SX)
      DOUBLE PRECISION W(LDW,2)
      INTEGER LDUTILDE
      DOUBLE PRECISION UTILDE(LDUTILDE,LDUTILDE)
      INTEGER LDVTTILDE
      DOUBLE PRECISION VTTILDE(LDVTTILDE,LDVTTILDE)
      INTEGER I
      IF (X(N) == 1) THEN
             DO I=1, K,1
                     W(I,1) = W(I,1) + UTILDE(N,I)/FACT
             END DO
      ELSE
              DO I=1,K,1
                     W(I,2) = W(I,2) + VTTILDE(I,N)/FACT
              END DO
      ENDIF
      END
************************
*      SUBROUTINE TILDE (MATRIX BY SCALAR AND CONVERSION TO THE NEAREST INT)
************************
      SUBROUTINE CALC_TILDE(K, A, M, N, LDA)
      IMPLICIT NONE
      INTEGER K, M, N, LDA
      DOUBLE PRECISION A(LDA,*)
      INTEGER I, J
      DOUBLE PRECISION FACT
      FACT = M*K**2d0
      DO I=1, N, 1
             CALL DSCAL(M, FACT,A(1,I), 1); 
      END DO
      DO I=1, N, 1
              DO J=1, M, 1
                     A(J,I) = NINT(A(J,I))
              END DO
      END DO
      END
************************
*      SUBROUTINE SVD
************************
      SUBROUTINE CALC_SVD(M,N,A,LDA, U, LDU, VT, LDVT, S, LDS)
      IMPLICIT NONE
*     Local integer      
      DOUBLE PRECISION A(LDA,LDA), U(LDU,LDU), VT(LDVT,LDVT), S(LDS)
      INTEGER LDA, LDU, LDVT, LDS
      INTEGER M,N
      INTEGER INFO, LWORK, IWORK (8*LDA)
      INTEGER LWMAX
      PARAMETER (LWMAX=1000)
      DOUBLE PRECISION WORK(LWMAX)
      LWORK = -1     
      CALL DGESDD ('Singular vectors', M,N, A, LDA, S, U, LDU, VT, LDVT,
     $        WORK, LWORK, IWORK, INFO)
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
*
*     Compute SVD.
*
      CALL DGESDD('Singular vectors', M, N, A, LDA, S, U, LDU, VT, LDVT,
     $             WORK, LWORK, IWORK, INFO )
*
*
*     Check for convergence.
*
      IF( INFO.GT.0 ) THEN
          WRITE(*,*)'The algorithm computing SVD failed to converge.'
          STOP
      END IF
      END
