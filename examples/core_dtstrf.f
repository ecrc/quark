      SUBROUTINE CORE_DTSTRF( M, N, IB, NB, U, LDU, A, LDA,
     $                        L, LDL, IPIV, INFO )

      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER            M, N, IB, NB, LDU, LDA, LDL, INFO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION            U( LDU, * ), A( LDA, * ), L( LDL, * )
      INTEGER            IPIV( * )
*     ..
*     Internal variables ..
      DOUBLE PRECISION            Ltmp( LDA, N )
      INTEGER            I, II, J, SB, IM, IP, IDAMAX, IINFO
      EXTERNAL           IDAMAX
*     ..
*     .. Parameters ..
      DOUBLE PRECISION            MDONE, DZERO
      PARAMETER          ( MDONE = -1.0D+0 )
      PARAMETER          ( DZERO = 0.0D+0 )
*     ..
*     Test the input parameters.
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( IB.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDU.LT.MAX( 1, M ) ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -8
      ELSE IF( LDL.LT.MAX( 1, IB ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CORE_DTSTRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. IB.EQ.0 )
     $   RETURN
*
      IP = 1
*
      DO 20 II = 1, N, IB
         SB = MIN( N-II+1, IB )
         DO 10 I = 1, SB
*
            IM = IDAMAX( M, A( 1, II+I-1 ), 1 )
            IPIV( IP ) = II+I-1
*
            IF( ABS( A( IM, II+I-1 ) ) 
     $          .GT. ABS( U( II+I-1, II+I-1 ) ) ) THEN
*
*              Swap behind.
*
               CALL DSWAP( I-1, L( I, II ), LDL,
     $                    Ltmp( IM, II ),LDA )
*
*              Swap ahead.
*
               CALL DSWAP( SB-I+1, U( II+I-1, II+I-1 ), LDU,
     $                    A( IM, II+I-1 ), LDA )
*
*              Set IPIV.
*
               IPIV( IP ) = NB + IM
*
               DO 50 J = 1, I-1
                  A( IM, II+J-1 ) = DZERO
 50            CONTINUE
*
            END IF
            IF( INFO.EQ.0 .AND. ABS(A(IM,II+I-1)).EQ.DZERO 
     $          .AND. ABS(U(II+I-1, II+I-1)).EQ.DZERO ) THEN
               INFO = II+I-1
            END IF
*
            CALL DSCAL( M, 1/U( II+I-1, II+I-1 ), A( 1, II+I-1 ), 1 )
            CALL DCOPY( M, A( 1, II+I-1 ), 1, Ltmp( 1, II+I-1 ), 1 )
            CALL DGER( M, SB-I, MDONE, A( 1, II+I-1 ), 1,
     $                U( II+I-1, II+I ), LDU,
     $                A( 1, II+I ), LDA )
*
            IP = IP+1
 10      CONTINUE
*
*        Apply the subpanel to the rest of the panel.
*
         IF( II+I-1.LE. N ) THEN
            DO 80 J = II, II+SB-1
               IF( IPIV( J ).LE.NB ) THEN
                  IPIV( J ) = IPIV( J )-II+1
               ENDIF
 80         CONTINUE
*
            CALL CORE_DSSSSM( NB, M, N-( II+SB-1 ), SB, SB,
     $                       U( II, II+SB ), LDU,
     $                       A( 1, II+SB ), LDA,
     $                       L( 1, II ), LDL,
     $                       Ltmp( 1, II ), LDA, IPIV( II ),
     $                       IINFO )
*
            DO 70 J = II, II+SB-1
               IF( IPIV( J ).LE.NB ) THEN
                  IPIV( J ) = IPIV( J )+II-1
               ENDIF
 70         CONTINUE
         END IF
 20   CONTINUE
*
      RETURN
*
*     End of CORE_DTSTRF.
*
      END
