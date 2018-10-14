      SUBROUTINE CORE_DGESSM( M, N, K, IB, IPIV, L, LDL, A, LDA, INFO )

      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER            M, N, K, IB, LDL, LDA, INFO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION            L(LDL,*), A(LDA,*)
      INTEGER            IPIV( * )
*     ..
*     Internal variables ..
      INTEGER            J, SB
*     ..
*     .. Parameters ..
      DOUBLE PRECISION   DONE, MDONE
      PARAMETER          ( DONE = 1.0D+0 )
      PARAMETER          ( MDONE = -1.0D+0 )
*     ..
*     Test the input parameters.
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( K.LT.0 ) THEN
         INFO = -3
      ELSE IF( IB.LT.0 ) THEN
         INFO = -4      
      ELSE IF( LDL.LT.MAX( 1, M ) ) THEN
         INFO = -7
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CORE_DGESSM', -INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 .OR. IB.EQ.0  )
     $   RETURN
*
      DO 10 J = 1, K, IB
         SB = MIN( K-J+1, IB )
*
*        Apply interchanges to columns J*IB+1:IB*( J+1 )+1.
*
         CALL DLASWP( N, A, LDA, J, J+SB-1, IPIV, 1 )
*
*        Compute block row of U.
*
         CALL DTRSM( 'Left', 'Lower', 'Notranspose', 'Unit', SB,
     $              N, DONE, L( J, J ), LDL,
     $              A( J, 1 ), LDA )
*
         IF( J+SB.LE.M ) THEN
*
*           Update trailing submatrix.
*
            CALL DGEMM( 'Notranspose', 'Notranspose', M-( J+SB-1 ),
     $                 N, SB, MDONE,
     $                 L( J+SB, J ), LDL,
     $                 A( J, 1 ), LDA, DONE,
     $                 A( J+SB, 1 ), LDA )
         END IF
   10 CONTINUE
*
      RETURN
*
*     End of CORE_DGESSM.
*
      END
