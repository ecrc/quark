      SUBROUTINE CORE_DSSSSM( M1, M2, NN, IB, K, A0, LDA0, A1, LDA1,
     $                       L0, LDL0, L1, LDL1, IPIV, INFO )

      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER            M1, M2, IB, NN, K
      INTEGER            LDL0, LDL1, LDA0, LDA1, INFO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION            A0( LDA0, * ), A1( LDA1, * )
      DOUBLE PRECISION            L0( LDL0, * ), L1( LDL1, * )
      INTEGER            IPIV( * )
*     ..
*     .. Internal variables ..
      INTEGER            II, I, IP, IM, SB
*     .. Parameters ..
      DOUBLE PRECISION            DONE, MDONE
      PARAMETER          ( DONE = 1.0D+0 )
      PARAMETER          ( MDONE = -1.0D+0 )
*     ..
*     Test the input parameters.
      INFO = 0
      IF( M1.LT.0 ) THEN
         INFO = -1
      ELSE IF( M2.LT.0 ) THEN
         INFO = -2
      ELSE IF( NN.LT.0 ) THEN
         INFO = -3
      ELSE IF( IB.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CORE_DTSTRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( M1.EQ.0 .OR. M2.EQ.0 .OR. NN.EQ.0 .OR. IB.EQ.0  .OR. K.EQ.0 )
     $   RETURN
*
      IP = 1
*
      DO 10 II = 1, K, IB
         SB = MIN( K-II+1, IB )
         DO 20 I = 1, IB
*
            IM = IPIV( IP )
            IF( IM.NE.II+I-1 ) THEN
               IM = IM - M1
               CALL DSWAP( NN, A0( II+I-1, 1 ), LDA0, 
     $                    A1( IM, 1 ), LDA1 )
            END IF
            IP = IP+1
*
 20      CONTINUE
         CALL DTRSM( 'Left', 'Lower', 'Notranspose', 'Unit',
     $              SB, NN, DONE,
     $              L0( 1, II ), LDL0,
     $              A0( II, 1 ), LDA0 )
*
         CALL DGEMM( 'Notranspose', 'Notranspose',
     $              M2, NN, SB, MDONE,
     $              L1( 1, II ), LDL1,
     $              A0( II, 1 ), LDA0,
     $              DONE, A1, LDA1 )
     $
 10   CONTINUE
*
      RETURN
*
*     End of CORE_DSSSSM.
*
      END
