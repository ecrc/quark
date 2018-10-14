      SUBROUTINE CORE_DGETRF( M, N, IB, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N, IB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION            A( LDA, * )
      INTEGER            IPIV( * )
*     ..
*
*  Purpose
*  =======
*
*  DGETF2 computes an LU factorization of a general m-by-n matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 2 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION               SFMIN 
      INTEGER            I, J, K, SB, IINFO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION               DLAMCH
      INTEGER            IDAMAX
      EXTERNAL           DLAMCH, IDAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGETF2, CORE_DGESSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     Test the input parameters.
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( IB.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CORE_DGETRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. IB.EQ.0 )
     $   RETURN
*
*     Compute machine safe minimum.
*
      SFMIN = DLAMCH('S')
*
      K = MIN( M, N )
*
      DO 10 I = 1, K, IB
            SB = MIN( K-I+1, IB )
*
*           Factor diagonal and subdiagonal blocks and test for exact
*           singularity.
*
            CALL DGETF2( M-I+1, SB, A( I, I ), LDA,
     $                   IPIV( I ), IINFO )
*
*           Adjust INFO and the pivot indices.
*
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $         INFO = IINFO + I-1
*
            IF( I+SB.LE.N ) THEN
               CALL CORE_DGESSM( M-I+1, N-( I+SB-1 ), SB, SB,
     $                           IPIV( I ), A( I, I ),
     $                           LDA, A( I, I+SB ), LDA,
     $                           IINFO)
            END IF
*
            DO 20 J = I, I+SB-1
               IPIV( J ) = I + IPIV( J ) -1
   20       CONTINUE
   10 CONTINUE
*
      RETURN
*
*     End of CORE_DGETRF.
*
      END
