///////////////////////////////////////////////////////////////////////////////////////////////////
/*                                                                                               */
/* INPUT:  ./tileChInv (MATRIX_SIZE) (TILE_SIZE) (NB_THREADS)                                    */
/*                                                                                               */
/* OUTPUT:                                                                                       */
/*    when (#undef  SEQCHECK) (#define NUMCHECK)                                                 */
/*    QUARK CHINV (MATRIX_SIZE) (TILE_SIZE) (NB_THREADS) (PERF) (CHECK)                           */
/*                                                                                               */
/*    (PERF)  is in GFlops/sec, # of flops is n^3, time is taken with get_current_time()         */
/*    (CHECK) should be less than 1.00, check is computed by LAPACK DPOT03 (sequential) and      */
/*            computes norm( I - A * Ainv, 1 ) / ( N * norm( A, 1) * norm( Ainv, 1 ) * EPS )     */
/*                                                                                               */
/* LIMITATION:                                                                                   */
/*    (MATRIX_SIZE) needs to be a multiple of (TILE_SIZE)                                        */
/*    ONLY LOWER (no UPPER)                                                                      */
/*                                                                                               */
/* NOTE                                                                                          */
/*     When you compare QUARK with different thread numbers of against sequential version, the    */
/*     output matrices are exactly the same. This is expected. Same for comparing                */
/*     INPLACE/OUTOFPLACE                                                                        */
/*                                                                                               */
///////////////////////////////////////////////////////////////////////////////////////////////////

// (replace-string "quark new(" "quark_new(" nil nil nil)
// (replace-string "delete(quark)" "quark_delete(quark)" nil nil nil)
// (replace-string "quark->Insert_Task(" "quark_insert_task(quark, " nil nil nil)
// (replace-string "quark->rank()" "quark->rank(quark)" nil nil nil)


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>
#include <fcntl.h>
#include <sys/mman.h>

//#include "Quark.hh"
#include "quark.h"

#undef  DEBUG
#undef  ALLINOUT

#define NUMCHECK           // NUMERICAL CHECK
#undef SEQCHECK           // CHECK AGAINST SEQUENTIAL EXECUTION (You need to have NUMCHECK on for this to happen!)
#undef OUTOFPLACE         // (if "undef OUTOFPLACE" then you are INPLACE)
#define PIPELINE           // to break the pipeline of quark taks: undef the PIPELINE variable

#define INNERLOOP_POTRF_UP
#undef  INNERLOOP_TRTRI_UP
#define INNERLOOP_LAUUM_UP

#if defined ALLINOUT // EA2ALL: FIX for handling antidependecies bug.
#undef INPUT
#undef OUTPUT
#define INPUT  0x03
#define OUTPUT 0x03
#endif // ALLINOUT

///////////////////////////////////////////////////////////////////////////////////////////////////
#include <time.h>
#define MAX_THREADS 48
#define MAX_EVENTS 16384
int    event_num        [MAX_THREADS];
double event_start_time [MAX_THREADS];
double event_end_time   [MAX_THREADS];
double event_log        [MAX_THREADS][MAX_EVENTS];
int log_events = 1;
#define core_event_start(my_core_id)                    \
  event_start_time[my_core_id] = get_current_time();
#define core_event_end(my_core_id)                      \
  event_end_time[my_core_id] = get_current_time();
#define core_log_event(event, my_core_id, tag)                          \
  event_log[my_core_id][event_num[my_core_id]+0] = tag;                 \
  event_log[my_core_id][event_num[my_core_id]+1] = event_start_time[my_core_id]; \
  event_log[my_core_id][event_num[my_core_id]+2] = event_end_time[my_core_id]; \
  event_log[my_core_id][event_num[my_core_id]+3] = (event);             \
  event_num[my_core_id] += (log_events << 2);                           \
  if (event_num[my_core_id] >= MAX_EVENTS-5) event_num[my_core_id] -= (log_events << 2);
//  event_num[my_core_id] &= (MAX_EVENTS-1);
/* Return current time */
double get_current_time(void) {
  struct timeval time_val;
  gettimeofday(&time_val, NULL);
  return (double) (time_val.tv_sec) + (double) (time_val.tv_usec) / 1000000.0;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
void diff_matrix(double *A, double *B, int NB, int BBM, int BBN, int M, int N)
{
    int X, Y, x, y, i;

    printf("\n");
    for (Y = 0; Y < BBM; Y++) {
      for (y = 0; y < NB; y++) {
        for (X = 0; X < BBN; X++) {
          for (x = 0; x < NB; x++) {
            if (Y*NB + y < M && X*NB + x < N) {

              double a, b, c, d, e;
              a = fabs(A[(Y*NB+y) + (X*NB+x)*M]);
              b = fabs(B[(Y*NB+y) + (X*NB+x)*M]);
              c = max(a, b);
              d = min(a, b);
              e = (c - d) / d;

              if (c == 0.0 && d == 0.0)
                printf(".");
              else
                printf("%c", e < 0.00000000001 ? '.' : '#');

              if (x == 2) x = NB-2;
            }
            else
            {
              printf("%c", '+');
            }
          }
          printf("|");
        }
        printf("\n");
        if (y == 2) y = NB-2;
      }
      if (Y < BBM-1)
        for (i = 0; i < BBM*5; i++) printf("=");
      printf("\n");
    }
    printf("\n");
}

void dump_trace(int cores_num)
{
    char trace_file_name[32];
    FILE *trace_file;
    int event;
    int core;
//  double scale = 500000.0;
//  double scale = 300.0;
    double scale = 150000.0;

    //sprintf(trace_file_name, "trace_%d.svg", (int)(time(NULL)));
    sprintf(trace_file_name, "trace.svg");
    trace_file = fopen(trace_file_name, "w");
    assert(trace_file != NULL);

    fprintf(trace_file,
        "<?xml version=\"1.0\" standalone=\"no\"?>\n"
        "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n"
        "<svg width=\"200mm\" height=\"40mm\" viewBox=\"0 0 20000 4000\" version=\"1.1\" \n"
        "xmlns=\"http://www.w3.org/2000/svg\">\n"
        "  <g font-size=\"20\">\n");

    for (core = 0; core < cores_num; core++)
        for (event = 0; event < event_num[core]; event += 4)
        {
            int    tag   = (int)event_log[core][event+0];
            double start =      event_log[core][event+1];
            double end   =      event_log[core][event+2];
            int    color = (int)event_log[core][event+3];

            start -= event_log[0][1];
            end   -= event_log[0][1];

            fprintf(trace_file,
                "    "
                "<rect x=\"%.2lf\" y=\"%.0lf\" width=\"%.2lf\" height=\"%.0lf\" "
                "fill=\"#%06x\" stroke=\"#000000\" stroke-width=\"1\"/>\n",
                start * scale,
                core * 100.0,
                (end - start) * scale,
                90.0,
                color);

            fprintf(trace_file,
                "    "
                "<text x=\"%.2lf\" y=\"%.0lf\" font-size=\"20\" fill=\"black\">"
                "%d"
                "</text>\n",
                start * scale + 10,
                core * 100.0 + 20,
                (int)tag);
        }

    fprintf(trace_file,
        "  </g>\n"
        "</svg>\n");

    fclose(trace_file);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif
void dpotrf_( char *uplo, int *n, double *a, int *lda, int *info );
void dsyrk_( char *uplo, char *trans, int *n, int *k, double *alpha, double *a, int *lda, double *beta, double *c, int *ldc);
void dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
void dtrsm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n, double *alpha, double *a, int *lda, double *b, int *ldb);
double ddot_( int *n, double *x, int *incx, double *dy, int *incy);

void dtrtri_( char *uplo, char *diag, int *n, double *a, int *lda, int *info );
void dtrmm_( char *side, char *uplo, char *transa, char *diag, int *m, int *n, double *alpha, double *a, int *lda, double *b, int *ldb);
void dlauum_( char *uplo, int *n, double *a, int *lda, int *info );
#if defined OUTOFPLACE
void dcopy_( int *n,double *dx,int *incx,double *dy,int *incy);
#endif // OUTOFPLACE
void dpot03_( char *uplo, int *n, double *a, int *lda, double *ainv, int *ldainv, double *work, int *ldwork, double *rwork, double *rcond, double *resid );
double dlansy_( char *norm, char *uplo, int *n, double *a, int *lda, double *work );

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////
void SCHED_dgemm(Quark *quark)
{

  char transa;
  char transb;
  int m;
  int n;
  int k; 
  double alpha;
  double *A; 
  int lda;
  double *B;
  int ldb;
  double beta;
  double *C;
  int ldc;

  quark_unpack_args_13(quark, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);

  core_event_start(QUARK_Thread_Rank(quark));
  dgemm_(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
  core_event_end(QUARK_Thread_Rank(quark));
  core_log_event(0xD0F040, QUARK_Thread_Rank(quark), n);
}

void SCHED_dsyrk(Quark *quark)
{

  char uplo;
  char trans;
  int n;
  int k;
  double alpha;
  double *A;
  int lda;
  double beta;
  double *C;
  int ldc;


  quark_unpack_args_10(quark, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);

  core_event_start(QUARK_Thread_Rank(quark));
  dsyrk_(&uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc);
  core_event_end(QUARK_Thread_Rank(quark));
  core_log_event(0x00C0F0, QUARK_Thread_Rank(quark), n);
}

void SCHED_dpotrf(Quark *quark) {
  char uplo;
  int n;
  double *A;
  int lda;
  int info;

  quark_unpack_args_5(quark, uplo, n, A, lda, info);

  core_event_start(QUARK_Thread_Rank(quark));
  dpotrf_(&uplo, &n, A, &lda, &info);
  core_event_end(QUARK_Thread_Rank(quark));
  core_log_event(0x006680, QUARK_Thread_Rank(quark), n);
}

void SCHED_dtrsm(Quark *quark) {
  char side;
  char uplo;
  char transa;
  char diag;
  int m;
  int n;
  double alpha;
  double *A;
  int lda;
  double *B;
  int ldb;

  quark_unpack_args_11(quark, side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb);

  core_event_start(QUARK_Thread_Rank(quark));
  dtrsm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb);
  core_event_end(QUARK_Thread_Rank(quark));
  core_log_event(0xB060D0, QUARK_Thread_Rank(quark), n);
}

void SCHED_dtrtri(Quark *quark) {
  char uplo;
  char diag;
  int n;
  double *A;
  int lda;
  int info;
  
  quark_unpack_args_6(quark, uplo, diag, n, A, lda, info); //TODO_long_term: get info out of the call

  core_event_start(QUARK_Thread_Rank(quark));
  dtrtri_(&uplo, &diag, &n, A, &lda, &info);
  core_event_end(QUARK_Thread_Rank(quark));
  core_log_event(0xB060D0, QUARK_Thread_Rank(quark), n); //TODO:change color
}

void SCHED_dtrmm(Quark *quark) {
  char side;
  char uplo;
  char transa;
  char diag;
  int m;
  int n;
  double alpha;
  double *A;
  int lda; 
  double *B;
  int ldb;
  quark_unpack_args_11(quark, side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb);

  core_event_start(QUARK_Thread_Rank(quark));
  dtrmm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb);
  core_event_end(QUARK_Thread_Rank(quark));
  core_log_event(0xB060D0, QUARK_Thread_Rank(quark), n);

}

void SCHED_dlauum(Quark *quark) {
  char uplo;
  int n; 
  double *A;
  int lda;
  int info;

  quark_unpack_args_5(quark, uplo, n, A, lda, info);

  core_event_start(QUARK_Thread_Rank(quark));
  dlauum_(&uplo, &n, A, &lda, &info);
  core_event_end(QUARK_Thread_Rank(quark));
  core_log_event(0xB060D0, QUARK_Thread_Rank(quark), n);

}

#if defined OUTOFPLACE
void SCHED_dcopy(Quark *quark)
{

        int n;
        double *X;
        int incx;
        double *Y;
        int incy;

        quark_unpack_args_5( quark, n, X, incx, Y, incy );

        core_event_start(QUARK_Thread_Rank(quark));
        dcopy_ ( &n, X, &incx, Y, &incy );
        core_event_end(QUARK_Thread_Rank(quark));
        core_log_event(0xD0F040, QUARK_Thread_Rank(quark), n);

}
#endif // OUTOFPLACE


///////////////////////////////////////////////////////////////////////////////////////////////////
#define A(m,n) &A[NBNBSIZE*(m)+NBNBSIZE*MT*(n)]
#if defined OUTOFPLACE
#define Ac(m,n) &Ac[NBNBSIZE*(m)+NBNBSIZE*MT*(n)]
#define Ai(m,n) &Ai[NBNBSIZE*(m)+NBNBSIZE*MT*(n)]
#else // INPLACE
#define Ac(m,n) &A[NBNBSIZE*(m)+NBNBSIZE*MT*(n)]
#define Ai(m,n) &A[NBNBSIZE*(m)+NBNBSIZE*MT*(n)]
#endif // INPLACE or OUTOFPLACE

#if defined SEQCHECK
void tile_ch_sequential(double *A, int n, int nb_default)
{

        int i, j, k, p;
        int *nb;
        int info;
        double mone=-1.0e+00;
        double pone=+1.0e+00;
        char c_lower = 'L';
        char c_left = 'L';
        char c_right = 'R';
        char c_notranspose = 'N';
        char c_transpose = 'T';
        char c_nonunit = 'N';

        struct timeval tp;
        double t1,t2,elapsed;
        int rtn;

        p = n / nb_default; if( p*nb_default != n ) p++;
        nb = (int *)malloc(p*sizeof(int)) ;
        for (i=0;i<p-1;i++) nb[i] = nb_default; nb[p-1] = n - nb_default*(p-1);

        int NBNBSIZE = nb_default*nb_default;
        int MT = p;

    /************************************************************************************************/
    /* Cholesky factorization (LAPACK POTRF)                                                        */
    /************************************************************************************************/

        rtn=gettimeofday(&tp, NULL);
        t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;

        for ( j = 0 ; j < p ; j++ ){

          for ( k = 0 ; k < j ; k++ ) {
#ifdef DEBUG        
            printf("dsyrk %d %d\n", j, k);
#endif // DEBUG
            dsyrk_( &c_lower, &c_notranspose, &nb[j], &nb[k], &mone, A(j,k), &nb[j], &pone, A(j,j), &nb[j]);
          }

#ifdef DEBUG        
          printf("dpotrf %d\n", j);
#endif // DEBUG
          dpotrf_( &c_lower, &nb[j], A(j,j), &nb[j], &info );

          for ( i = j+1 ; i < p ; i++ ) {

#ifdef INNERLOOP_POTRF_UP
            for ( k = 0 ; k < j ; k++ ) {
#else // INNERLOOP_POTRF_UP
            for ( k = j-1 ; k  > -1 ; k-- ) {
#endif // INNERLOOP_POTRF_UP
              
#ifdef DEBUG        
              printf("dgemm %d %d %d\n", j, i, k);
#endif // DEBUG
              dgemm_( &c_notranspose, &c_transpose, &nb[i], &nb[j], &nb[k], &mone, A(i,k), &nb[i], A(j,k), &nb[j], &pone, A(i,j), &nb[i]);
            }
          }

          for ( i = j+1 ; i < p ; i++ ) {
#ifdef DEBUG        
            printf("dtrsm %d %d\n", j, i);
#endif // DEBUG
            dtrsm_( &c_right, &c_lower, &c_transpose, &c_nonunit, &nb[i], &nb[j], &pone, A(j,j), &nb[j], A(i,j), &nb[i]);
          }
          
        }

        rtn=gettimeofday(&tp, NULL);
        t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
        elapsed=t2-t1;
        /* printf("%fsec\n",elapsed); */

        /************************************************************************************************/
        /* Inplace Triangular Inversion (LAPACK TRTRI)                                                  */
        /************************************************************************************************/

        rtn=gettimeofday(&tp, NULL);
        t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
        
        for ( j = p-1 ; j > -1 ; j-- ){
          
          //dtrtri_( &c_lower, &c_nonunit, &nb[j], A(j,j), &nb[j],  &info );

          for ( i = p-1 ; i > j ; i-- ){
            
            dtrmm_( &c_left, &c_lower, &c_notranspose, &c_nonunit, &nb[i], &nb[j], &pone, A(i,i), &nb[i], A(i,j), &nb[i] );
            
#ifdef INNERLOOP_TRTRI_UP
            for ( k = j+1 ; k < i ; k++ ){
#else // INNERLOOP_TRTRI_UP
            for ( k = i-1 ; k > j ; k-- ){
#endif // INNERLOOP_TRTRI_UP
              
              dgemm_( &c_notranspose, &c_notranspose, &nb[i], &nb[j], &nb[k], &pone, A(i,k), &nb[i], A(k,j), &nb[k], &pone, A(i,j), &nb[i]);
              
            }

            //dtrmm_( &c_right, &c_lower, &c_notranspose, &c_nonunit, &nb[i], &nb[j], &mone, A(j,j), &nb[j], A(i,j), &nb[i] );
            dtrsm_( &c_right, &c_lower, &c_notranspose, &c_nonunit, &nb[i], &nb[j], &mone, A(j,j), &nb[j], A(i,j), &nb[i] );

          }
          
          dtrtri_( &c_lower, &c_nonunit, &nb[j], A(j,j), &nb[j],  &info );

        }

        rtn=gettimeofday(&tp, NULL);
        t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
        elapsed=t2-t1;
        /* printf("%fsec\n",elapsed); */

        /************************************************************************************************/
        /* Inplace xLAUUM (LAPACK LAUUM)                                                                */
        /************************************************************************************************/

        rtn=gettimeofday(&tp, NULL);
        t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;

        for ( i = 0 ; i < p ; i++ ){
          
          for ( j = 0 ; j < i ; j++ ){
            dtrmm_( &c_left, &c_lower, &c_transpose, &c_nonunit, &nb[i], &nb[j], &pone, A(i,i), &nb[i], A(i,j), &nb[i] );
          }
                
          dlauum_( &c_lower, &nb[i], A(i,i), &nb[i], &info );
                
          for ( j = 0 ; j < i ; j++ ) {
            
#ifdef INNERLOOP_LAUUM_UP
            for ( k = i+1 ; k < p ; k++ ) {
#else // INNERLOOP_LAUUM_UP
            for ( k = p-1 ; k > i ; k-- ) {
#endif // INNERLOOP_LAUUM_UP
              
              dgemm_( &c_transpose, &c_notranspose, &nb[i], &nb[j], &nb[k], &pone, A(k,i), &nb[k], A(k,j), &nb[k], &pone, A(i,j), &nb[i]);
            }
          }
                
          for ( k = i+1 ; k < p ; k++ ){
            dsyrk_( &c_lower, &c_transpose, &nb[i], &nb[k], &pone, A(k,i), &nb[k], &pone, A(i,i), &nb[i]);
          }

        }

        rtn=gettimeofday(&tp, NULL);
        t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
        elapsed=t2-t1;
        /* printf("%fsec\n",elapsed); */
}
#endif // SEQCHECK
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
#if defined OUTOFPLACE
void tile_potrf_parallel(Quark *quark, double *A, double *Ac, double *Ai, int n, int nb_default)
#else // INPLACE
void tile_potrf_parallel(Quark *quark, double *A, int n, int nb_default)
#endif // INPLACE or OUTOFPLACE
{
        int i, j, k, p;
        int *nb;
        int info;
        double mone=-1.0e+00;
        double pone=+1.0e+00;

        //struct timeval tp;
        //double t1,t2,elapsed;
        //int rtn;

        p = n / nb_default; if( p*nb_default != n ) p++;
        nb = (int *)malloc(p*sizeof(int)) ;
        for (i=0;i<p-1;i++) nb[i] = nb_default; nb[p-1] = n - nb_default*(p-1);

        int NBNBSIZE = nb_default*nb_default;
        int MT = p;
        int NB = nb_default;
        char *tasklabel = calloc(200,sizeof(char));
        char *taskcolor = calloc(200,sizeof(char));

    /************************************************************************************************/
    /* Cholesky factorization (LAPACK POTRF)                                                        */
    /************************************************************************************************/

        //rtn=gettimeofday(&tp, NULL);
        //t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;


        for ( j = 0 ; j < p ; j++ ){
          
          for ( k = 0 ; k < j ; k++ ) { 
            core_event_start(0);
#ifdef DEBUG        
            printf("dsyrk %d %d\n", j, k);
#endif // DEBUG
            snprintf( tasklabel, 200, "dsyrk %d %d", j, k);
            QUARK_Insert_Task(quark, SCHED_dsyrk, 0, 
                               sizeof(char),         "L",        VALUE,
                               sizeof(char),         "N",        VALUE,
                               sizeof(int),          &nb[j],     VALUE,
                               sizeof(int),          &nb[k],     VALUE,
                               sizeof(double),       &mone,      VALUE, 
                               sizeof(double)*NB*NB, A(j,k),     INPUT,
                               sizeof(int),          &nb[j],     VALUE,
                               sizeof(double),       &pone,      VALUE, 
                               sizeof(double)*NB*NB, A(j,j),     INOUT | LOCALITY,
                               sizeof(int),          &nb[j],     VALUE,
                               strlen(tasklabel)+1,  tasklabel,  VALUE | TASKLABEL,
                               0);
            core_event_end(0);
            core_log_event(0xA0A0A0, 0, 0);
          }
          
          core_event_start(0);
#ifdef DEBUG        
          printf("dpotrf %d\n", j);
#endif // DEBUG
          snprintf( tasklabel, 200, "dpotrf %d", j);
          QUARK_Insert_Task(quark, SCHED_dpotrf, 0, 
                             sizeof(char),         "L",        VALUE,
                             sizeof(int),          &nb[j],     VALUE,
                             sizeof(double)*NB*NB, A(j,j),     INOUT | LOCALITY,
                             sizeof(int),          &nb[j],     VALUE,
                             sizeof(int),          &info,      VALUE,
                             strlen(tasklabel)+1,  tasklabel,  VALUE | TASKLABEL,
                             0);
          core_event_end(0);
          core_log_event(0xA0A0A0, 0, 0);



          for ( i = j+1 ; i < p ; i++ ) {

#ifdef INNERLOOP_POTRF_UP
            for ( k = 0 ; k < j ; k++ ) {
#else // INNERLOOP_POTRF_UP
            for ( k = j-1 ; k  > -1 ; k-- ) {
#endif // INNERLOOP_POTRF_UP
              
              core_event_start(0);
#ifdef DEBUG        
              printf("dgemm %d %d %d\n", j, i, k);
#endif // DEBUG
              snprintf( tasklabel, 200, "dgemm %d %d %d", j,i,k);
              QUARK_Insert_Task(quark, SCHED_dgemm, 0, 
                                 sizeof(char),         "N",        VALUE,
                                 sizeof(char),         "T",        VALUE,
                                 sizeof(int),          &nb[i],     VALUE,
                                 sizeof(int),          &nb[j],     VALUE,
                                 sizeof(int),          &nb[k],     VALUE,
                                 sizeof(double),       &mone,      VALUE, 
                                 sizeof(double)*NB*NB, A(i,k),     INPUT,
                                 sizeof(int),          &nb[i],     VALUE,
                                 sizeof(double)*NB*NB, A(j,k),     INPUT,
                                 sizeof(int),          &nb[j],     VALUE,
                                 sizeof(double),       &pone,      VALUE, 
                                 sizeof(double)*NB*NB, A(i,j),     INOUT | LOCALITY | NOACCUMULATOR,
                                 sizeof(int),          &nb[i],     VALUE,
                                 strlen(tasklabel)+1,  tasklabel,  VALUE | TASKLABEL,
                                 0);
              core_event_end(0);
              core_log_event(0xA0A0A0, 0, 0);
            }
          }
          
          for ( i = j+1 ; i < p ; i++ ) {
#ifdef DEBUG        
            printf("dtrsm %d %d\n", j, i);
#endif // DEBUG
            // dtrsm_("R", "L", "T", "N", &nb[i], &nb[j], &pone, A(j,j), &nb[j], A(i,j), &nb[i]);
            snprintf( tasklabel, 200, "dtrsm %d %d", j,i);
            core_event_start(0);
            QUARK_Insert_Task(quark, SCHED_dtrsm, 0, 
                               sizeof(char),         "R",        VALUE,
                               sizeof(char),         "L",        VALUE,
                               sizeof(char),         "T",        VALUE,
                               sizeof(char),         "N",        VALUE,
                               sizeof(int),          &nb[i],     VALUE,
                               sizeof(int),          &nb[j],     VALUE,
                               sizeof(double),       &pone,      VALUE,
                               sizeof(double)*NB*NB, A(j,j),     INPUT,
                               sizeof(int),          &nb[j],     VALUE,
                               sizeof(double)*NB*NB, A(i,j),     INOUT | LOCALITY,
                               sizeof(int),          &nb[i],     VALUE,
                               strlen(tasklabel)+1,  tasklabel,  VALUE | TASKLABEL,
                               0);
            core_event_end(0);
            core_log_event(0xA0A0A0, 0, 0);

          }
          
        }
}
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
#if defined OUTOFPLACE
void tile_trtri_parallel(Quark *quark, double *A, double *Ac, double *Ai, int n, int nb_default)
#else // INPLACE
void tile_trtri_parallel(Quark *quark, double *A, int n, int nb_default)
#endif // INPLACE or OUTOFPLACE
{
        int i, j, k, p;
        int *nb;
        int info;
        double mone=-1.0e+00;
        double pone=+1.0e+00;
#if defined OUTOFPLACE
        int ione=1;
#endif

        p = n / nb_default; if( p*nb_default != n ) p++;
        nb = (int *)malloc(p*sizeof(int)) ;
        for (i=0;i<p-1;i++) nb[i] = nb_default; nb[p-1] = n - nb_default*(p-1);

        int NBNBSIZE = nb_default*nb_default;
        int MT = p;
        int NB = nb_default;
        char *tasklabel = calloc(200,sizeof(char));
        char *taskcolor = calloc(200,sizeof(char));
        unsigned long long int taskid = 0;
        int retval = 0;

        /************************************************************************************************/
        /* Inplace Triangular Inversion (LAPACK TRTRI)                                                  */
        /************************************************************************************************/

#if defined OUTOFPLACE
        // JL2ALL: think about the smart to do the looping j=0:p-1, i=0:p-1 ... other order?
        // JL2ALL: only copy the symmetric part of diagonal tiles
        for ( j = 0 ; j < p ; j++ ){
                for ( i = j ; i < p ; i++ ){
                        core_event_start(0);
#ifdef DEBUG        
                        printf("dcopy  %d %d\n", i, j);
#endif // DEBUG
                        snprintf( tasklabel, 200, "dcopy %d %d", i, j);
                        k = nb[j]*nb[i];
                        QUARK_Insert_Task(quark, SCHED_dcopy, 0, 
                                           sizeof(int),          &k,         VALUE,
                                           sizeof(double)*NB*NB, A(i,j),     INPUT,
                                           sizeof(int),          &pone,      VALUE,
                                           sizeof(double)*NB*NB, Ac(i,j),    OUTPUT,
                                           sizeof(int),          &pone,      VALUE,
                                           strlen(tasklabel)+1,  tasklabel,  VALUE | TASKLABEL,
                                           0);
                        core_event_end(0);
                        core_log_event(0x000000, 0, 0);
                }
        }
#endif // OUTOFPLACE

        for ( j = p-1 ; j > -1 ; j-- ){
          
/*
          // dtrtri_( "L", "N", &nb[j], A(j,j), &nb[j],  &info );
           core_event_start(0);
#ifdef DEBUG        
           printf("dtrtri %d\n", j);
#endif // DEBUG
           QUARK_Insert_Task(quark, SCHED_dtrtri, 0, 
           sizeof(char),         "L",        VALUE,
           sizeof(char),         "N",        VALUE,
           sizeof(int),          &nb[j],     VALUE,
           sizeof(double)*NB*NB, Ac(j,j),     INOUT | NOCLOCALITY,
           sizeof(int),          &nb[j],     VALUE,
           sizeof(int),          &info,      VALUE,
           0);
              core_event_end(0);
              core_log_event(0xA0A0A0, 0, 0);
*/

          for ( i = p-1 ; i > j ; i-- ){
            
            //dtrmm_( "L", "L", "N", "N", &nb[i], &nb[j], &pone, A(i,i), &nb[i], A(i,j), &nb[i] );
            core_event_start(0);
#ifdef DEBUG        
            printf("dtrmm %d %d\n", i, j);
#endif // DEBUG
            snprintf( tasklabel, 200, "dtrmm %d %d", i, j);
            QUARK_Insert_Task(quark, SCHED_dtrmm, 0, 
                               sizeof(char),         "L",        VALUE,
                               sizeof(char),         "L",        VALUE,
                               sizeof(char),         "N",        VALUE,
                               sizeof(char),         "N",        VALUE,
                               sizeof(int),          &nb[i],     VALUE,
                               sizeof(int),          &nb[j],     VALUE,
                               sizeof(double),       &pone,      VALUE,
                               sizeof(double)*NB*NB, Ac(i,i),     INPUT,
                               sizeof(int),          &nb[i],     VALUE,
                               sizeof(double)*NB*NB, Ac(i,j),     INOUT | LOCALITY,
                               sizeof(int),          &nb[i],     VALUE,
                               strlen(tasklabel)+1,  tasklabel,  VALUE | TASKLABEL,
                               0);
            core_event_end(0);
            core_log_event(0xA0A0A0, 0, 0);
            
#ifdef INNERLOOP_TRTRI_UP
            for ( k = j+1 ; k < i ; k++ ){
#else // INNERLOOP_TRTRI_UP
            for ( k = i-1 ; k > j ; k-- ){
#endif // INNERLOOP_TRTRI_UP
              
              
              //dgemm_( "N", "N", &nb[i], &nb[j], &nb[k], &pone, A(i,k), &nb[i], A(k,j), &nb[k], &pone, A(i,j), &nb[i]);
              core_event_start(0);
#ifdef DEBUG        
              printf("dgemm %d %d %d\n", i, j, k);
#endif // DEBUG
              snprintf( tasklabel, 200, "dgemm %d %d %d", i, j, k);
              snprintf( taskcolor, 200, "red");
              QUARK_Insert_Task(quark, SCHED_dgemm, 0, 
                                 sizeof(char),         "N",        VALUE,
                                 sizeof(char),         "N",        VALUE,
                                 sizeof(int),          &nb[i],     VALUE,
                                 sizeof(int),          &nb[j],     VALUE,
                                 sizeof(int),          &nb[k],     VALUE,
                                 sizeof(double),       &pone,      VALUE, 
                                 sizeof(double)*NB*NB, Ac(i,k),    INPUT,
                                 sizeof(int),          &nb[i],     VALUE,
                                 sizeof(double)*NB*NB, A(k,j),     INPUT,
                                 sizeof(int),          &nb[k],     VALUE,
                                 sizeof(double),       &pone,      VALUE, 
                                 sizeof(double)*NB*NB, Ac(i,j),   INOUT | LOCALITY | NOACCUMULATOR,
                                 sizeof(int),          &nb[i],     VALUE,
                                 strlen(tasklabel)+1,  tasklabel,  VALUE | TASKLABEL,
                                 strlen(tasklabel)+1,  taskcolor,  VALUE | TASKCOLOR,
                                 0);
              core_event_end(0);
              core_log_event(0xA0A0A0, 0, 0);
              
            }

            //dtrmm_( "R", "L", "N", "N", &nb[i], &nb[j], &mone, A(j,j), &nb[j], A(i,j), &nb[i] );
            core_event_start(0);
#ifdef DEBUG        
            printf("dtrsm %d %d\n", i, j);
#endif // DEBUG
            snprintf( tasklabel, 200, "dtrsm %d %d", i, j);
            QUARK_Insert_Task(quark, SCHED_dtrsm, 0,
                               sizeof(char),         "R",        VALUE,
                               sizeof(char),         "L",        VALUE,
                               sizeof(char),         "N",        VALUE,
                               sizeof(char),         "N",        VALUE,
                               sizeof(int),          &nb[i],     VALUE,
                               sizeof(int),          &nb[j],     VALUE,
                               sizeof(double),       &mone,      VALUE,
                               sizeof(double)*NB*NB, A(j,j),     INPUT,
                               sizeof(int),          &nb[j],     VALUE,
                               sizeof(double)*NB*NB, Ac(i,j),    INOUT | LOCALITY,
                               sizeof(int),          &nb[i],     VALUE,
                               strlen(tasklabel)+1,  tasklabel,  VALUE | TASKLABEL,
                               0);
            core_event_end(0);
            core_log_event(0xA0A0A0, 0, 0);

          }

          // dtrtri_( "L", "N", &nb[j], A(j,j), &nb[j],  &info );
           core_event_start(0);
#ifdef DEBUG        
           printf("dtrtri %d\n", j);
#endif // DEBUG
           snprintf( tasklabel, 200, "dtrtri %d", j);
           QUARK_Insert_Task(quark, SCHED_dtrtri, 0, 
                              sizeof(char),         "L",        VALUE,
                              sizeof(char),         "N",        VALUE,
                              sizeof(int),          &nb[j],     VALUE,
                              sizeof(double)*NB*NB, Ac(j,j),     INOUT | LOCALITY,
                              sizeof(int),          &nb[j],     VALUE,
                              sizeof(int),          &info,      VALUE,
                              strlen(tasklabel)+1,  tasklabel,  VALUE | TASKLABEL,
                              0);
              core_event_end(0);
              core_log_event(0xA0A0A0, 0, 0);

        }
}
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
#if defined OUTOFPLACE
void tile_lauum_parallel(Quark *quark, double *A, double *Ac, double *Ai, int n, int nb_default)
#else // INPLACE
void tile_lauum_parallel(Quark *quark, double *A, int n, int nb_default)
#endif // INPLACE or OUTOFPLACE
{
        int i, j, k, p;
        int *nb;
        int info;
        double pone=+1.0e+00;
#if defined OUTOFPLACE
        int ione=1;
#endif

        p = n / nb_default; if( p*nb_default != n ) p++;
        nb = (int *)malloc(p*sizeof(int)) ;
        for (i=0;i<p-1;i++) nb[i] = nb_default; nb[p-1] = n - nb_default*(p-1);

        int NBNBSIZE = nb_default*nb_default;
        int MT = p;
        int NB = nb_default;
        char *tasklabel = calloc(200,sizeof(char));

        /************************************************************************************************/
        /* Inplace xLAUUM (LAPACK LAUUM)                                                                */
        /************************************************************************************************/

#if defined OUTOFPLACE
        for ( j = 0 ; j < p ; j++ ){
                for ( i = j ; i < p ; i++ ){
                        core_event_start(0);
#ifdef DEBUG        
                        printf("dcopy  %d %d\n", i, j);
#endif // DEBUG
                        k = nb[j]*nb[i];
                        snprintf( tasklabel, 200, "dcopy %d %d", i, j);
                        QUARK_Insert_Task(quark, SCHED_dcopy, 0,
                                           sizeof(int),          &k,         VALUE,
                                           sizeof(double)*NB*NB, Ac(i,j),    INPUT,
                                           sizeof(int),          &pone,      VALUE,
                                           sizeof(double)*NB*NB, Ai(i,j),    OUTPUT,
                                           sizeof(int),          &pone,      VALUE,
                                           strlen(tasklabel)+1,  tasklabel,  VALUE | TASKLABEL,
                                           0);
                        core_event_end(0);
                        core_log_event(0x000000, 0, 0);
                }
        }
#endif // OUTOFPLACE

        
        for ( i = 0 ; i < p ; i++ ){
          
          for ( j = 0 ; j < i ; j++ ){
            //dtrmm_( "L", "L", "T", "N", &nb[i], &nb[j], &pone, A(i,i), &nb[i], A(i,j), &nb[i] );
            core_event_start(0);
#ifdef DEBUG        
            printf("dtrmm %d %d\n", i, j);
#endif // DEBUG
            snprintf( tasklabel, 200, "dtrmm %d %d", i, j);
            QUARK_Insert_Task(quark, SCHED_dtrmm, 0, 
                               sizeof(char),         "L",        VALUE,
                               sizeof(char),         "L",        VALUE,
                               sizeof(char),         "T",        VALUE,
                               sizeof(char),         "N",        VALUE,
                               sizeof(int),          &nb[i],     VALUE,
                               sizeof(int),          &nb[j],     VALUE,
                               sizeof(double),       &pone,      VALUE,
                               sizeof(double)*NB*NB, Ac(i,i),     INPUT,
                               sizeof(int),          &nb[i],     VALUE,
                               sizeof(double)*NB*NB, Ai(i,j),     INOUT | LOCALITY,
                               sizeof(int),          &nb[i],     VALUE,
                               strlen(tasklabel)+1,  tasklabel,  VALUE | TASKLABEL,
                               0);
            core_event_end(0);
            core_log_event(0xA0A0A0, 0, 0);

          }
                
          //dlauum_( "L", &nb[i], A(i,i), &nb[i], &info );
          core_event_start(0);
#ifdef DEBUG        
          printf("dlauum %d\n", i);
#endif // DEBUG
          snprintf( tasklabel, 200, "dlauum %d", i);
          QUARK_Insert_Task(quark, SCHED_dlauum, 0, 
                             sizeof(char),         "L",        VALUE,
                             sizeof(int),          &nb[i],     VALUE,
                             sizeof(double)*NB*NB, Ai(i,i),     INOUT | LOCALITY,
                             sizeof(int),          &nb[i],     VALUE,
                             sizeof(int),          &info,      VALUE,
                             strlen(tasklabel)+1,  tasklabel,  VALUE | TASKLABEL,
                             0);
          
          for ( j = 0 ; j < i ; j++ ) {
            
#ifdef INNERLOOP_LAUUM_UP
            for ( k = i+1 ; k < p ; k++ ) {
#else // INNERLOOP_LAUUM_UP
            for ( k = p-1 ; k > i ; k-- ) {
#endif // INNERLOOP_LAUUM_UP
              
              //dgemm_( "T", "N", &nb[i], &nb[j], &nb[k], &pone, A(k,i), &nb[k], A(k,j), &nb[k], &pone, A(i,j), &nb[i]);
              core_event_start(0);
#ifdef DEBUG        
              printf("dgemm %d %d %d\n", i, j, k);
#endif // DEBUG
              snprintf( tasklabel, 200, "dgemm %d %d %d", i, j, k);
              QUARK_Insert_Task(quark, SCHED_dgemm, 0, 
                                 sizeof(char),         "T",        VALUE,
                                 sizeof(char),         "N",        VALUE,
                                 sizeof(int),          &nb[i],     VALUE,
                                 sizeof(int),          &nb[j],     VALUE,
                                 sizeof(int),          &nb[k],     VALUE,
                                 sizeof(double),       &pone,      VALUE, 
                                 sizeof(double)*NB*NB, Ac(k,i),     INPUT,
                                 sizeof(int),          &nb[k],     VALUE,
                                 sizeof(double)*NB*NB, Ac(k,j),     INPUT,
                                 sizeof(int),          &nb[k],     VALUE,
                                 sizeof(double),       &pone,      VALUE, 
                                 sizeof(double)*NB*NB, Ai(i,j),     INOUT | LOCALITY | NOACCUMULATOR,
                                 sizeof(int),          &nb[i],     VALUE,
                                 strlen(tasklabel)+1,  tasklabel,  VALUE | TASKLABEL,
                                 0);
              core_event_end(0);
              core_log_event(0xA0A0A0, 0, 0);

            }
          }
                
          for ( k = i+1 ; k < p ; k++ ){
            //dsyrk_( "L", "T", &nb[i], &nb[k], &pone, A(k,i), &nb[k], &pone, A(i,i), &nb[i]);
            core_event_start(0);
#ifdef DEBUG        
            printf("dsyrk %d %d\n", i, k);
#endif // DEBUG
            snprintf( tasklabel, 200, "dsyrk %d %d", i, k);
            QUARK_Insert_Task(quark, SCHED_dsyrk, 0, 
                               sizeof(char),         "L",        VALUE,
                               sizeof(char),         "T",        VALUE,
                               sizeof(int),          &nb[i],     VALUE,
                               sizeof(int),          &nb[k],     VALUE,
                               sizeof(double),       &pone,      VALUE, 
                               sizeof(double)*NB*NB, Ac(k,i),     INPUT,
                               sizeof(int),          &nb[k],     VALUE,
                               sizeof(double),       &pone,      VALUE, 
                               sizeof(double)*NB*NB, Ai(i,i),     INOUT | LOCALITY,
                               sizeof(int),          &nb[i],     VALUE,
                               strlen(tasklabel)+1,  tasklabel,  VALUE | TASKLABEL,
                               0);
            core_event_end(0);
            core_log_event(0xA0A0A0, 0, 0);
          }
        }

}
///////////////////////////////////////////////////////////////////////////////////////////////////
#undef A
#undef Ac
#undef Ai
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <pthread.h>

#define HUGE_PAGE_SIZE 2048*1024
int main (int argc, char **argv)
{
    assert(argc == 4);
    int N = atoi(argv[1]);
    int NB = atoi(argv[2]);
    int BB = N / NB; if( BB*NB != N ) BB++;
    int NxN = N*N;
    int THREADS = atoi(argv[3]); assert(THREADS <= MAX_THREADS);

    double *Ablk  = (double*)malloc(NxN*sizeof(double));
    double *Ablk2 = (double*)malloc(NxN*sizeof(double));
    double *A     = (double*)malloc(NxN*sizeof(double));
    double *A2    = (double*)malloc(NxN*sizeof(double));

    double elapsed;
    double GFLOPS;

    Quark *quark;

#if defined PIPELINE
#else // PIPELINE
    double elapsed_potrf;
    double elapsed_trtri;
    double elapsed_lauum;
#endif // PIPELINE

#if defined OUTOFPLACE
    double *Ac = (double*)malloc(NxN*sizeof(double));
    double *Ai = (double*)malloc(NxN*sizeof(double));
#endif

#if defined NUMCHECK
    char c_lower='L';
#if defined SEQCHECK
    double GFLOPSseq;
    double residseq;
#endif // SEQCHECK
    double *Asv   = (double*)malloc(NxN*sizeof(double));
    double rcond, resid;
    double *work, *rwork;
    work = (double*)malloc(N*N*sizeof(double)); // will not work if N=0 !!! needs to be one.
    rwork = (double*)malloc(N*sizeof(double));
#endif // NUMCHECK

    // Initialize A and A2
    int i;
    for (i = 0; i < NxN; i++)
        A[i] = A2[i] = 0.5 - (double)rand() / RAND_MAX; 
    for (i = 0; i < N; i++)
        A[i*N+i] = A2[i*N+i] = N;

#if defined NUMCHECK
    // Make copies to use for checks
    for (i = 0; i < NxN; i++)
        Asv[i] = A[i]; 
#endif // NUMCHECK

    // Move from F77 to BDL
    {int X, Y, x, y;
    for (X = 0; X < BB; X++)
      for (Y = 0; Y < BB; Y++)
        for (x = 0; x < NB; x++)
          for (y = 0; y < NB; y++)
            Ablk2[Y*NB*NB + y + X*NB*NB*BB + x*NB] = A2[Y*NB + y + X*NB*N + x*N];}

    {int X, Y, x, y;
    for (X = 0; X < BB; X++)
      for (Y = 0; Y < BB; Y++)
        for (x = 0; x < NB; x++)
          for (y = 0; y < NB; y++)
            Ablk[Y*NB*NB + y + X*NB*NB*BB + x*NB] = A[Y*NB + y + X*NB*N + x*N];}

#if defined NUMCHECK
#if defined SEQCHECK
    // Sequential tile Cholesky in BDL
    //printf("Enter Sequential Cholesky\n");
    elapsed =- get_current_time();
    tile_ch_sequential(Ablk, N, NB);
    //printf("Leave Sequential Cholesky\n");
    elapsed += get_current_time();
    GFLOPSseq = ((double)N)*((double)N)*((double)N) / elapsed / 1.0e+9;
#endif // SEQCHECK
#endif // NUMCHECK


#if defined OUTOFPLACE
#if defined PIPELINE

    elapsed =- get_current_time();
    quark = QUARK_New(THREADS); // quark_new(THREADS);  
    tile_potrf_parallel(quark, Ablk2, Ac, Ai, N, NB);
    tile_trtri_parallel(quark, Ablk2, Ac, Ai, N, NB);
    tile_lauum_parallel(quark, Ablk2, Ac, Ai, N, NB);
    QUARK_Delete(quark); // quark_delete(quark);
    Ablk2 = Ai;
    elapsed += get_current_time();

#else // PIPELINE

    elapsed =- get_current_time();
    elapsed_potrf = elapsed;
    quark = QUARK_New(THREADS);
    tile_potrf_parallel(quark, Ablk2, Ac, Ai, N, NB);
    QUARK_Delete(quark);
    elapsed_potrf += get_current_time();
    elapsed_trtri =- get_current_time();
    quark = QUARK_New(THREADS);
    tile_trtri_parallel(quark, Ablk2, Ac, Ai, N, NB);
    QUARK_Delete(quark);
    elapsed_trtri += get_current_time();
    elapsed_lauum =- get_current_time();
    quark = QUARK_New(THREADS);
    tile_lauum_parallel(quark, Ablk2, Ac, Ai, N, NB);
    QUARK_Delete(quark);
    elapsed_lauum += get_current_time();
    Ablk2 = Ai;
    elapsed += get_current_time();

#endif // PIPELINE
#else // INPLACE
#if defined PIPELINE

    elapsed =- get_current_time();
    quark = QUARK_New(THREADS);
    tile_potrf_parallel(quark, Ablk2, N, NB);
    tile_trtri_parallel(quark, Ablk2, N, NB);
    tile_lauum_parallel(quark, Ablk2, N, NB);
    QUARK_Delete(quark);
    elapsed += get_current_time();

#else // PIPELINE

    elapsed =- get_current_time();
    elapsed_potrf = elapsed;
    quark = QUARK_New(THREADS);
    tile_potrf_parallel(quark, Ablk2, N, NB);
    QUARK_Delete(quark);
    elapsed_potrf += get_current_time();
    elapsed_trtri =- get_current_time();
    quark = QUARK_New(THREADS);
    tile_trtri_parallel(quark, Ablk2, N, NB);
    QUARK_Delete(quark);
    elapsed_trtri += get_current_time();
    elapsed_lauum =- get_current_time();
    quark = QUARK_New(THREADS);
    tile_lauum_parallel(quark, Ablk2, N, NB);
    QUARK_Delete(quark);
    elapsed_lauum += get_current_time();
    elapsed += get_current_time();

#endif // PIPELINE
#endif // INPLACE or OUTOFPLACE

    GFLOPS = ((double)N)*((double)N)*((double)N) / elapsed / 1.0e+9;

    /* dump_trace(THREADS); */

#if defined SEQCHECK
    // Move from BDL to F77
    {int X, Y, x, y;
    for (X = 0; X < BB; X++)
      for (Y = 0; Y < BB; Y++)
        for (x = 0; x < NB; x++)
          for (y = 0; y < NB; y++)
            A[Y*NB + y + X*NB*N + x*N] = Ablk[Y*NB*NB + y + X*NB*NB*BB + x*NB];}
#endif // SEQCHECK

    {int X, Y, x, y;
    for (X = 0; X < BB; X++)
      for (Y = 0; Y < BB; Y++)
        for (x = 0; x < NB; x++)
          for (y = 0; y < NB; y++) {
            A2[Y*NB + y + X*NB*N + x*N] = Ablk2[Y*NB*NB + y + X*NB*NB*BB + x*NB];
          }
    }

#if defined PIPELINE

#if defined NUMCHECK

#if defined SEQCHECK
    dpot03_( &c_lower, &N, Asv, &N, A, &N, work, &N, rwork, &rcond, &residseq );
#if defined OUTOFPLACE
    printf("SEQ  CHINV OUTOFPLACE %6d %4d %3d %6.2lf %6.2e\n", N, NB, 1, GFLOPSseq, residseq);
#else // INPLACE
    printf("SEQ  CHINV INPLACE    %6d %4d %3d %6.2lf %6.2e\n", N, NB, 1, GFLOPSseq, residseq);
#endif // INPLACE or OUTOFPLACE
#endif // SEQCHECK

    dpot03_( &c_lower, &N, Asv, &N, A2, &N, work, &N, rwork, &rcond, &resid );

#if defined OUTOFPLACE
    printf("QUARK CHINV OUTOFPLACE %6d %4d %3d %6.2lf %6.2e\n", N, NB, THREADS, GFLOPS, resid);
#else // INPLACE
    printf("QUARK CHINV INPLACE    %6d %4d %3d %6.2lf %6.2e\n", N, NB, THREADS, GFLOPS, resid);
#endif // INPLACE or OUTOFPLACE

    free(work);
    free(rwork);
#else // (no NUMCHECK)

#if defined OUTOFPLACE
    printf("QUARK CHINV OUTOFPLACE %6d %4d %3d %8.2lf\n", N, NB, THREADS, GFLOPS);
#else // INPLACE
    printf("QUARK CHINV INPLACE    %6d %4d %3d %8.2lf\n", N, NB, THREADS, GFLOPS);
#endif // INPLACE or OUTOFPLACE

#endif // NUMCHECK

#else // PIPELINE

#if defined NUMCHECK

#if defined SEQCHECK
    dpot03_( &c_lower, &N, Asv, &N, A, &N, work, &N, rwork, &rcond, &residseq );
#if defined OUTOFPLACE
    printf("SEQ  CHINV OUTOFPLACE %6d %4d %3d %6.2lf %6.2e %10.5f %10.5f %10.5f %10.5f\n", N, NB, 1, GFLOPSseq, residseq, elapsed, elapsed_potrf, elapsed_trtri, elapsed_lauum);
#else // INPLACE
    printf("SEQ  CHINV INPLACE    %6d %4d %3d %6.2lf %6.2e %10.5f %10.5f %10.5f %10.5f\n", N, NB, 1, GFLOPSseq, residseq, elapsed, elapsed_potrf, elapsed_trtri, elapsed_lauum);
#endif // INPLACE or OUTOFPLACE
#endif // SEQCHECK

    dpot03_( &c_lower, &N, Asv, &N, A2, &N, work, &N, rwork, &rcond, &resid );

#if defined OUTOFPLACE
    printf("QUARK CHINV OUTOFPLACE %6d %4d %3d %6.2lf %6.2e %10.5f %10.5f %10.5f %10.5f\n", N, NB, THREADS, GFLOPS, resid, elapsed, elapsed_potrf, elapsed_trtri, elapsed_lauum);
#else // INPLACE
    printf("QUARK CHINV INPLACE    %6d %4d %3d %6.2lf %6.2e %10.5f %10.5f %10.5f %10.5f\n", N, NB, THREADS, GFLOPS, resid, elapsed, elapsed_potrf, elapsed_trtri, elapsed_lauum);
#endif // INPLACE or OUTOFPLACE

    free(work);
    free(rwork);
#else // (no NUMCHECK)

#if defined OUTOFPLACE
    printf("QUARK CHINV OUTOFPLACE %6d %4d %3d %8.2lf %10.5f %10.5f %10.5f %10.5f\n", N, NB, THREADS, GFLOPS, elapsed, elapsed_potrf, elapsed_trtri, elapsed_lauum);
#else // INPLACE
    printf("QUARK CHINV INPLACE    %6d %4d %3d %8.2lf %10.5f %10.5f %10.5f %10.5f\n", N, NB, THREADS, GFLOPS, elapsed, elapsed_potrf, elapsed_trtri, elapsed_lauum);
#endif // INPLACE or OUTOFPLACE

#endif // NUMCHECK

#endif // PIPELINE


}
