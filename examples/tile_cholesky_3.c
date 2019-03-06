//##################################################################################################

#include <stdio.h>    // Standard stuff
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include <fcntl.h>    // Huge TLB pages
#include <sys/mman.h>

#include <math.h>     // Math and MKL

typedef enum {CblasRowMajor=101, CblasColMajor=102} CBLAS_ORDER; 
typedef enum {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113} CBLAS_TRANSPOSE;
typedef enum {CblasUpper=121, CblasLower=122} CBLAS_UPLO; 
typedef enum {CblasNonUnit=131, CblasUnit=132} CBLAS_DIAG;
typedef enum {CblasLeft=141, CblasRight=142} CBLAS_SIDE;                                                           
void cblas_ssyrk(const  CBLAS_ORDER Order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const float alpha, const float *A, const int lda,
                 const float beta, float *C, const int ldc);
void cblas_sgemm(const  CBLAS_ORDER Order, const  CBLAS_TRANSPOSE TransA,
                 const  CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const float alpha, const float *A,
                 const int lda, const float *B, const int ldb,
                 const float beta, float *C, const int ldc);
void cblas_strsm(const  CBLAS_ORDER Order, const  CBLAS_SIDE Side,
                 const  CBLAS_UPLO Uplo, const  CBLAS_TRANSPOSE TransA,
                 const  CBLAS_DIAG Diag, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 float *B, const int ldb);

#include <sys/time.h> // Time

#include "quark.h"
#include "event_trace.h"


#include <pthread.h>

#define MAX_NB 1024
#define MAX_BB 128
#define HUGE_PAGE_SIZE 2048*1024

#define BSIZE 200
typedef float block_t[BSIZE][BSIZE];
typedef block_t **matrix_t;

#define LL
//#define RL

// #define CHECK
// #define LOG

#define SPOTRF_COLOR "red"
#define SSYRK_COLOR "#B060D0"
#define SGEMM_COLOR "#D0F040"
#define STRSM_COLOR "#00C0F0"

////////////////////////////////////////////////////////////////////////////////////////////////////
int    event_num        [MAX_THREADS];
double event_start_time [MAX_THREADS];
double event_end_time   [MAX_THREADS];
double event_log        [MAX_THREADS][MAX_EVENTS];
int log_events = 1;




struct {
    int cores_num;
    float *A;
    int NB;
    int BB;
} core_in_all;

double GFLOPS;

char Left = 'L', Transpose = 'T', Forward = 'F', Columnwise = 'C', Upper = 'U', Lower = 'L';

void dump_trace(int);
void diff_matrix(float *A, float *B, int NB, int BB);
//extern "C" void slarnv(int*, int*, int*, float*);
//extern "C" void spotrf(char*, int*, float*, int*, int*);
void slarnv_(int*, int*, int*, float*);
void spotrf_(char*, int*, float*, int*, int*);

void SCHED_tile_ssyrk(Quark *);
void SCHED_tile_spotrf(Quark *);
void SCHED_tile_sgemm(Quark *);
void SCHED_tile_strsm(Quark *);
////////////////////////////////////////////////////////////////////////////////////////////////////

double get_current_time(void)
{
    struct timeval  time_val;
//    struct timezone time_zone;

    gettimeofday(&time_val, NULL);

    return (double)(time_val.tv_sec) + (double)(time_val.tv_usec) / 1000000.0;
}

/*
double get_current_time(void) { 
  struct timeval tp; 
  gettimeofday( &tp, NULL );
  return tp.tv_sec + 1e-6 * tp.tv_usec; 
}
*/
////////////////////////////////////////////////////////////////////////////////////////////////////

int my_thread_id()
{
    static int num_threads = 0;
    static pthread_t threads[MAX_THREADS];
    pthread_t self;
    int i;

    self = pthread_self();
    for (i = 0; i < num_threads; i++)
      if (pthread_equal(self, threads[i]))
        return (i);

    threads[num_threads++] = self;
    return (num_threads);
}
////////////////////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char **argv)
{
    assert(argc == 3);
    int BB = atoi(argv[1]); assert(BB <= MAX_BB);
    int THREADS = atoi(argv[2]); assert(THREADS <= MAX_THREADS);
    int NB = BSIZE;
    int N = BB*BSIZE;
    int NxN = N*N;

    double start_time, end_time, elapsed_time;
    int m, n;
    int i, j;
    int step;
    int X, Y, x, y;
    matrix_t a;
    int INFO;
    int ONE = 1;
    int ISEED[4] = {0,0,0,1};

    float *A1   = (float*)malloc(NxN*sizeof(float));

    a = (matrix_t)malloc(BB * sizeof(block_t*));

    for (n = 0; n < BB; n++){
        a[n] = (block_t*)malloc(BB * sizeof(block_t));}

    slarnv_(&ONE, ISEED, &NxN, A1);
    // Make the matrix SYM
    for (i=0;i<N;i++) {
        for (j=i;j<N;j++) {
            *(A1+j*N+i)=*(A1+i*N+j);
            if (i==j)
               *(A1+j*N+i)=*(A1+i*N+j)+100*N;
        }
    }

    // Move from F77 to BDL
    for (X = 0; X < BB; X++)
      for (Y = 0; Y < BB; Y++)
        for (x = 0; x < BSIZE; x++)
          for (y = 0; y < BSIZE; y++)
              a[Y][X][x][y] = A1[Y*BSIZE + X*BSIZE*N + y + x*N];

    Quark *quark = QUARK_New(THREADS);
    assert(quark != NULL);
    start_time = get_current_time();
    char *tasklabel = calloc(200,sizeof(char));
    char *taskcolor = calloc(200,sizeof(char));

#if defined LL

    // Left-looking tile Cholesky
    for (step = 0; step < BB; step++)
    {
        for (n = 0; n < step; n++){
            //tile_ssyrk((float*)(a[step][n]), (float*)(a[step][step]), BSIZE);
            snprintf( tasklabel, 200, "SSYRK");
            snprintf( taskcolor, 200, SSYRK_COLOR);
            QUARK_Insert_Task(quark,SCHED_tile_ssyrk, 0, 
                               sizeof(int),          &NB,                     QUARK_VALUE,
                               sizeof(float)*NB*NB,  (float*)(a[step][n]),    QUARK_INPUT,
                               sizeof(float)*NB*NB,  (float*)(a[step][step]), QUARK_INOUT,
                               sizeof(int),          &step,                   QUARK_VALUE,
                               strlen(tasklabel)+1,  tasklabel,  QUARK_VALUE | TASKLABEL,
                               strlen(taskcolor)+1,  taskcolor,  QUARK_VALUE | TASKCOLOR,
                               0);
        }
        //tile_spotrf((float*)(a[step][step]), BSIZE);
        snprintf( tasklabel, 200, "SPOTRF");
        snprintf( taskcolor, 200, SPOTRF_COLOR);
        QUARK_Insert_Task(quark,SCHED_tile_spotrf, 0, 
                           sizeof(int),             &NB,                      QUARK_VALUE,
                           sizeof(float)*NB*NB,     (float*)(a[step][step]),   QUARK_INOUT | LOCALITY,
                           sizeof(int),              &INFO,                    QUARK_INOUT,
                           sizeof(int),             &step,                     QUARK_VALUE,
                           strlen(tasklabel)+1,  tasklabel,  QUARK_VALUE | TASKLABEL,
                           strlen(taskcolor)+1,  taskcolor,  QUARK_VALUE | TASKCOLOR,
                           0);

        for (m = step+1; m < BB; m++)
        {
            for (n = 0; n < step; n++){
                //tile_sgemm((float*)(a[step][n]), (float*)(a[m][n]), (float*)(a[m][step]), BSIZE);
                snprintf( tasklabel, 200, "SGEMM");
                snprintf( taskcolor, 200, SGEMM_COLOR);
                QUARK_Insert_Task(quark,SCHED_tile_sgemm, 0,
                                   sizeof(int),          &NB,                     QUARK_VALUE,
                                   sizeof(float)*NB*NB,  (float*)(a[step][n]),    QUARK_INPUT,
                                   sizeof(float)*NB*NB,  (float*)(a[m][n]),       QUARK_INPUT,
                                   sizeof(float)*NB*NB,  (float*)(a[m][step]),    QUARK_INOUT,
                                   sizeof(int),          &step,                   QUARK_VALUE,
                                   strlen(tasklabel)+1,  tasklabel,  QUARK_VALUE | TASKLABEL,
                                   strlen(taskcolor)+1,  taskcolor,  QUARK_VALUE | TASKCOLOR,
                                   0);
            }

            //tile_strsm((float*)(a[step][step]), (float*)(a[m][step]), BSIZE);
            snprintf( tasklabel, 200, "STRSM");
            snprintf( taskcolor, 200, STRSM_COLOR);
            QUARK_Insert_Task(quark,SCHED_tile_strsm, 0, 
                               sizeof(int),          &NB,                     QUARK_VALUE,
                               sizeof(float)*NB*NB,  (float*)(a[step][step]), QUARK_INPUT,
                               sizeof(float)*NB*NB,  (float*)(a[m][step]),    QUARK_INOUT,
                               sizeof(int),          &step,                   QUARK_VALUE,
                               strlen(tasklabel)+1,  tasklabel,  QUARK_VALUE | TASKLABEL,
                               strlen(taskcolor)+1,  taskcolor,  QUARK_VALUE | TASKCOLOR,
                               0);
        }
        
    }

#else 

    // Right-looking tile Cholesky
    for (step = 0; step < BB; step++)
    {
        //tile_spotrf((float*)(a[step][step]), BSIZE);
        snprintf( tasklabel, 200, "SPOTRF");
        snprintf( taskcolor, 200, SPOTRF_COLOR);
        QUARK_Insert_Task(quark,SCHED_tile_spotrf, HIGH_PRIORITY, 
                           sizeof(int),          &NB,                     QUARK_VALUE,
                           sizeof(float)*NB*NB,  (float*)(a[step][step]), QUARK_INOUT | LOCALITY,
                           sizeof(int),          &INFO,                      QUARK_INOUT,
                           sizeof(int),          &step,                   QUARK_VALUE,
                           strlen(tasklabel)+1,  tasklabel,  QUARK_VALUE | TASKLABEL,
                           strlen(taskcolor)+1,  taskcolor,  QUARK_VALUE | TASKCOLOR,
                           0);
        for (n = step+1; n < BB; n++){
            //tile_strsm((float*)(a[step][step]), (float*)(a[n][step]), BSIZE);
            snprintf( tasklabel, 200, "STRSM");
            snprintf( taskcolor, 200, STRSM_COLOR);
            QUARK_Insert_Task(quark,SCHED_tile_strsm, 0, 
                               sizeof(int),          &NB,                    QUARK_VALUE,
                               sizeof(float)*NB*NB, (float*)(a[step][step]), QUARK_INPUT,
                               sizeof(float)*NB*NB,  (float*)(a[n][step]),    QUARK_INOUT,
                               sizeof(int),          &step,                   QUARK_VALUE,
                               strlen(tasklabel)+1,  tasklabel,  QUARK_VALUE | TASKLABEL,
                               strlen(taskcolor)+1,  taskcolor,  QUARK_VALUE | TASKCOLOR,
                               0);
            //tile_ssyrk((float*)(a[n][step]), (float*)(a[n][n]), BSIZE);
            snprintf( tasklabel, 200, "SSYRK");
            snprintf( taskcolor, 200, SSYRK_COLOR);
            QUARK_Insert_Task(quark,SCHED_tile_ssyrk, 0, 
                               sizeof(int),          &NB,                     QUARK_VALUE,
                               sizeof(float)*NB*NB,  (float*)(a[n][step]),    QUARK_INPUT,
                               sizeof(float)*NB*NB,  (float*)(a[n][n]),       QUARK_INOUT,
                               sizeof(int),          &step,                   QUARK_VALUE,
                               strlen(tasklabel)+1,  tasklabel,  QUARK_VALUE | TASKLABEL,
                               strlen(taskcolor)+1,  taskcolor,  QUARK_VALUE | TASKCOLOR,
                               0);
        }
        for (m = step+2; m < BB; m++)
            for (n = step+1; n < m; n++) {
                //tile_sgemm((float*)(a[n][step]), (float*)(a[m][step]), (float*)(a[m][n]), BSIZE);
                snprintf( tasklabel, 200, "SGEMM");
                snprintf( taskcolor, 200, SGEMM_COLOR);
                QUARK_Insert_Task(quark,SCHED_tile_sgemm, 0, 
                                   sizeof(int),          &NB,                     QUARK_VALUE,
                                   sizeof(float)*NB*NB,  (float*)(a[n][step]),    QUARK_INPUT,
                                   sizeof(float)*NB*NB,  (float*)(a[m][step]),    QUARK_INPUT,
                                   sizeof(float)*NB*NB,  (float*)(a[m][n]),       QUARK_INOUT,
                                   sizeof(int),          &step,                   QUARK_VALUE,
                                   strlen(tasklabel)+1,  tasklabel,  QUARK_VALUE | TASKLABEL,
                                   strlen(taskcolor)+1,  taskcolor,  QUARK_VALUE | TASKCOLOR,
                                   0);
            }
    }

#endif

    QUARK_Delete(quark);


    end_time = get_current_time();
    elapsed_time = end_time - start_time;
    GFLOPS = 1.0*N*N*N/3.0 / elapsed_time / 1000000000;

#if defined CHECK
    {
    int INFO;
    spotrf_(&Lower, &N, A1, &N, &INFO);
    assert(INFO == 0);
    }
    float *A2   = (float*)malloc(NxN*sizeof(float));
    // Move from BDL to F77
    {int X, Y, x, y;
    for (X = 0; X < BB; X++)
      for (Y = 0; Y < BB; Y++)
        for (x = 0; x < BSIZE; x++)
          for (y = 0; y < BSIZE; y++)
            A2[Y*BSIZE + X*BSIZE*N + y + x*N] = a[Y][X][x][y];}
    diff_matrix(A1, A2, BSIZE, BB);
#endif

#if defined LOG
    dump_trace(THREADS);
#endif


#if defined LL
    printf("%s %s \t\t%d \t\t%d \t\t%d \t\t%d \t\t%.2lf \t\t%.2lf\n","azz","Left--Looking",NB,BB,N,THREADS,GFLOPS,GFLOPS / (2.393895 * 8 * 16) * 100.0);
#else
    printf("%s %s \t\t%d \t\t%d \t\t%d \t\t%d \t\t%.2lf \t\t%.2lf\n","azz","Right-Looking",NB,BB,N,THREADS,GFLOPS,GFLOPS / (2.393895 * 8 * 16) * 100.0);
#endif

}




////////////////////////////////////////////////////////////////////////////////////////////////////

void SCHED_tile_ssyrk(Quark *quark)
{
    int NB;
    float *A;
    float *T;
    int k;
    quark_unpack_args_4(quark, NB, A, T, k);

    //printf("SYRK\n");

#if defined LOG
    event_start(QUARK_Thread_Rank(quark));
#endif

    cblas_ssyrk(
        CblasColMajor,
        CblasLower, CblasNoTrans,
        NB, NB,
       -1.0, A, NB,
        1.0, T, NB);

#if defined LOG
    event_end(QUARK_Thread_Rank(quark));
    log_event(SSYRK_COLOR, QUARK_Thread_Rank(quark), k);
#endif
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SCHED_tile_spotrf(Quark *quark)
{
    int NB;
    float *A;
    int INFO;
    int k;
    quark_unpack_args_4(quark, NB, A, INFO, k);

    //printf("POTRF\n");

#if defined LOG
    event_start(QUARK_Thread_Rank(quark));
#endif

    spotrf_("L", &NB, A, &NB, &INFO);
    assert(INFO == 0);

#if defined LOG
    event_end(QUARK_Thread_Rank(quark));
    log_event(SPOTRF_COLOR, QUARK_Thread_Rank(quark), k);
#endif
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void SCHED_tile_sgemm(Quark *quark)
{
    int NB;
    float *A;
    float *B;
    float *C;
    int k;
    quark_unpack_args_5(quark, NB, A, B, C, k);

    //printf("GEMM\n");

#if defined LOG
    event_start(QUARK_Thread_Rank(quark));
#endif

    cblas_sgemm(
        CblasColMajor,
        CblasNoTrans, CblasTrans,
        NB, NB, NB,
       -1.0, B, NB,
             A, NB,
        1.0, C, NB);

#if defined LOG
    event_end(QUARK_Thread_Rank(quark));
    log_event(SGEMM_COLOR, QUARK_Thread_Rank(quark), k);
#endif
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SCHED_tile_strsm(Quark *quark)
{
    int NB;
    float *T;
    float *B;
    int k;
    quark_unpack_args_4(quark, NB, T, B, k);

    //printf("TRSM\n");

#if defined LOG
    event_start(QUARK_Thread_Rank(quark));
#endif

    cblas_strsm(
        CblasColMajor,
        CblasRight, CblasLower,
        CblasTrans, CblasNonUnit,
        NB, NB,
        1.0, T, NB,
             B, NB);

#if defined LOG
    event_end(QUARK_Thread_Rank(quark));
    log_event(STRSM_COLOR, QUARK_Thread_Rank(quark), k);
#endif
}

//##################################################################################################

void dump_trace(int cores_num)
{
    char trace_file_name[32];
    FILE *trace_file;
    int event;
    int core;
//  double scale = 500000.0;
//  double scale = 300.0;
    double scale = 150000.0;

    sprintf(trace_file_name, "trace_%d.svg", (int)(time(NULL)));
    trace_file = fopen(trace_file_name, "w");
    assert(trace_file != NULL);

    fprintf(trace_file,
        "<svg width=\"200mm\" height=\"40mm\" viewBox=\"0 0 20000 4000\">\n"
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


////////////////////////////////////////////////////////////////////////////////////////////////////

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

void diff_matrix(float *A, float *B, int NB, int BB)
{
    int X, Y, x, y;
    int N = NB*BB;

    printf("\n");
    for (Y = 0; Y < BB; Y++) {
      for (y = 0; y < NB; y++) {
        for (X = 0; X < BB; X++) {
          for (x = 0; x < NB; x++) {

            float a, b, c, d, e;
            a = fabs(A[(Y*NB+y) + (X*NB+x)*N]);
            b = fabs(B[(Y*NB+y) + (X*NB+x)*N]);
            c = max(a, b);
            d = min(a, b);
            e = (c - d) / d;

            //printf("%c", e < 0.001 ? '.' : '#');
                if(e>0.001) {
                        printf("%s%20.17lf %20.17f %20.17f %s %5d%5d%5d%5d\n"," ERROR for value  ",a,b,e," at Y y X x  ", Y, y, X, x);
                //      printf("%s%f%f%s%d%d%d%d\n"," ERROR for value  ",a,b, " at Y y X x ", Y,y,X,x);
                }
            //if (x == 3) x = NB-5;
            //if (x == 7) x = NB-1;
          }
          //printf("  |");
        }
       // printf("\n");
        //if (y == 3) y = NB-5;
        //if (y == 7) y = NB-1;
      }
      //if (Y < BB-1) for (i = 0; i < BB*12; i++) printf("=");  
      //printf("\n");
    }
    printf("\n");
}

////////////////////////////////////////////////////////////////////////////////////////////////////


