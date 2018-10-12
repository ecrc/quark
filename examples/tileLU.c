///////////////////////////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <string.h>
/* char *strdup(const char *s1); */

//#include "Quark.hh"
#include "quark.h"

#include "event_trace.h"
int    event_num        [MAX_THREADS];
double event_start_time [MAX_THREADS];
double event_end_time   [MAX_THREADS];
double event_log        [MAX_THREADS][MAX_EVENTS];
char *event_label      [MAX_THREADS][MAX_EVENTS];
int log_events = 1;
char tasklabel[200];
char taskcolor[200];

/* struct timeval {time_t tv_sec; suseconds_t tv_usec;}; */
double get_current_time(void) { 
  struct timeval tp; 
  gettimeofday( &tp, NULL );
  return tp.tv_sec + 1e-6 * tp.tv_usec; 
}



void diff_matrix(double *A, double *B, int NB, int BBM, int BBN, int M, int N)
{
    int X, Y, x, y;

    printf("Checking ...\n");
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
                  ; //printf(".");
              else if (e < 0.00000000001 )
                  ; //printf("%c", e < 0.00000000001 ? '.' : '#');
              else
                  printf("%c", e < 0.00000000001 ? '.' : '#');

              if (x == 2) x = NB-2;
            }
            else
            {
              printf("%c", '+');
            }
          }
          //printf("|");
        }
        //printf("\n");
        if (y == 2) y = NB-2;
      }
      if (Y < BBM-1)
          ; //int i;
          ; //for (i = 0; i < BBM*5; i++) printf("=");
      //printf("\n");
    }
    printf("\n");
}
///////////////////////////////////////////////////////////////////////////////////////////////////
#include <time.h>
void dump_trace(int cores_num)
{
    char trace_file_name[32];
    FILE *trace_file;
    double end;
    int event, core;
//    double scale = 500000.0;
    double scale = 30000.0;
//    double scale = 150000.0;

    //sprintf(trace_file_name, "trace_%d.svg", (int)(time(NULL)));
    sprintf(trace_file_name, "trace.svg");
    trace_file = fopen(trace_file_name, "w");
    assert(trace_file != NULL);

//     fprintf(trace_file,
//         "<svg width=\"200mm\" height=\"40mm\" viewBox=\"0 0 3000 400\">\n"
//         "  <g font-size=\"20\">\n");
    int x_end_max = -1;
    for (core = 0; core < cores_num; core++) {
      event = event_num[core]-4;
      end = event_log[core][event+2] - event_log[0][1];
      x_end_max = max(x_end_max, (int)(end * scale) + 20);
    }
/*     int y_max = (cores_num -1) * 100 + 40; */
/*     fprintf(trace_file, */
/*          "<?xml version=\"1.0\" standalone=\"no\"?>" */
/*          "<svg version=\"1.1\" baseProfile=\"full\" xmlns=\"http://www.w3.org/2000/svg\" " */
/*          "xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns:ev=\"http://www.w3.org/2001/xml-events\"  " */
/*          "width=\"200mm\" height=\"40mm\" viewBox=\"0 0 %d %d\">\n" */
/*          "  <g font-size=\"20\">\n", x_end_max, y_max ); */
    fprintf(trace_file,
            "<?xml version=\"1.0\" standalone=\"no\"?>"
            "<svg version=\"1.1\" baseProfile=\"full\" xmlns=\"http://www.w3.org/2000/svg\" "
            "xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns:ev=\"http://www.w3.org/2001/xml-events\"  "
            ">\n"
            "  <g font-size=\"20\">\n");
    for (core = 0; core < cores_num; core++)
        for (event = 0; event < event_num[core]; event += 4)
        {
            int    tag   = (int)event_log[core][event+0];
            double start =      event_log[core][event+1];
            double end   =      event_log[core][event+2];
            int    color = (int)event_log[core][event+3];
            char   *label = event_label[core][event+0];

            start -= event_log[0][1];
            end   -= event_log[0][1];

            fprintf(trace_file,
                "    "
                "<rect x=\"%.2lf\" y=\"%.0lf\" width=\"%.2lf\" height=\"%.0lf\" "
                "fill=\"#%06x\" stroke=\"#000000\" stroke-width=\"1\"/>\n",
                start * scale,  // x
                core * 100.0,   // y
                (end - start) * scale, // width
                90.0,           // height
                color);

/*             fprintf(trace_file, */
/*                 "    " */
/*                 "<text x=\"%.2lf\" y=\"%.0lf\" font-size=\"20\" fill=\"black\">" */
/*                 "%d" */
/*                 "</text>\n", */
/*                 start * scale + 10, // x */
/*                 core * 100.0 + 20, // y */
/*                 (int)tag); */

            fprintf(trace_file,
                "    "
                "<text x=\"%.2lf\" y=\"%.0lf\" font-size=\"20\" fill=\"black\">"
                "%s"
                "</text>\n",
                start * scale + 10, // x
                core * 100.0 + 20, // y
                label);
        }

    fprintf(trace_file,
        "  </g>\n"
        "</svg>\n");

    fclose(trace_file);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
#define core_dgetrf core_dgetrf_
#define core_dtstrf core_dtstrf_
#define core_dgessm core_dgessm_
#define core_dssssm core_dssssm_

#ifdef __cplusplus
extern "C" {
#endif

void core_dgetrf(int*, int*, int*, double*, int*, int*, int*);
void core_dtstrf(int*, int*, int*, int*, double*, int*, double*, int*, double*, int*, int*, int*);
void core_dgessm(int*, int*, int*, int*, int*, double*, int*, double*, int*, int*);
void core_dssssm(int*, int*, int*, int*, int*, double*, int*, double*, int*, double*, int*,
                 double*, int*, int*, int*);

#ifdef __cplusplus
}
#endif


///////////////////////////////////////////////////////////////////////////////////////////////////
void CORE_dgetrf(int M, int N, int IB, double *A, int LDA, int *IPIV, int *INFO)
{
    core_dgetrf(&M, &N, &IB, A, &LDA, IPIV, INFO);
}

void CORE_dtstrf(int M, int N, int IB, int NB, double *U, int LDU, double *A, int LDA, double *L,
                 int LDL, int *IPIV, int *INFO)
{
  core_dtstrf(  &M, &N, &IB, &NB, U, &LDU, A, &LDA, L, &LDL, IPIV, INFO);
}

void CORE_dgessm(int M, int N, int K, int IB, int *IPIV, double *L, int LDL, double *A, int LDA)
{
    int INFO;
    core_dgessm( &M, &N, &K, &IB, IPIV, L, &LDL, A, &LDA, &INFO);
}

void CORE_dssssm(int M1, int M2, int NN, int IB, int K, double *A1, int LDA1, double *A2, int LDA2,
                 double *L1, int LDL1, double *L2, int LDL2, int *IPIV)
{
    int INFO;
    core_dssssm(&M1, &M2, &NN, &IB, &K, A1, &LDA1, A2, &LDA2, L1, &LDL1, L2, &LDL2, IPIV, &INFO);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void SCHED_dgetrf(Quark *quark)
{
    int M;
    int N;
    int IB;
    double *A;
    int LDA;
    int *IPIV;
    int INFO;
    int k;
    quark_unpack_args_7(quark, M, N, IB, A, LDA, IPIV, k);
/* printf("DGETRF\n"); */
    core_event_start(QUARK_Thread_Rank(quark));
    core_dgetrf(
        &M, &N, &IB,
        A, &LDA,
        IPIV, &INFO);
    core_event_end(QUARK_Thread_Rank(quark));
    core_log_event_label(DGETRF_COLOR, QUARK_Thread_Rank(quark), k, QUARK_Get_Task_Label(quark) );
}

//////////////////////////////////////////////////////////////////////////////////////////
void SCHED_dtstrf(Quark *quark)
{
    int M;
    int N;
    int IB;
    int NB;
    double *U;
    int LDU;
    double *A;
    int LDA;
    double *L;
    int LDL;
    int *IPIV;
    int INFO;
    int k;
    quark_unpack_args_12(quark, M, N, IB, NB, U, LDU, A, LDA, L, LDL, IPIV, k);
/* printf("DTSTRF\n"); */
    core_event_start(QUARK_Thread_Rank(quark));
    core_dtstrf(
        &M, &N, &IB, &NB,
        U, &LDU,
        A, &LDA,
        L, &LDL,
        IPIV, &INFO);
    core_event_end(QUARK_Thread_Rank(quark));
    core_log_event_label(DTSTRF_COLOR, QUARK_Thread_Rank(quark), k, QUARK_Get_Task_Label(quark) );
}

//////////////////////////////////////////////////////////////////////////////////////////
void SCHED_dgessm(Quark *quark)
{
    int M;
    int N;
    int K;
    int IB;
    int *IPIV;
    double *L;
    int LDL;
    double *A;
    int LDA;
    int INFO;
    int n;
    quark_unpack_args_10(quark, M, N, K, IB, IPIV, L, LDL, A, LDA, n);
/* printf("DGESSM\n"); */
    core_event_start(QUARK_Thread_Rank(quark));
    core_dgessm(
        &M, &N, &K, &IB,
        IPIV,
        L, &LDL,
        A, &LDA,
        &INFO);
    core_event_end(QUARK_Thread_Rank(quark));
    core_log_event_label(DGESSM_COLOR, QUARK_Thread_Rank(quark), n, QUARK_Get_Task_Label(quark) );
}

//////////////////////////////////////////////////////////////////////////////////////////
void SCHED_dssssm(Quark *quark)
{
    int M1;
    int M2;
    int NN;
    int IB;
    int K;
    double *A1;
    int LDA1;
    double *A2;
    int LDA2;
    double *L1;
    int LDL1;
    double *L2;
    int LDL2;
    int *IPIV;
    int INFO;
    int n;
    quark_unpack_args_15(quark, M1, M2, NN, IB, K, A1, LDA1, A2, LDA2, L1, LDL1, L2, LDL2, IPIV, n);
/* printf("DSSSSM\n"); */
    core_event_start(QUARK_Thread_Rank(quark));
    core_dssssm(
        &M1, &M2, &NN, &IB, &K,
        A1, &LDA1,
        A2, &LDA2,
        L1, &LDL1,
        L2, &LDL2,
        IPIV,
        &INFO);
    core_event_end(QUARK_Thread_Rank(quark));
    core_log_event_label(DSSSSM_COLOR, QUARK_Thread_Rank(quark), n, QUARK_Get_Task_Label(quark) );
}
///////////////////////////////////////////////////////////////////////////////////////////////////
#define A(m,n) &A[NBNBSIZE*(m)+NBNBSIZE*MT*(n)]
#define L(m,n) &L[IBNBSIZE*(m)+IBNBSIZE*MT*(n)]
#define IPIV(m,n) &IPIV[NB*(m)+NB*MT*(n)]

void tile_el_you_sequential(double *A, double *L, int *IPIV, int IB, int NB, int BB)
{
    int NBNBSIZE = NB*NB;
    int IBNBSIZE = NB*IB;
    int MT = BB;

    int k, m, n;
    int IINFO;

    for (k = 0; k < BB; k++)
    {
        CORE_dgetrf(
            NB, NB, IB,
            A(k, k), NB,
            IPIV(k, k), &IINFO);

        for (n = k+1; n < BB; n++)
            CORE_dgessm(
                NB, NB, NB, IB,
                IPIV(k, k),
                A(k, k), NB,
                A(k, n), NB);

        for (m = k+1; m < BB; m++)
        {
            CORE_dtstrf(
                NB, NB, IB, NB,
                A(k, k), NB,
                A(m, k), NB,
                L(m, k), IB,
                IPIV(m, k), &IINFO);

            for (n = k+1; n < BB; n++)
                CORE_dssssm(
                    NB, NB, NB, IB, NB,
                    A(k, n), NB,
                    A(m, n), NB,
                    L(m, k), IB,
                    A(m, k), NB,
                    IPIV(m, k));
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void tile_el_you_parallel(Quark *quark, double *A, double *L, int *IPIV, int IB, int NB, int BB)
{
    int NBNBSIZE = NB*NB;
    int IBNBSIZE = NB*IB;
    int MT = BB;
    int zero=0;
    void *dummy = &zero;
    int k, m, n;
    
    Quark_sequence_t *seq1 = QUARK_Sequence_Create( quark );
    void *req1 = NULL;
    int priority = 0;

    for (k = 0; k < BB; k++) {

        core_event_start(0);
        snprintf( taskcolor, 200, "red");
        snprintf( tasklabel, 200, "%d %d", k, k);
        priority = 100;
        QUARK_Insert_Task(quark, SCHED_dgetrf, 0,
                           sizeof(int),          &NB,        VALUE,
                           sizeof(int),          &NB,        VALUE,
                           sizeof(int),          &IB,        VALUE,
                           sizeof(double)*NB*NB,        A(k, k),    INOUT | LOCALITY,
                           sizeof(int),          &NB,        VALUE,
                           sizeof(double)*NB,           IPIV(k, k), OUTPUT,
//            sizeof(double)*NB,   IPIV(k, k),         INPUT,
                           sizeof(int),          &k,         VALUE,
                           strlen(taskcolor)+1,  taskcolor,  VALUE | TASKCOLOR,
                           strlen(tasklabel)+1,  tasklabel,  VALUE | TASKLABEL,
                           sizeof(int),       &priority,  VALUE | TASK_PRIORITY,
/*                            sizeof(double)*NB*NB,        A(k, k),    INPUT, /\* Double one of the pointers TESTING TODO FIXME *\/ */
/*                            sizeof(double)*NB*NB,        A(k, k),    OUTPUT, /\* Double one of the pointers TESTING TODO FIXME *\/ */
/*                            sizeof(double)*NB*NB,        A(k, k),    INOUT, /\* Double one of the pointers TESTING TODO FIXME *\/ */
                           0);
        core_event_end(0);
        core_log_event(QUARK_COLOR, 0, 0);

        for (n = k+1; n < BB; n++) {
            snprintf( taskcolor, 200, "cyan");
            snprintf( tasklabel, 200, "%d %d", k, n);
/*             core_event_start(0); */
            priority = 50;
        priority = 100;
            QUARK_Insert_Task_In_Sequence(quark, SCHED_dgessm, seq1, req1, 0,
                               sizeof(int),          &NB,        VALUE,
                               sizeof(int),          &NB,        VALUE,
                               sizeof(int),          &NB,        VALUE,
                               sizeof(int),          &IB,        VALUE,
                               sizeof(double)*NB,           IPIV(k, k), INPUT,
                               sizeof(double)*NB*NB,        A(k, k),    NODEP | LOCALITY,
                               sizeof(int),          &NB,        VALUE,
                               sizeof(double)*NB*NB,        A(k, n),    INOUT,
//                sizeof(double)*NB*NB,   A(k, n),         INOUT ,
//                sizeof(double)*NB*NB,   A(k, n),         INPUT ,
                               sizeof(int),          &NB,        VALUE,
                               sizeof(int),          &n,         VALUE,
                               strlen(taskcolor)+1,  taskcolor,  VALUE | TASKCOLOR,
                               strlen(tasklabel)+1,  tasklabel,  VALUE | TASKLABEL,
                           sizeof(int),       &priority,  VALUE | TASK_PRIORITY,
                               0);
/*             core_event_end(0); */
/*             core_log_event(QUARK_COLOR, 0, 0); */
            
        }
        //QUARK_Sequence_Wait(quark, seq1);
        
        for (m = k+1; m < BB; m++) {
            snprintf( taskcolor, 200, "plum");
            snprintf( tasklabel, 200, "%d %d", m, k);
/*             core_event_start(0); */
            priority = 75;
        priority = 100;
            QUARK_Insert_Task(quark, SCHED_dtstrf, 0,
                               sizeof(int),          &NB,        VALUE,
                               sizeof(int),          &NB,        VALUE,
                               sizeof(int),          &IB,        VALUE,
                               sizeof(int),          &NB,        VALUE,
                               sizeof(double)*NB*NB,        A(k, k),    INOUT | LOCALITY,
                               sizeof(int),          &NB,        VALUE,
                               sizeof(double)*NB*NB,        A(m, k),    INOUT,
                               sizeof(int),          &NB,        VALUE,
                               sizeof(double)*NB*IB,        L(m, k),    OUTPUT,
                               sizeof(int),          &IB,        VALUE,
                               sizeof(int)*NB,           IPIV(m, k), OUTPUT,
                               sizeof(int),          &k,         VALUE,
                               2000,          NULL,         SCRATCH,
                               strlen(taskcolor)+1,  taskcolor,  VALUE | TASKCOLOR,
                               strlen(tasklabel)+1,  tasklabel,  VALUE | TASKLABEL,
                             sizeof(int),       &priority,  VALUE | TASK_PRIORITY,
                               0);
/*             core_event_end(0); */
/*             core_log_event(QUARK_COLOR, 0, 0); */

            for (n = k+1; n < BB; n++) {
                snprintf( taskcolor, 200, "green");
                snprintf( tasklabel, 200, "%d %d", m, n);
/*                 core_event_start(0); */
                priority = 10;
        priority = 100;
                QUARK_Insert_Task(quark, SCHED_dssssm, 0,
                                   sizeof(int),          &NB,        VALUE,
                                   sizeof(int),          &NB,        VALUE,
                                   sizeof(int),          &NB,        VALUE,
                                   sizeof(int),          &IB,        VALUE,
                                   sizeof(int),          &NB,        VALUE,
                                   sizeof(double)*NB*NB,        A(k, n),    INOUT | LOCALITY,
                                   sizeof(int),          &NB,        VALUE,
                                   sizeof(double)*NB*NB,        A(m, n),    INOUT,
                                   sizeof(int),          &NB,        VALUE,
                                   sizeof(double)*NB*IB,        L(m, k),    INPUT,
                                   sizeof(int),          &IB,        VALUE,
                                   sizeof(double)*NB*NB,        A(m, k),    INPUT,
                                   sizeof(int),          &NB,        VALUE,
                                   sizeof(int)*NB,           IPIV(m, k), INPUT,
                                   sizeof(int),          &n,         VALUE,
                                   strlen(taskcolor)+1,  taskcolor,  VALUE | TASKCOLOR,
                                   strlen(tasklabel)+1,  tasklabel,  VALUE | TASKLABEL,
                                  sizeof(int),       &priority,  VALUE | TASK_PRIORITY,
                                   0);
/*                 core_event_end(0); */
/*                 core_log_event(QUARK_COLOR, 0, 0); */
            }
        }
    }


}
/* #undef A */
/* #undef L */
/* #undef IPIV */

///////////////////////////////////////////////////////////////////////////////////////////////////
/* void tile_el_you_parallel(Quark *quark, double *A, double *L, int *IPIV, int IB, int NB, int BB) */
/* { */
/*     int NBNBSIZE = NB*NB; */
/*     int IBNBSIZE = NB*IB; */
/*     int MT = BB; */
/*     int zero=0; */
/*     void *dummy = &zero; */
/*     int k, m, n; */
    
/*     Quark_sequence_t *seq1 = QUARK_Sequence_Create( quark ); */
/*     void *req1 = NULL; */
/*     int priority = 0; */

/*     for (k = 0; k < BB; k++) { */
/*         core_event_start(0); */
/*         snprintf( taskcolor, 200, "red"); */
/*         snprintf( tasklabel, 200, "%d %d", k, k); */
/*         priority = 100; */
/*         QUARK_Insert_Task(quark, SCHED_dgetrf, 0,  */
/*                            sizeof(int),          &NB,        VALUE, */
/*                            sizeof(int),          &NB,        VALUE, */
/*                            sizeof(int),          &IB,        VALUE, */
/*                            sizeof(double)*NB*NB,        A(k, k),    INOUT | LOCALITY, */
/*                            sizeof(int),          &NB,        VALUE, */
/*                            sizeof(double)*NB,           IPIV(k, k), OUTPUT, */
/* //            sizeof(double)*NB,   IPIV(k, k),         INPUT, */
/*                            sizeof(int),          &k,         VALUE, */
/*                            strlen(taskcolor)+1,  taskcolor,  VALUE | TASKCOLOR, */
/*                            strlen(tasklabel)+1,  tasklabel,  VALUE | TASKLABEL, */
/*                            sizeof(int),       &priority,  VALUE | TASK_PRIORITY, */
/* /\*                            sizeof(double)*NB*NB,        A(k, k),    INPUT, /\\* Double one of the pointers TESTING TODO FIXME *\\/ *\/ */
/* /\*                            sizeof(double)*NB*NB,        A(k, k),    OUTPUT, /\\* Double one of the pointers TESTING TODO FIXME *\\/ *\/ */
/* /\*                            sizeof(double)*NB*NB,        A(k, k),    INOUT, /\\* Double one of the pointers TESTING TODO FIXME *\\/ *\/ */
/*                            0); */
/*         core_event_end(0); */
/*         core_log_event(QUARK_COLOR, 0, 0); */

/*         for (m = k+1; m < BB; m++) { */
/*             snprintf( taskcolor, 200, "plum"); */
/*             snprintf( tasklabel, 200, "%d %d", m, k); */
/* /\*             core_event_start(0); *\/ */
/*             priority = 75; */
/*             QUARK_Insert_Task(quark, SCHED_dtstrf, 0, */
/*                                sizeof(int),          &NB,        VALUE, */
/*                                sizeof(int),          &NB,        VALUE, */
/*                                sizeof(int),          &IB,        VALUE, */
/*                                sizeof(int),          &NB,        VALUE, */
/*                                sizeof(double)*NB*NB,        A(k, k),    INOUT | LOCALITY, */
/*                                sizeof(int),          &NB,        VALUE, */
/*                                sizeof(double)*NB*NB,        A(m, k),    INOUT, */
/*                                sizeof(int),          &NB,        VALUE, */
/*                                sizeof(double)*NB*IB,        L(m, k),    OUTPUT, */
/*                                sizeof(int),          &IB,        VALUE, */
/*                                sizeof(int)*NB,           IPIV(m, k), OUTPUT, */
/*                                sizeof(int),          &k,         VALUE, */
/*                                strlen(taskcolor)+1,  taskcolor,  VALUE | TASKCOLOR, */
/*                                strlen(tasklabel)+1,  tasklabel,  VALUE | TASKLABEL, */
/*                                2000,          NULL,         SCRATCH, */
/*                                sizeof(int),       &priority,  VALUE | TASK_PRIORITY, */
/*                                0); */
/* /\*             core_event_end(0); *\/ */
/* /\*             core_log_event(QUARK_COLOR, 0, 0); *\/ */
/*         } */

/*         for (n = k+1; n < BB; n++) { */
/*             snprintf( taskcolor, 200, "cyan"); */
/*             snprintf( tasklabel, 200, "%d %d", k, n); */
/* /\*             core_event_start(0); *\/ */
/*             priority = 50; */
/*             QUARK_Insert_Task_In_Sequence(quark, SCHED_dgessm, seq1, req1, 0,  */
/*                                sizeof(int),          &NB,        VALUE, */
/*                                sizeof(int),          &NB,        VALUE, */
/*                                sizeof(int),          &NB,        VALUE, */
/*                                sizeof(int),          &IB,        VALUE, */
/*                                sizeof(double)*NB,           IPIV(k, k), INPUT, */
/*                                sizeof(double)*NB*NB,        A(k, k),    NODEP, */
/*                                sizeof(int),          &NB,        VALUE, */
/*                                sizeof(double)*NB*NB,        A(k, n),    INOUT | LOCALITY, */
/* //                sizeof(double)*NB*NB,   A(k, n),         INOUT , */
/* //                sizeof(double)*NB*NB,   A(k, n),         INPUT , */
/*                                sizeof(int),          &NB,        VALUE, */
/*                                sizeof(int),          &n,         VALUE, */
/*                                strlen(taskcolor)+1,  taskcolor,  VALUE | TASKCOLOR, */
/*                                strlen(tasklabel)+1,  tasklabel,  VALUE | TASKLABEL, */
/*                            sizeof(int),       &priority,  VALUE | TASK_PRIORITY, */
/*                                0); */
/* /\*             core_event_end(0); *\/ */
/* /\*             core_log_event(QUARK_COLOR, 0, 0); *\/ */


/*             for (m = k+1; m < BB; m++) { */
/*                 snprintf( taskcolor, 200, "green"); */
/*                 snprintf( tasklabel, 200, "%d %d", m, n); */
/* /\*                 core_event_start(0); *\/ */
/*                 priority = 10; */
/*                 QUARK_Insert_Task(quark, SCHED_dssssm, 0, */
/*                                    sizeof(int),          &NB,        VALUE, */
/*                                    sizeof(int),          &NB,        VALUE, */
/*                                    sizeof(int),          &NB,        VALUE, */
/*                                    sizeof(int),          &IB,        VALUE, */
/*                                    sizeof(int),          &NB,        VALUE, */
/*                                    sizeof(double)*NB*NB,        A(k, n),    INOUT | LOCALITY, */
/*                                    sizeof(int),          &NB,        VALUE, */
/*                                    sizeof(double)*NB*NB,        A(m, n),    INOUT, */
/*                                    sizeof(int),          &NB,        VALUE, */
/*                                    sizeof(double)*NB*IB,        L(m, k),    INPUT, */
/*                                    sizeof(int),          &IB,        VALUE, */
/*                                    sizeof(double)*NB*NB,        A(m, k),    INPUT, */
/*                                    sizeof(int),          &NB,        VALUE, */
/*                                    sizeof(int)*NB,           IPIV(m, k), INPUT, */
/*                                    sizeof(int),          &n,         VALUE, */
/*                                    strlen(taskcolor)+1,  taskcolor,  VALUE | TASKCOLOR, */
/*                                    strlen(tasklabel)+1,  tasklabel,  VALUE | TASKLABEL, */
/*                                    sizeof(int),       &priority,  VALUE | TASK_PRIORITY, */
/*                                    0); */
/* /\*                 core_event_end(0); *\/ */
/* /\*                 core_log_event(QUARK_COLOR, 0, 0); *\/ */
/*             } */
/*         } */
/*     } */
/* } */
/* #undef A */
/* #undef L */
/* #undef IPIV */
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <pthread.h>

#define HUGE_PAGE_SIZE 2048*1024
int main (int argc, char **argv)
{
    assert(argc == 5);
    int IB = atoi(argv[1]);     // eg 20
    int NB = atoi(argv[2]); assert(NB%IB == 0); // eg 40, 80, 120
    int BB = atoi(argv[3]);     // eg 5
    int N = BB*NB;
    int M = N;
    //int NxN = N*N;
    int MxN = M*N;
    //int MxNB = M*NB;
    int THREADS = atoi(argv[4]); assert(THREADS <= MAX_THREADS); // eg 2, 4
    //int thread;

    //char mem_file_name[32];
    //unsigned long huge_size;
    //int  fmem;
    //char *mem_block = 0;
/*
    // Allocate memory in huge TLB pages
    double *Ablk  = (double*)mem_block; mem_block += MxN*sizeof(double);
    double *Ablk2 = (double*)mem_block; mem_block += MxN*sizeof(double);
    double *A     = (double*)mem_block; mem_block += MxN*sizeof(double);
    double *A2    = (double*)mem_block; mem_block += MxN*sizeof(double);
    double *L     = (double*)mem_block; mem_block += MxN*sizeof(double);
    int    *IPIV  = (int*)   mem_block; mem_block += MxN*sizeof(int);

    huge_size = (unsigned long)mem_block;
    huge_size = (huge_size + HUGE_PAGE_SIZE-1) & ~(HUGE_PAGE_SIZE-1);
    sprintf(mem_file_name, "/huge/huge_tlb_page.bin");
    assert((fmem = open(mem_file_name, O_CREAT | O_RDWR, 0755)) != -1);
    remove(mem_file_name);
    mem_block = (char*)mmap(0, huge_size, PROT_READ | PROT_WRITE, MAP_PRIVATE, fmem, 0);
    assert(mem_block != MAP_FAILED);

    Ablk  = (double*)(mem_block + (unsigned long)Ablk);
    Ablk2 = (double*)(mem_block + (unsigned long)Ablk2);
    A     = (double*)(mem_block + (unsigned long)A);
    A2    = (double*)(mem_block + (unsigned long)A2);
    L     = (double*)(mem_block + (unsigned long)L);
    IPIV  =    (int*)(mem_block + (unsigned long)IPIV);
*/
    double *Ablk  = (double*)malloc(MxN*sizeof(double));
    double *Ablk2 = (double*)malloc(MxN*sizeof(double));
    double *A     = (double*)malloc(MxN*sizeof(double));
    double *A2    = (double*)malloc(MxN*sizeof(double));
    double *L     = (double*)malloc(MxN*sizeof(double));
    int    *IPIV  = (int*)   malloc(MxN*sizeof(int));

    // Initialize A and A2
    int i;
    for (i = 0; i < MxN; i++)
        A[i] = A2[i] = 0.5 - (double)rand() / RAND_MAX;

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

    // Sequential tile LU in BDL
    memset(L, 0, MxN*sizeof(double));
    memset(IPIV, 0, MxN*sizeof(int));
    //tile_el_you_sequential(Ablk, L, IPIV, IB, NB, BB);

    memset(L, 0, MxN*sizeof(double));
    memset(IPIV, 0, MxN*sizeof(int));

    core_event_start(0);
    Quark *quark = QUARK_New(THREADS);
    core_event_end(0);
    core_log_event(QUARK_COLOR, 0, 0);
    assert(quark != NULL);
    double elapsed = -get_current_time();
    tile_el_you_parallel(quark, Ablk2, L, IPIV, IB, NB, BB);
    QUARK_Delete(quark);
    elapsed += get_current_time();

/*     core_event_start(0); */
/*     quark = quark_new(THREADS); */
/*     core_event_end(0); */
/*     core_log_event(QUARK_COLOR, 0, 0); */
/*     assert(quark != NULL); */
/*     elapsed = -get_current_time(); */
/*     tile_el_you_parallel(quark, Ablk2, L, IPIV, IB, NB, BB); */
/*     quark_delete(quark); */
/*     elapsed += get_current_time(); */

    double GFLOPS = 2.0*N*N*N/3.0 / elapsed / 1000000000;
    DBGPRINTF("Called tile_el_you_parallel\n");

    // Move from BDL to F77
    {int X, Y, x, y;
    for (X = 0; X < BB; X++)
      for (Y = 0; Y < BB; Y++)
        for (x = 0; x < NB; x++)
          for (y = 0; y < NB; y++)
            A[Y*NB + y + X*NB*N + x*N] = Ablk[Y*NB*NB + y + X*NB*NB*BB + x*NB];}

    {int X, Y, x, y;
    for (X = 0; X < BB; X++)
      for (Y = 0; Y < BB; Y++)
        for (x = 0; x < NB; x++)
          for (y = 0; y < NB; y++)
            A2[Y*NB + y + X*NB*N + x*N] = Ablk2[Y*NB*NB + y + X*NB*NB*BB + x*NB];}

    /* diff_matrix(A, A2, NB, BB, BB, M, N); */
    dump_trace(THREADS);
    printf("%s %d %d %d %d %d %.2lf\n", argv[0], IB, NB, BB, THREADS, N, GFLOPS);

    free(Ablk);
    free(Ablk2);
    free(A);
    free(A2);
    free(L);
    free(IPIV);

    return 0;
}
