/**
 * @file eztrace_convert_quark.c
 *
 * QUARK (QUeuing And Runtime for Kernels) provides a runtime
 * enviornment for the dynamic execution of precedence-constrained
 * tasks.
 *
 * QUARK is a software package provided by Univ. of Tennessee,
 * Univ. of California Berkeley and Univ. of Colorado Denver.
 *
 * @version 1.0.0
 * @author Vijay Joshi
 * @date 2013-01-23
 *
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <GTG.h>
#include <ev_codes.h>
#include <eztrace_list.h>
#include <eztrace_convert.h>
#include "quark.h"
#include "quark_trace.h"

#ifndef min
#define min( a, b ) ( (a) < (b) ? (a) : (b) )
#endif
#ifndef max
#define max( a, b ) ( (a) > (b) ? (a) : (b) )
#endif

#define QUARK_STATE      "ST_Thread"
#define QUARK_TASK_NAME  "Submitted Tasks counter"
#define QUARK_TASK_ALIAS "STasks"

#define QUARK_TASKR_NAME  "Global Ready Tasks counter"
#define QUARK_TASKR_ALIAS "GRTasks"

#define QUARK_TASKWR_NAME  "Local Ready Tasks counter"
#define QUARK_TASKWR_ALIAS "LRTasks"

#define QUARK_THREADS_MAX 4096


#define QUARK_NAME(ev) quark_array[ (int)(QUARK_GET_CODE(ev)) ].name
#define GTG_RANDOM gtg_get_random_color()

#define QUARK_INIT_EVENT( _idx_, _name_, _color_ )     \
    quark_array[_idx_].name  = _name_;                 \
    quark_array[_idx_].color = _color_;                \
    quark_array[_idx_].nb    = 0;                      \


/*
 * Check EZTrace version for compatibility
 */
#if !defined(EZTRACE_API_VERSION) || !(EZTRACE_API_VERSION > 0x00000400)
#error "EZTrace 0.7 or greater is required"
#endif

/*
 * Data for each event
 *   @name: name of the kernel (several kernels can use the same name, for
 *          example: gemm & gemm_f1, or laswp & laswpc
 *   @color: color assigne to the kernel
 *   @nb: number of calls to the kernel (with eztrace_stats)
 *   @sum: Total time spent executing this kernel (with eztrace_stats)
 *   @min: Minimal execution time of this kernel (with eztrace_stats)
 *   @max: Maximal execution time of this kernel (with eztrace_stats)
 */
typedef struct quark_s {
    char *name;
    gtg_color_t color;
    int  nb;
    double sum;
    double min;
    double max;
} quark_t;

static int quark_array_initialized = 0;
static quark_t quark_array[QUARK_NBMAX_EVENTS];
static gtg_color_t colors_array[20];

static int nbtasks_submitted = 0;
static int nbtasks_ready = 0;

/*
 * Keep information on thread status.
 *   @tid: Thread Id
 *   @active: Thread is active/inactive
 *   @lasttime: Start time of this task
 */
typedef struct quark_thrdstate_s {
    unsigned int tid;
    int          active;
    double       lasttime;
    int          nbtasks;
} quark_thrdstate_t;

static quark_thrdstate_t *thrdstate = NULL;
static int nbtrhd = 0;

static inline gtg_color_t gtg_get_random_color(){
    static int i = -1;
    i = (i+1)%20;
    return colors_array[i];
}


static int get_short_tid( unsigned int thread_id )
{
    int tid;

    for (tid=0; tid<nbtrhd; tid++) {
        if ( thrdstate[tid].tid == thread_id )
            break;
    }

    /* Thread not found, we add it */
    if ( tid == nbtrhd ) {
        if ( tid < QUARK_THREADS_MAX ) {
            thrdstate[nbtrhd].tid = thread_id;
            nbtrhd++;
        }
        else {
            fprintf(stderr, "Too many threads, increase QUARK_THREADS_MAX and recompile\n");
            return -1;
        }
    }
    return tid;
}

void handle_quark_start(eztrace_event_t *ev)
{
    FUNC_NAME;
    DECLARE_THREAD_ID_STR(_threadstr, CUR_INDEX, CUR_THREAD_ID);
    CHANGE() setState (CURRENT, QUARK_STATE, _threadstr, QUARK_NAME(ev) );
}

void handle_quark_task (eztrace_event_t *ev)
{
    FUNC_NAME;
    DECLARE_PROCESS_ID_STR( process_id, CUR_INDEX );
    //assert( GET_NBPARAMS(ev) == 1 );

    int value = (int)(GET_PARAM(ev, 0));
    fprintf(stderr, "task: %d %d\n", GET_NBPARAMS(ev), value );
    nbtasks_submitted += value;
    //assert( nbtasks_submitted >= 0 );
    if ( nbtasks_submitted < 0 ) {
        fprintf(stderr, "WARNING: Negative number of submitted tasks\n");
        /* Set the value to 0, to avoid bug in GTG (Eztrace <= 1.1) */
        nbtasks_submitted = 0;
    }
    CHANGE() setVar (CURRENT, QUARK_TASK_ALIAS, process_id, (varPrec)nbtasks_submitted);
}

void handle_quark_taskw (eztrace_event_t *ev)
{
    FUNC_NAME;
    //assert( GET_NBPARAMS(ev) == 2 );

    uint64_t thrdid = GET_PARAM(ev, 0);
    int value = GET_PARAM(ev, 1);
    int tid   = get_short_tid( thrdid );
    assert( tid != -1 );

    fprintf(stderr, "taskW: %d \n", GET_NBPARAMS(ev) );

    nbtasks_ready += value;
    thrdstate[tid].nbtasks += value;

    if ( nbtasks_ready < 0 ) {
        fprintf(stderr, "WARNING: Negative number of ready tasks\n");
        /* Set the value to 0, to avoid bug in GTG (Eztrace <= 1.1) */
        nbtasks_ready = 0;
    }
    if ( thrdstate[tid].nbtasks < 0 ) {
        fprintf(stderr, "WARNING: Negative number of ready tasks ready on thread: %lu\n", thrdid );
        /* Set the value to 0, to avoid bug in GTG (Eztrace <= 1.1) */
        thrdstate[tid].nbtasks = 0;
    }
    DECLARE_PROCESS_ID_STR( process_id, CUR_INDEX );
    DECLARE_THREAD_ID_STR(  thread_id,  CUR_INDEX, thrdid );
    CHANGE() setVar (CURRENT, QUARK_TASKR_ALIAS,  process_id, nbtasks_ready);
    CHANGE() setVar (CURRENT, QUARK_TASKWR_ALIAS, thread_id,  thrdstate[tid].nbtasks);
}

void
handle_quark_stop ()
{
    FUNC_NAME;
    DECLARE_THREAD_ID_STR(_threadstr, CUR_INDEX, CUR_THREAD_ID);
    CHANGE() setState (CURRENT, QUARK_STATE, _threadstr, "wait");
}

int
eztrace_convert_quark_init()
{
    int i;

    if ( quark_array_initialized == 0 ) {

        /* Initialize the colors_array */
        colors_array[ 0] = GTG_RED;
        colors_array[ 1] = GTG_GREEN;
        colors_array[ 2] = GTG_BLUE;
        colors_array[ 3] = GTG_WHITE;
        colors_array[ 4] = GTG_TEAL;
        colors_array[ 5] = GTG_DARKGREY;
        colors_array[ 6] = GTG_YELLOW;
        colors_array[ 7] = GTG_PURPLE;
        colors_array[ 8] = GTG_LIGHTBROWN;
        colors_array[ 9] = GTG_DARKBLUE;
        colors_array[10] = GTG_PINK;
        colors_array[11] = GTG_DARKPINK;
        colors_array[12] = GTG_SEABLUE;
        colors_array[13] = GTG_KAKI;
        colors_array[14] = GTG_REDBLOOD;
        colors_array[15] = GTG_BROWN;
        colors_array[16] = GTG_GRENAT;
        colors_array[17] = GTG_ORANGE;
        colors_array[18] = GTG_MAUVE;
        colors_array[19] = GTG_LIGHTPINK;

        /* First initialization to fill in the gap */
        for(i=0; i<QUARK_NBMAX_EVENTS; i++) {
            quark_array[i].name  = "";
            quark_array[i].color = GTG_RANDOM;
            quark_array[i].nb    = -1;
            quark_array[i].sum   = 0.;
            quark_array[i].min   = 999999999999.;
            quark_array[i].max   = 0.;
        }

        QUARK_INIT_EVENT(QUARK_WORKING, "Working",  GTG_BLUE);
        QUARK_INIT_EVENT(QUARK_INSERT_TASK, "QUARK_Insert_Task",  GTG_WHITE);
        QUARK_INIT_EVENT(QUARK_INSERT_TASK_PACKED, "QUARK_Insert_Task_Packed",  GTG_GREEN);
        QUARK_INIT_EVENT(QUARK_PROCESS_COMPLETED_TASKS, "QUARK_Process_Completed_Tasks",  GTG_BROWN);

        if ( thrdstate == NULL ) {
            thrdstate = (quark_thrdstate_t*)malloc(QUARK_THREADS_MAX * sizeof(quark_thrdstate_t));
            memset( thrdstate, 0, QUARK_THREADS_MAX * sizeof(quark_thrdstate_t));
        }

        quark_array_initialized = 1;
    }

    /*
     * Init only for trace conversion, not for statistics
     */
    if(get_mode() == EZTRACE_CONVERT) {
        addVarType( QUARK_TASK_ALIAS,   QUARK_TASK_NAME,   "CT_Process" );
        addVarType( QUARK_TASKR_ALIAS,  QUARK_TASKR_NAME,  "CT_Process" );
        addVarType( QUARK_TASKWR_ALIAS, QUARK_TASKWR_NAME, "CT_Thread"  );

        for(i=0; i<QUARK_NBMAX_EVENTS; i++) {
            if ( quark_array[i].nb == -1 )
                continue;

            addEntityValue( quark_array[i].name,
                            QUARK_STATE,
                            quark_array[i].name,
                            quark_array[i].color );
        }

        /* plasma quark */
        addEntityValue ("wait", QUARK_STATE, "wait", GTG_BLACK );
    }
    return 0;
}

int
handle_quark_events(eztrace_event_t *ev)
{
    switch (LITL_READ_GET_CODE(ev))
        {
        case FUT_QUARK(STOP)  : handle_quark_stop();
                                break;
        case FUT_QUARK(TASK)  : handle_quark_task(ev);
                                break;
        case FUT_QUARK(TASKW) : handle_quark_taskw(ev);
                                break;
        default:
            if ( IS_A_QUARK_EV(ev) ) {
                handle_quark_start(ev);
            }
            else
            {
                return 0;
            }
        }

    return 1;

}

void
eztrace_convert_quark_finalize()
{

}

int
handle_quark_stats(eztrace_event_t *ev)
{
    double time;

    if ( IS_A_QUARK_EV(ev) ) {
        switch (LITL_READ_GET_CODE(ev)) {
        case FUT_QUARK(STOP)  :
        {
            int tid = get_short_tid( CUR_THREAD_ID );
            assert(tid != -1);

            if ( thrdstate[tid].active == 0 ) {
                fprintf(stderr, "WARNING: The end of a state appears before the beginning\n");
                return 0;
            }

            time = ( CURRENT - thrdstate[tid].lasttime );

            /* Check that we have an existing state */
            assert(  quark_array[ thrdstate[tid].active ].nb >= 0 );

            if( quark_array[ thrdstate[tid].active ].nb == 0 ) {
                quark_array[ thrdstate[tid].active ].sum = 0.;
                quark_array[ thrdstate[tid].active ].max = 0.;
                quark_array[ thrdstate[tid].active ].min = 999999999999.;
            }
            quark_array[ thrdstate[tid].active ].nb++;
            quark_array[ thrdstate[tid].active ].sum += time;
            quark_array[ thrdstate[tid].active ].max = max( quark_array[ thrdstate[tid].active ].max, time );
            quark_array[ thrdstate[tid].active ].min = min( quark_array[ thrdstate[tid].active ].min, time );

            thrdstate[tid].active = 0;
            thrdstate[tid].lasttime = 0;
            return 1;
        }
        break;

        case FUT_QUARK(TASK):
        {
            int value = (int)(GET_PARAM(ev, 0));
            nbtasks_submitted += value;
            return 1;
        }
        break;

        case FUT_QUARK(TASKW):
        {
            int tid   = get_short_tid( GET_PARAM(ev, 0) );
            int value = GET_PARAM(ev, 1);
            assert(tid != -1);
            nbtasks_ready += value;
            thrdstate[tid].nbtasks += value;
            return 1;
        }
        break;

        default: /* All the different states */
        {
            int ev_code = (int)(QUARK_GET_CODE(ev));
            int tid = get_short_tid( CUR_THREAD_ID );

            assert(tid != -1);
            if ( thrdstate[tid].active != 0 ) {
                fprintf(stderr, "QUARK WARNING: thread %d change to state %d before to stop previous state %d\n",
                        (int)CUR_THREAD_ID, ev_code, thrdstate[tid].active );
            }

            thrdstate[tid].active = ev_code;
            thrdstate[tid].lasttime = CURRENT;
            return 1;
        }
        }
    }
    else
        return 0;
}

/*
 * Print the results of statistics.
 */
void print_quark_stats() {
    int i;

    printf ( "\nquark Module:\n");
    printf (   "-----------\n");

    for(i=0; i<QUARK_NBMAX_EVENTS; i++) {
        if ( quark_array[ i ].nb > 0 ) {
            printf ( "%s : %d calls\n"
                     "\tAverage time: %.3f ms\n"
                     "\tMaximun time: %.3f ms\n"
                     "\tMinimun time: %.3f ms\n",
                     quark_array[ i ].name,
                     quark_array[ i ].nb,
                     quark_array[ i ].sum / (double)(quark_array[ i ].nb),
                     quark_array[ i ].max,
                     quark_array[ i ].min);
        }
    }
}

struct eztrace_convert_module quark_module;

void libinit(void) __attribute__ ((constructor));
void libinit(void)
{
  quark_module.api_version = EZTRACE_API_VERSION;

  /* Specify the initialization function.
   * This function will be called once all the plugins are loaded
   * and the trace is started.
   * This function usually declared StateTypes, LinkTypes, etc.
   */
  quark_module.init = eztrace_convert_quark_init;

  /* Specify the function to call for handling an event
   */
  quark_module.handle = handle_quark_events;

  /* Specify the function to call for handling an event when eztrace_stats is called
   */
  quark_module.handle_stats = handle_quark_stats;

  /* Print the results of statistics
   */
  quark_module.print_stats = print_quark_stats;

  /* Specify the module prefix */
  quark_module.module_prefix = QUARK_EVENTS_ID;

  if ( asprintf(&quark_module.name, "quark") < 0 ) {
      fprintf(stderr, "Failed to create module name\n");
      exit(-1);
  }
  if ( asprintf(&quark_module.description, "Module for tracing the QUARK runtime") < 0 ) {
      fprintf(stderr, "Failed to create module description\n");
      exit(-1);
  }

  quark_module.token.data = &quark_module;

  /* Register the module to eztrace_convert */
  eztrace_convert_register_module(&quark_module);
}

void libfinalize(void) __attribute__ ((destructor));
void libfinalize(void)
{
}
