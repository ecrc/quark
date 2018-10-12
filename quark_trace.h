/**
 * @file quark_trace.h
 *
 * QUARK (QUeuing And Runtime for Kernels) provides a runtime
 * enviornment for the dynamic execution of precedence-constrained
 * tasks.
 *
 * QUARK is a software package provided by Univ. of Tennessee,
 * Univ. of California Berkeley and Univ. of Colorado Denver.
 *
 * @version 2.8.0
 * @author Mathieu Faverge (2010-11-15)
 * @author Vijay Joshi (2013-05-15)
 * @author Asim YarKhan (2013-05-15)
 * @date 2010-11-15
 *
 */
#ifndef _QUARK_TRACE_H_
#define _QUARK_TRACE_H_

enum quark_ev_code_e {
    QUARK_STOP,
    QUARK_WORKING,
    QUARK_TASK,
    QUARK_TASKW,
    QUARK_INSERT_TASK,
    QUARK_INSERT_TASK_PACKED,
    QUARK_PROCESS_COMPLETED_TASKS,
    QUARK_NBMAX_EVENTS,
};

#ifdef TRACE_QUARK

#include <eztrace.h>
#include <ev_codes.h>

#define QUARK_EVENTS_ID    USER_MODULE_ID(0x02)
#define QUARK_PREFIX       GENERATE_USER_MODULE_PREFIX(QUARK_EVENTS_ID)

#define PREFIX_MASK  (((1 << NB_BITS_PREFIX) -1) << NB_BITS_EVENTS)
#define EVENTS_MASK  (~PREFIX_MASK)

#define IS_A_QUARK_EV(ev)  ((LITL_READ_GET_CODE(ev) & PREFIX_MASK) == QUARK_PREFIX)
#define QUARK_GET_CODE(ev) (LITL_READ_GET_CODE(ev) & EVENTS_MASK)
#define FUT_QUARK(event)   (QUARK_PREFIX | QUARK_##event)

#define quark_trace_addtask()             EZTRACE_EVENT1(FUT_QUARK(TASK),  1)
#define quark_trace_deltask()             EZTRACE_EVENT1(FUT_QUARK(TASK), -1)
#define quark_trace_addtask2worker(__tid) EZTRACE_EVENT2(FUT_QUARK(TASKW), __tid,  1)
#define quark_trace_deltask2worker(__tid) EZTRACE_EVENT2(FUT_QUARK(TASKW), __tid, -1)

#define quark_trace_event_start(ev_code) EZTRACE_EVENT0(FUT_QUARK(ev_code));
#define quark_trace_event_end()          EZTRACE_EVENT0(FUT_QUARK(STOP));

#define quark_trace_function_start(task)                                \
    EZTRACE_EVENT2( task->start_code, task->taskid, task->sequence )
#define quark_trace_function_stop(task)         \
    EZTRACE_EVENT0( task->stop_code )

#define QUARK_TASK_DEFAULT_START_CODE FUT_QUARK(WORKING)
#define QUARK_TASK_DEFAULT_STOP_CODE  FUT_QUARK(STOP)

#else

#define quark_trace_addtask()
#define quark_trace_deltask()
#define quark_trace_addtask2worker(__tid)
#define quark_trace_deltask2worker(__tid)

#define quark_trace_event_start(ev_code)
#define quark_trace_event_end()

#define quark_trace_function_start(task) do {} while (0)
#define quark_trace_function_stop(task)  do {} while (0)

#define QUARK_TASK_DEFAULT_START_CODE 0
#define QUARK_TASK_DEFAULT_STOP_CODE  0

#endif

#endif /* _QUARK_TRACE_H_ */
