/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * DARTS - A fine-grain dataflow-inspired runtime system.                          *
 * Copyright (C) 2011-2014  University of Delaware                                 *
 *                                                                                 *
 * This library is free software; you can redistribute it and/or                   *
 * modify it under the terms of the GNU Lesser General Public                      *
 * License as published by the Free Software Foundation; either                    *
 * version 2.1 of the License, or (at your option) any later version.              *
 *                                                                                 *
 * This library is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of                  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU               *
 * Lesser General Public License for more details.                                 *
 *                                                                                 *
 * You should have received a copy of the GNU Lesser General Public                *
 * License along with this library; if not, write to the Free Software             *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA  *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef SIMPLIFYING_DARTS_H
#define SIMPLIFYING_DARTS_H

#include <stdint.h>
#include <darts.h>
#include "getClock.h"
//#include <cuda.h>
//#include <cuda_runtime_api.h>

using namespace darts;

#define TPPUSHFULL 0
#define TPROUNDROBIN 1
#define MCSTANDARD 0

//extern size_t darts::g_nCU, darts::g_nSU;
extern size_t g_nCU, g_nSU;

#define TOTAL_NUM_CORES (g_nCU+1)*g_nSU
#define TOTAL_NUM_CU    (g_nCU+0)*g_nSU

//#include <ctime>
//inline uint64_t getTime(void)
//{
//    timespec time;
//    clock_gettime(CLOCK_REALTIME, &time);
//    return (static_cast<uint64_t>(time.tv_sec)*1000000000 + static_cast<uint64_t>(time.tv_nsec));
//}


using namespace darts;

#ifdef TRACE
	#define DEF_TP(TPName) struct TPName : public ThreadedProcedureTime
#else
	#define DEF_TP(TPName) struct TPName : public ThreadedProcedure
#endif


#define DEF_CODELET_ITER(name,deps,wait)               \
class name : public darts::Codelet                     \
{                                                      \
private:                                               \
    uint64_t _id;                                      \
public:                                                \
    name(uint32_t           dep,   uint32_t reset,     \
         ThreadedProcedure *frame, uint64_t meta,      \
         uint64_t id)                                  \
    : Codelet(dep,reset,frame,meta,id)                 \
    {                                                  \
    }                                                  \
    name(darts::ThreadedProcedure *frame=0)            \
    : Codelet(deps,deps,frame,wait)                    \
    , _id(0)                                           \
    {                                                  \
    }                                                  \
    virtual void fire();                               \
}                                                      

#define DEF_CODELET(name,deps,wait)                    \
struct name : public darts::Codelet                    \
{                                                      \
    name(uint32_t           dep=0,   uint32_t reset=0, \
         ThreadedProcedure *frame=0, uint64_t meta=0)  \
    : Codelet(dep,reset,frame,meta)                    \
    {                                                  \
    }                                                  \
    name(darts::ThreadedProcedure *frame)              \
    : Codelet(deps,deps,frame,wait)                    \
    {                                                  \
    }                                                  \
    virtual void fire();                               \
} 

#define LOAD_FRAME(TPName) TPName* frame = static_cast<TPName*>(myTP_)
#define GET_FRAME()        (*frame)
#define FRAME(field)       (*frame).field
#define INVOKE(TPName,...) invoke<TPName>(frame,__VA_ARGS__)
#define PLACE(clusterID, TPName,...) invoke<TPName>(frame,__VA_ARGS__)
#define SIGNAL(field)      (*frame).field->decDep()
#define SYNC(field)        (*frame).field.decDep()
#define INCR(field)		   (*frame).field.incDep()
//#define ADD(cd_name)       frame->cd_name.add()
//#define ADD(cd_name)       frame->add(frame->cd_name)
#define ADD(cd_name)       add(frame->cd_name)
#define RESET(field)        (*frame).field.resetCodelet()

#define SETSYNC(field,dep,reset)        (*frame).field.setSync(dep,reset)

#define CDSYNC(field) field.decDep()


#define DARTS_EXIT() Runtime::finalSignal.decDep()
#define EXIT_TP() return




#endif // SIMPLIFYING_DARTS_H




