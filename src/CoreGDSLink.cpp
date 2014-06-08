// ===========================================================
//     _/_/_/   _/_/_/  _/_/_/_/    _/_/_/_/  _/_/_/   _/_/_/
//      _/    _/       _/             _/    _/    _/   _/   _/
//     _/    _/       _/_/_/_/       _/    _/    _/   _/_/_/
//    _/    _/       _/             _/    _/    _/   _/
// _/_/_/   _/_/_/  _/_/_/_/_/     _/     _/_/_/   _/_/
// ===========================================================
//
// CoreGDSLink.cpp: C interface for CoreArray's dynamic library
//
// Copyright (C) 2011	Xiuwen Zheng
//
// This file is part of CoreArray.
//
// CoreArray is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License Version 3 as
// published by the Free Software Foundation.
//
// CoreArray is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with CoreArray.
// If not, see <http://www.gnu.org/licenses/>.

#include <CoreGDSLink.h>

#if defined(COREARRAY_UNIX)
#  include <dlfcn.h>
#elif defined(COREARRAY_WINDOWS)
#  include <windows.h>
#else
#  error "not supported"
#endif

using namespace GDSInterface;

// ErrCoreArray

void ErrCoreArray::Init(const char *fmt, va_list arglist)
{
	char buf[1024];
	vsnprintf(buf, sizeof(buf), fmt, arglist);
	fMessage = buf;
}


/// get the degree of dimension
typedef int (*TAttrNameIndex)(TdGDSObj obj, const char *Name);
static TAttrNameIndex _AttrNameIndex = NULL;
int GDSInterface::gds_AttrNameIndex(TdGDSObj obj, const char *Name)
	{ return (*_AttrNameIndex)(obj, Name); }


/// get the degree of dimension
typedef int (*TSeqDimCnt)(TdSequenceX obj);
static TSeqDimCnt _SeqDimCnt = NULL;
int GDSInterface::gds_SeqDimCnt(TdSequenceX obj)
	{ return (*_SeqDimCnt)(obj); }


/// get the dimensions
typedef int (*TSeqGetDim)(TdSequenceX obj, int *OutBuf);
static TSeqGetDim _SeqGetDim = NULL;
bool GDSInterface::gds_SeqGetDim(TdSequenceX obj, int *OutBuf)
	{ return (*_SeqGetDim)(obj, OutBuf); }


/// read the data
typedef bool (*TSeqRData)(TdSequenceX obj, Int32 const* Start,
		Int32 const* Length, void *OutBuf, TSVType OutSV);
static TSeqRData _SeqrData = NULL;
bool GDSInterface::gds_rData(TdSequenceX obj, Int32 const* Start,
		Int32 const* Length, void *OutBuf, TSVType OutSV)
	{ return (*_SeqrData)(obj, Start, Length, OutBuf, OutSV); }

/// read the data
typedef bool (*TSeqRDataEx)(TdSequenceX obj, Int32 const* Start,
		Int32 const* Length, CBOOL *Selection[], void *OutBuf, TSVType OutSV);
static TSeqRDataEx _SeqrDataEx = NULL;
bool GDSInterface::gds_rDataEx(TdSequenceX obj, Int32 const* Start,
		Int32 const* Length, CBOOL *Selection[], void *OutBuf, TSVType OutSV)
	{ return (*_SeqrDataEx)(obj, Start, Length, Selection, OutBuf, OutSV); }

/// write the data
typedef bool (*TSeqWData)(TdSequenceX obj, Int32 const* Start,
		Int32 const* Length, const void *InBuf, TSVType InSV);
static TSeqWData _SeqwData = NULL;
bool GDSInterface::gds_wData(TdSequenceX obj, Int32 const* Start,
		Int32 const* Length, const void *InBuf, TSVType InSV)
	{ return (*_SeqwData)(obj, Start, Length, InBuf, InSV); }

/// append the data
typedef bool (*TSeqAppendData)(TdSequenceX obj, int Cnt, const void *InBuf,
	TSVType OutSV);
static TSeqAppendData _SeqAppendData = NULL;
bool GDSInterface::gds_AppendData(TdSequenceX obj, int Cnt, const void *InBuf, TSVType InSV)
	{ return (*_SeqAppendData)(obj, Cnt, InBuf, InSV); }



// ******************************************************************
// ****  the functions for machine configuration
//

// NaN
typedef float (*TF32)();
static TF32 _F32_NaN = NULL;
float GDSInterface::conf_F32_NaN() { return (*_F32_NaN)(); }
typedef double (*TF64)();
static TF64 _F64_NaN = NULL;
double GDSInterface::conf_F64_NaN() { return (*_F64_NaN)(); }

// Return infinity
static TF32 _F32_Inf = NULL;
float GDSInterface::conf_F32_Inf() { return (*_F32_Inf)(); }
static TF64 _F64_Inf = NULL;
double GDSInterface::conf_F64_Inf() { return (*_F64_Inf)(); }
// Return negative infinity
static TF32 _F32_NegInf = NULL;
float GDSInterface::conf_F32_NegInf() { return (*_F32_NegInf)(); }
static TF64 _F64_NegInf = NULL;
double GDSInterface::conf_F64_NegInf() { return (*_F64_NegInf)(); }

// If it is finite
typedef bool (*TIsFinite32)(float val);
static TIsFinite32 _IsFinite32 = NULL;
bool GDSInterface::conf_IsFinite32(float val) { return (*_IsFinite32)(val); }
typedef bool (*TIsFinite64)(double val);
static TIsFinite64 _IsFinite64 = NULL;
bool GDSInterface::conf_IsFinite64(double val) { return (*_IsFinite64)(val); }

typedef int (*TSingleFunc)();
static TSingleFunc _GetNumberOfCPU = NULL;
static TSingleFunc _GetL1CacheMemory = NULL;
static TSingleFunc _GetL2CacheMemory = NULL;

int GDSInterface::conf_GetNumberOfCPU() { return (*_GetNumberOfCPU)(); }
int GDSInterface::conf_GetL1CacheMemory() { return (*_GetL1CacheMemory)(); }
int GDSInterface::conf_GetL2CacheMemory() { return (*_GetL2CacheMemory)(); }



// ******************************************************************
// ****  the functions for parellel computing
//

// create a mutex object
typedef TdMutex (*TInitMutex)();
static TInitMutex _InitMutex = NULL;
TdMutex GDSInterface::plc_InitMutex() { return (*_InitMutex)(); }

// destroy the mutex object
typedef bool (*TObjMutex)(TdMutex obj);
static TObjMutex _DoneMutex = NULL;
bool GDSInterface::plc_DoneMutex(TdMutex obj) { return (*_DoneMutex)(obj); }

// lock the mutex object
static TObjMutex _LockMutex = NULL;
bool GDSInterface::plc_LockMutex(TdMutex obj) { return (*_LockMutex)(obj); }

// unlock the mutex object
static TObjMutex _UnlockMutex = NULL;
bool GDSInterface::plc_UnlockMutex(TdMutex obj) { return (*_UnlockMutex)(obj); }

// parallel functions
typedef bool (*TDoBaseThread)(void (*Proc)(TdThread, int, void*), void *param, int nThread);
static TDoBaseThread _DoBaseThread = NULL;
bool GDSInterface::plc_DoBaseThread(void (*Proc)(TdThread, int, void*), void *param, int nThread)
	{ return _DoBaseThread(Proc, param, nThread); }


// initialize a thread suspending object
typedef TdThreadsSuspending (*TInitSuspend)();
static TInitSuspend _InitSuspend = NULL;
TdThreadsSuspending GDSInterface::plc_InitSuspend() { return (*_InitSuspend)(); }

// destroy the thread suspending object
typedef bool (*TObjSuspend)(TdThreadsSuspending obj);
static TObjSuspend _DoneSuspend = NULL;
bool GDSInterface::plc_DoneSuspend(TdThreadsSuspending obj) { return (*_DoneSuspend)(obj); }

// suspend the thread suspending object
static TObjSuspend _Suspend = NULL;
bool GDSInterface::plc_Suspend(TdThreadsSuspending obj) { return (*_Suspend)(obj); }

// wakeup the thread suspending object
static TObjSuspend _WakeUp = NULL;
bool GDSInterface::plc_WakeUp(TdThreadsSuspending obj) { return (*_WakeUp)(obj); }



// ******************************************************************
// ****  the functions for error messages  ****
//

typedef std::string &(*TLastError)();
static TLastError _LastError = NULL;
std::string & GDSInterface::gds_LastError()
	{ return (*_LastError)(); }




// ******************************************************************

template <typename TO, typename FROM> TO nasty_cast(FROM f)
{
	union {
		FROM f; TO t;
	} u;
	u.f = f;
	return u.t;
}

#if defined(COREARRAY_UNIX)

	static void* GDS_Handle = NULL;
	#define LOAD(var, type, name)	\
		var = nasty_cast<type>(dlsym(GDS_Handle, name)); \
		if ((err = dlerror()) != NULL) \
			{ throw ErrCoreArray(err); DoneGDSInterface(); }

#elif defined(COREARRAY_WINDOWS)

	static HMODULE GDS_Handle = NULL;
	#define LOAD(var, type, name)	\
		var = nasty_cast<type>(GetProcAddress(GDS_Handle, name)); \
		if (var == NULL) \
			{ throw ErrCoreArray("No %s function.", name); }

#else
#  error "not supported"
#endif



void GDSInterface::InitGDSInterface(const char *lib_fn)
{
	#if defined(COREARRAY_UNIX)

		// open dll
		GDS_Handle = dlopen(lib_fn, RTLD_LAZY);
		if (!GDS_Handle) throw ErrCoreArray(dlerror());
		dlerror();  // Clear any existing error
		const char *err;

	#elif defined(COREARRAY_WINDOWS)

		GDS_Handle = LoadLibrary(lib_fn);
		if (GDS_Handle == NULL)
		{
			char buf[1024];
			memset((void*)buf, 0, sizeof(buf));
			FormatMessage(FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS |
				FORMAT_MESSAGE_ARGUMENT_ARRAY, NULL,
				GetLastError(), 0, buf, sizeof(buf), NULL);
			throw ErrCoreArray(buf);
		}

	#else
	#  error "not supported"
	#endif

	// load functions
	LOAD(_AttrNameIndex, TAttrNameIndex, "gds_AttrNameIndex");
	LOAD(_SeqDimCnt, TSeqDimCnt, "gds_SeqDimCnt");
	LOAD(_SeqGetDim, TSeqGetDim, "gds_SeqGetDim");
	LOAD(_SeqrData, TSeqRData, "gds_rData");
	LOAD(_SeqrDataEx, TSeqRDataEx, "gds_rDataEx");
	LOAD(_SeqwData, TSeqWData, "gds_wData");
	LOAD(_SeqAppendData, TSeqAppendData, "gds_AppendData");

	// ****  the functions for parellel computing  ****
	LOAD(_InitMutex, TInitMutex, "plc_InitMutex");
	LOAD(_DoneMutex, TObjMutex, "plc_DoneMutex");
	LOAD(_LockMutex, TObjMutex, "plc_LockMutex");
	LOAD(_UnlockMutex, TObjMutex, "plc_UnlockMutex");
	LOAD(_InitSuspend, TInitSuspend, "plc_InitSuspend");
	LOAD(_DoneSuspend, TObjSuspend, "plc_DoneSuspend");
	LOAD(_Suspend, TObjSuspend, "plc_Suspend");
	LOAD(_WakeUp, TObjSuspend, "plc_WakeUp");
	LOAD(_DoBaseThread, TDoBaseThread, "plc_DoBaseThread");

	// ****  the functions for error messages  ****
	LOAD(_LastError, TLastError, "gds_LastError");

	// ****  the functions for machine configuration  ****
	LOAD(_F32_NaN, TF32, "conf_F32_NaN");
	LOAD(_F64_NaN, TF64, "conf_F64_NaN");
	LOAD(_F32_Inf, TF32, "conf_F32_Inf");
	LOAD(_F64_Inf, TF64, "conf_F64_Inf");
	LOAD(_F32_NegInf, TF32, "conf_F32_NegInf");
	LOAD(_F64_NegInf, TF64, "conf_F64_NegInf");
	LOAD(_IsFinite32, TIsFinite32, "conf_IsFinite32");
	LOAD(_IsFinite64, TIsFinite64, "conf_IsFinite64");
	LOAD(_GetNumberOfCPU, TSingleFunc, "conf_GetNumberOfCPU");
	LOAD(_GetL1CacheMemory, TSingleFunc, "conf_GetL1CacheMemory");
	LOAD(_GetL2CacheMemory, TSingleFunc, "conf_GetL2CacheMemory");
}

void GDSInterface::DoneGDSInterface()
{
	#if defined(COREARRAY_UNIX)
		if (GDS_Handle)
		{
			dlclose(GDS_Handle);
			GDS_Handle = NULL;
		}
	#elif defined(COREARRAY_WINDOWS)
		if (GDS_Handle != NULL)
			FreeLibrary(GDS_Handle);
	#else
	#  error "not supported"
	#endif
}
