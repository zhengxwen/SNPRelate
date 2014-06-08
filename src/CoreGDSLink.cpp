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
// Copyright (C) 2013	Xiuwen Zheng
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


// Functions for CdContainer - TdIterator

typedef int (*TIterGet)(PdContainer Node, TdIterator &Out);
static TIterGet _IterGetStart = NULL;
bool GDSInterface::gds_IterGetStart(PdContainer Node, TdIterator &Out)
	{ return (*_IterGetStart)(Node, Out); }

static TIterGet _IterGetEnd = NULL;
bool GDSInterface::gds_IterGetEnd(PdContainer Node, TdIterator &Out)
	{ return (*_IterGetEnd)(Node, Out); }

typedef bool (*TIterStep)(TdIterator &Iter);
static TIterStep _IterAdv = NULL;
bool GDSInterface::gds_IterAdv(TdIterator &Iter)
	{ return (*_IterAdv)(Iter); }

typedef bool (*TIterStepEx)(TdIterator &Iter, const ssize_t offset);
static TIterStepEx _IterAdvEx = NULL;
bool GDSInterface::gds_IterAdvEx(TdIterator &Iter, const ssize_t offset)
	{ return (*_IterAdvEx)(Iter, offset); }

static TIterStep _IterPrev = NULL;
bool GDSInterface::gds_IterPrev(TdIterator &Iter)
	{ return (*_IterPrev)(Iter); }

static TIterStepEx _IterPrevEx = NULL;
bool GDSInterface::gds_IterPrevEx(TdIterator &Iter, const ssize_t offset)
	{ return (*_IterPrevEx)(Iter, offset); }

typedef bool (*TIterRData)(TdIterator &Iter, void *OutBuf, size_t Cnt, TSVType OutSV);
static TIterRData _IterRData = NULL;
size_t GDSInterface::gds_IterRData(TdIterator &Iter, void *OutBuf, size_t Cnt, TSVType OutSV)
	{ return (*_IterRData)(Iter, OutBuf, Cnt, OutSV); }

typedef bool (*TIterWData)(TdIterator &Iter, const void *InBuf, size_t Cnt, TSVType InSV);
static TIterWData _IterWData = NULL;
size_t GDSInterface::gds_IterWData(TdIterator &Iter, const void *InBuf, size_t Cnt, TSVType InSV)
	{ return (*_IterWData)(Iter, InBuf, Cnt, InSV); }



// get the degree of dimension
typedef int (*TAttrNameIndex)(PdGDSObj obj, const char *Name);
static TAttrNameIndex _AttrNameIndex = NULL;
int GDSInterface::gds_AttrNameIndex(PdGDSObj obj, const char *Name)
{
	return (*_AttrNameIndex)(obj, Name);
}


// get an object of GDS variable by a path
typedef PdGDSObj* (*TNodePath)(PdGDSObj Obj, const char *Path);
static TNodePath _NodePath = NULL;
PdGDSObj GDSInterface::gds_NodePath(PdGDSObj Obj, const char *Path)
{
	return (*_NodePath)(Obj, Path);
}

// get the class name of a node
typedef int (*TNodeClassName)(PdGDSObj Obj, char *OutStr, int OutBufLen);
static TNodeClassName _NodeClassName = NULL;
int GDSInterface::gds_NodeClassName(PdGDSObj Obj, char *OutStr, int OutBufLen)
{
	return (*_NodeClassName)(Obj, OutStr, OutBufLen);
}

// remove a GDS node
typedef bool (*TNodeDelete)(PdGDSObj Node);
static TNodeDelete _NodeDelete = NULL;
bool GDSInterface::gds_NodeDelete(PdGDSObj Node)
{
	return (*_NodeDelete)(Node);
}



/// get the degree of dimension
typedef int (*TSeqDimCnt)(PdSequenceX);
static TSeqDimCnt _SeqDimCnt = NULL;
int GDSInterface::gds_SeqDimCnt(PdSequenceX obj)
	{ return (*_SeqDimCnt)(obj); }

/// get the dimensions
typedef int (*TSeqGetDim)(PdSequenceX, int *);
static TSeqGetDim _SeqGetDim = NULL;
bool GDSInterface::gds_SeqGetDim(PdSequenceX obj, int *OutBuf)
	{ return (*_SeqGetDim)(obj, OutBuf); }

/// get the total count of elements
typedef Int64 (*TSeqGetCount)(PdSequenceX);
static TSeqGetCount _SeqGetCount = NULL;
Int64 GDSInterface::gds_SeqGetCount(PdSequenceX obj)
	{ return (*_SeqGetCount)(obj); }

/// get SVType
typedef int (*TSeqSVType)(PdSequenceX);
static TSeqSVType _SeqSVType = NULL;
int GDSInterface::gds_SeqSVType(PdSequenceX obj)
	{ return (*_SeqSVType)(obj); }


/// read the data
typedef bool (*TSeqRData)(PdSequenceX, CoreArray::Int32 const*,
		CoreArray::Int32 const*, void *, TSVType);
static TSeqRData _SeqrData = NULL;
bool GDSInterface::gds_rData(PdSequenceX obj, CoreArray::Int32 const* Start,
		CoreArray::Int32 const* Length, void *OutBuf, TSVType OutSV)
	{ return (*_SeqrData)(obj, Start, Length, OutBuf, OutSV); }

/// read the data
typedef bool (*TSeqRDataEx)(PdSequenceX, CoreArray::Int32 const*,
		CoreArray::Int32 const*, const CBOOL *const[], void *, TSVType);
static TSeqRDataEx _SeqrDataEx = NULL;
bool GDSInterface::gds_rDataEx(PdSequenceX obj, CoreArray::Int32 const* Start,
		CoreArray::Int32 const* Length, const CBOOL *const Selection[], void *OutBuf,
		TSVType OutSV)
	{ return (*_SeqrDataEx)(obj, Start, Length, Selection, OutBuf, OutSV); }


/// write the data
typedef bool (*TSeqWData)(PdSequenceX, CoreArray::Int32 const*,
		CoreArray::Int32 const*, const void *, TSVType);
static TSeqWData _SeqwData = NULL;
bool GDSInterface::gds_wData(PdSequenceX obj, CoreArray::Int32 const* Start,
		CoreArray::Int32 const* Length, const void *InBuf, TSVType InSV)
	{ return (*_SeqwData)(obj, Start, Length, InBuf, InSV); }


/// append the data
typedef bool (*TSeqAppendData)(PdSequenceX, int, const void *, TSVType);
static TSeqAppendData _SeqAppendData = NULL;
bool GDSInterface::gds_AppendData(PdSequenceX obj, int Cnt, const void *InBuf, TSVType InSV)
	{ return (*_SeqAppendData)(obj, Cnt, InBuf, InSV); }


typedef bool (*TSeqAppendString)(PdSequenceX obj, int Cnt, const char *buffer[]);
static TSeqAppendString _SeqAppendString = NULL;
bool GDSInterface::gds_AppendString(PdSequenceX obj, int Cnt, const char *buffer[])
	{ return (*_SeqAppendString)(obj, Cnt, buffer); }

bool GDSInterface::gds_AppendString(PdSequenceX obj, const char *text)
	{ return gds_AppendString(obj, 1, &text); }


typedef bool (*TSeqAssign)(PdSequenceX, PdSequenceX, bool);
static TSeqAssign _SeqAssign = NULL;
bool GDSInterface::gds_Assign(PdSequenceX dest_obj, PdSequenceX src_obj, bool append)
	{ return (*_SeqAssign)(dest_obj, src_obj, append); }



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
// ****	 the functions for block read
//

/// read an array-oriented object margin by margin
typedef PdArrayRead (*T_ArrayRead_Init)(PdSequenceX, int, TSVType,
	const CBOOL *const [], bool);
static T_ArrayRead_Init _ArrayRead_Init = NULL;
PdArrayRead GDSInterface::gds_ArrayRead_Init(PdSequenceX Obj,
		int Margin, TSVType SVType, const CBOOL *const Selection[],
		bool buf_if_need)
	{ return (*_ArrayRead_Init)(Obj, Margin, SVType, Selection, buf_if_need); }

/// free a 'CdArrayRead' object
typedef bool (*T_ArrayRead_Free)(PdArrayRead);
static T_ArrayRead_Free _ArrayRead_Free = NULL;
bool GDSInterface::gds_ArrayRead_Free(PdArrayRead Obj)
	{ return (*_ArrayRead_Free)(Obj); }

/// read data
typedef bool (*T_ArrayRead_Read)(PdArrayRead, void *);
static T_ArrayRead_Read _ArrayRead_Read = NULL;
bool GDSInterface::gds_ArrayRead_Read(PdArrayRead Obj, void *Buffer)
	{ return (*_ArrayRead_Read)(Obj, Buffer); }

/// return true, if it is of the end
typedef bool (*T_ArrayRead_Eof)(PdArrayRead);
static T_ArrayRead_Eof _ArrayRead_Eof = NULL;
bool GDSInterface::gds_ArrayRead_Eof(PdArrayRead Obj)
	{ return (*_ArrayRead_Eof)(Obj); }

/// reallocate the buffer with specified size with respect to array
typedef bool (*T_Balance_ArrayRead_Buffer)(PdArrayRead [], int, Int64);
static T_Balance_ArrayRead_Buffer _Balance_ArrayRead_Buffer = NULL;
bool GDSInterface::gds_Balance_ArrayRead_Buffer(PdArrayRead array[],
		int n, Int64 buffer_size)
	{ return (*_Balance_ArrayRead_Buffer)(array, n, buffer_size); }



// ******************************************************************
// ****  the functions for error messages  ****
//

typedef std::string &(*TLastError)();
static TLastError _LastError = NULL;
std::string & GDSInterface::gds_LastError()
	{ return (*_LastError)(); }



// ******************************************************************
// ****  the functions for R
//

typedef bool (*T_Is_R_Logical)(PdGDSObj Obj);
static T_Is_R_Logical _Is_R_Logical = NULL;
bool GDSInterface::gds_Is_R_Logical(PdGDSObj Obj)
	{ return (*_Is_R_Logical)(Obj); }

typedef int (*T_Set_If_R_Factor)(PdGDSObj, _SEXP_);
static T_Set_If_R_Factor _Set_If_R_Factor = NULL;
int GDSInterface::gds_Set_If_R_Factor(PdGDSObj Obj, _SEXP_ val)
	{ return (*_Set_If_R_Factor)(Obj, val); }


typedef _SEXP_ (*T_Read_SEXP)(PdSequenceX, CoreArray::Int32 const*,
		CoreArray::Int32 const*, const CBOOL *const []);
static T_Read_SEXP _Read_SEXP = NULL;
_SEXP_ GDSInterface::gds_Read_SEXP(PdSequenceX Obj, CoreArray::Int32 const* Start,
		CoreArray::Int32 const* Length, const CBOOL *const Selection[])
	{ return (*_Read_SEXP)(Obj, Start, Length, Selection); }




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

	LOAD(_IterGetStart, TIterGet, "gds_IterGetStart");
	LOAD(_IterGetEnd, TIterGet, "gds_IterGetEnd");
	LOAD(_IterAdv, TIterStep, "gds_IterAdv");
	LOAD(_IterAdvEx, TIterStepEx, "gds_IterAdvEx");
	LOAD(_IterPrev, TIterStep, "gds_IterPrev");
	LOAD(_IterPrevEx, TIterStepEx, "gds_IterPrevEx");
	LOAD(_IterRData, TIterRData, "gds_IterRData");
	LOAD(_IterWData, TIterWData, "gds_IterWData");


	LOAD(_AttrNameIndex, TAttrNameIndex, "gds_AttrNameIndex");

	LOAD(_NodePath, TNodePath, "gds_NodePath");
	LOAD(_NodeClassName, TNodeClassName, "gds_NodeClassName");
	LOAD(_NodeDelete, TNodeDelete, "gds_NodeDelete");

	LOAD(_SeqDimCnt, TSeqDimCnt, "gds_SeqDimCnt");
	LOAD(_SeqGetDim, TSeqGetDim, "gds_SeqGetDim");
	LOAD(_SeqGetCount, TSeqGetCount, "gds_SeqGetCount");
	LOAD(_SeqSVType, TSeqSVType, "gds_SeqSVType");

	LOAD(_SeqrData, TSeqRData, "gds_rData");
	LOAD(_SeqrDataEx, TSeqRDataEx, "gds_rDataEx");
	LOAD(_Read_SEXP, T_Read_SEXP, "gds_Read_SEXP");
	LOAD(_SeqwData, TSeqWData, "gds_wData");
	LOAD(_SeqAppendData, TSeqAppendData, "gds_AppendData");
	LOAD(_SeqAppendString, TSeqAppendString, "gds_AppendString");
	LOAD(_SeqAssign, TSeqAssign, "gds_Assign");

	// ****	 the functions for block read
	LOAD(_ArrayRead_Init, T_ArrayRead_Init, "gds_ArrayRead_Init");
	LOAD(_ArrayRead_Free, T_ArrayRead_Free, "gds_ArrayRead_Free");
	LOAD(_ArrayRead_Read, T_ArrayRead_Read, "gds_ArrayRead_Read");
	LOAD(_ArrayRead_Eof, T_ArrayRead_Eof, "gds_ArrayRead_Eof");
	LOAD(_Balance_ArrayRead_Buffer, T_Balance_ArrayRead_Buffer, "gds_Balance_ArrayRead_Buffer");

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

	// **** the R functions
	LOAD(_Is_R_Logical, T_Is_R_Logical, "gds_Is_R_Logical");
	LOAD(_Set_If_R_Factor, T_Set_If_R_Factor, "gds_Set_If_R_Factor");
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
