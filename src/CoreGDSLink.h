// ===========================================================
//     _/_/_/   _/_/_/  _/_/_/_/    _/_/_/_/  _/_/_/   _/_/_/
//      _/    _/       _/             _/    _/    _/   _/   _/
//     _/    _/       _/_/_/_/       _/    _/    _/   _/_/_/
//    _/    _/       _/             _/    _/    _/   _/
// _/_/_/   _/_/_/  _/_/_/_/_/     _/     _/_/_/   _/_/
// ===========================================================
//
// CoreGDSLink.h: Link C interface from the CoreArray's dynamic library "CoreGDS"
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

/**
 *	\file     CoreGDSLink.h
 *	\author   Xiuwen Zheng
 *	\version  1.0
 *	\date     2007 - 2013
 *	\brief    Link C interface from the CoreArray's dynamic library "CoreGDS"
 *	\details
**/


#ifndef _CoreGDSLink_H_
#define _CoreGDSLink_H_

#include <dType.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>

#include <limits>
#include <string>


namespace GDSInterface
{
	using namespace CoreArray;

	/// Error Macro
	#define _COREARRAY_ERRMACRO_(x) { \
		va_list args; va_start(args, x); \
		Init(x, args); \
		va_end(args); \
	}

	/// The root class of exception for GWAS library
	class ErrCoreArray: public std::exception
	{
	public:
		ErrCoreArray() {};
		ErrCoreArray(const char *fmt, ...) { _COREARRAY_ERRMACRO_(fmt); }
		ErrCoreArray(const std::string &msg) { fMessage = msg; }
		virtual const char *what() const throw() { return fMessage.c_str(); }
		virtual ~ErrCoreArray() throw() {};
	protected:
		std::string fMessage;
		void Init(const char *fmt, va_list arglist);
	};


	/// initialize
	void InitGDSInterface(const char *lib_fn);
	/// finalize
	void DoneGDSInterface();


	// ******************************************************************
	// ****  the functions for a GDS file
	//

	/// the interface for CoreArray library

	/// the pointer to a GDS file
	typedef void* PdGDSFile;
	/// the pionter to a GDS node
	typedef void* PdGDSObj;
	/// the pointer to a sequence object
	typedef void* PdSequenceX;



	/// get the degree of dimension
	int gds_AttrNameIndex(PdGDSObj obj, const char *Name);



	PdGDSObj gds_NodePath(PdGDSObj Obj, const char *Path);



	/// get the degree of dimension
	int gds_SeqDimCnt(PdSequenceX obj);
	/// get the dimensions
	bool gds_SeqGetDim(PdSequenceX obj, int *OutBuf);

	/// get the total count of elements
	Int64 gds_SeqGetCount(PdSequenceX obj);

	/// read the data
	bool gds_rData(PdSequenceX obj, CoreArray::Int32 const* Start,
		CoreArray::Int32 const* Length, void *OutBuf, TSVType OutSV);

	/// read the data
	bool gds_rDataEx(PdSequenceX obj, CoreArray::Int32 const* Start,
		CoreArray::Int32 const* Length, CBOOL *Selection[], void *OutBuf, TSVType OutSV);

	/// write the data
	bool gds_wData(PdSequenceX obj, CoreArray::Int32 const* Start,
		CoreArray::Int32 const* Length, const void *InBuf, TSVType InSV);

	/// append the data
	bool gds_AppendData(PdSequenceX obj, int Cnt, const void *InBuf, TSVType InSV);
	bool gds_AppendString(PdSequenceX obj, int Cnt, const char *buffer[]);
	bool gds_AppendString(PdSequenceX obj, const char *text);




	// ******************************************************************
	// ****  the functions for machine configuration
	//

	/// Return NaN
	double conf_F64_NaN();
	float conf_F32_NaN();
	/// Return infinity
	float conf_F32_Inf();
	double conf_F64_Inf();
	/// Return negative infinity
	float conf_F32_NegInf();
	double conf_F64_NegInf();
	/// If it is finite
	bool conf_IsFinite32(float val);
	bool conf_IsFinite64(double val);

	int conf_GetNumberOfCPU();
	int conf_GetL1CacheMemory();
	int conf_GetL2CacheMemory();



	// ******************************************************************
	// ****  the functions for parellel computing
	//

	/// the class of mutex object
	typedef void* TdMutex;

	/// create a mutex object
	TdMutex plc_InitMutex();
	/// destroy the mutex object
	bool plc_DoneMutex(TdMutex obj);
	/// lock the mutex object
	bool plc_LockMutex(TdMutex obj);
	/// unlock the mutex object
	bool plc_UnlockMutex(TdMutex obj);

	/// automatic mutex object
	struct TdAutoMutex
	{
		TdMutex mutex;
		TdAutoMutex(TdMutex m) { mutex = m; if (m) plc_LockMutex(m); };
		~TdAutoMutex() { if (mutex) plc_UnlockMutex(mutex); };
		inline void Reset(TdMutex m)
		{
			if (m != mutex) {
				if (mutex) plc_UnlockMutex(mutex);
				mutex = m;
				if (m) plc_LockMutex(m);
			}
		}
	};

	/// the class of suspending object
	typedef void* TdThreadsSuspending;

	/// initialize a thread suspending object
	TdThreadsSuspending plc_InitSuspend();
	/// destroy the thread suspending object
	bool plc_DoneSuspend(TdThreadsSuspending obj);
	/// suspend the thread suspending object
	bool plc_Suspend(TdThreadsSuspending obj);
	/// wakeup the thread suspending object
	bool plc_WakeUp(TdThreadsSuspending obj);





	/// the class of thread object
	typedef void* TdThread;

	bool plc_DoBaseThread(void (*Proc)(TdThread, int, void*), void *param, int nThread);



	// ******************************************************************
	// ****  the functions for error messages
	//

	/// get the last error message
	std::string & gds_LastError();
}

#endif /* _CoreGDSLink_H_ */

