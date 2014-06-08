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
	/// the pointer to a container object
	typedef void* PdContainer;
	/// the pointer to a sequence object
	typedef void* PdSequenceX;

	/// Iterator for CoreArray array-oriented container
	/** sizeof(TdIterator) = 32 **/
	struct TdIterator
	{
        /// the handler of this iterator
		PdContainer* Handler;
		unsigned char VByte[32 - sizeof(PdContainer)];
	};


	// Functions for CdContainer - TdIterator

	bool gds_IterGetStart(PdContainer Node, TdIterator &Out);
	bool gds_IterGetEnd(PdContainer Node, TdIterator &Out);
	bool gds_IterAdv(TdIterator &Iter);
	// TODO: offset
	bool gds_IterAdvEx(TdIterator &Iter, const ssize_t offset);
	bool gds_IterPrev(TdIterator &Iter);
	bool gds_IterPrevEx(TdIterator &Iter, const ssize_t offset);
	size_t gds_IterRData(TdIterator &Iter, void *OutBuf, size_t Cnt, TSVType OutSV);
	size_t gds_IterWData(TdIterator &Iter, const void *InBuf, size_t Cnt, TSVType InSV);



	/// get the degree of dimension
	int gds_AttrNameIndex(PdGDSObj obj, const char *Name);


	PdGDSObj gds_NodePath(PdGDSObj Obj, const char *Path);

	int gds_NodeClassName(PdGDSObj Obj, char *OutStr, int OutBufLen);


	/// remove a GDS node
	bool gds_NodeDelete(PdGDSObj Node);


	/// get the degree of dimension
	int gds_SeqDimCnt(PdSequenceX obj);
	/// get the dimensions
	bool gds_SeqGetDim(PdSequenceX obj, int *OutBuf);

	/// get the total count of elements
	Int64 gds_SeqGetCount(PdSequenceX obj);

	/// get SVType
	int gds_SeqSVType(PdSequenceX obj);


	/// read the data
	bool gds_rData(PdSequenceX obj, CoreArray::Int32 const* Start,
		CoreArray::Int32 const* Length, void *OutBuf, TSVType OutSV);

	/// read the data
	bool gds_rDataEx(PdSequenceX obj, CoreArray::Int32 const* Start,
		CoreArray::Int32 const* Length, const CBOOL *const Selection[],
		void *OutBuf, TSVType OutSV);


	/// write the data
	bool gds_wData(PdSequenceX obj, CoreArray::Int32 const* Start,
		CoreArray::Int32 const* Length, const void *InBuf, TSVType InSV);

	/// append the data
	bool gds_AppendData(PdSequenceX obj, int Cnt, const void *InBuf, TSVType InSV);
	bool gds_AppendString(PdSequenceX obj, int Cnt, const char *buffer[]);
	bool gds_AppendString(PdSequenceX obj, const char *text);

	/// assign
	bool gds_Assign(PdSequenceX dest_obj, PdSequenceX src_obj, bool append);



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
	// ****	 the functions for block read
	//

	typedef void* PdArrayRead;

	/// read an array-oriented object margin by margin
	PdArrayRead gds_ArrayRead_Init(PdSequenceX Obj,
		int Margin, TSVType SVType, const CBOOL *const Selection[],
		bool buf_if_need=true);

	/// free a 'CdArrayRead' object
	bool gds_ArrayRead_Free(PdArrayRead Obj);

	/// read data
	bool gds_ArrayRead_Read(PdArrayRead Obj, void *Buffer);

	/// return true, if it is of the end
	bool gds_ArrayRead_Eof(PdArrayRead Obj);

	/// reallocate the buffer with specified size with respect to array
	bool gds_Balance_ArrayRead_Buffer(PdArrayRead array[],
		int n, Int64 buffer_size=-1);


	/// the class of read array
	class CArrayRead
	{
	public:
		CArrayRead(PdSequenceX Obj, int Margin, TSVType SVType,
			const CBOOL *const Selection[], bool buf_if_need=true)
		{
			_Obj = gds_ArrayRead_Init(Obj, Margin, SVType, Selection, buf_if_need);
			if (!_Obj)
				throw ErrCoreArray("Error 'initialize CArrayRead'.");
		}
		~CArrayRead() { gds_ArrayRead_Free(_Obj); }

		/// read data
		bool Read(void *Buffer)
		{
			return gds_ArrayRead_Read(_Obj, Buffer);
		}

		/// return true, if it is of the end
		bool Eof()
		{
			return gds_ArrayRead_Eof(_Obj);
		}

	protected:
		PdArrayRead _Obj;
	};



	// ******************************************************************
	// ****  the functions for R
	//

	// define _SEXP to avoid include <R.h> and <Rdefineh.h>
	typedef void* _SEXP_;

	bool gds_Is_R_Logical(PdGDSObj Obj);
	int gds_Set_If_R_Factor(PdGDSObj Obj, _SEXP_ val);
	_SEXP_ gds_Read_SEXP(PdSequenceX Obj, CoreArray::Int32 const* Start,
		CoreArray::Int32 const* Length, const CBOOL *const Selection[]);


	// ******************************************************************
	// ****  the functions for error messages
	//

	/// get the last error message
	std::string & gds_LastError();



	// ******************************************************************
	// ****  the R function
	//


}

#endif /* _CoreGDSLink_H_ */
