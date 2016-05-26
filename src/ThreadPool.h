// ===========================================================
//
// ThreadPool.h: the C++ implementation of thread pool
//
// Copyright (C) 2016    Xiuwen Zheng
//
// This file is part of SNPRelate.
//
// SNPRelate is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License Version 3 as
// published by the Free Software Foundation.
//
// SNPRelate is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU General Public
// License along with SNPRelate.
// If not, see <http://www.gnu.org/licenses/>.

/**
 *	\file     ThreadPool.h
 *	\author   Xiuwen Zheng [zhengxwen@gmail.com]
 *	\version  1.0
 *	\date     2016
 *	\brief    the C++ implementation of thread pool
 *	\details
**/

#ifndef _HEADER_THREAD_POOL_
#define _HEADER_THREAD_POOL_

#include <CoreDEF.h>
#include <R_GDS_CPP.h>
#include <vector>
#include <queue>

#ifdef COREARRAY_PLATFORM_WINDOWS
#   include <windows.h>
#endif

#ifdef COREARRAY_POSIX_THREAD
#   include <pthread.h>
#endif


namespace CoreArray
{
	using namespace std;

	// =====================================================================

	/// Mutex object
	class COREARRAY_DLL_DEFAULT CMutex
	{
	public:
		friend class CCondition;

		/// constructor
		CMutex();
		/// destructor
		~CMutex();

		/// lock the mutex object
		void Lock();
		/// unlock the mutex object
		void Unlock();
		/// attempt to lock the mutex object, return true if succeed
		bool TryLock();

	protected:

	#if defined(COREARRAY_PLATFORM_WINDOWS)
		#ifdef COREARRAY_CC_GNU_MINGW32
			typedef CRITICAL_SECTION TdMutex;
		#else
			typedef RTL_CRITICAL_SECTION TdMutex;
		#endif
	#elif defined(COREARRAY_POSIX_THREAD)
		typedef pthread_mutex_t TdMutex;
	#endif

		TdMutex mutex;
	};


	/// Auto object for locking and unlocking a mutex object
	struct COREARRAY_DLL_DEFAULT CAutoMutex
	{
	public:
		CAutoMutex(CMutex *m) { mutex = m; if (m) m->Lock(); }
		CAutoMutex(CMutex &m) { mutex = &m; m.Lock(); }
		~CAutoMutex() { if (mutex) mutex->Unlock(); mutex = NULL; }

		/// reset the mutex object
		void Reset(CMutex *m)
		{
			if (m != mutex)
			{
				if (mutex) mutex->Unlock();
				mutex = m;
				if (m) m->Lock();
			}
		}
	private:
		CMutex *mutex;
	};


	// =====================================================================

	/// Condition object
	class COREARRAY_DLL_DEFAULT CCondition
	{
	public:
		friend class CMutex;

		/// constructor
		CCondition();
		/// destructor
		~CCondition();

		/// signal a condition object
		void Signal();
		/// broadcast a condition object
		void Broadcast();
		/// wait for a condition object to become signaled
		void Wait(CMutex &mutex);

	protected:

	#if defined(COREARRAY_PLATFORM_WINDOWS)
		typedef struct {
			/// signal and broadcast event HANDLEs
			HANDLE events[2];
			/// count of the number of waiters
			size_t waiter_count;
			/// serialize access to waiter_count
			#ifdef COREARRAY_CC_GNU_MINGW32
				CRITICAL_SECTION waiter_count_mutex;
			#else
				RTL_CRITICAL_SECTION waiter_count_mutex;
			#endif
		} TCondition;
	#elif defined(COREARRAY_POSIX_THREAD)
		typedef pthread_cond_t TCondition;
	#endif

		TCondition cond;
	};


	// =====================================================================

    /// Thread class
	class COREARRAY_DLL_DEFAULT CThread
	{
	public:
	#if defined(COREARRAY_POSIX_THREAD)
		typedef pthread_t TThread;
	#elif defined(COREARRAY_PLATFORM_WINDOWS)
		typedef struct {
			HANDLE Handle;
			DWORD ThreadID;
		} TThread;
	#endif

		/// constructor
		CThread();
		/// destructor
		virtual ~CThread();

		/// start the threads
		void BeginThread();

		/// run the threads safely, calling RunThread()
		int RunThreadSafe();
		/// terminate the threads
		int EndThread();
		/// terminate the threads
		void Terminate();

		COREARRAY_INLINE bool Terminated() const { return terminated; }
		COREARRAY_INLINE TThread &Thread() { return thread; }
		COREARRAY_INLINE int &ExitCode() { return fExitCode; }
        COREARRAY_INLINE string &ErrorInfo() { return fErrorInfo; }

	protected:
		TThread thread;
		bool terminated;
		int fExitCode;
		string fErrorInfo;

		/// the thread procedure
		virtual int RunThread();

		void _BeginThread(void *ptr); // need vData
		void _Done();
	};


	// =====================================================================

	typedef void (*TProc)(size_t i, void *ptr);

	/// Thread pool
	class COREARRAY_DLL_DEFAULT CThreadPool
	{
	public:	
		CThreadPool(int num_threads);
		~CThreadPool();

		void AddWork(TProc proc, size_t i, void *ptr);
		void Wait();

	protected:

		class COREARRAY_DLL_DEFAULT CThread_in_Pool: public CThread
		{
		public:
			CThreadPool *Owner;
			CThread_in_Pool();
		protected:
			virtual int RunThread();
		};

		struct TProcData
		{
			TProc proc;
			size_t data;
			void *ptr;
			TProcData() { }
			TProcData(TProc p, size_t d, void *s) { proc=p; data=d; ptr=s; }
		};

		/// a collection of threads
		vector<CThread_in_Pool> workers;
		/// the task queue
		queue<TProcData> tasks;
		/// the number of working threads
		size_t num_threads_working;

		// synchronization
		CMutex mutex;
		CCondition thread_wait_cond;
		CCondition main_wait_cond;
		bool stop;
	};



	template<typename T1, typename T2> struct PARAM2
	{
		T1 p1;
		T2 p2;
		PARAM2() { }
		PARAM2(T1 v1, T2 v2) { p1=v1; p2=v2; }
	};



	// =====================================================================

	/// error exception for the gdsfmt package
	class COREARRAY_DLL_EXPORT ErrThread: public ErrCoreArray
	{
	public:
		ErrThread() {}
		ErrThread(const char *fmt, ...) { _COREARRAY_ERRMACRO_(fmt); }
		ErrThread(const std::string &msg) { fMessage = msg; }
	};

}

#endif /* _HEADER_THREAD_POOL_ */
