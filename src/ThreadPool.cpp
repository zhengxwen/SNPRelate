// ===========================================================
//
// ThreadPool.cpp: the implementation of thread pool
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
 *	\file     dVect.h
 *	\author   Xiuwen Zheng [zhengxwen@gmail.com]
 *	\version  1.0
 *	\date     2007 - 2015
 *	\brief    Classess and functions for vectorization
 *	\details
**/

#include "ThreadPool.h"


namespace CoreArray
{

// =========================================================================
// CMutex

#ifdef COREARRAY_POSIX_THREAD
static const char *ERR_PTHREAD = "'%s' returns an error code (%d).";
#endif

CMutex::CMutex()
{
#if defined(COREARRAY_POSIX_THREAD)
   	int v = pthread_mutex_init(&mutex, NULL);
   	if (v != 0)
   		throw ErrThread(ERR_PTHREAD, "pthread_mutex_init", v);
#elif defined(COREARRAY_PLATFORM_WINDOWS)
	InitializeCriticalSection(&mutex);
#endif
}

CMutex::~CMutex()
{
#if defined(COREARRAY_POSIX_THREAD)
	int v = pthread_mutex_destroy(&mutex);
	if (v != 0)
		throw ErrThread(ERR_PTHREAD, "pthread_mutex_destroy", v);
#elif defined(COREARRAY_PLATFORM_WINDOWS)
	DeleteCriticalSection(&mutex);
#endif
}

void CMutex::Lock()
{
#if defined(COREARRAY_POSIX_THREAD)
	int v = pthread_mutex_lock(&mutex);
	if (v != 0)
		throw ErrThread(ERR_PTHREAD, "pthread_mutex_lock", v);
#elif defined(COREARRAY_PLATFORM_WINDOWS)
	EnterCriticalSection(&mutex);
#endif
}

void CMutex::Unlock()
{
#if defined(COREARRAY_POSIX_THREAD)
	int v = pthread_mutex_unlock(&mutex);
	if (v != 0)
		throw ErrThread(ERR_PTHREAD, "pthread_mutex_unlock", v);
#elif defined(COREARRAY_PLATFORM_WINDOWS)
	LeaveCriticalSection(&mutex);
#endif
}

bool CMutex::TryLock()
{
#if defined(COREARRAY_POSIX_THREAD)
	return (pthread_mutex_trylock(&mutex) == 0);
#elif defined(COREARRAY_PLATFORM_WINDOWS)
	return (TryEnterCriticalSection(&mutex) != 0);
#endif
}


// =========================================================================
// CCondition

CCondition::CCondition()
{
#if defined(COREARRAY_POSIX_THREAD)

	int v = pthread_cond_init(&cond, NULL);
	if (v != 0)
		throw ErrThread(ERR_PTHREAD, "pthread_cond_init", v);

#elif defined(COREARRAY_PLATFORM_WINDOWS)

	static const char *ERR_CONDITION =
		"Internal error when initializing CCondition.";

	cond.waiter_count = 0;
	InitializeCriticalSection(&cond.waiter_count_mutex);

	// initialize the events
	cond.events[0] = CreateEvent(NULL, FALSE, FALSE, NULL);
	if (cond.events[0] == NULL)
		throw ErrThread(ERR_CONDITION);
	cond.events[1] = CreateEvent(NULL, TRUE, FALSE, NULL);
	if (cond.events[1] == NULL)
	{
		CloseHandle(cond.events[0]);
		throw ErrThread(ERR_CONDITION);
	}

#endif
}

CCondition::~CCondition()
{
#if defined(COREARRAY_POSIX_THREAD)
	pthread_cond_destroy(&cond);
#elif defined(COREARRAY_PLATFORM_WINDOWS)
	if (cond.events[0] != NULL)
		CloseHandle(cond.events[0]);
	if (cond.events[1] != NULL)
		CloseHandle(cond.events[1]);
	DeleteCriticalSection(&cond.waiter_count_mutex);
#endif
}

void CCondition::Signal()
{
#if defined(COREARRAY_POSIX_THREAD)
	int v = pthread_cond_signal(&cond);
	if (v != 0)
		throw ErrThread(ERR_PTHREAD, "pthread_cond_signal", v);
#elif defined(COREARRAY_PLATFORM_WINDOWS)
	EnterCriticalSection(&cond.waiter_count_mutex);
	bool has_any_waiter = (cond.waiter_count > 0);
	LeaveCriticalSection(&cond.waiter_count_mutex);
	if (has_any_waiter)
	{
		if (SetEvent(cond.events[0]) == 0)
			RaiseLastOSError<ErrThread>();
	}
#endif
}

void CCondition::Broadcast()
{
#if defined(COREARRAY_POSIX_THREAD)
	int v = pthread_cond_broadcast(&cond);
	if (v != 0)
		throw ErrThread(ERR_PTHREAD, "pthread_cond_broadcast", v);
#elif defined(COREARRAY_PLATFORM_WINDOWS)
	EnterCriticalSection(&cond.waiter_count_mutex);
	bool has_any_waiter = (cond.waiter_count > 0);
	LeaveCriticalSection(&cond.waiter_count_mutex);
	if (has_any_waiter)
	{
		if (SetEvent(cond.events[1]) == 0)
			RaiseLastOSError<ErrThread>();
	}
#endif
}

void CCondition::Wait(CMutex &mutex)
{
#if defined(COREARRAY_POSIX_THREAD)

	int v = pthread_cond_wait(&cond, &mutex.mutex);
	if (v != 0)
		throw ErrThread(ERR_PTHREAD, "pthread_cond_wait", v);

#elif defined(COREARRAY_PLATFORM_WINDOWS)
	// increase the number of waiters
	EnterCriticalSection(&cond.waiter_count_mutex);
	cond.waiter_count ++;
	LeaveCriticalSection(&cond.waiter_count_mutex);

	// release the mutex object
	mutex.Unlock();

	// wait for either event to become signaled
	DWORD rv = WaitForMultipleObjects(2, cond.events, FALSE, INFINITE);
	if (rv == WAIT_TIMEOUT)
		throw ErrThread("condition object wait time out.");
	else if (rv == WAIT_FAILED)
		RaiseLastOSError<ErrThread>();

	// check if we are the last waiter
	EnterCriticalSection(&cond.waiter_count_mutex);
	cond.waiter_count --;
	bool IsLast = (rv == (WAIT_OBJECT_0 + 1)) && (cond.waiter_count == 0);
	LeaveCriticalSection(&cond.waiter_count_mutex);

	// reset the event if last waiter
	if (IsLast)
	{
		if (ResetEvent(cond.events[1]) == 0)
			RaiseLastOSError<ErrThread>();
	}

	// Re-acquire the mutex
	mutex.Lock();
#endif
}

}


// =========================================================================
// CThread

extern "C"
{
#if defined(COREARRAY_POSIX_THREAD)
	void* COREARRAY_CALL_ALIGN ThreadWrap(void *lpParam)
	{
		CoreArray::CThread *p = (CoreArray::CThread *)lpParam;
		ssize_t rv = p->RunThreadSafe();
		return (void*)rv;
	}
#elif defined(COREARRAY_PLATFORM_WINDOWS)
	DWORD WINAPI COREARRAY_CALL_ALIGN ThreadWrap(LPVOID lpParam)
	{
		CoreArray::CThread *p = (CoreArray::CThread *)lpParam;
		return p->RunThreadSafe();
	}
#endif
}


namespace CoreArray
{

CThread::CThread()
{
	terminated = false;
	fExitCode = 0;
	memset(&thread, 0, sizeof(thread));
}

CThread::~CThread()
{
	try {
		Terminate();
		EndThread();
	} catch (...) {
		_Done(); throw;
	}
	_Done();
}

void CThread::_Done()
{
#if defined(COREARRAY_POSIX_THREAD)
	if (thread)
	{
		pthread_detach(thread);
		thread = 0;
	}
#elif defined(COREARRAY_PLATFORM_WINDOWS)
	if (thread.Handle != NULL)
	{
		CloseHandle(thread.Handle);
		thread.Handle = NULL;
	}
#endif
}

void CThread::BeginThread()
{
#if defined(COREARRAY_POSIX_THREAD)
	if (!thread)
	{
		int v = pthread_create(&thread, NULL, ThreadWrap, (void*)this);
		if (v != 0)
			throw ErrThread(ERR_PTHREAD, "pthread_create", v);
	} else
    	throw ErrThread("BeginThread");
#elif defined(COREARRAY_PLATFORM_WINDOWS)
	if (thread.Handle == NULL)
	{
		SECURITY_ATTRIBUTES attr;
			attr.nLength = sizeof(attr);
			attr.lpSecurityDescriptor = NULL;
			attr.bInheritHandle = true;
		thread.Handle = CreateThread(&attr, 0, ThreadWrap, (void*)this,
			0, &thread.ThreadID);
		if (thread.Handle == NULL)
			RaiseLastOSError<ErrThread>();
	} else
		throw ErrThread("BeginThread");
#endif
}

int CThread::RunThreadSafe()
{
	try {
		fExitCode = RunThread();
	}
	catch (exception &E) {
    	fErrorInfo = E.what(); fExitCode = -1;
	}
	catch (const char *E) {
		fErrorInfo = E; fExitCode = -1;
	}
	catch (...) {
        fExitCode = -1;
    }
	return fExitCode;
}

int CThread::RunThread()
{
	return 0;
}

void CThread::Terminate()
{
	terminated = true;
}

int CThread::EndThread()
{
#if defined(COREARRAY_POSIX_THREAD)
	if (thread)
	{
		int v = pthread_join(thread, NULL);
		if (v != 0)
			throw ErrThread(ERR_PTHREAD, "pthread_join", v);
		_Done();
	}
#elif defined(COREARRAY_PLATFORM_WINDOWS)
	if (thread.Handle != NULL)
	{
		if (WaitForSingleObject(thread.Handle, INFINITE) == WAIT_FAILED)
			RaiseLastOSError<ErrThread>();
		_Done();
	}
#endif
	return fExitCode;
}


// =====================================================================

CThreadPool::CThread_in_Pool::CThread_in_Pool(): CThread()
{
	Owner = NULL;
}

int CThreadPool::CThread_in_Pool::RunThread()
{
	while (Owner)
	{
		CThreadPool::TProcData pc;
		{
			CAutoMutex lck(Owner->mutex);
			while (!Owner->stop && Owner->tasks.empty())
				Owner->thread_wait_cond.Wait(Owner->mutex);
			if (Owner->stop && Owner->tasks.empty())
				return 0;
			pc = Owner->tasks.front();
			Owner->tasks.pop();
			Owner->num_threads_working ++;
		}
		(*pc.proc)(pc.data, pc.ptr);
		{
			CAutoMutex lck(Owner->mutex);
			Owner->num_threads_working --;
		}
		Owner->main_wait_cond.Signal();
	}

	return 0;
}


CThreadPool::CThreadPool(int num_threads)
{
	stop = false;
	num_threads_working = 0;
	if (num_threads > 1)
	{
		workers.resize(num_threads);
		for(int i=0; i < num_threads; i++)
		{
			workers[i].Owner = this;
			workers[i].BeginThread();
		}
	}
}

CThreadPool::~CThreadPool()
{
	{
		CAutoMutex lck(mutex);
		stop = true;
	}
	thread_wait_cond.Broadcast();
	main_wait_cond.Broadcast();
	workers.clear();
}

void CThreadPool::AddWork(TProc proc, size_t i, void *ptr)
{
	if (workers.empty())
	{
		num_threads_working ++;
		(*proc)(i, ptr);
		num_threads_working --;
	} else {
		{
			CAutoMutex lck(mutex);
			if(stop)
				throw "AddWork on stopped CThreadPool";
			tasks.push(TProcData(proc, i, ptr));
		}
		thread_wait_cond.Signal();
	}
}

void CThreadPool::Wait()
{
	if (!workers.empty())
	{
		CAutoMutex lck(mutex);
		while (!stop && (num_threads_working>0 || !tasks.empty()))
			main_wait_cond.Wait(mutex);
	}
}

}
