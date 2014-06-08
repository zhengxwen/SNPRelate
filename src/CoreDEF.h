// ===========================================================
//     _/_/_/   _/_/_/  _/_/_/_/    _/_/_/_/  _/_/_/   _/_/_/
//      _/    _/       _/             _/    _/    _/   _/   _/
//     _/    _/       _/_/_/_/       _/    _/    _/   _/_/_/
//    _/    _/       _/             _/    _/    _/   _/
// _/_/_/   _/_/_/  _/_/_/_/_/     _/     _/_/_/   _/_/
// ===========================================================
//
// CoreDEF.h: CoreArray global macro
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
 *	\file     CoreDEF.h
 *	\author   Xiuwen Zheng
 *	\version  1.0
 *	\date     2007 - 2013
 *	\brief    CoreArray global macro
 *	\details
**/

#ifndef _COREDEF_H_
#define _COREDEF_H_


// C++ Compiler MACRO

//  GNU C++:
#if defined __GNUC__
#   define COREARRAY_GNUG
#   if defined __MINGW32__
#      define COREARRAY_MINGW32
#   endif

//  Kai C++
#elif defined __KCC
#   define COREARRAY_KCC

//  SGI MIPSpro C++
#elif defined __sgi
#   define COREARRAY_SGICC

//  Compaq Tru64 Unix cxx
#elif defined __DECCXX
#   define COREARRAY_DECCXX

//  Greenhills C++
#elif defined __ghs
#   define COREARRAY_GHS

//  Borland
#elif defined __BORLANDC__
#   define COREARRAY_BORLANDC

//  Intel
#elif defined(__ICL) || defined(__ICC)
#   define COREARRAY_ICC

//  Metrowerks CodeWarrior
#elif defined  __MWERKS__
#   define COREARRAY_MWERKS

//  Sun Workshop Compiler C++
#elif defined  __SUNPRO_CC
#   define COREARRAY_SUNPROCC

//  HP aCC
#elif defined __HP_aCC
#   define COREARRAY_HPACC

//  MPW MrCpp or SCpp
#elif defined(__MRC__) || defined(__SC__)
#   define COREARRAY_MRC

//  IBM Visual Age
#elif defined(__IBMCPP__)
#   define COREARRAY_IBMCPP

//  Comeau C++
# elif defined(__COMO__)
#   define COREARRAY_COMO

//  Microsoft Visual C++
#elif defined(_MSC_VER)
//
//  Must remain the last #elif since some other vendors (Metrowerks, for
//  example) also #define _MSC_VER
#  define COREARRAY_MSC

#else

// not recognise the compiler:
#  error "Unknown compiler."

#endif


// Platform MACRO

// linux:
#if defined(linux) || defined(__linux) || defined(__linux__)
#  define COREARRAY_LINUX

// BSD:
#elif defined(__FreeBSD__) || defined(__NetBSD__) || defined(__OpenBSD__)
#  define COREARRAY_BSD

// solaris:
#elif defined(sun) || defined(__sun)
#  define COREARRAY_SUN
#  if defined(__SVR4) || defined(__svr4__)
#    define COREARRAY_SOLARIS
#  else
#    define COREARRAY_SUNOS
#  endif

// SGI Irix:
#elif defined(__sgi)
#  define COREARRAY_SGI

// hp unix:
#elif defined(__hpux)
#  define COREARRAY_HPUNIX

// cygwin is not win32:
#elif defined(__CYGWIN__)
#  define COREARRAY_CYGWIN

// win32:
#elif defined(_WIN32) || defined(__WIN32__) || defined(WIN32) || defined(_WIN32_WCE)
#  define COREARRAY_WIN32
#  define COREARRAY_WINDOWS

// BeOS
#elif defined(__BEOS__)
#  define COREARRAY_BEOS

// MacOS
#elif defined(macintosh) || defined(__APPLE__) || defined(__APPLE_CC__)
#  define COREARRAY_MACOS

// IBM
#elif defined(__IBMCPP__)
#  define COREARRAY_IBM

// not recognise the platform:
#elif defined(COREARRAY_ASSERT_CONFIG)
#  error "Unknown platform!"

#endif


#if defined(unix) || defined(__unix) || defined(__unix__) || defined(COREARRAY_MACOS)
#  define COREARRAY_UNIX
#endif




// ******************************************************************************
// MACRO for POSIX thread

#if defined(COREARRAY_UNIX) || defined(posix) || defined(_posix) || defined(__posix)
#  define COREARRAY_POSIX_THREAD
#endif




// ******************************************************************************
// Streaming SIMD Extensions

#if defined(__SSE__)
#    define COREARRAY_SIMD_SSE
#    ifndef COREARRAY_VT_SIMD
#        define COREARRAY_VT_SIMD
#    endif
#endif

#if defined(__SSE2__)
#    define COREARRAY_SIMD_SSE2
#    ifndef COREARRAY_VT_SIMD
#        define COREARRAY_VT_SIMD
#    endif
#endif

#if defined(__SSE3__)
#    define COREARRAY_SIMD_SSE3
#    ifndef COREARRAY_VT_SIMD
#        define COREARRAY_VT_SIMD
#    endif
#endif

#if defined(__SSE4__)
#    define COREARRAY_SIMD_SSE4
#    ifndef COREARRAY_VT_SIMD
#        define COREARRAY_VT_SIMD
#    endif
#endif

#ifdef COREARRAY_DONT_SIMD
#  ifdef COREARRAY_VT_SIMD
#    undef COREARRAY_VT_SIMD
#    undef COREARRAY_SIMD_SSE
#    undef COREARRAY_SIMD_SSE2
#    undef COREARRAY_SIMD_SSE3
#    undef COREARRAY_SIMD_SSE4
#  endif
#endif




// ******************************************************************************
// Detecting the endianness (byte order)

#if (!defined(COREARRAY_LITTLE_ENDIAN)) && (!defined(COREARRAY_BIG_ENDIAN))
#
#  if defined(COREARRAY_UNIX)
#    if defined(COREARRAY_MACOS)
#      if defined(__i386__) || defined(__x86_64__)
#        define COREARRAY_LITTLE_ENDIAN
#      elif defined(__ppc__) || defined(__ppc64__) || defined(__arm__)
#        define COREARRAY_BIG_ENDIAN
#      else
#         error "Unsupported Apply Mac architecture!"
#      endif
#    elif defined(COREARRAY_SUN)
#      if defined(__i386) || defined(__x86_64) || defined(__amd64)
#        define COREARRAY_LITTLE_ENDIAN
#      else
#        define COREARRAY_BIG_ENDIAN
#      endif
#    else
#      include <endian.h>
#      if __BYTE_ORDER == __LITTLE_ENDIAN
#        define COREARRAY_LITTLE_ENDIAN
#      elif __BYTE_ORDER == __BIG_ENDIAN
#        define COREARRAY_BIG_ENDIAN
#      endif
#    endif
#  elif defined(COREARRAY_WINDOWS)
#    define COREARRAY_LITTLE_ENDIAN
#  else
#    error "Unsupported architecture!"
#  endif
#
#endif





// ******************************************************************************
// Floating point data type

#ifdef COREARRAY_HAVE_FLOAT128
#  undef COREARRAY_HAVE_FLOAT128
#endif

#ifdef COREARRAY_LONGFLOAT_IS_DOUBLE
#  undef COREARRAY_LONGFLOAT_IS_DOUBLE
#endif

#include <float.h>
#if (LDBL_MANT_DIG == 113)
#  define COREARRAY_HAVE_FLOAT128
#elif (LDBL_MANT_DIG == 53)
#  define COREARRAY_LONGFLOAT_IS_DOUBLE
#elif (LDBL_MANT_DIG != 64)
//#  error "Unable to determine long double."
#endif




// ******************************************************************************
// The stack is forced to be 16-byte aligned

#if defined(__GNUC__) && (__GNUC__ > 4 || __GNUC__ == 4 && __GNUC_MINOR__>1) && defined(__MINGW32__)
#  define COREARRAY_CALL_ALIGN	__attribute__((force_align_arg_pointer))
#  ifndef COREARRAY_CALL_ALIGN_NEED
#    define COREARRAY_CALL_ALIGN_NEED
#  endif
#else
#  define COREARRAY_CALL_ALIGN
#  undef COREARRAY_CALL_ALIGN_NEED
#endif




// ******************************************************************************
// Function attributes

/// fastcall keyword
#define COREARRAY_SUPPORT_FASTCALL

#if defined(COREARRAY_BORLANDC) || defined(COREARRAY_MSC)
#   define COREARRAY_FASTCALL	__fastcall
#elif defined(COREARRAY_SUN)
#   define COREARRAY_FASTCALL
#elif (defined(__GNUC__) && ((__GNUC__ > 3) || (__GNUC__ == 3 && __GNUC_MINOR__ >= 3))) && (defined(__386__) || defined(__i386__))
#   define COREARRAY_FASTCALL	__attribute__((fastcall))
#else
#   define COREARRAY_FASTCALL
#   undef COREARRAY_SUPPORT_FASTCALL
#endif


/// inline keyword
#define COREARRAY_SUPPORT_INLINE

#define COREARRAY_INLINE    inline


/// force inline keyword
#ifndef COREARRAY_NO_FORCEINLINE

#  define COREARRAY_SUPPORT_FORCE_INLINE

#  if defined(COREARRAY_BORLANDC) || defined(COREARRAY_MSC)
#     define COREARRAY_FORCE_INLINE	__forceinline
#  elif defined(COREARRAY_SUN)
#     define COREARRAY_FORCE_INLINE	COREARRAY_INLINE
#  elif (defined(__GNUC__) && ((__GNUC__ > 3) || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1)))
#     define COREARRAY_FORCE_INLINE	__attribute__((always_inline))
#  else
#     define COREARRAY_FORCE_INLINE	COREARRAY_INLINE
#     undef COREARRAY_SUPPORT_FORCE_INLINE
#  endif

#else

#  ifdef COREARRAY_SUPPORT_FORCE_INLINE
#    undef COREARRAY_SUPPORT_FORCE_INLINE
#  endif
#  define COREARRAY_FORCE_INLINE	COREARRAY_INLINE

#endif


/// dllexport keyword
#define COREARRAY_SUPPORT_DLL_EXPORT

#if defined(COREARRAY_BORLANDC) || defined(COREARRAY_MSC)
#   define COREARRAY_DLLEXPORT	__declspec(dllexport)
#elif (defined(__GNUC__) && ((__GNUC__ > 2) || (__GNUC__ == 2 && __GNUC_MINOR__ >= 8)))
#   define COREARRAY_DLLEXPORT	__attribute__((dllexport))
#else
#   define COREARRAY_DLLEXPORT
#   undef COREARRAY_SUPPORT_DLL_EXPORT
#endif






//

#define COREARRAY_DEBUG_CODE
#define COREARRAY_WARN_CODE
#define COREARRAY_HINT_CODE


#endif /* _COREDEF_H_ */
