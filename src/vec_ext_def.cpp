// ===========================================================
//
// vec_ext_def.cpp: scalar reference kernels + runtime dispatcher
//
// Copyright (C) 2026    Xiuwen Zheng
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
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with SNPRelate.
// If not, see <http://www.gnu.org/licenses/>.

#include "vec_ext.h"
#include <cstdlib>
#include <cstring>


namespace SNPvec
{

// ---------------------------------------------------------------------------
// Scalar reference kernels: 64-bit bitwise + popcount, identical logic to the
// original genIBS/genKING inner loops.  These are both the portable fallback
// and the correctness oracle the SIMD paths are validated against.

void ibs_count_def(const uint8_t *gi, const uint8_t *gj, size_t npack,
	uint32_t &ibs0, uint32_t &ibs1, uint32_t &ibs2)
{
	const uint8_t *p1a=gi, *p1b=gi+npack, *p2a=gj, *p2b=gj+npack;
	uint64_t c0=0, c2=0, cm=0;
	size_t m = npack;
	for (; m >= 8; m -= 8)
	{
		uint64_t g1_1,g1_2,g2_1,g2_2;
		__builtin_memcpy(&g1_1,p1a,8); __builtin_memcpy(&g1_2,p1b,8);
		__builtin_memcpy(&g2_1,p2a,8); __builtin_memcpy(&g2_2,p2b,8);
		p1a+=8; p1b+=8; p2a+=8; p2b+=8;
		const uint64_t mask = (g1_1 | ~g1_2) & (g2_1 | ~g2_2);
		const uint64_t v0 = (~((g1_1 ^ ~g2_1) | (g1_2 ^ ~g2_2))) & mask;
		const uint64_t v2 = (~((g1_1 ^  g2_1) | (g1_2 ^  g2_2))) & mask;
		c0 += __builtin_popcountll(v0);
		c2 += __builtin_popcountll(v2);
		cm += __builtin_popcountll(mask);
	}
	ibs_tail(p1a, p1b, p2a, p2b, m, c0, c2, cm);
	ibs0 += (uint32_t)c0; ibs2 += (uint32_t)c2; ibs1 += (uint32_t)(cm-c0-c2);
}

void king_robust_def(const uint8_t *gi, const uint8_t *gj, size_t npack,
	TKINGRobust &out)
{
	const uint8_t *p1a=gi, *p1b=gi+npack, *p2a=gj, *p2b=gj+npack;
	uint64_t c_ibs0=0, c_mask=0, c_het=0, c_aa1=0, c_aa2=0;
	size_t m = npack;
	for (; m >= 8; m -= 8)
	{
		uint64_t g1_1,g1_2,g2_1,g2_2;
		__builtin_memcpy(&g1_1,p1a,8); __builtin_memcpy(&g1_2,p1b,8);
		__builtin_memcpy(&g2_1,p2a,8); __builtin_memcpy(&g2_2,p2b,8);
		p1a+=8; p1b+=8; p2a+=8; p2b+=8;
		const uint64_t mask = (g1_1 | ~g1_2) & (g2_1 | ~g2_2);
		const uint64_t ibs0 = (~((g1_1 ^ ~g2_1) | (g1_2 ^ ~g2_2))) & mask;
		const uint64_t het  = ((g1_1 ^ g1_2) ^ (g2_1 ^ g2_2)) & mask;
		const uint64_t aa1  = g1_1 & ~g1_2 & mask;
		const uint64_t aa2  = g2_1 & ~g2_2 & mask;
		c_ibs0 += __builtin_popcountll(ibs0);
		c_mask += __builtin_popcountll(mask);
		c_het  += __builtin_popcountll(het);
		c_aa1  += __builtin_popcountll(aa1);
		c_aa2  += __builtin_popcountll(aa2);
	}
	king_tail(p1a, p1b, p2a, p2b, m, c_ibs0, c_mask, c_het, c_aa1, c_aa2);
	out.ibs0  += (uint32_t)c_ibs0;
	out.nloci += (uint32_t)c_mask;
	out.sumsq += (uint32_t)(c_het + 4*c_ibs0);
	out.n1_aa += (uint32_t)c_aa1;
	out.n2_aa += (uint32_t)c_aa2;
}

// Strict-IEEE reduction (4-unrolled), identical to the original MulAdd scalar
// path -- so a forced-scalar build is bit-identical to <= v1.46.
double dot_f64_def(const double *a, const double *b, size_t n)
{
	double rv = 0; size_t k = 0;
	for (; k+4 <= n; k += 4)
	{
		rv += a[k]*b[k]; rv += a[k+1]*b[k+1];
		rv += a[k+2]*b[k+2]; rv += a[k+3]*b[k+3];
	}
	for (; k < n; k++) rv += a[k]*b[k];
	return rv;
}


// ---------------------------------------------------------------------------
// Runtime dispatcher

static inline bool eq(const char *a, const char *b)
{
	return a && !std::strcmp(a, b);
}

static kernels make_def()
{
	return kernels{ &ibs_count_def, &king_robust_def, &dot_f64_def, "scalar" };
}

// pick the best ISA the running CPU supports (no name constraint)
static kernels select_auto()
{
#if defined(COREARRAY_SIMD_NEON)
	return kernels{ &ibs_count_neon, &king_robust_neon, &dot_f64_neon, "neon" };
#elif defined(SNP_VEC_X86)
	__builtin_cpu_init();
#  ifdef SNP_VEC_X86_TARGET
	if (__builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512bw"))
		return kernels{ &ibs_count_avx512, &king_robust_avx512, &dot_f64_avx512, "avx512" };
	if (__builtin_cpu_supports("avx2"))
		return kernels{ &ibs_count_avx2, &king_robust_avx2, &dot_f64_avx2, "avx2" };
#  endif
	if (__builtin_cpu_supports("sse2"))
		return kernels{ &ibs_count_sse2, &king_robust_sse2, &dot_f64_sse2, "sse2" };
	return make_def();
#else
	return make_def();
#endif
}

// pick a specific ISA by name (for SNPRELATE_FORCE_ISA / set_isa); if the named
// ISA is unavailable here, fall back to select_auto() (so the returned name will
// differ from the request and callers/tests can detect "not available").
static kernels select_named(const char *name)
{
	if (!name || !*name) return select_auto();
	if (eq(name,"def") || eq(name,"scalar")) return make_def();
#if defined(COREARRAY_SIMD_NEON)
	if (eq(name,"neon"))
		return kernels{ &ibs_count_neon, &king_robust_neon, &dot_f64_neon, "neon" };
#endif
#if defined(SNP_VEC_X86)
	__builtin_cpu_init();
	if (eq(name,"sse2") && __builtin_cpu_supports("sse2"))
		return kernels{ &ibs_count_sse2, &king_robust_sse2, &dot_f64_sse2, "sse2" };
#  ifdef SNP_VEC_X86_TARGET
	if (eq(name,"avx2") && __builtin_cpu_supports("avx2"))
		return kernels{ &ibs_count_avx2, &king_robust_avx2, &dot_f64_avx2, "avx2" };
	if (eq(name,"avx512") &&
		__builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512bw"))
		return kernels{ &ibs_count_avx512, &king_robust_avx512, &dot_f64_avx512, "avx512" };
#  endif
#endif
	return select_auto();
}

static kernels select_init()
{
	const char *force = std::getenv("SNPRELATE_FORCE_ISA");
	return (force && *force) ? select_named(force) : select_auto();
}

// The active table. Initialized once, lazily and thread-safely (C++11 magic
// static). Mutable only via set_isa() (tests), never during a computation.
static kernels& table()
{
	static kernels K = select_init();
	return K;
}

const kernels& cpu()        { return table(); }
const char* active_isa()    { return table().name; }
const char* set_isa(const char *name)
{
	table() = select_named(name);
	return table().name;
}

}  // namespace SNPvec
