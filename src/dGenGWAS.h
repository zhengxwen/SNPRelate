// ===========================================================
//
// dGenGWAS.h: Workspace of Genome-Wide Association Studies
//
// Copyright (C) 2011-2015    Xiuwen Zheng
//
// This file is part of SNPRelate.
//
// SNPRelate is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License Version 3 as published
// by the Free Software Foundation.
//
// SNPRelate is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with SNPRelate.
// If not, see <http://www.gnu.org/licenses/>.

#ifndef _HEADER_GWAS_
#define _HEADER_GWAS_

#include <R_GDS_CPP.h>
#include <ctime>
#include <vector>
#include <string>
#include <algorithm>
#include <memory>

#ifndef NO_COREARRAY_VECTORIZATION
#   include <dVect.h>
#endif


namespace GWAS
{
	using namespace std;
	using namespace CoreArray;


	// ===================================================================== //

	/** return the number of valid SNP genotypes
	 *  \param pGeno    the pointer to SNP genotype array
	 *  \param NumGeno  the number of SNPs
	 */
	COREARRAY_DLL_LOCAL long GENO_Get_ValidNumSNP(C_UInt8 *pGeno,
		long NumGeno);

	/** get the numbers of genotypes AA, AB, BB
	 *  \param pGeno    the pointer to SNP genotype array
	 *  \param NumGeno  the number of SNPs
	 *  \param NumAA    output the number of AA
	 *  \param NumAB    output the number of AB
	 *  \param NumBB    output the number of BB
	 */
	COREARRAY_DLL_LOCAL void GENO_Get_Num(C_UInt8 *pGeno, long NumGeno,
		long &NumAA, long &NumAB, long &NumBB);

	/** get the sum of SNP genotypes (AA -- 2, AB -- 1, BB -- 0)
	 *  \param pGeno           the pointer to SNP genotype array
	 *  \param NumGeno         the number of SNPs
	 *  \param OutValidNumSNP  output the number of SNP genotypes
	 */
	COREARRAY_DLL_LOCAL long GENO_Get_Sum_ValidNumSNP(
		C_UInt8 *pGeno, long NumGeno, long *OutValidNumSNP=NULL);



	// ===================================================================== //

	/// the base class for working space
	class COREARRAY_DLL_LOCAL CdBaseWorkSpace
	{
	public:
		enum TTypeGenoDim
		{
			RDim_Sample_X_SNP,    ///< genotype matrix: sample X snp
			RDim_SNP_X_Sample     ///< genotype matrix: snp X sample
		};

		CdBaseWorkSpace();
		virtual ~CdBaseWorkSpace();

		virtual void InitSelection() = 0;
		virtual void InitSelectionSampOnly() = 0;
		virtual void InitSelectionSNPOnly() = 0;

		virtual void snpRead(C_Int32 SnpStart, C_Int32 SnpCount,
			C_UInt8 *OutBuf, bool SnpOrder) = 0;
		virtual void sampleRead(C_Int32 SampStart, C_Int32 SampCount,
			C_UInt8 *OutBuf, bool SnpOrder) = 0;

		/// set the selection among the selected SNPs
		virtual void Set_SNPSelection(C_BOOL flag[]) = 0;
		/// set the selection among the selected samples
		virtual void Set_SampSelection(C_BOOL flag[]) = 0;

		inline TTypeGenoDim GenoDimType() const { return fGenoDimType; }
		inline C_Int32 TotalSampleNum() const { return fTotalSampleNum; };
		inline C_Int32 TotalSNPNum() const { return fTotalSNPNum; };
		inline C_Int32 SampleNum() const { return fSampleNum; };
		inline C_Int32 SNPNum() const { return fSNPNum; };
		inline C_BOOL *SampleSelection() { return &fSampleSelection[0]; };
		inline C_BOOL *SNPSelection() { return &fSNPSelection[0]; };

		/// the sum of working genotypes (the count of reference alleles)
		C_Int64 SumOfGenotype();

		void GetMissingRates(double OutRate[]);
		void GetSampValidNum(int OutNum[]);
		void GetSampMissingRates(double OutRate[]);
		void GetAlleleFreqs(double OutFreq[]);
		void GetMinorAlleleFreqs(double OutFreq[]);
		void GetABNumPerSNP(int AA[], int AB[], int BB[]);

		/** To select SNPs based on criteria, and return # of SNPs deleted
		 *  \param remove_mono    whether remove monomorphic snps or not
		 *  \param maf            the threshold of minor allele frequencies, keeping ">= maf"
		 *  \param missrate       the threshold of missing rates, keeping "<= missing.rate"
		**/
		int Select_SNP_Base(bool remove_mono, double maf, double missrate,
			C_BOOL *out_sel=NULL);

		/** To select SNPs based on criteria, and return # of SNPs deleted, specifying
		 *    the allele frequencies
		 *  \param afreq          the allele frequencies
		 *  \param remove_mono    whether remove monomorphic snps or not
		 *  \param maf            the threshold of minor allele frequencies, keeping ">= maf"
		 *  \param missrate       the threshold of missing rates, keeping "<= missing.rate"
		**/
		int Select_SNP_Base_Ex(const double afreq[], bool remove_mono,
			double maf, double missrate, C_BOOL *out_sel=NULL);

	protected:
		TTypeGenoDim fGenoDimType;
		C_Int32 fTotalSampleNum, fTotalSNPNum;
		C_Int32 fSampleNum, fSNPNum;

		vector<C_BOOL> fSampleSelection;
		vector<C_BOOL> fSNPSelection;
	};


	/// the working space for SNP genotypes in the SNP GDS file
	class COREARRAY_DLL_LOCAL CdGenoWorkSpace: public CdBaseWorkSpace
	{
	public:
		CdGenoWorkSpace();
		virtual ~CdGenoWorkSpace();

		/// set the pointer to snp genotypes
		void SetGeno(PdAbstractArray vGeno, bool _InitSelection=true);

		virtual void InitSelection();
		virtual void InitSelectionSampOnly();
		virtual void InitSelectionSNPOnly();

		virtual void snpRead(C_Int32 SnpStart, C_Int32 SnpCount,
			C_UInt8 *OutBuf, bool SnpOrder);
		virtual void sampleRead(C_Int32 SampStart, C_Int32 SampCount,
			C_UInt8 *OutBuf, bool SnpOrder);

		virtual void Set_SNPSelection(C_BOOL flag[]);
		virtual void Set_SampSelection(C_BOOL flag[]);

		void ExtractSNPs(long Start, long Length);
		void ExtractSamples(long Start, long Length);

		inline PdAbstractArray Geno() { return fGeno; };

	protected:
		PdAbstractArray fGeno;

	private:
		vector<C_Int32> vSampleIndex, vSNPIndex;
		vector<C_UInt8> vBuf;
		size_t vBufSize;

		void _NeedBuffer(size_t NewSize);
	};



	// ===================================================================== //

	/// the buffer object for SNP genotypes
	class COREARRAY_DLL_LOCAL CdBufSpace
	{
	public:
		enum TAccessFlag { acDec=0, acInc=1, acRandom=2 };

		CdBufSpace(CdBaseWorkSpace &space, bool SNPorSamp, TAccessFlag AF,
			long _bufsize=0);
		~CdBufSpace();

		C_UInt8 *ReadGeno(long idx);
		C_UInt8 *ReadPackedGeno(long idx, C_UInt8 *out_buf);
		C_UInt8 *ReadPackedGeno4b(long idx, C_UInt8 *out_buf);

		inline CdBaseWorkSpace &Space() { return *fSpace; }
		inline bool ifSNP() const { return fSNPorSamp; }
		inline bool ifSamp() const { return !fSNPorSamp; }
		inline TAccessFlag AccessFlag() const { return fAccessFlag; }
		inline void SetAccessFlag(TAccessFlag AF) { fAccessFlag = AF; }
		inline long BufSize() const { return fBufSize; }
		inline long BufElmSize() const { return fBufElmSize; }
		inline long IdxStart() const { return fIdxStart; }
		inline long IdxEnd() const { return fIdxEnd; }
		inline long IdxCnt() const { return fIdxCnt; }
	protected:
		CdBaseWorkSpace *fSpace;
		bool fSNPorSamp;
		TAccessFlag fAccessFlag;
		long fBufSize, fBufElmSize;
		C_UInt8 *_buf;
		long fIdxCnt, fIdxStart, fIdxEnd;

		void _RequireIdx(long idx);
	};



	// ===================================================================== //

	/// Four SNP genotypes are packed into one byte
	/** (s7,s6,s5,s4,s3,s2,s1,s0):
	 *    genotype 1: (s1,s0), genotype 2: (s3,s2),
	 *    genotype 3: (s5,s4), genotype 4: (s7,s6)
	**/
	COREARRAY_DLL_LOCAL C_UInt8 *PackGeno2b(const C_UInt8 *src,
		size_t cnt, C_UInt8 *dest);

	/// Two SNP genotypes are packed into one byte
	/** (s7,s6,s5,s4,s3,s2,s1,s0):
	 *    genotype 1: (s3,s2,s1,s0), s2 = s3 = 0,
	 *    genotype 2: (s7,s6,s5,s4), s6 = s7 = 0
	**/
	COREARRAY_DLL_LOCAL C_UInt8 *PackGeno4b(const C_UInt8 *src,
		size_t cnt, C_UInt8 *dest);


	// ===================================================================== //

	/// The basic class for progress object
	class COREARRAY_DLL_LOCAL CdProgression
	{
	public:
		/// The associated information
		string Info;

		/// Constructor
		CdProgression(int type=0, bool show=true);
		~CdProgression();

		void Init(C_Int64 TotalCnt, bool ShowInit=true);
		bool Forward(C_Int64 step = 1, bool Show=true);
		void ShowProgress();

        /// Return the current percentile
		inline int Percent() const { return fPercent; }
		/// Return the total number
		inline C_Int64 Total() const { return fTotal; }
		/// Return the current position
		inline C_Int64 Current() const { return fCurrent; }
		/// Whether show information
		inline bool &Show() { return fShow; }

	protected:
		int fType;
		C_Int64 fTotal, fCurrent;
		int fPercent;
		bool fShow;
		clock_t TimeInterval, OldTime;
	};


	// ===================================================================== //

	/// get a string of current time
	COREARRAY_DLL_LOCAL string NowDateToStr();



	// Matrix index

	/// Exceptions for conversion
	class COREARRAY_DLL_EXPORT ErrMatIndex: public ErrCoreArray
	{
	public:
		ErrMatIndex() {};
		ErrMatIndex(const char *fmt, ...) { _COREARRAY_ERRMACRO_(fmt); }
		ErrMatIndex(const string &msg) { fMessage = msg; }
	};


	struct COREARRAY_DLL_LOCAL IdMat
	{
	public:
		IdMat() {};

		IdMat(int n, int m); // n-by-m matrix

		IdMat & operator+=(C_Int64 val);
		IdMat & operator-=(C_Int64 val);
		IdMat & operator++ ();
		IdMat & operator-- ();
		IdMat & operator= (C_Int64 val);
		bool operator== (const IdMat &val) const;
		bool operator!= (const IdMat &val) const;

		inline int Row() const { return fOffset / fM; }
		inline int Column() const { return fOffset % fM; }
		inline C_Int64 Offset() const { return fOffset; }
	private:
		int fN, fM;
		C_Int64 fCnt, fOffset;
	};

	struct COREARRAY_DLL_LOCAL IdMatTri
	{
	public:
		IdMatTri() {};

		IdMatTri(int n); // n-by-n matrix

		IdMatTri & operator+=(int val);
		IdMatTri & operator-=(int val);
		IdMatTri & operator++ ();
		IdMatTri & operator-- ();
		IdMatTri & operator= (C_Int64 val);

		inline int Row() const { return fRow; }
		inline int Column() const { return fColumn; }
		inline C_Int64 Offset() const { return fOffset; }
	private:
		int fN, fRow, fColumn;
		C_Int64 fOffset;
	};

	struct COREARRAY_DLL_LOCAL IdMatTriD
	{
	public:
		IdMatTriD() {};

		IdMatTriD(int n); // n-by-n matrix

		IdMatTriD & operator+=(int val);
		IdMatTriD & operator-=(int val);
		IdMatTriD & operator++ ();
		IdMatTriD & operator-- ();
		IdMatTriD & operator= (int val);

		void reset(int n);
		inline int Row() const { return fRow; }
		inline int Column() const { return fColumn; }
		inline ssize_t Offset() const { return fOffset; }
	private:
		int fN, fRow, fColumn;
		ssize_t fOffset;
	};


	// matrix class

#ifndef NO_COREARRAY_VECTORIZATION

	template<typename Tx, size_t vAlign = Vectorization::_SIMD_ALIGN_>
	class COREARRAY_DLL_LOCAL CdMatTri
	{
	public:
		CdMatTri()
			{ fN = 0; }
		CdMatTri(size_t n)
			{ fN = 0; Reset(n); }
		CdMatTri(size_t n, const Tx InitVal)
			{ fN = 0; Reset(n); Clear(InitVal); }

		void Reset(size_t n)
		{
			if (n != fN)
			{
				ptr.Reset(0);
				if (n != 0) ptr.Reset(n*(n+1)/2);
				fN = n;
			}
		}
		void Clear(const Tx val)
		{
        	Tx IVAL = val, *p = ptr.get();
			for (size_t n = fN*(fN+1)/2; n > 0; n--)
				*p++ = IVAL;
		}
		void GetRow(Tx *outbuf, size_t i)
		{
			for (size_t j = 0; j < i; j++)
				outbuf[j] = ptr.get()[i + j*(2*fN-j-1)/2];
			for (size_t j = i; j < fN; j++)
				outbuf[j] = ptr.get()[j + i*(2*fN-i-1)/2];
		}
		Tx Trace()
		{
			Tx rv = 0, *p = ptr.get();
			size_t n = fN;
			while (n > 0) { rv += *p; p += n; n--; }
			return rv;
		}
		Tx Sum()
		{
			Tx rv = 0, *p = ptr.get();
			size_t n = Size();
			while (n--) rv += *p++;
			return rv;
		}
		size_t Index(size_t row, size_t col)
		{
			if (row > col)
			{ size_t t = row; row = col; col = t; }
			return col + row*(2*fN-row-1)/2;
		}

		template<typename OUTTYPE> void SaveTo(OUTTYPE *n_n_buffer)
		{
			vector<Tx> buf(fN);
			for (size_t i=0; i < fN; i++)
			{
				GetRow(&buf[0], i);
				for (size_t j=0; j < fN; j++)
					*n_n_buffer++ = (OUTTYPE)(buf[j]);
			}
		}

		inline Tx *get() { return ptr.get(); }
		inline size_t N() const { return fN; }
		inline size_t Size() const { return fN*(fN+1)/2; }

	protected:
		CoreArray::Vectorization::TdAlignPtr<Tx, vAlign> ptr;
		size_t fN;
	};


	template<typename Tx, size_t vAlign = Vectorization::_SIMD_ALIGN_>
	class COREARRAY_DLL_LOCAL CdMatTriDiag
	{
	public:
		CdMatTriDiag()
			{ fN = 0; }
		CdMatTriDiag(const Tx vDiag)
			{ fN = 0; fDiag = vDiag; }
		CdMatTriDiag(const Tx vDiag, size_t n)
			{ fN = 0; fDiag = vDiag; Reset(n); }
		CdMatTriDiag(const Tx vDiag, size_t n, const Tx InitVal)
			{ fN = 0; fDiag = vDiag; Reset(n); Clear(InitVal); }

		void Reset(size_t n)
		{
			if (n != fN)
			{
				ptr.Reset(0);
				if (n != 0) ptr.Reset((n-1)*n/2);
				fN = n;
			}
		}
		void Clear(const Tx val)
		{
			Tx IVAL = val, *p = ptr.get();
			fDiag = IVAL;
			for (size_t n = fN*(fN-1)/2; n > 0; n--)
				*p++ = IVAL;
		}
		void GetRow(Tx *outbuf, size_t i)
		{
			for (size_t j = 0; j < i; j++)
				outbuf[j] = ptr.get()[TriPtr(i-1, j, fN-1)];
			outbuf[i] = fDiag;
			for (size_t j = i+1; j < fN; j++)
				outbuf[j] = ptr.get()[TriPtr(j-1, i, fN-1)];
		}
		size_t Index(size_t row, size_t col)
		{
			if (row == col) throw "row should not be col.";
			if (row > col)
			{ size_t t = row; row = col; col = t; }
			return TriPtr(col-1, row, fN-1);
		}

		inline Tx *get() { return ptr.get(); }
		inline size_t N() const { return fN; }
		inline size_t Size() const { return fN*(fN-1)/2; }
		inline Tx &Diag() { return fDiag; }

	protected:
		CoreArray::Vectorization::TdAlignPtr<Tx, vAlign> ptr;
		size_t fN;
		Tx fDiag;
		inline size_t TriPtr(size_t i1, size_t i2, size_t n)
		{
			return i1 + i2*(2*n-i2-1)/2;
		}
	};
#endif



	// ===================================================================== //

	/// The number of SNPs in a block, optimized for memory cache
	extern long BlockNumSNP;
	/// The number of samples in a block
	extern long BlockSamp;
	/// The mutex object for the variable "Progress" and the function "RequireWork"
	extern PdThreadMutex _Mutex;
	/// The starting point of SNP, used in the function "RequireWork"
	extern long SNPStart;
	/// The starting point of samples, used in the function "RequireWorkSamp"
	extern long SampStart;

	/// The genotypes will be filled in the buffer, one genotype per byte.
	/** 0 -- BB, 1 -- AB, 2 -- AA, otherwise missing
	 *  \param buf        the output buffer
	 *  \param _SNPstart  output parameter, the starting SNP of block
	 *  \param _SNPlen    output parameter, the SNP length of block
	 *  \param SNPOrder   genotypes are stored in the SNP-major order if true,
	 *                    genotypes are stored in the sample-major order if false
	**/
	COREARRAY_DLL_LOCAL bool RequireWork(C_UInt8 *buf, long &_SNPstart,
		long &_SNPlen, bool SNPOrder);
	COREARRAY_DLL_LOCAL bool RequireWork_NoMutex(C_UInt8 *buf,
		long &_SNPstart, long &_SNPlen, bool SNPOrder);

	/// The genotypes will be filled in the buffer, one genotype per byte.
	/** 0 -- BB, 1 -- AB, 2 -- AA, otherwise missing
	 *  \param buf        the output buffer
	 *  \param _SampStart  output parameter, the starting SNP of block
	 *  \param _SampLen    output parameter, the SNP length of block
	 *  \param SNPOrder   genotypes are stored in the SNP-major order if true,
	 *                    genotypes are stored in the sample-major order if false
	**/
	COREARRAY_DLL_LOCAL bool RequireWorkSamp(C_UInt8 *buf, long &_SampStart,
		long &_SampLen, bool SNPOrder);
	COREARRAY_DLL_LOCAL bool RequireWorkSamp_NoMutex(C_UInt8 *buf,
		long &_SampStart, long &_SampLen, bool SNPOrder);



	// ===================================================================== //

	class COREARRAY_DLL_LOCAL CMultiCoreWorkingGeno
	{
	public:
		typedef void (*TDoBlockRead)(C_UInt8 *GenoBuf, long Start,
			long Cnt, void* Param);
		typedef void (*TDoEachThread)(int ThreadIndex, long Start,
			long Cnt, void* Param);

		/// The progression information
		CdProgression Progress;

		CMultiCoreWorkingGeno();
		~CMultiCoreWorkingGeno();

		/// Initialize the working space according to the SNP GDS file
		void InitSNPGDSFile(PdAbstractArray vGeno, bool _InitSelection);

		/// Initialize the parameters for multiple cores
		void InitParam(bool snp_direction, bool read_snp_order,
			long block_size);

		/// Run the user-defined functions
		void Run(int nThread, TDoBlockRead do_read, TDoEachThread do_thread,
			void *Param);

		static void SplitJobs(int nJob, int MatSize, IdMatTri outMatIdx[],
			C_Int64 outMatCnt[]);
		static void SplitJobs(int nJob, int MatSize, IdMatTriD outMatIdx[],
			C_Int64 outMatCnt[]);

		// internal uses
		void _DoThread_WorkingGeno(PdThread Thread, int ThreadIndex);

		inline CdBaseWorkSpace &Space() { return _Space; }

	protected:
		/// The working genotypes
		CdGenoWorkSpace _Space;

		/// The genotypes will be filled in the buffer, one genotype per byte.
		/** 0 -- BB, 1 -- AB, 2 -- AA, otherwise missing
		 *  \param buf        the output buffer
		 *  \param _SNPstart  output parameter, the starting SNP of block
		 *  \param _SNPlen    output parameter, the SNP length of block
		 *  \param SNPOrder   genotypes are stored in the SNP-major order if true,
		 *                    genotypes are stored in the sample-major order if false
		**/

		/// if TRUE perform computing SNP by SNP, otherwise sample by sample
		bool _SNP_Direction;
		/// if TRUE read genotypes in the major of SNP order, otherwise in the major of sample order
		bool _Read_SNP_Order;
		/// The number of SNPs or samples in a block
		long _Block_Size;
		/// The starting point of SNP or sample
		long _Start_Position;
		/// The temparory genotype buffer
		vector<C_UInt8> _Geno_Block;

		/// The mutex object
		PdThreadMutex _Mutex;
		PdThreadsSuspending _Suspend;

		// The internal parameter
		void *_Param;            /// the internal parameter
		int _Num_Thread;         /// the number of threads
		TDoBlockRead _DoRead;
		TDoEachThread _DoThread;
		int _Num_Use;
		bool _If_End;
		long _StepCnt, _StepStart;
	};

	extern CMultiCoreWorkingGeno MCWorkingGeno;

	/// Thread variables
	const int N_MAX_THREAD = 256;

	extern IdMatTri  Array_Thread_MatIdx[N_MAX_THREAD];
	extern IdMatTriD Array_Thread_MatIdxD[N_MAX_THREAD];
	extern C_Int64   Array_Thread_MatCnt[N_MAX_THREAD];



	// ===================================================================== //

	class COREARRAY_DLL_LOCAL CSummary_AvgSD
	{
	public:
		double Sum;    /// the sum: \sum_{x_i}
		double SqSum;  /// the sum: \sum_{x_i^2}
		int Num;       /// the number of valid elements
		double Avg;
		double SD;

		CSummary_AvgSD();
		void Add(double Elm);
		void Add(const double array[], size_t n);
		void CalcAvgSD();
	};



	// ===================================================================== //

	/// Detect the argument 'verbose'
	COREARRAY_DLL_LOCAL bool SEXP_Verbose(SEXP Verbose);

	///
	COREARRAY_DLL_LOCAL void CachingSNPData(const char *Msg, bool Verbose);


	/** Detect the optimized number of SNPs in a block according to
	 *    L2 and L3 cache memory, and the value is assigned to 'BlockNumSNP'
	**/
	void DetectOptimizedNumOfSNP(int nSamp, size_t atleast);


	/// The packed genotype buffer
	extern vector<C_UInt8> Array_PackedGeno;
	/// The allele frequencies
	extern vector<double> Array_AlleleFreq;
}

#endif /* _HEADER_GWAS_ */
