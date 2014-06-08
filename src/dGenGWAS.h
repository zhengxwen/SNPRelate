// ===========================================================
//     _/_/_/   _/_/_/  _/_/_/_/    _/_/_/_/  _/_/_/   _/_/_/
//      _/    _/       _/             _/    _/    _/   _/   _/
//     _/    _/       _/_/_/_/       _/    _/    _/   _/_/_/
//    _/    _/       _/             _/    _/    _/   _/
// _/_/_/   _/_/_/  _/_/_/_/_/     _/     _/_/_/   _/_/
// ===========================================================
// Name        : dGenGWAS
// Author      : Xiuwen Zheng
// Version     : 1.0.0.0
// Copyright   : Xiuwen Zheng (GPL v3.0)
// Created     : 04/22/2012
// Description :
// ===========================================================

#ifndef _dGenGWAS_H_
#define _dGenGWAS_H_

#include <dType.h>
#include <CoreGDSLink.h>
#include <ctime>
#include <vector>
#include <string>
#include <algorithm>
#include <memory>

#ifndef NO_COREARRAY_VECTORIZATION
#  include <dVect.h>
#endif


namespace GWAS
{
	using namespace CoreArray;
	using namespace GDSInterface;


	// ===========================================================

	/**
	 *  return the number of valid SNP genotypes
	 *
	 *  \param pGeno    the pointer to SNP genotype array
	 *  \param NumGeno  the number of SNPs
	 */
	long GENO_Get_ValidNumSNP(UInt8 *pGeno, long NumGeno);

	/**
	 *  get the numbers of genotypes AA, AB, BB
	 *
	 *  \param pGeno    the pointer to SNP genotype array
	 *  \param NumGeno  the number of SNPs
	 *  \param NumAA    output the number of AA
	 *  \param NumAB    output the number of AB
	 *  \param NumBB    output the number of BB
	 */
	void GENO_Get_Num(UInt8 *pGeno, long NumGeno, long &NumAA, long &NumAB, long &NumBB);

	/**
	 *  get the sum of SNP genotypes (AA -- 2, AB -- 1, BB -- 0)
	 *
	 *  \param pGeno           the pointer to SNP genotype array
	 *  \param NumGeno         the number of SNPs
	 *  \param OutValidNumSNP  output the number of SNP genotypes
	 */
	long GENO_Get_Sum_ValidNumSNP(UInt8 *pGeno, long NumGeno, long *OutValidNumSNP=NULL);


	/// the working space for SNP genotypes
	class CdGenoWorkSpace
	{
	public:
		CdGenoWorkSpace();
		virtual ~CdGenoWorkSpace();

		/// set the pointer to snp genotypes
		void SetGeno(PdSequenceX vGeno, bool _InitSelection=true);
		void InitSelection();

		Int64 GenoSum();

		void snpRead(CoreArray::Int32 SnpStart, CoreArray::Int32 SnpCount,
			UInt8 *OutBuf, bool SnpOrder);
		void sampleRead(CoreArray::Int32 SampStart, CoreArray::Int32 SampCount,
			UInt8 *OutBuf, bool SnpOrder);

		void ExtractSNPs(long Start, long Length);
		void ExtractSamples(long Start, long Length);

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
		int Select_SNP_Base(bool remove_mono, double maf, double missrate, CBOOL *out_sel=NULL);

		/** To select SNPs based on criteria, and return # of SNPs deleted, specifying
		 *    the allele frequencies
		 *  \param afreq          the allele frequencies
		 *  \param remove_mono    whether remove monomorphic snps or not
		 *  \param maf            the threshold of minor allele frequencies, keeping ">= maf"
		 *  \param missrate       the threshold of missing rates, keeping "<= missing.rate"
		**/
		int Select_SNP_Base_Ex(const double afreq[], bool remove_mono, double maf,
			double missrate, CBOOL *out_sel=NULL);

		void Set_SNPSelection(CBOOL flag[]);
		void Set_SampSelection(CBOOL flag[]);

		inline PdSequenceX Geno() { return fGeno; };
		inline bool SNPOrder() const { return fSNPOrder; };
		inline CoreArray::Int32 TotalSampleNum() const { return fTotalSampleNum; };
		inline CoreArray::Int32 TotalSNPNum() const { return fTotalSNPNum; };
		inline CoreArray::Int32 SampleNum() const { return fSampleNum; };
		inline CoreArray::Int32 SNPNum() const { return fSNPNum; };
		inline CBOOL *SampleSelection() { return &fSampleSelection[0]; };
		inline CBOOL *SNPSelection() { return &fSNPSelection[0]; };

	protected:
		PdSequenceX fGeno;
		bool fSNPOrder;
		CoreArray::Int32 fTotalSampleNum, fTotalSNPNum;
		CoreArray::Int32 fSampleNum, fSNPNum;
		std::vector<CBOOL> fSampleSelection, fSNPSelection;
	private:
		std::vector<CoreArray::Int32> vSampleIndex, vSNPIndex;
		std::auto_ptr<UInt8> vBuf;
		size_t vBufSize;

		void _NeedBuffer(size_t NewSize);
		void _CheckGeno();
	};

	/// the buffer object for SNP genotypes
	class CdBufSpace
	{
	public:
		enum TAccessFlag { acDec=0, acInc=1, acRandom=2 };

		CdBufSpace(CdGenoWorkSpace &space, bool SNPorSamp, TAccessFlag AF,
			long _bufsize=0);
		~CdBufSpace();

		UInt8 * ReadGeno(long idx);
		UInt8 * ReadPackedGeno(long idx, UInt8 *out_buf);

		inline CdGenoWorkSpace &Space() { return *fSpace; }
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
		CdGenoWorkSpace *fSpace;
		bool fSNPorSamp;
		TAccessFlag fAccessFlag;
		long fBufSize, fBufElmSize;
		UInt8 *_buf;
		long fIdxCnt, fIdxStart, fIdxEnd;

		void _RequireIdx(long idx);
	};


	/// the memory object for SNP genotypes
	class CdBaseGenoMem
	{
	public:
		CdBaseGenoMem();
		CdBaseGenoMem(CdGenoWorkSpace &space);
		~CdBaseGenoMem();
	protected:
		CdGenoWorkSpace *fSpace;
		UInt8 *fMemory;
		int fRow, fColumn, fElmSize;
	};

	/// the memory object for samples
	class CdPackSampGenoMem: public CdBaseGenoMem
	{
	public:
		CdPackSampGenoMem(CdGenoWorkSpace &space);

		UInt8 at(int iSamp, int iSNP) const;
		inline int SampleNum() { return fRow; }
		inline UInt8 *PackedGeno(int iSamp) { return(fMemory + iSamp*fElmSize); }
	};

	/// the memory object for samples
	class CdSampGenoMem: public CdBaseGenoMem
	{
	public:
		CdSampGenoMem();
		CdSampGenoMem(CdGenoWorkSpace &space);

		void SetGeno(CdGenoWorkSpace &space);
		void SetGeno(int n_snp, int n_samp);

		UInt8 at(int iSamp, int iSNP) const;
		inline int SampleNum() { return fRow; }
		inline UInt8 *PtrGeno(int iSamp) { return(fMemory + iSamp*fElmSize); }
	};


	// ===========================================================

	/// four genotypes are packed into one byte
	UInt8 *PackGenotypes(const UInt8 *src, long cnt, UInt8 *dest);


	// ===========================================================

	/// The basic class for progress object
	class CdProgression
	{
	public:
		/// The associated information
		std::string Info;

		/// Constructor
		CdProgression();

		void Init(Int64 TotalCnt, bool ShowInit=true);
		bool Forward(Int64 step = 1, bool Show=true);
		void ShowProgress();

        /// Return the current percentile
		inline int Percent() const { return fPercent; }
		/// Return the total number
		inline Int64 Total() const { return fTotal; }
		/// Return the current position
		inline Int64 Current() const { return fCurrent; }
		/// Whether show information
		inline bool &Show() { return fShow; }
	protected:
		Int64 fTotal, fCurrent;
		int fPercent;
		bool fShow;
		clock_t OldTime;
	};


	// ===========================================================

	std::string NowDateToStr();



	// Matrix index

	/// Exceptions for conversion
	class ErrMatIndex: public ErrCoreArray
	{
	public:
		ErrMatIndex() {};
		ErrMatIndex(const char *fmt, ...) { _COREARRAY_ERRMACRO_(fmt); }
		ErrMatIndex(const std::string &msg) { fMessage = msg; }
	};


	struct IdMat
	{
	public:
		IdMat() {};

		IdMat(int n, int m); // n-by-m matrix

		IdMat & operator+=(Int64 val);
		IdMat & operator-=(Int64 val);
		IdMat & operator++ ();
		IdMat & operator-- ();
		IdMat & operator= (Int64 val);
		bool operator== (const IdMat &val) const;
		bool operator!= (const IdMat &val) const;

		inline int Row() const { return fOffset / fM; }
		inline int Column() const { return fOffset % fM; }
		inline Int64 Offset() const { return fOffset; }
	private:
		int fN, fM;
		Int64 fCnt, fOffset;
	};

	struct IdMatTri
	{
	public:
		IdMatTri() {};

		IdMatTri(int n); // n-by-n matrix

		IdMatTri & operator+=(int val);
		IdMatTri & operator-=(int val);
		IdMatTri & operator++ ();
		IdMatTri & operator-- ();
		IdMatTri & operator= (Int64 val);

		inline int Row() const { return fRow; }
		inline int Column() const { return fColumn; }
		inline Int64 Offset() const { return fOffset; }
	private:
		int fN, fRow, fColumn;
		Int64 fOffset;
	};

	struct IdMatTriD
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

	template<typename Tx, size_t vAlign = Vectorization::_SSE_16_Align_>
	class CdMatTri
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
			std::auto_ptr<Tx> buf(new Tx[fN]);
			for (size_t i=0; i < fN; i++)
			{
				GetRow(buf.get(), i);
				for (size_t j=0; j < fN; j++)
					*n_n_buffer++ = (OUTTYPE)(buf.get()[j]);
			}
		}

		inline Tx *get() { return ptr.get(); }
		inline size_t N() const { return fN; }
		inline size_t Size() const { return fN*(fN+1)/2; }

	protected:
		CoreArray::Vectorization::TdAlignPtr<Tx, vAlign> ptr;
		size_t fN;
	};


	template<typename Tx, size_t vAlign = Vectorization::_SSE_16_Align_>
	class CdMatTriDiag
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


	// ===========================================================

	/// The number of SNPs in a block, the number of samples in a block
	extern long BlockSNP, BlockSamp;
	/// The mutex object for the variable "Progress" and the function "RequireWork"
	extern TdMutex _Mutex;
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
	bool RequireWork(UInt8 *buf, long &_SNPstart, long &_SNPlen, bool SNPOrder);
	bool RequireWork_NoMutex(UInt8 *buf, long &_SNPstart, long &_SNPlen, bool SNPOrder);

	/// The genotypes will be filled in the buffer, one genotype per byte.
	/** 0 -- BB, 1 -- AB, 2 -- AA, otherwise missing
	 *  \param buf        the output buffer
	 *  \param _SampStart  output parameter, the starting SNP of block
	 *  \param _SampLen    output parameter, the SNP length of block
	 *  \param SNPOrder   genotypes are stored in the SNP-major order if true,
	 *                    genotypes are stored in the sample-major order if false
	**/
	bool RequireWorkSamp(UInt8 *buf, long &_SampStart, long &_SampLen, bool SNPOrder);
	bool RequireWorkSamp_NoMutex(UInt8 *buf, long &_SampStart, long &_SampLen, bool SNPOrder);


	// ===========================================================

	class CMultiCoreWorkingGeno
	{
	public:
		typedef void (*TDoBlockRead)(UInt8 *GenoBuf, long Start, long Cnt, void* Param);
		typedef void (*TDoEachThread)(int ThreadIndex, long Start, long Cnt, void* Param);

		/// The working genotypes
		CdGenoWorkSpace Space;
		/// The progression information
		CdProgression Progress;

		CMultiCoreWorkingGeno();
		~CMultiCoreWorkingGeno();

		void InitParam(bool snp_direction, bool read_snp_order, long block_size);

		void Run(int nThread, TDoBlockRead do_read, TDoEachThread do_thread, void *Param);

		static void SplitJobs(int nJob, int MatSize, IdMatTri outMatIdx[], Int64 outMatCnt[]);
		static void SplitJobs(int nJob, int MatSize, IdMatTriD outMatIdx[], Int64 outMatCnt[]);

		// internal uses
		void _DoThread_WorkingGeno(TdThread Thread, int ThreadIndex);

	protected:
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
		std::auto_ptr<UInt8> _Geno_Block;

		/// The mutex object
		TdMutex _Mutex;
		TdThreadsSuspending _Suspend;

		// The internal parameter
		void *_Param;            /// The internal parameter
		int _Num_Thread;         /// The number of threads
		TDoBlockRead _DoRead;
		TDoEachThread _DoThread;
		int _Num_Use;
		bool _If_End;
		long _StepCnt, _StepStart;
	};

	extern CMultiCoreWorkingGeno MCWorkingGeno;
}


#ifdef  __cplusplus
extern "C" {
#endif
void Rprintf(const char *, ...);
#ifdef  __cplusplus
}
#endif

#endif /* _dGenGWAS_H_ */
