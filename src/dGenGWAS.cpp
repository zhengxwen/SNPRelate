// ===========================================================
//
// dGenGWAS.cpp: Workspace of Genome-Wide Association Studies
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

#include "dGenGWAS.h"
#include <R_ext/Rdynload.h>

using namespace std;
using namespace GWAS;


// ===================================================================== //

long GWAS::GENO_Get_ValidNumSNP(C_UInt8 *pGeno, long NumGeno)
{
	long rv = 0;
	for (; NumGeno > 0; NumGeno--)
		if (*pGeno++ < 3) rv++;
	return rv;
}

void GWAS::GENO_Get_Num(C_UInt8 *pGeno, long NumGeno,
	long &NumAA, long &NumAB, long &NumBB)
{
	NumAA = NumAB = NumBB = 0;
	for (; NumGeno > 0; NumGeno--, pGeno++)
	{
		switch (*pGeno)
		{
			case 0: NumBB++; break;
			case 1: NumAB++; break;
			case 2: NumAA++; break;
		}
	}
}

long GWAS::GENO_Get_Sum_ValidNumSNP(C_UInt8 *pGeno, long NumGeno,
	long *OutValidNumSNP)
{
	long OutSum=0, ValidNumSNP=0;
	for (; NumGeno > 0; NumGeno--, pGeno++)
		if (*pGeno < 3) { OutSum+= *pGeno; ValidNumSNP++; }
	if (OutValidNumSNP)
		*OutValidNumSNP = ValidNumSNP;
	return OutSum;
}



// ===================================================================== //
// CdBaseWorkSpace

CdBaseWorkSpace::CdBaseWorkSpace()
{
	fGenoDimType = RDim_Sample_X_SNP;
	fTotalSampleNum = fTotalSNPNum = 0;
	fSampleNum = fSNPNum = 0;
}

CdBaseWorkSpace::~CdBaseWorkSpace()
{ }

void CdBaseWorkSpace::InitSelection()
{
	InitSelectionSampOnly();
	InitSelectionSNPOnly();
}

C_Int64 CdBaseWorkSpace::SumOfGenotype()
{
	C_Int64 rv = 0;

	if (fGenoDimType == RDim_Sample_X_SNP)
	{
		vector<C_UInt8> buf(fSampleNum);
		for (int i=0; i < fSNPNum; i++)
		{
			snpRead(i, 1, &buf[0], RDim_Sample_X_SNP);
			C_UInt8 *p = &buf[0];
			for (int j=fSampleNum; j > 0; j--)
			{
				if (*p <= 2)
					rv += *p;
				p ++;
			}
		}
	} else if (fGenoDimType == RDim_SNP_X_Sample)
	{
		vector<C_UInt8> buf(fSNPNum);
		for (int i=0; i < fSampleNum; i++)
		{
			sampleRead(i, 1, &buf[0], RDim_SNP_X_Sample);
			C_UInt8 *p = &buf[0];
			for (int j=fSNPNum; j > 0; j--)
			{
				if (*p <= 2)
					rv += *p;
				p ++;
			}
		}
	}

	return rv;
}

void CdBaseWorkSpace::GetMissingRates(double OutRate[])
{
	if (fGenoDimType == RDim_SNP_X_Sample)
	{
		// initialize
		for (int i=0; i < fSNPNum; i++)
			OutRate[i] = 0;
		vector<C_UInt8> buf(fSNPNum);

		// for-loop for each sample
		for (int iSamp=0; iSamp < fSampleNum; iSamp++)
		{
			sampleRead(iSamp, 1, &buf[0], RDim_SNP_X_Sample);
			for (int i=0; i < fSNPNum; i++)
			{
				if (buf[i] > 2)
					OutRate[i] ++;
			}
		}

		// average
		for (int i=0; i < fSNPNum; i++)
			OutRate[i] /= fSampleNum;

	} else {
		// initialize
		vector<C_UInt8> buf(fSampleNum);

		// for-loop for each snp
		for (int isnp=0; isnp < fSNPNum; isnp++)
		{
			double &val = OutRate[isnp];
			val = 0;
			snpRead(isnp, 1, &buf[0], RDim_Sample_X_SNP);
			for (int i=0; i < fSampleNum; i++)
			{
				if (buf[i] > 2)
					val++;
			}
			val /= fSampleNum;
		}
	}
}

void CdBaseWorkSpace::GetSampValidNum(int OutNum[])
{
	if (fGenoDimType == RDim_SNP_X_Sample)
	{
		// initialize
		vector<C_UInt8> buf(fSNPNum);

		// for-loop for each sample
		for (int iSamp=0; iSamp < fSampleNum; iSamp++)
		{
			sampleRead(iSamp, 1, &buf[0], RDim_SNP_X_Sample);
			int &cnt = OutNum[iSamp];
			cnt = 0;
			for (int i=0; i < fSNPNum; i++)
			{
				if (buf[i] <= 2)
					cnt ++;
			}
		}
	} else {
		// initialize
		vector<C_UInt8> buf(fSampleNum);
		for (int i=0; i < fSampleNum; i++)
			OutNum[i] = 0;

		// for-loop for each snp
		for (int isnp=0; isnp < fSNPNum; isnp++)
		{
			snpRead(isnp, 1, &buf[0], RDim_Sample_X_SNP);
			for (int i=0; i < fSampleNum; i++)
			{
				if (buf[i] <= 2)
					OutNum[isnp] ++;
			}
		}
	}
}

void CdBaseWorkSpace::GetSampMissingRates(double OutRate[])
{
	if (fGenoDimType == RDim_SNP_X_Sample)
	{
		// initialize
		vector<C_UInt8> buf(fSNPNum);

		// for-loop for each sample
		for (int iSamp=0; iSamp < fSampleNum; iSamp++)
		{
			sampleRead(iSamp, 1, &buf[0], RDim_SNP_X_Sample);
			double &cnt = OutRate[iSamp];
			cnt = 0;
			for (int i=0; i < fSNPNum; i++)
			{
				if (buf[i] > 2)
					cnt ++;
			}
			cnt /= fSNPNum;
		}
	} else {
		// initialize
		vector<C_UInt8> buf(fSampleNum);
		for (int i=0; i < fSampleNum; i++)
			OutRate[i] = 0;

		// for-loop for each snp
		for (int isnp=0; isnp < fSNPNum; isnp++)
		{
			snpRead(isnp, 1, &buf[0], RDim_Sample_X_SNP);
			for (int i=0; i < fSampleNum; i++)
			{
				if (buf[i] > 2)
					OutRate[i] ++;
			}
		}

		// average
		for (int i=0; i < fSampleNum; i++)
			OutRate[i] /= fSNPNum;
	}
}

void CdBaseWorkSpace::GetAlleleFreqs(double OutFreq[])
{
	if (fGenoDimType == RDim_SNP_X_Sample)
	{
		// initialize
		vector<C_UInt8> buf(fSNPNum);
		vector<int> n(fSNPNum);
		for (int i=0; i < fSNPNum; i++) n[i] = 0;
		for (int i=0; i < fSNPNum; i++) OutFreq[i] = 0;

		// for-loop for each sample
		for (int iSamp=0; iSamp < fSampleNum; iSamp++)
		{
			sampleRead(iSamp, 1, &buf[0], RDim_SNP_X_Sample);
			for (int i=0; i < fSNPNum; i++)
			{
				C_UInt8 &v = buf[i];
				if (v <= 2)
				{
					OutFreq[i] += v;
					n[i] += 2;
				}
			}
		}

		// average
		for (int i=0; i < fSNPNum; i++)
			OutFreq[i] /= n[i];

	} else {
		// initialize
		vector<C_UInt8> buf(fSampleNum);

		// for-loop for each snp
		for (int isnp=0; isnp < fSNPNum; isnp++)
		{
			int n = 0;
			double &val = OutFreq[isnp];
			val = 0;
			snpRead(isnp, 1, &buf[0], RDim_Sample_X_SNP);
			for (int i=0; i < fSampleNum; i++)
			{
				C_UInt8 &v = buf[i];
				if (v <= 2)
				{
					val += v;
					n += 2;
				}
			}
			val /= n;
		}
	}
}

void CdBaseWorkSpace::GetMinorAlleleFreqs(double OutFreq[])
{
	GetAlleleFreqs(OutFreq);
	for (int i=0; i < fSNPNum; i++)
	{
		double &freq = OutFreq[i];
		freq = min(freq, 1-freq);
	}
}

void CdBaseWorkSpace::GetABNumPerSNP(int AA[], int AB[], int BB[])
{
	// initialize the outputs
	memset(AA, 0, sizeof(int)*fSNPNum);
	memset(AB, 0, sizeof(int)*fSNPNum);
	memset(BB, 0, sizeof(int)*fSNPNum);

	if (fGenoDimType == RDim_SNP_X_Sample)
	{
		// initialize
		vector<C_UInt8> buf(fSNPNum);

		// for-loop for each sample
		for (int iSamp=0; iSamp < fSampleNum; iSamp++)
		{
			sampleRead(iSamp, 1, &buf[0], RDim_SNP_X_Sample);
			for (int i=0; i < fSNPNum; i++)
			{
				switch (buf[i])
				{
					case 0: BB[i] ++; break;
					case 1: AB[i] ++; break;
					case 2: AA[i] ++; break;
				}
			}
		}
	} else {
		// initialize
		vector<C_UInt8> buf(fSampleNum);

		// for-loop for each snp
		for (int isnp=0; isnp < fSNPNum; isnp++)
		{
			snpRead(isnp, 1, &buf[0], RDim_Sample_X_SNP);
			for (int i=0; i < fSampleNum; i++)
			{
				switch (buf[i])
				{
					case 0: BB[isnp] ++; break;
					case 1: AB[isnp] ++; break;
					case 2: AA[isnp] ++; break;
				}
			}
		}
	}
}

int CdBaseWorkSpace::Select_SNP_Base(bool remove_mono, double maf,
	double missrate, C_BOOL *out_sel)
{
	// initial variables
	vector<double> AFreq(fSNPNum);
	vector<double> MissRate(fSNPNum);

	if (fGenoDimType == RDim_SNP_X_Sample)
	{
		// initialize
		vector<C_UInt8> buf(fSNPNum);
		vector<int> n(fSNPNum);
		for (int i=0; i < fSNPNum; i++) n[i] = 0;
		for (int i=0; i < fSNPNum; i++) AFreq[i] = 0;

		// for-loop for each sample
		for (int iSamp=0; iSamp < fSampleNum; iSamp++)
		{
			sampleRead(iSamp, 1, &buf[0], RDim_SNP_X_Sample);
			for (int i=0; i < fSNPNum; i++)
			{
				C_UInt8 &v = buf[i];
				if (v <= 2)
				{
					AFreq[i] += v;
					n[i] += 2;
				}
			}
		}

		// average
		for (int i=0; i < fSNPNum; i++)
			AFreq[i] = (n[i] > 0) ? (AFreq[i]/n[i]) : R_NaN;
		for (int i=0; i < fSNPNum; i++)
			MissRate[i] = 1 - (0.5*n[i]) / fSampleNum;

	} else {
		// initialize
		vector<C_UInt8> buf(fSampleNum);

		// for-loop for each snp
		for (int isnp=0; isnp < fSNPNum; isnp++)
		{
			int n = 0;
			double &val = AFreq[isnp];
			double &miss = MissRate[isnp];
			val = 0;
			snpRead(isnp, 1, &buf[0], RDim_Sample_X_SNP);
			for (int i=0; i < fSampleNum; i++)
			{
				C_UInt8 &v = buf[i];
				if (v <= 2)
				{
					val += v;
					n += 2;
				}
			}
			val = (n > 0) ? (val/n) : R_NaN;
			miss = 1.0 - (0.5*n) / fSampleNum;
		}
	}

	// SNP selections
	vector<C_BOOL> sel(fSNPNum);
	for (int i=0; i < fSNPNum; i++)
	{
		bool flag = true;
		if (R_FINITE(AFreq[i]))
		{
			double MF = min(AFreq[i], 1-AFreq[i]);
			double MR = MissRate[i];
			if (remove_mono && (MF<=0)) flag = false;
			if (flag && (MF<maf)) flag = false;
			if (flag && (MR>missrate)) flag = false;
		} else
			flag = false;
		sel[i] = flag;
	}
	if (out_sel)
		memmove(out_sel, &sel[0], sizeof(C_BOOL)*fSNPNum);

	int cnt = 0;
	for (int i=0; i < fSNPNum; i++)
		if (!sel[i]) cnt ++;
	Set_SNPSelection(&sel[0]);

	// result
	return cnt;
}

int CdBaseWorkSpace::Select_SNP_Base_Ex(const double afreq[],
	bool remove_mono, double maf, double missrate, C_BOOL *out_sel)
{
	// initial variables
	vector<double> MissRate(fSNPNum);
	if (fGenoDimType == RDim_SNP_X_Sample)
	{
		// initialize
		vector<C_UInt8> buf(fSNPNum);
		vector<int> n(fSNPNum);
		for (int i=0; i < fSNPNum; i++) n[i] = 0;

		// for-loop for each sample
		for (int iSamp=0; iSamp < fSampleNum; iSamp++)
		{
			sampleRead(iSamp, 1, &buf[0], RDim_SNP_X_Sample);
			for (int i=0; i < fSNPNum; i++)
			{
				C_UInt8 &v = buf[i];
				if (v <= 2) n[i] ++;
			}
		}

		// average
		for (int i=0; i < fSNPNum; i++)
			MissRate[i] = 1 - double(n[i]) / fSampleNum;

	} else {
		// initialize
		vector<C_UInt8> buf(fSampleNum);

		// for-loop for each snp
		for (int isnp=0; isnp < fSNPNum; isnp++)
		{
			int n = 0;
			snpRead(isnp, 1, &buf[0], RDim_Sample_X_SNP);
			for (int i=0; i < fSampleNum; i++)
			{
				if (buf[i] <= 2)
					n ++;
			}
			MissRate[isnp] = 1.0 - double(n) / fSampleNum;
		}
	}

	// SNP selections
	vector<C_BOOL> sel(fSNPNum);
	for (int i=0; i < fSNPNum; i++)
	{
		bool flag = true;
		if (R_FINITE(afreq[i]))
		{
			double MF = min(afreq[i], 1-afreq[i]);
			double MR = MissRate[i];
			if (remove_mono && (MF<=0)) flag = false;
			if (flag && (MF<maf)) flag = false;
			if (flag && (MR>missrate)) flag = false;
		} else
			flag = false;
		sel[i] = flag;
	}
	if (out_sel)
		memmove(out_sel, &sel[0], sizeof(C_BOOL)*fSNPNum);

	int cnt = 0;
	for (int i=0; i < fSNPNum; i++)
		if (!sel[i]) cnt ++;
	Set_SNPSelection(&sel[0]);

	// result
	return cnt;
}



// ===================================================================== //
// CdSNPWorkSpace

CdSNPWorkSpace::CdSNPWorkSpace(): CdBaseWorkSpace()
{
	fGeno = NULL;
	vBufSize = 0;
}

CdSNPWorkSpace::~CdSNPWorkSpace()
{ }

void CdSNPWorkSpace::SetSNPGeno(PdAbstractArray vGeno, bool _InitSelection)
{
	if (vGeno)
	{
		// checking
		if (GDS_Array_DimCnt(vGeno) != 2)
			throw ErrCoreArray("Invalid dimension of genotype dataset.");

		// determine sample or snp order
		bool sample = (GDS_Attr_Name2Index(vGeno, "sample.order")>=0);
		bool snp = (GDS_Attr_Name2Index(vGeno, "snp.order")>=0);
		if (sample && snp)
		{
			throw ErrCoreArray(
				"Unable to determine the dimension of genotype dataset.");
		}
		if (snp)
			fGenoDimType = RDim_SNP_X_Sample;
		else if (sample)
			fGenoDimType = RDim_Sample_X_SNP;
		else
			fGenoDimType = RDim_SNP_X_Sample;

		// determine numbers of samples and snps
		C_Int32 DLen[2];
		GDS_Array_GetDim(vGeno, DLen, 2);
		if (fGenoDimType == RDim_SNP_X_Sample)
		{
			fTotalSampleNum = DLen[0]; fTotalSNPNum = DLen[1];
		} else {
			fTotalSampleNum = DLen[1]; fTotalSNPNum = DLen[0];
		}
		fSampleNum = fSNPNum = 0;

		// selection of sample
		if (fTotalSampleNum > 0)
		{
			fSampleSelection.resize(fTotalSampleNum);
			memset(&fSampleSelection[0], TRUE, fTotalSampleNum);
		} else
			fSampleSelection.clear();

		// selection of snp
		if (fTotalSNPNum > 0)
		{
			fSNPSelection.resize(fTotalSNPNum);
			memset(&fSNPSelection[0], TRUE, fTotalSNPNum);
		} else
			fSNPSelection.clear();
	} else {
		throw ErrCoreArray("'genotype' does not exist in the GDS file.");
	}

	fGeno = vGeno;
	if (_InitSelection) InitSelection();
}

void CdSNPWorkSpace::InitSelectionSampOnly()
{
	// samples
	if (fTotalSampleNum > 0)
	{
		C_BOOL *s = &fSampleSelection[0];
		fSampleNum = 0;
		for (int L=fTotalSampleNum; L > 0; L--)
			if (*s++) fSampleNum++;
		if (fSampleNum > 0)
		{
			vSampleIndex.resize(fSampleNum);
			C_Int32 *p = &vSampleIndex[0];
			s = &fSampleSelection[0];
			for (int i=0; i < fTotalSampleNum; i++)
				if (*s++) *p++ = i;
		} else {
			fSampleNum = 0;
			vSampleIndex.clear();
        }
	} else {
		fSampleNum = 0;
		vSampleIndex.clear();
	}
}

void CdSNPWorkSpace::InitSelectionSNPOnly()
{
	// snps
	if (fTotalSNPNum > 0)
	{
		C_BOOL *s = &fSNPSelection[0];
		fSNPNum = 0;
		for (int L=fTotalSNPNum; L > 0; L--)
			if (*s++) fSNPNum++;
		if (fSNPNum > 0)
		{
			vSNPIndex.resize(fSNPNum);
			C_Int32 *p = &vSNPIndex[0];
			s = &fSNPSelection[0];
			for (int i=0; i < fTotalSNPNum; i++)
				if (*s++) *p++ = i;
		} else {
			fSNPNum = 0;
			vSNPIndex.clear();
		}
	} else {
		fSNPNum = 0;
		vSNPIndex.clear();
	}
}

void CdSNPWorkSpace::snpRead(C_Int32 SnpStart, C_Int32 SnpCount,
	C_UInt8 *OutBuf, TTypeGenoDim OutDim)
{
	if ((SnpStart < 0) || (SnpStart >= fSNPNum) || (SnpCount < 0) ||
		(SnpStart+SnpCount > fSNPNum) || (fSampleNum <= 0))
	{
		throw ErrCoreArray("Invalid SnpStart and SnpCount.");
	}

	if (SnpCount > 0)
	{
		if (fGenoDimType == RDim_SNP_X_Sample)
		{
			C_Int32 st[2] =
				{ vSampleIndex[0], vSNPIndex[SnpStart] };
			C_Int32 cnt[2] =
				{ vSampleIndex[fSampleNum-1] - st[0] + 1,
				  vSNPIndex[SnpStart+SnpCount-1] - st[1] + 1 };
			C_BOOL *Sel[2] =
				{ &fSampleSelection[st[0]], &fSNPSelection[ st[1] ] };
			if ((OutDim==RDim_SNP_X_Sample) || (SnpCount==1))
			{
				GDS_Array_ReadDataEx(fGeno, st, cnt, Sel, (void*)OutBuf, svUInt8);
			} else {
				_NeedBuffer(fSampleNum*SnpCount);
				GDS_Array_ReadDataEx(fGeno, st, cnt, Sel, &vBuf[0], svUInt8);
				// transpose
				for (int i1=0; i1 < SnpCount; i1++)
				{
					for (int i0=0; i0 < fSampleNum; i0++)
						*OutBuf++ = vBuf[i0*SnpCount + i1];
                }
			}
		} else {
			C_Int32 st[2] =
				{ vSNPIndex[SnpStart], vSampleIndex[0] };
			C_Int32 cnt[2] =
				{ vSNPIndex[SnpStart+SnpCount-1]-st[0]+1,
				  vSampleIndex[fSampleNum-1]-st[1]+1 };
			C_BOOL *Sel[2] =
				{ &fSNPSelection[st[0]], &fSampleSelection[ st[1] ] };
			if ((OutDim==RDim_SNP_X_Sample) && (SnpCount>1))
			{
				_NeedBuffer(fSampleNum*SnpCount);
				GDS_Array_ReadDataEx(fGeno, st, cnt, Sel, &vBuf[0], svUInt8);
				// transpose
				for (int i1=0; i1 < fSampleNum; i1++)
				{
					for (int i0=0; i0 < SnpCount; i0++)
						*OutBuf++ = vBuf[i0*fSampleNum+i1];
                }
			} else {
				GDS_Array_ReadDataEx(fGeno, st, cnt, Sel, OutBuf, svUInt8);
            }
		}
	}
}

void CdSNPWorkSpace::sampleRead(C_Int32 SampStart, C_Int32 SampCount,
	C_UInt8 *OutBuf, TTypeGenoDim OutDim)
{
	if ((SampStart < 0) || (SampStart >= fSampleNum) || (SampCount < 0) ||
		(SampStart+SampCount > fSampleNum) || (fSNPNum <= 0))
	{
		throw ErrCoreArray("Invalid SnpStart and SnpCount.");
	}

	if (SampCount > 0)
	{
		if (fGenoDimType == RDim_SNP_X_Sample)
		{
			C_Int32 st[2] =
				{ vSampleIndex[SampStart], vSNPIndex[0] };
			C_Int32 cnt[2] =
				{ vSampleIndex[SampStart+SampCount-1] - st[0] + 1,
					vSNPIndex[fSNPNum-1] - st[1] + 1 };
			C_BOOL *Sel[2] =
				{ &fSampleSelection[st[0]], &fSNPSelection[st[1]] };
			if ((OutDim==RDim_SNP_X_Sample) || (SampCount==1))
				GDS_Array_ReadDataEx(fGeno, st, cnt, Sel, OutBuf, svUInt8);
			else {
                _NeedBuffer(SampCount*fSNPNum);
				GDS_Array_ReadDataEx(fGeno, st, cnt, Sel, &vBuf[0], svUInt8);
				// transpose
				for (int i1=0; i1 < fSNPNum; i1++)
					for (int i0=0; i0 < SampCount; i0++)
						*OutBuf++ = vBuf[i0*fSNPNum + i1];
			}
		} else {
			C_Int32 st[2] =
				{ vSNPIndex[0], vSampleIndex[SampStart] };
			C_Int32 cnt[2] =
				{ vSNPIndex[fSNPNum-1]-st[0]+1,
					vSampleIndex[SampStart+SampCount-1]-st[1]+1 };
			C_BOOL *Sel[2] =
				{ &fSNPSelection[st[0]], &fSampleSelection[st[1]] };
			if ((OutDim==RDim_SNP_X_Sample) && (SampCount>1))
			{
				_NeedBuffer(SampCount*fSNPNum);
				GDS_Array_ReadDataEx(fGeno, st, cnt, Sel, &vBuf[0], svUInt8);
				// transpose
				for (int i1=0; i1 < SampCount; i1++)
					for (int i0=0; i0 < fSNPNum; i0++)
						*OutBuf++ = vBuf[i0*SampCount + i1];
			} else
				GDS_Array_ReadDataEx(fGeno, st, cnt, Sel, OutBuf, svUInt8);
		}
	}
}

void CdSNPWorkSpace::ExtractSNPs(long Start, long Length)
{
	for (long i=0; i < Start; i++)
		fSNPSelection[vSNPIndex[i]] = FALSE;
	for (long i=Start+Length; i < fSNPNum; i++)
		fSNPSelection[vSNPIndex[i]] = FALSE;
	InitSelectionSNPOnly();
}

void CdSNPWorkSpace::ExtractSamples(long Start, long Length)
{
	for (long i=0; i < Start; i++)
		fSampleSelection[vSampleIndex[i]] = FALSE;
	for (long i=Start+Length; i < fSampleNum; i++)
		fSampleSelection[vSampleIndex[i]] = FALSE;
	InitSelectionSampOnly();
}

void CdSNPWorkSpace::Set_SNPSelection(C_BOOL flag[])
{
	for (int i=0; i < fSNPNum; i++)
		fSNPSelection[vSNPIndex[i]] = flag[i];
	InitSelectionSNPOnly();
}

void CdSNPWorkSpace::Set_SampSelection(C_BOOL flag[])
{
	for (int i=0; i < fSampleNum; i++)
		fSampleSelection[vSampleIndex[i]] = flag[i];
	InitSelectionSampOnly();
}

void CdSNPWorkSpace::_NeedBuffer(size_t NewSize)
{
	if (NewSize > vBufSize)
	{
		vBuf.resize(NewSize);
		vBufSize = NewSize;
	}
}


// ===================================================================== //
// CdSeqWorkSpace

static bool HasInitSeqArrayProc = false;

typedef void (*Type_Seq_Func)(void*);
static Type_Seq_Func fn_seq_InitSeqArray = NULL;
static Type_Seq_Func fn_seq_DoneSeqArray = NULL;

typedef void (*Type_Seq_InitSel)(C_BOOL*, void*);
static Type_Seq_InitSel fn_seq_InitSelSampOnly = NULL;
static Type_Seq_InitSel fn_seq_InitSelSNPOnly = NULL;

typedef void (*Type_Seq_Read)(C_Int32, C_Int32, C_UInt8*, TTypeGenoDim, void*);
static Type_Seq_Read fn_seq_SnpRead = NULL;
static Type_Seq_Read fn_seq_SampleRead = NULL;

typedef void (*Type_Seq_SetSel)(C_BOOL*, void*);
static Type_Seq_SetSel fn_seq_SNPSelection = NULL;
static Type_Seq_SetSel fn_seq_SampleSelection = NULL;

static void InitSeqArrayProc()
{
	if (!HasInitSeqArrayProc)
	{
		static const char *SeqArray_pkg_name = "SeqArray";

		#define LOAD(var, name)    \
			*((DL_FUNC*)&var) = R_GetCCallable(SeqArray_pkg_name, name)

		LOAD(fn_seq_InitSeqArray, "SNPRelate_InitSeqArray");
		LOAD(fn_seq_DoneSeqArray, "SNPRelate_DoneSeqArray");
		LOAD(fn_seq_InitSelSampOnly, "SNPRelate_InitSelSampOnly");
		LOAD(fn_seq_InitSelSNPOnly,  "SNPRelate_InitSelSNPOnly");
		LOAD(fn_seq_SnpRead, "SNPRelate_SnpRead");
		LOAD(fn_seq_SampleRead, "SNPRelate_SampleRead");
		LOAD(fn_seq_SNPSelection, "SNPRelate_SetSnpSelection");
		LOAD(fn_seq_SampleSelection, "SNPRelate_SetSampSelection");

		#undef LOAD
		
		HasInitSeqArrayProc = true;
	}
}


/// get the list element named str, or return NULL
inline static SEXP GetListElement(SEXP list, const char *str)
{
	SEXP elmt = R_NilValue;
	SEXP names = getAttrib(list, R_NamesSymbol);
	for (R_len_t i = 0; i < XLENGTH(list); i++)
	{
		if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0)
		{
			elmt = VECTOR_ELT(list, i);
			break;
		}
	}
	return elmt;
}

CdSeqWorkSpace::CdSeqWorkSpace(): CdBaseWorkSpace()
{
	InitSeqArrayProc();

	fParam.pGenoDimType = &fGenoDimType;
	fParam.pTotalSampleNum = &fTotalSampleNum;
	fParam.pTotalSNPNum = &fTotalSNPNum;
	fParam.pSampleNum = &fSampleNum;
	fParam.pSNPNum = &fSNPNum;

	fParam.SeqGDSFile = NULL;
	fParam.Object = NULL;
	fParam.GenoBuffer = NULL;
	fParam.Index = 0;
}

CdSeqWorkSpace::~CdSeqWorkSpace()
{
	if (fn_seq_DoneSeqArray)
	{
		if (fParam.SeqGDSFile)
		{
			(*fn_seq_DoneSeqArray)(&fParam);
			fParam.SeqGDSFile = NULL;
		}
	}
}

void CdSeqWorkSpace::SetSeqArray(SEXP SeqFile, bool _InitSelection)
{
	if (fParam.SeqGDSFile)
	{
		(*fn_seq_DoneSeqArray)(&fParam);
		fParam.SeqGDSFile = NULL;
	}
	fParam.SeqGDSFile = SeqFile;
	(*fn_seq_InitSeqArray)(&fParam);

	// selection of sample
	if (fTotalSampleNum > 0)
	{
		fSampleSelection.resize(fTotalSampleNum);
		memset(&fSampleSelection[0], TRUE, fTotalSampleNum);
	} else
		fSampleSelection.clear();

	// selection of snp
	if (fTotalSNPNum > 0)
	{
		fSNPSelection.resize(fTotalSNPNum);
		memset(&fSNPSelection[0], TRUE, fTotalSNPNum);
	} else
		fSNPSelection.clear();

	if (_InitSelection) InitSelection();
}

void CdSeqWorkSpace::InitSelectionSampOnly()
{
	(*fn_seq_InitSelSampOnly)(&fSampleSelection[0], (void*)&fParam);
}

void CdSeqWorkSpace::InitSelectionSNPOnly()
{
	(*fn_seq_InitSelSNPOnly)(&fSNPSelection[0], (void*)&fParam);
}

void CdSeqWorkSpace::snpRead(C_Int32 SnpStart, C_Int32 SnpCount,
	C_UInt8 *OutBuf, TTypeGenoDim OutDim)
{
	if ((SnpStart < 0) || (SnpStart >= fSNPNum) || (SnpCount < 0) ||
		(SnpStart+SnpCount > fSNPNum) || (fSampleNum <= 0))
	{
		throw ErrCoreArray("Invalid SnpStart and SnpCount.");
	}

	(*fn_seq_SnpRead)(SnpStart, SnpCount, OutBuf, OutDim, (void*)&fParam);
}

void CdSeqWorkSpace::sampleRead(C_Int32 SampStart, C_Int32 SampCount,
	C_UInt8 *OutBuf, TTypeGenoDim OutDim)
{
	if ((SampStart < 0) || (SampStart >= fSampleNum) || (SampCount < 0) ||
		(SampStart+SampCount > fSampleNum) || (fSNPNum <= 0))
	{
		throw ErrCoreArray("Invalid SnpStart and SnpCount.");
	}

	(*fn_seq_SampleRead)(SampStart, SampCount, OutBuf, OutDim, (void*)&fParam);
}

void CdSeqWorkSpace::Set_SNPSelection(C_BOOL flag[])
{
	(*fn_seq_SNPSelection)(flag, (void*)&fParam);
}

void CdSeqWorkSpace::Set_SampSelection(C_BOOL flag[])
{
	(*fn_seq_SampleSelection)(flag, (void*)&fParam);
}



// ===================================================================== //

// CdBufSpace

CdBufSpace::CdBufSpace(CdBaseWorkSpace &space, bool SNPorSamp, TAccessFlag AF,
	long _bufsize)
{
	fSpace = &space;
	fSNPorSamp = SNPorSamp;
	fAccessFlag = AF;
	if (_bufsize <= 0)
		fBufSize = SNPorSamp ? 128 : 32;
	else
		fBufSize = _bufsize;

	if (SNPorSamp)
	{
		// per SNP
		fBufElmSize = space.SampleNum();
		_buf = new C_UInt8[fBufSize * fBufElmSize];
		fIdxCnt = space.SNPNum();
	} else {
		// per sample
		fBufElmSize = space.SNPNum();
		_buf = new C_UInt8[fBufSize * fBufElmSize];
		fIdxCnt = space.SampleNum();
	}

	fIdxStart = fIdxEnd = 0;
}

CdBufSpace::~CdBufSpace()
{
	if (_buf)
		delete [] _buf;
	_buf = NULL;
}

C_UInt8 *CdBufSpace::ReadGeno(long idx)
{
	_RequireIdx(idx);
	return _buf + (idx-fIdxStart)*fBufElmSize;
}

C_UInt8 *CdBufSpace::ReadPackedGeno(long idx, C_UInt8 *out_buf)
{
	_RequireIdx(idx);
	return PackGeno2b(_buf + (idx-fIdxStart)*fBufElmSize, fBufElmSize, out_buf);
}

C_UInt8 *CdBufSpace::ReadPackedGeno4b(long idx, C_UInt8 *out_buf)
{
	_RequireIdx(idx);
	return PackGeno2b(_buf + (idx-fIdxStart)*fBufElmSize, fBufElmSize, out_buf);
}

void CdBufSpace::_RequireIdx(long idx)
{
	if ((idx < 0) || (idx >= fIdxCnt))
		throw ErrCoreArray("Invalid index %d in the buffer object.", idx);
	if ((idx < fIdxStart) || (idx >= fIdxEnd))
	{
		// determine the starting and ending positions
		switch (fAccessFlag)
		{
			case acDec:
				fIdxEnd = idx + 1;fIdxStart = fIdxEnd - fBufSize;
				if (fIdxStart < 0)
				{
					fIdxStart = 0;
					fIdxEnd = fBufSize;
					if (fIdxEnd > fIdxCnt) fIdxEnd = fIdxCnt;
				}
				break;
			case acInc:
				fIdxStart = idx; fIdxEnd = fIdxStart + fBufSize;
				if (fIdxEnd > fIdxCnt)
				{
					fIdxEnd = fIdxCnt;
					fIdxStart = fIdxCnt - fBufSize;
					if (fIdxStart < 0) fIdxStart = 0;
				}
				break;
			case acRandom:
				fIdxStart = idx - fBufSize/2; fIdxEnd = fIdxStart + fBufSize;
				if (fIdxStart < 0) fIdxStart = 0;
				if (fIdxEnd > fIdxCnt) fIdxEnd = fIdxCnt;
				break;
		}
		// load data
		if (fSNPorSamp)
		{
			fSpace->snpRead(fIdxStart, fIdxEnd-fIdxStart, _buf,
				RDim_Sample_X_SNP);
		} else {
			fSpace->sampleRead(fIdxStart, fIdxEnd-fIdxStart, _buf,
				RDim_SNP_X_Sample);
		}
	}
}



// ===================================================================== //

C_UInt8 *GWAS::PackGeno2b(const C_UInt8 *src, size_t cnt, C_UInt8 *dest)
{
	for (size_t n=(cnt >> 2); n > 0; n--)
	{
		*dest++ =
			(src[0] & 0x03) | ((src[1] & 0x03) << 2) |
			((src[2] & 0x03) << 4) | ((src[3] & 0x03) << 6);
		src += 4;
	}

	switch (cnt & 0x03)
	{
		case 1:
			*dest++ = (src[0] & 0x03) | 0xFC;
			break;
		case 2:
			*dest++ = (src[0] & 0x03) | ((src[1] & 0x03) << 2) | 0xF0;
			break;
		case 3:
			*dest++ = (src[0] & 0x03) | ((src[1] & 0x03) << 2) |
				((src[2] & 0x03) << 4) | 0xC0;
			break;
	}

	return dest;
}

C_UInt8 *GWAS::PackGeno4b(const C_UInt8 *src, size_t cnt, C_UInt8 *dest)
{
	for (size_t n=(cnt >> 1); n > 0; n--)
	{
		*dest++ = (src[0] & 0x03) | ((src[1] & 0x03) << 4);
		src += 2;
	}

	if (cnt & 0x01)
		*dest++ = (src[0] & 0x03) | 0x30;

	return dest;
}




// ===================================================================== //

// CdProgression

CdProgression::CdProgression(int type, bool show)
{
	fShow = show;
	Init(0, false);

	switch (fType = type)
	{
	case 0:
		TimeInterval = 30 * CLOCKS_PER_SEC;
		break;
	case 1:
		TimeInterval = 5 * CLOCKS_PER_SEC;
		break;
	}
}

CdProgression::~CdProgression()
{
	if (fType == 1)
	{
		string s(64, '=');
		Rprintf("\r%s\n", s.c_str());
	}
}

void CdProgression::Init(C_Int64 TotalCnt, bool ShowInit)
{
	if (TotalCnt < 0) TotalCnt = 0;
	fTotal = TotalCnt;
	fCurrent = fPercent = 0;
	OldTime = clock();
	if (ShowInit) ShowProgress();
}

bool CdProgression::Forward(C_Int64 step, bool Show)
{
	fCurrent += step;
	int p = int((100.0*fCurrent) / fTotal);
	if ((p != fPercent) || (p == 100))
	{
		clock_t Now = clock();
		if (((Now - OldTime) >= TimeInterval) || (p == 100))
		{
			fPercent = p;
			if (Show) ShowProgress();
			OldTime = Now;
			return true;
		}
	}
	return false;
}

void CdProgression::ShowProgress()
{
	if (fShow)
	{
		switch (fType)
		{
		case 0:
			{
				time_t tm; time(&tm);
				string s(ctime(&tm));
				s.erase(s.size()-1, 1);
				if (Info.empty())
					Rprintf("%s\t%d%%\n", s.c_str(), fPercent);
				else
					Rprintf("%s\t%s\t%d%%\n", Info.c_str(), s.c_str(), fPercent);
				break;
			}
		case 1:
			{
				int n = (int)round(fPercent * 0.64);
				string s1(n, '>');
				string s2(64-n, ' ');
				Rprintf("\r%s%s", s1.c_str(), s2.c_str());
				break;
			}
		}
	}
}



// ===================================================================== //

string GWAS::NowDateToStr()
{
	time_t tm;
	time(&tm);
	string rv(ctime(&tm));
	rv.erase(rv.size()-1, 1);
	return rv;
}




// IdMat

IdMat::IdMat(int n, int m)
{
	if ((n <= 0) || (m <= 0))
        throw ErrMatIndex("Invalid n and m: %d x %d", n, m);
	fN = n; fM = m; fOffset = 0;
    fCnt = C_Int64(n) * m;
}

IdMat & IdMat::operator+=(C_Int64 val)
{
	C_Int64 p = fOffset + val;
	if ((p < 0) || (p > fCnt))
		throw ErrMatIndex("Invalid operator += for IdMat");
	else
		fOffset = p;
	return *this;
}

IdMat & IdMat::operator-=(C_Int64 val)
{
	C_Int64 p = fOffset - val;
	if ((p < 0) || (p > fCnt))
		throw ErrMatIndex("Invalid operator -= for IdMat");
	else
		fOffset = p;
	return *this;
}

IdMat & IdMat::operator++ ()
{
	C_Int64 p = fOffset + 1;
	if (p > fCnt)
		throw ErrMatIndex("Invalid operator ++ for IdMat");
	else
		fOffset = p;
	return *this;
}

IdMat & IdMat::operator-- ()
{
	C_Int64 p = fOffset - 1;
	if (p < 0)
		throw ErrMatIndex("Invalid operator -- for IdMat");
	else
		fOffset = p;
	return *this;
}

IdMat & IdMat::operator= (C_Int64 val)
{
	if ((val < 0) || (val > fCnt))
		throw ErrMatIndex("Invalid operator = for IdMat");
	fOffset = val;
	return *this;
}

bool IdMat::operator== (const IdMat &val) const
{
    return (fN==val.fN) && (fM==val.fM) && (fOffset==val.fOffset);
}

bool IdMat::operator!= (const IdMat &val) const
{
    return (fN!=val.fN) || (fM!=val.fM) || (fOffset!=val.fOffset);
}


// IdMatTri

IdMatTri::IdMatTri(int n)
{
	if (n <= 0)
		throw ErrMatIndex("Invalid n: %d", n);
	fN = n; fRow = fColumn = fOffset = 0;
}

IdMatTri & IdMatTri::operator+= (int val)
{
	if (val > 0)
	{
		while (val > 0)
		{
			int L = val;
			if (L > (fN-fColumn)) L = fN - fColumn;
			val -= L; fOffset += L; fColumn += L;
			if (fColumn >= fN)
			{
				fColumn = ++fRow;
				if (fRow >= fN)
				{
					if (val > 0)
						throw ErrMatIndex("Invalid operator += in IdMatTri");
					break;
				}
			}
		}
	} else if (val < 0)
	{
		val = -val;
		while (val > 0)
		{
			int L = val;
			if (L > (fColumn-fRow+1)) L = fColumn-fRow+1;
			val -= L; fOffset -= L; fColumn -= L;
			if (fColumn < fRow)
			{
				--fRow; fColumn = fN-1;
				if (fRow < 0)
				{
					if (val > 0)
						throw ErrMatIndex("Invalid operator += in IdMatTri");
					break;
                }
			}
		}
	}
	return *this;
}

IdMatTri & IdMatTri::operator-= (int val)
{
	*this += (-val);
	return *this;
}

IdMatTri & IdMatTri::operator++ ()
{
	fOffset ++; fColumn ++;
	if (fColumn >= fN)
	{
		fColumn = ++fRow;
	}
	return *this;
}

IdMatTri & IdMatTri::operator-- ()
{
	fOffset --; fColumn --;
	if (fColumn < fRow)
	{
		--fRow; fColumn = fN-1;
	}
	return *this;
}

IdMatTri & IdMatTri::operator= (C_Int64 val)
{
	if ((val < 0) || (val > fN*(fN+1)/2))
		throw ErrMatIndex("Invalid operator = in IdMatTri");
	fRow = fColumn = 0;
	fOffset = val;
	while (val > 0)
	{
		if (val >= (fN-fColumn))
		{
            fColumn = ++fRow;
            val -= fN-fColumn;
		} else {
			fColumn += val;
			break;
        }
    }
	return *this;
}


// IdMatTriD

IdMatTriD::IdMatTriD(int n)
{
	fN = n; fRow = 0; fColumn = 1; fOffset = 0;
}

IdMatTriD & IdMatTriD::operator+= (int val)
{
	if (val > 0)
	{
		while (val > 0)
		{
			int L = val;
			if (L > (fN-fColumn)) L = fN-fColumn;
			val -= L; fOffset += L; fColumn += L;
			if (fColumn >= fN)
			{
				fColumn = (++fRow) + 1;
				if (fRow >= fN)
				{
					if (val > 0)
						throw ErrMatIndex("Invalid operator += in IdMatTri");
					break;
				}
			}
		}
	} else if (val < 0)
	{
		val = -val;
		while (val > 0)
		{
			int L = val;
			if (L > (fColumn-fRow)) L = fColumn-fRow;
			val -= L; fOffset -= L; fColumn -= L;
			if (fColumn <= fRow)
			{
				--fRow; fColumn = fN-1;
				if (fRow < 0)
				{
					if (val > 0)
						throw ErrMatIndex("Invalid operator += in IdMatTri");
					break;
                }
			}
		}
	}
	return *this;
}

IdMatTriD & IdMatTriD::operator-= (int val)
{
	*this += (-val);
	return *this;
}

IdMatTriD & IdMatTriD::operator++ ()
{
	*this += 1;
	return *this;
}

IdMatTriD & IdMatTriD::operator-- ()
{
	*this -= 1;
	return *this;
}

IdMatTriD & IdMatTriD::operator= (int val)
{
	if ((val < 0) || (val > fN*(fN-1)/2))
		throw ErrMatIndex("Invalid operator = in IdMatTriD");
	fRow = 0; fColumn = 1;
	fOffset = val;
	while (val > 0)
	{
		if (val >= (fN-fColumn))
		{
			val -= fN-fColumn;
			fColumn = (++fRow) + 1;
		} else {
			fColumn += val;
			break;
        }
    }
	return *this;
}

void IdMatTriD::reset(int n)
{
	if (n <= 1)
		throw ErrMatIndex("Invalid n: %d", n);
	fN = n; fRow = 0; fColumn = 1; fOffset = 0;
}



// ===================================================================== //

/// The number of SNPs in a block
long GWAS::BlockNumSNP = 256;
/// The number of samples in a block
long GWAS::BlockSamp = 32;
/// The mutex object for the variable "Progress" and the function "RequireWork"
PdThreadMutex GWAS::_Mutex = NULL;
/// The starting point of SNP, used in the function "RequireWork"
long GWAS::SNPStart = 0;
/// The starting point of sample, used in the function "RequireWorkSamp"
long GWAS::SampStart = 0;

bool GWAS::RequireWork(C_UInt8 *buf, long &_SNPstart, long &_SNPlen,
	TTypeGenoDim DimOrder)
{
	// auto Lock and Unlock
	TdAutoMutex _m(_Mutex);

	long Cnt = MCWorkingGeno.Space().SNPNum() - SNPStart;
	if (Cnt <= 0) return false;
	if (Cnt > BlockNumSNP) Cnt = BlockNumSNP;

	MCWorkingGeno.Space().snpRead(SNPStart, Cnt, buf, DimOrder);
	_SNPstart = SNPStart; _SNPlen = Cnt;
	SNPStart += Cnt;
	return true;
}

bool GWAS::RequireWork_NoMutex(C_UInt8 *buf, long &_SNPstart, long &_SNPlen,
	TTypeGenoDim DimOrder)
{
	long Cnt = MCWorkingGeno.Space().SNPNum() - SNPStart;
	if (Cnt <= 0) return false;
	if (Cnt > BlockNumSNP) Cnt = BlockNumSNP;

	MCWorkingGeno.Space().snpRead(SNPStart, Cnt, buf, DimOrder);
	_SNPstart = SNPStart; _SNPlen = Cnt;
	SNPStart += Cnt;
	return true;
}

bool GWAS::RequireWorkSamp(C_UInt8 *buf, long &_SampStart, long &_SampLen,
	TTypeGenoDim DimOrder)
{
	// auto Lock and Unlock
	TdAutoMutex _m(_Mutex);

	long Cnt = MCWorkingGeno.Space().SampleNum() - SampStart;
	if (Cnt <= 0) return false;
	if (Cnt > BlockSamp) Cnt = BlockSamp;

	MCWorkingGeno.Space().sampleRead(SampStart, Cnt, buf, DimOrder);
	_SampStart = SampStart; _SampLen = Cnt;
	SampStart += Cnt;
	return true;
}

bool GWAS::RequireWorkSamp_NoMutex(C_UInt8 *buf, long &_SampStart,
	long &_SampLen, TTypeGenoDim DimOrder)
{
	long Cnt = MCWorkingGeno.Space().SampleNum() - SampStart;
	if (Cnt <= 0) return false;
	if (Cnt > BlockSamp) Cnt = BlockSamp;

	MCWorkingGeno.Space().sampleRead(SampStart, Cnt, buf, DimOrder);
	_SampStart = SampStart; _SampLen = Cnt;
	SampStart += Cnt;
	return true;
}


// CMultiCoreWorkingGeno

CMultiCoreWorkingGeno GWAS::MCWorkingGeno;

CMultiCoreWorkingGeno::CMultiCoreWorkingGeno()
{
	_Space = NULL;

	_SNP_Direction = true;
	_Read_Order = RDim_SNP_X_Sample;
	_Block_Size = _Start_Position = 0;

	_Param = NULL;
	_Num_Thread = 0;
	_Num_Use = 0;
	_If_End = false;
	_StepCnt = _StepStart = 0;

	_Mutex = NULL;
	_Suspend = NULL;
}

CMultiCoreWorkingGeno::~CMultiCoreWorkingGeno()
{
	if (_Space)
	{
		delete _Space;
		_Space = NULL;
	}
	if (_Mutex)
		GDS_Parallel_DoneMutex(_Mutex);
	if (_Suspend)
		GDS_Parallel_DoneSuspend(_Suspend);
}

void CMultiCoreWorkingGeno::InitSNPGDSFile(PdAbstractArray vGeno,
	bool _InitSelection)
{
	if (!dynamic_cast<CdSNPWorkSpace*>(_Space))
	{
		if (_Space)
			delete _Space;
		_Space = new CdSNPWorkSpace;
	}

	static_cast<CdSNPWorkSpace*>(_Space)->SetSNPGeno(vGeno, _InitSelection);
}

void CMultiCoreWorkingGeno::InitSeqGDSFile(SEXP GDSFile, bool _InitSelection)
{
	if (!dynamic_cast<CdSeqWorkSpace*>(_Space))
	{
		if (_Space)
			delete _Space;
		_Space = new CdSeqWorkSpace;
	}

	static_cast<CdSeqWorkSpace*>(_Space)->SetSeqArray(GDSFile, _InitSelection);
}

void CMultiCoreWorkingGeno::InitParam(bool snp_direction,
	TTypeGenoDim read_order, long block_size)
{
	if (_Mutex == NULL) _Mutex = GDS_Parallel_InitMutex();
	if (_Suspend == NULL) _Suspend = GDS_Parallel_InitSuspend();

	_SNP_Direction = snp_direction;
	_Read_Order = read_order;
	_Block_Size = block_size;
	if (snp_direction)
	{
		_Geno_Block.resize(block_size * Space().SampleNum());
		Progress.Init(Space().SNPNum());
	} else {
		_Geno_Block.resize(block_size * Space().SNPNum());
		Progress.Init(Space().SampleNum());
	}

	// init the internal variables
	_Start_Position = 0;
}

static void __DoThread_WorkingGeno(PdThread Thread, int ThreadIndex,
	void* Param)
{
	CMultiCoreWorkingGeno *obj = (CMultiCoreWorkingGeno*)Param;
	obj->_DoThread_WorkingGeno(Thread, ThreadIndex);
}

void CMultiCoreWorkingGeno::Run(int nThread, TDoBlockRead do_read,
	TDoEachThread do_thread, void *Param)
{
	_Num_Thread = nThread;
	_DoRead = do_read; _DoThread = do_thread;
	_Param = Param;
	_If_End = false;
	_StepCnt = 0; _StepStart = 0;

	_Num_Use = _Num_Thread;
	GDS_Parallel_RunThreads(__DoThread_WorkingGeno, this, nThread);
}

void CMultiCoreWorkingGeno::_DoThread_WorkingGeno(PdThread Thread,
	int ThreadIndex)
{
	// initialize ...
	GDS_Parallel_LockMutex(_Mutex);
	_Num_Use --;
	GDS_Parallel_UnlockMutex(_Mutex);

	// if not main thread, suspend
	if (ThreadIndex == 0)
	{
		// main thread, wait until the other threads are all suspended
		int I;
		do {
			GDS_Parallel_LockMutex(_Mutex);
			I = _Num_Use;
			GDS_Parallel_UnlockMutex(_Mutex);
		} while (I > 0);
	} else {
		GDS_Parallel_Suspend(_Suspend);
	}

	// for loop
	while (!_If_End)
	{
		if (ThreadIndex == 0)
		{
			// check no other thread is using ...
			int I;
			do {
				GDS_Parallel_LockMutex(_Mutex);
				I = _Num_Use;
				GDS_Parallel_UnlockMutex(_Mutex);
			} while (I > 0);

			// progression information
			{
				long L = _Start_Position - _StepStart;
				if (L > 0) Progress.Forward(L);
			}

			// reading ...
			if (_SNP_Direction)
			{
				_StepCnt = Space().SNPNum() - _Start_Position;
				if (_StepCnt <= 0)
					{ _If_End = true; break; }
				if (_StepCnt > _Block_Size) _StepCnt = _Block_Size;

				Space().snpRead(_Start_Position, _StepCnt,
					&_Geno_Block[0], _Read_Order);
			} else {
				_StepCnt = Space().SampleNum() - _Start_Position;
				if (_StepCnt <= 0)
					{ _If_End = true; break; }
				if (_StepCnt > _Block_Size) _StepCnt = _Block_Size;

				Space().sampleRead(_Start_Position, _StepCnt,
					&_Geno_Block[0], _Read_Order);
			}
			_StepStart = _Start_Position;
			_Start_Position += _StepCnt;

			// handle reading ...
			_DoRead(&_Geno_Block[0], _StepStart, _StepCnt, _Param);

			// Wake up other threads
			_Num_Use = _Num_Thread;
			GDS_Parallel_WakeUp(_Suspend);

			// handle each thread
			_DoThread(ThreadIndex, _StepStart, _StepCnt, _Param);

			GDS_Parallel_LockMutex(_Mutex);
			_Num_Use --;
			GDS_Parallel_UnlockMutex(_Mutex);
		} else {
			// handle each thread
			_DoThread(ThreadIndex, _StepStart, _StepCnt, _Param);
			//
			GDS_Parallel_LockMutex(_Mutex);
			_Num_Use --;
			GDS_Parallel_UnlockMutex(_Mutex);
			// Suspend immediately
			GDS_Parallel_Suspend(_Suspend);
		}
	}

	// end ...
	// if (ThreadIndex == 0)
	GDS_Parallel_WakeUp(_Suspend);
}

void CMultiCoreWorkingGeno::SplitJobs(int nJob, int MatSize, IdMatTri outMatIdx[],
	C_Int64 outMatCnt[])
{
	#ifdef COREARRAY_NO_MULTICORE
		nJob = 1; 
	#endif

	if (nJob <= 0) nJob = 1;
	IdMatTri pt(MatSize);
	double ratio = 0.5*(MatSize+1)*MatSize / nJob, st = 0;
	C_Int64 s = 0;
	for (int i=0; i < nJob; i++)
	{
		st += ratio;
		C_Int64 p = (C_Int64)(st + 0.5);
		outMatIdx[i] = pt; outMatCnt[i] = p - s;
		pt += p - s; s = p;
	}
}

void CMultiCoreWorkingGeno::SplitJobs(int nJob, int MatSize, IdMatTriD outMatIdx[],
	C_Int64 outMatCnt[])
{
	#ifdef COREARRAY_NO_MULTICORE
		nJob = 1; 
	#endif

	if (nJob <= 0) nJob = 1;
	IdMatTriD pt(MatSize);
	double ratio = 0.5*(MatSize-1)*MatSize / nJob, st = 0;
	C_Int64 s = 0;
	for (int i=0; i < nJob; i++)
	{
		st += ratio;
		C_Int64 p = (C_Int64)(st + 0.5);
		outMatIdx[i] = pt; outMatCnt[i] = p - s;
		pt += p - s; s = p;
	}
}



bool GWAS::SEXP_Verbose(SEXP Verbose)
{
	int flag = Rf_asLogical(Verbose);
	if (flag == NA_LOGICAL)
		error("'verbose' must be TRUE or FALSE.");
	return (flag == TRUE);
}

void GWAS::CachingSNPData(const char *Msg, bool Verbose)
{
	double SumOfGenotype = MCWorkingGeno.Space().SumOfGenotype();
	if (Verbose)
	{
		Rprintf(
			"%s:\tthe sum of all working genotypes (0, 1 and 2) = %.0f\n",
				Msg, SumOfGenotype);
	}
}

void GWAS::DetectOptimizedNumOfSNP(int nSamp, size_t atleast)
{
	C_UInt64 L2Cache = GDS_Mach_GetCPULevelCache(2);
	C_UInt64 L3Cache = GDS_Mach_GetCPULevelCache(3);
	C_UInt64 Cache = (L2Cache > L3Cache) ? L2Cache : L3Cache;
	if ((C_Int64)Cache <= 0)
		Cache = 1024*1024; // 1MiB

	BlockNumSNP = (Cache - atleast - 8*1024) / nSamp * 2;
	BlockNumSNP = (BlockNumSNP / 8) * 8;
	if (BlockNumSNP < 16) BlockNumSNP = 16;
}


IdMatTri  GWAS::Array_Thread_MatIdx[N_MAX_THREAD];
IdMatTriD GWAS::Array_Thread_MatIdxD[N_MAX_THREAD];
C_Int64   GWAS::Array_Thread_MatCnt[N_MAX_THREAD];

vector<C_UInt8> GWAS::Array_PackedGeno;
vector<double>  GWAS::Array_AlleleFreq;


// =====================================================================

CSummary_AvgSD::CSummary_AvgSD()
{
	Sum = SqSum = 0;
	Num = 0;
}

void CSummary_AvgSD::Add(double Elm)
{
	if (R_FINITE(Elm))
	{
		Sum += Elm; SqSum += Elm*Elm;
		Num ++;
	}
}

void CSummary_AvgSD::Add(const double array[], size_t n)
{
	for (size_t i=0; i < n; i++)
	{
		const double v = array[i];
		if (R_FINITE(v))
		{
			Sum += v; SqSum += v*v;
			Num ++;
		}
	}
}

void CSummary_AvgSD::CalcAvgSD()
{
	if (Num > 0)
	{
		if (Num > 1)
		{
			Avg = Sum / Num;
			SD  = sqrt((SqSum - Num*Avg*Avg) / (Num - 1));
		} else {
			Avg = Sum; SD = R_NaN;
		}
	} else {
		Avg = SD = R_NaN;
	}
}
