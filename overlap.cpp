#include "overlap.h"
#include "utility.h"
#include <malloc.h>
#include <mutex>
#include <iostream>
#include <algorithm>
#include <random>
#include <cstdlib>
#include <fstream>
#include "boost/graph/copy.hpp"
#include "boost/graph/depth_first_search.hpp"
#include "boost/format.hpp"

using namespace std;
using namespace seqan;
using namespace cwd;
void compSeqInRange(cwd::seqData_t& seq, uint r1, uint r2, uint start1, uint start2, uint end1, uint end2, uint len, bool orient = true);

mutex fileMutex;
mutex overlapMutex;

uint KMER_STEP = 1;
int KMER_LIMIT = 51;
int thread_i = 1;
int CHAIN_LEN = 2;
double DETECT_RATIO = 0.2;
uint KMER_LEN = 31;

vector<assemblyInfo_t> overlap;
seqData_t assemblySeq;

kmerHashTable_t cwd::createKmerHashTable(const seqData_t& seq, bool isFull)
{
	kmerHashTable_t kmerHashTable;
	uint readID = 0;
	std::default_random_engine dre;
	std::uniform_int_distribution<int> di(1, KMER_LIMIT);
	cerr << "Creating HashTable...\n";
	memset(complement, 'N' , 256);
	complement['A'] = 'T';
	complement['T'] = 'A';
	complement['C'] = 'G';
	complement['G'] = 'C';

	for (auto& read : seq)
	{
		int dlen = length(read) * DETECT_RATIO;
		uint pos = 0;
		if (!isFull)
		{
			auto headEnd = begin(read) + dlen - KMER_LEN;
			auto tailBegin = end(read) - dlen - KMER_LEN;
			for (auto i = begin(read); i < headEnd; i += KMER_STEP, pos += KMER_STEP)
			{
				KMER_STEP = di(dre);
				kmer_t kmer = { i, i + KMER_LEN }; // *kmer = string(left:b, right:s)
				hashValue_t kmerInfo{ readID, pos }; //, pos + KMER_LEN };
				kmerHashTable.emplace(kmer, kmerInfo);
			}

			uint pos2 = distance(begin(read), tailBegin);
			for (auto i = tailBegin; i < end(read) - KMER_LEN; i += KMER_STEP, pos2 += KMER_STEP)
			{
				KMER_STEP = di(dre);
				kmer_t kmer = { i, i + KMER_LEN }; // *kmer = string(left:b, right:s)
				hashValue_t kmerInfo{ readID, pos2 }; //, pos + KMER_LEN };
				kmerHashTable.emplace(kmer, kmerInfo);
			}
		}
		else
		{
			for (auto i = begin(read); i < end(read); i += KMER_STEP, pos += KMER_STEP)
			{
				KMER_STEP = di(dre);
				kmer_t kmer = { i, i + KMER_LEN }; // *kmer = string(left:b, right:s)
				hashValue_t kmerInfo{ readID, pos }; //, pos + KMER_LEN };
				kmerHashTable.emplace(kmer, kmerInfo);
			}
		}
		readID++;
	}
	cerr << "HashTable has been created!\n";
	return kmerHashTable; // depends return value optimization
}

vector<unique_ptr<vector<alignInfo_t>>> cwd::chainFromStart(const seqData_t& seq, vector<alignInfo_t>& cks, int k, int ks, int alpha, int beta, double gamma, int r, int t)
{
	auto chain = make_unique<vector<alignInfo_t>>();
	chain->reserve(cks.size() / 2);
	vector<decltype(chain)> chain_v;
	chain_v.reserve(10);
	sort(cks.begin(), cks.end(), [](const auto& a, const auto& b) { return a.SP1 < b.SP1; });
	chain->push_back(*cks.begin());
	for (auto ix = cks.begin(), nextx = next(ix); ix != cks.end() && nextx != cks.end();)
	{
		uint d1 = 0;
		uint d2 = 0;
		if (nextx->SP1 - ix->SP1 < 40 || ix->SP2 < nextx->SP2 && ix->orient && nextx->orient || ix->SP2 > nextx->SP2 && !ix->orient && !nextx->orient)
		{
			if (ix->SP2 < nextx->SP2)
			{
				d1 = nextx->SP1 - ix->SP1;
				d2 = nextx->SP2 - ix->SP2;
			}
			else
			{
				d1 = nextx->SP1 - ix->SP1;
				d2 = ix->SP2 - nextx->SP2;
			}
			if ((d1 < alpha && d2 < alpha))// and (double(max(d1, d2) - min(d1, d2)) / max(d1, d2) < gamma))
			{
				chain->push_back(*nextx);
				ix = nextx;
				nextx++;
			}
			else if (d1 < beta && d2 < beta)
			{
				if ((double)d1 / d2 > 0.9 || (double)d2 / d1 > 0.9)
				{
					int s = min(ix->SP1, nextx->SP1), s2 = min(ix->SP2, nextx->SP2);
					int min_d = min(d1, d2);
					if (findSmallerSameKmer(seq, r, t, ks, s, s2, min_d, ix->orient && nextx->orient))
					{
						chain->push_back(*nextx);
						ix = nextx;
						nextx++;
					}
					else
					{
						if (chain->size() > CHAIN_LEN)
						{
							chain_v.push_back(std::move(chain));
							chain = make_unique<vector<alignInfo_t>>();
							chain->reserve(cks.size() / 2);
							chain->push_back(*nextx);
						}
						else
						{
							chain->clear();
						}
						ix = nextx;
						nextx++;
					}
				}
			}
			else
			{
				if (chain->size() > CHAIN_LEN)
				{
					chain_v.push_back(std::move(chain));
					chain = make_unique<vector<alignInfo_t>>();
					chain->reserve(cks.size() / 2);
					chain->push_back(*nextx);
				}
				else
				{
					chain->clear();
				}
				ix = nextx;
				nextx++;
				continue;
			}
		}
		else
		{
			if (chain->size() > CHAIN_LEN)
				chain_v.push_back(std::move(chain));
			chain = make_unique<vector<alignInfo_t>>();
			chain->reserve(cks.size() / 2);
			ix = nextx;
			nextx++;
			continue;
		}
	}
	if (chain->size() > CHAIN_LEN)
		chain_v.push_back(std::move(chain));
	return chain_v;
}

uint cwd::maxKmerFrequency(ifstream& kmerFrequency)
{
	vector<uint> vFrequency;
	if (kmerFrequency.is_open())
	{
		kmer_t kmer;
		uint frequency;
		uint maxFrequency = 0;
		while (kmerFrequency >> kmer >> frequency)
		{
			vFrequency.push_back(frequency);
			if (maxFrequency < frequency)
			{
				maxFrequency = frequency;
			}
		}
		auto count = new uint[maxFrequency];
		memset(count, 0, maxFrequency);
		for (auto& x : vFrequency)
		{
			count[x]++;
		}
		double s = 0;
		for (uint i = 0; i < maxFrequency; i++)
		{
			// cout << count[i] << " ";
			s += count[i];
			if (s > 0.9 * maxFrequency)
			{
				if (i > 2) return i;
				else return 3;
			}
		}
	}
	return 0;
}

vector<assemblyInfo_t> cwd::finalOverlap(const vector<unique_ptr<vector<alignInfo_t>>>& chain_v, uint len1, uint len2, uint r, uint i, int chainLen, int ovLen)
{
	vector<assemblyInfo_t> res;
	res.reserve(chain_v.size());
	auto sumLen = 0;
	for (const auto& ch : chain_v)
	{
		if (ch->size() > chainLen)
		{
			uint P1 = ch->begin()->SP1; // P1
			uint Q1 = ch->begin()->SP2;
			uint Pnk = ch->rbegin()->SP1 + KMER_LEN; // Pn + k
			uint Qnk = ch->begin()->SP2;
			for (auto it = ch->begin(); it != ch->end(); ++it) {
				if (it->SP2 < Q1) Q1 = it->SP2;
				if (it->SP2 > Qnk) Qnk = it->SP2;
			}
			Qnk += KMER_LEN;
			uint ovl_str1 = P1, ovl_str2 = Q1, ovl_end1 = Pnk, ovl_end2 = Qnk;
			auto len = max(ovl_end1 - ovl_str1, ovl_end2 - ovl_str2);
			sumLen += len;
			auto ratio = 1.0 * sumLen / min(len1, len2);
			if (len > ovLen && ratio < 0.95)
			{
				res.emplace_back(r, i, ovl_str1, ovl_str2, ovl_end1, ovl_end2, ch->begin()->orient);
			}

			else if (ratio > 0.99)
			{
				//delReads.insert(len1 > len2 ? i : r);
			}
		}
	}
	return res;
}

unique_ptr<cwd::hash<uint, alignInfo_t>> cwd::findSameKmer(const kmerHashTable_t& kmerHashTable, const seqData_t& seq, uint r)
{
	//每一个读数一个表，用 ReadID 作为索引，记录 readx 与 readID 之间的相同的 kmer
	auto kmerSet = make_unique<cwd::hash<uint, alignInfo_t>>(); //[length(seq)];
	const auto& read1 = seq[r];
	kmer_t kmer1;
	for (auto i = begin(read1); i < end(read1) - KMER_LEN; i++)
	{
		kmer1 = { i, i + KMER_LEN };
		auto rangeP = kmerHashTable.equal_range(kmer1);
		auto rangeN = kmerHashTable.equal_range(revComp(kmer1));
		// kmer_t rkmer1 = revComp(kmer1);
		auto range = rangeP;
		bool orient = true;
		// int d = distance(range.first, range.second);

		if (distance(rangeP.first, rangeP.second) < distance(rangeN.first, rangeN.second))
		{
			orient = false;
			range = rangeN;
		}
		while (range.first != range.second)
		{
			uint readID = range.first->second.readID;
			uint startPos1 = distance(begin(read1), i);
			uint startPos2 = range.first->second.begin;
			//kmer_t search = { begin(seq[readID]) + startPos2, begin(seq[readID]) + startPos2 + KMER_LEN };
			//if (search == kmer1 || search == rkmer1)
			kmerSet->emplace(readID, alignInfo_t{ orient, startPos1, startPos2 });
			//else
			//{
			//	kmerSet->insert({ readID, {orient, startPos1, startPos2 } });
			//}
			++range.first; // Increment begin iterator
		}
	}
	return kmerSet;
}

void cwd::loadSeqData(const string& seqFileName, StringSet<CharString>& ID, seqData_t& seq)
{
	StringSet<CharString> id;
	SeqFileIn seqFileIn(seqFileName.c_str());
	cerr << "Reading seqFile...\n";
	readRecords(id, seq, seqFileIn);
	cerr << "seqFile has been read.\n";
	//ofstream dict("dict_dmel_trim10.txt", ios_base::out);
	//int i = 0;
	//for (auto& str : id)
	//{
	//	StringSet<CharString> split;
	//	strSplit(split, str);
	//	//erase(split[1], 0, 3);
	//	//appendValue(ID, split[1]);
	//	dict << split[0] << " " << i++ << endl;
	//}
	//dict.close();
}

void cwd::filterKmer(kmerHashTable_t& kmerHashTable, const string& kfFileName)
{
	ifstream kmerFrequency(kfFileName);
	if (kmerFrequency.is_open())
	{
		kmer_t kmer;
		uint frequency;
		uint maxFrequency = maxKmerFrequency(kmerFrequency);
		cout << "Kmer Frequncy Range : [ " << 2 << " , " << maxFrequency << " ]\n";
		kmerFrequency.clear();
		kmerFrequency.seekg(0);
		while (kmerFrequency >> kmer >> frequency)
		{
			if (frequency < 2 || frequency > maxFrequency)
			{
				kmerHashTable.erase(kmer);
			}
		}
	}
}

void cwd::outputOverlapInfo(uint r, uint i, const vector<unique_ptr<vector<alignInfo_t>>>& chain_v, const seqData_t& seq, const StringSet<CharString> & ID, ofstream& outFile, int minSize, int chainLen, int ovLen)
{
	auto v_ovl = finalOverlap(chain_v, length(seq[r]), length(seq[i]), r, i, chainLen, ovLen);
	overlap.insert(overlap.end(), v_ovl.begin(), v_ovl.end());
}

void cwd::mainProcess(const cwd::kmerHashTable_t& kmerHashTable, const seqData_t& seq, const StringSet<CharString>& ID, uint block1, uint block2, uint seqLen, ofstream& outFile, int chainLen, int ovLen)
{
	// 取出表中的一行 ，放到新的表 commonKmerSet 中，然后再去除重复的 kmer
	vector<assemblyInfo_t> ovls;
	ovls.reserve((block2 - block1) * 10);
	for (uint r = block1; r < block2; r++)
	{
		//每一个读数一个表，用 ReadID 作为索引，记录 readx 与 readID 之间的相同的 kmer
		auto kmerSet = findSameKmer(kmerHashTable, seq, r);
		for (uint i = r + 1; i < seqLen; i++)
		{
			auto range = kmerSet->equal_range(i);
			if (range.first == range.second)
			{
				continue;
			}
			else
			{
				auto commonKmerSet = getCommonKmerSet(range, seq[r], KMER_LEN);
				//kmerSet->erase(i);
				//malloc_trim(0);
				if (!commonKmerSet.empty())
				{
					auto chain_v = chainFromStart(seq, commonKmerSet, KMER_LEN, 15, 1000, 2000, 0.2, r, i);
					if (!chain_v.empty())
					{
						auto v_ovl = finalOverlap(chain_v, length(seq[r]), length(seq[i]), r, i, chainLen, ovLen);
						ovls.insert(ovls.end(), v_ovl.begin(), v_ovl.end());
						//outputOverlapInfo(r, i, chain_v, seq, ID, outFile, 600, chainLen, ovLen);
					}
				}
			}
		}
		//seqan::clear(seq[r]);
	}
	cerr << boost::format("mapped %d - %d reads.\n") % block1 % (block2 - 1);
	lock_guard<mutex> lock(fileMutex);
	if (outFile.is_open())
	{
		for (auto& ovl : ovls)
		{
	 		if (ovl.EP1 - ovl.SP1 > ovLen && ovl.EP2 - ovl.SP2 > ovLen)
	 		{
	 			outFile << boost::format("%u, %u, %u, %u, %u, %u, %u, %u, %u\n")
	 				% ovl.r1 % ovl.r2 % ovl.orient % ovl.SP1 % ovl.EP1 % ovl.SP2 % ovl.EP2 % length(seq[ovl.r1]) % length(seq[ovl.r2]);
	 		}
		}
	
	}
	overlap.insert(overlap.end(), ovls.begin(), ovls.end());
}

kmer_t cwd::revComp(const kmer_t& kmer)
{
	// 使用查表法优化：直接在一次迭代中完成反向互补	
	kmer_t rvcKmer(kmer.rbegin(), kmer.rend());
	for (auto& dNtp : rvcKmer)
	{
		dNtp = complement[static_cast<unsigned char>(dNtp)];
	}
	// cout << rvcKmer;
	return rvcKmer;
}

bool cwd::findSmallerSameKmer(const seqData_t& seq, uint r, uint t, uint SKMER_LEN, int s, int s2, int d, bool orient)
{
	//unique_ptr<cwd::hash<uint, alignInfo_t>> kmerSet(new cwd::hash<uint, alignInfo_t>);
	auto& read1 = seq[r];
	auto& read2 = seq[t];
	
	string gap1 = { begin(read1) + s, begin(read1) + s + d };
	string gap2 = { begin(read2) + s2, begin(read2) + s2 + d };
	if (!orient)
	{
		gap2 = cwd::revComp(gap2);
	}
	return hamming(gap1, gap2) < 0.5f;
}

void cwd::readPAF(const string & file, int minOverlapLen)
{
	uint r1, r2;
	ushort table[9] = { 0 };
	char op = 0;
	int orient = 0;
	FILE* fp = fopen(file.c_str(), "r");
	if (fp)
	{
		while (fscanf(fp, "%u, %u, %d, %hu, %hu, %hu, %hu, %hu, %hu, %hu\n",
			&r1, &r2, &orient, &table[0], &table[1], &table[2], &table[3], &table[4], &table[5], &table[6]) != EOF)
			//% ovl.r1 % ovl.r2 % ovl.orient % ovl.SP1 % ovl.EP1 % ovl.SP2 % ovl.EP2 % length(seq[ovl.r1]) % length(seq[ovl.r2]);
		{
			if (table[1] - table[2] > minOverlapLen)
			{
				overlap.emplace_back(assemblyInfo_t{r1, r2, table[0], table[1], table[2], table[3], (bool)orient});
			}
		}
		
	}
	else
	{
		throw invalid_argument("invalid file");
	}
	fclose(fp);
}