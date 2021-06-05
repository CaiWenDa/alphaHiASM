#include "proccess.h"

using namespace std;
using namespace seqan;
using namespace cwd;

kmerHashTable_t& cwd::createKmerHashTable(const StringSet<Dna5String>& seq)
{
	kmerHashTable_t* kmerHashTable = new kmerHashTable_t();
	size_t readID = 0;
	for (auto& read : seq)
	{
		size_t pos = 0;
		for (auto i = begin(read); i < end(read) - KMER_LEN; i += KMER_STEP, pos += KMER_STEP)
		{
			kmer_t kmer = { i, i + KMER_LEN }; // *kmer = string(left:b, right:s)
			hashValue_t kmerInfo{ readID, pos, pos + KMER_LEN };
			kmerHashTable->insert({ kmer, kmerInfo });
		}
		readID++;
	}
	cout << "Created HashTable!\n";
	return *kmerHashTable;
}

list<posPair_t>& cwd::chainFromStart(unordered_multimap<kmer_t, posPair_t>& CKS, int k, int ks, int alpha, double beta, double gamma)
{
	list<posPair_t>* chain = new list<posPair_t>();
	chain->push_back(CKS.begin()->second);
	vector<posPair_t> cks;
	// cout << CKS.size();
	for (auto& kmer : CKS)
	{
		cks.push_back(kmer.second);
	}
	sort(cks.begin(), cks.end(), [](posPair_t a, posPair_t b) { return a.first < b.first; });
	for (auto ix = cks.begin(); ix != prev(cks.end()); ix++)
	{
		size_t d1 = next(ix)->first - ix->first;
		size_t d2 = next(ix)->second - ix->second;
		if ((ix->first < ix->second) and ((ix + 1)->first < (ix + 1)->second) and (d1 < alpha && d2 < alpha) and ((max(d1, d2) - min(d1, d2)) / max(d1, d2) < gamma))
		{
			chain->push_back(*ix);
		}
	}
	return *chain;
}

size_t cwd::maxKmerFrequency(ifstream& kmerFrequency)
{
	vector<size_t> vFrequency;
	if (kmerFrequency.is_open())
	{
		kmer_t kmer;
		size_t frequency;
		size_t maxFrequency = 0;
		while (kmerFrequency >> kmer >> frequency)
		{
			vFrequency.push_back(frequency);
			if (maxFrequency < frequency)
			{
				maxFrequency = frequency;
			}
		}
		auto count = new size_t[maxFrequency];
		memset(count, 0, maxFrequency);
		for (auto& x : vFrequency)
		{
			count[x]++;
		}
		double s = 0;
		for (size_t i = 0; i < maxFrequency; i++)
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

overlapInfo_t cwd::finalOverlap(list<posPair_t>& chain, size_t len1, size_t len2)
{
	size_t P1 = chain.begin()->first; // P1
	size_t Q1 = chain.begin()->second; // Q1
	size_t Pnk = prev(chain.end())->first + KMER_LEN; // Pn + k
	size_t Qnk = prev(chain.end())->second + KMER_LEN; // Qn + k;

	size_t ovl_str1, ovl_str2, ovl_end1, ovl_end2;

	if (P1 > Q1 && len1 - P1 <= len2 - Qnk)
	{
		ovl_str1 = P1 - Q1 - 1;
		ovl_end1 = len1 - 1;
		ovl_str2 = 0;
		ovl_end2 = Qnk + len1 - Pnk - 1;
	}
	else if (P1 < Q1 && len1 - P1 <= len2 - ovl_end2)
	{
		ovl_str1 = 0;
		ovl_end1 = len1 - 1;
		ovl_str2 = Q1 - P1 - 1;
		ovl_end2 = Qnk + len1 - Pnk - 1;
	}
	else if (P1 > Q1 && len1 - Pnk > len2 - Qnk)
	{
		ovl_str1 = P1 - Q1;
		ovl_end1 = Pnk + len2 - Qnk - 1;
		ovl_str2 = 0;
		ovl_end2 = len2 - 1;
	}
	else
	{
		ovl_str1 = 0;
		ovl_end1 = Pnk + len2 - Qnk - 1;
		ovl_str2 = Q1 - P1 - 1;
		ovl_end2 = len2 - 1;
	}
	return { ovl_str1, ovl_end1, ovl_str2, ovl_end2 };
}

unordered_multimap<size_t, posPair_t>* cwd::findSameKmer(kmerHashTable_t& kmerHashTable, StringSet<Dna5String>& seq)
{
	//每一个读数一个表，用 ReadID 作为索引，记录 readx 与 readID 之间的相同的 kmer
	auto kmerSet = new unordered_multimap<size_t, posPair_t>[length(seq)];
	for (size_t ix = 0; ix < length(seq); ix++)
	{
		auto& read1 = seq[ix];
		kmer_t kmer1;
		for (auto i = begin(read1); i < end(read1) - KMER_LEN; i += KMER_STEP)
		{
			kmer1 = { i, i + KMER_LEN };
			auto range = kmerHashTable.equal_range(kmer1);
			while (range.first != range.second)
			{
				size_t readID = range.first->second.readID;
				size_t startPos1 = distance(begin(read1), i);
				size_t startPos2 = range.first->second.begin;
				if (readID != ix)
					kmerSet[ix].insert({ readID, {startPos1, startPos2 } });
				++range.first; // Increment begin iterator
			}
		}
		//if (kmers.size() > 1)
		//{
		//    cout << "seq: " << ix << " and seq: " << iy << " ";
		//    for (auto x : kmers)
		//    {
		//        cout << x << " ";
		//    }
		//    cout << endl;
		//}

		//Align<String<Dna> > ali;
		//resize(rows(ali), 2);
		//assignSource(row(ali, 0), *ix);
		//assignSource(row(ali, 1), *iy);
		//Score<int> scoring(2, -1, -2, 0);
		//LocalAlignmentEnumerator<Score<int>, Unbanded> enumerator(scoring, 5);
		//while (nextLocalAlignment(ali, enumerator))
		//{
		//    cout << "Score = " << getScore(enumerator) << endl;
		//    // cout << ali;
		//    cout << "Aligns Seq1[" << clippedBeginPosition(row(ali, 0)) << ":" << (clippedEndPosition(row(ali, 0)) - 1) << "]";
		//    cout << " and Seq2[" << clippedBeginPosition(row(ali, 1)) << ":" << (clippedEndPosition(row(ali, 1)) - 1) << "]" << endl << endl;
		//}
	//}
	}
	return kmerSet;
}

StringSet<Dna5String>* cwd::loadSeqData(const string& seqFileName)
{
	StringSet<CharString> id;
	// StringSet<Dna5String> seq;
	auto seq = new StringSet<Dna5String>();
	SeqFileIn seqFileIn(seqFileName.c_str());
	readRecords(id, *seq, seqFileIn);
	cout << "Read seqFile\n";
	return seq;
}

void cwd::filterKmer(kmerHashTable_t& kmerHashTable, const string& kfFileName)
{
	ifstream kmerFrequency(kfFileName);
	if (kmerFrequency.is_open())
	{
		kmer_t kmer;
		size_t frequency;
		size_t maxFrequency = maxKmerFrequency(kmerFrequency);
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

void cwd::outputOverlapInfo(const size_t& r, const size_t& i, std::__cxx11::list<cwd::posPair_t>& chain, seqan::StringSet<seqan::Dna5String>& seq)
{
	cout << "seq: " << r << " and: " << i << " ";
	cout << " has OVERLAP! " << endl;
	auto ovl = finalOverlap(chain, length(seq[r]), length(seq[i]));
	cout << ovl.SP1 << " , ";
	cout << ovl.EP1 << " , ";
	cout << ovl.SP2 << " , ";
	cout << ovl.EP2 << endl;
}