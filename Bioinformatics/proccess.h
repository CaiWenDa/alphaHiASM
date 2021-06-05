#pragma once
#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <seqan/seq_io.h>
#include <seqan/align.h>
#include <unordered_map>
#include <tuple>
#include <set>
#include <list>
#include <algorithm>

namespace cwd {

	const size_t KMER_LEN = 15;
	const size_t KMER_STEP = 1;

	struct kmer_hash;
	typedef struct {
		size_t readID;
		size_t begin;
		size_t end;
	} hashValue_t;
	typedef std::string kmer_t;
	// typedef std::tuple<size_t, size_t, size_t> hashValue_t; // (readId, begin, end)
	typedef std::unordered_multimap<kmer_t, hashValue_t> kmerHashTable_t;
	typedef std::pair<size_t, size_t> numPair_t;
	typedef numPair_t posPair_t;

	typedef struct {
		size_t SP1;
		size_t SP2;
		size_t EP1;
		size_t EP2;
	} overlapInfo_t;

	struct kmer_hash {
		size_t operator()(const kmer_t& kmer) const
		{
			size_t hash = 0;
			for (auto dntp : kmer)
			{
				switch (char(dntp))
				{
				case 'A': hash <<= 2;break;
				case 'T': hash = (hash << 2) + 1;break;
				case 'C': hash = (hash << 2) + 2;break;
				case 'G': hash = (hash << 2) + 3;break;
				default:
					break;
				}
			}
			return hash;
		}
	};

	kmerHashTable_t& createKmerHashTable(const seqan::StringSet<seqan::Dna5String>& seq);
	std::list<posPair_t>& chainFromStart(std::unordered_multimap<kmer_t, posPair_t>& CKS, int k, int ks, int alpha, double beta, double gamma);
	overlapInfo_t finalOverlap(std::list<posPair_t>& chain, size_t len1, size_t len2);
	size_t maxKmerFrequency(std::ifstream& kmerFrequency);
	std::unordered_multimap<size_t, posPair_t>* findSameKmer(kmerHashTable_t& kmerHashTable, seqan::StringSet<seqan::Dna5String>& seq);

	template<typename T, typename R>
	std::unordered_multimap<kmer_t, posPair_t> getCommonKmerSet(T range, R read)
	{

		std::unordered_multimap<kmer_t, posPair_t> commonKmerSet;
		std::set<kmer_t> keySet;
		while (range.first != range.second)
		{
			size_t startPos1 = std::get<0>(range.first->second);
			size_t startPos2 = std::get<1>(range.first->second);
			kmer_t kmer = { seqan::begin(read) + startPos1, seqan::begin(read) + startPos1 + KMER_LEN };
			commonKmerSet.insert({ kmer, {startPos1, startPos2} });
			keySet.insert(kmer);
			++range.first;
		}
		for (auto kmer : keySet)
		{
			if (commonKmerSet.count(kmer) > 1)
			{
				commonKmerSet.erase(kmer);
			}
		}
		return commonKmerSet;
	}

	seqan::StringSet<seqan::Dna5String>* loadSeqData(const std::string& seqFileName);
	void filterKmer(kmerHashTable_t& kmerHashTable, const std::string& kfFileName);
	void outputOverlapInfo(const size_t& r, const size_t& i, std::__cxx11::list<cwd::posPair_t>& chain, seqan::StringSet<seqan::Dna5String>& seq);

}