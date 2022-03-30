﻿#pragma once
#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "seqan/seq_io.h"
#include "seqan/align.h"
#include <unordered_map>
#include <tuple>
#include <set>
#include <list>
#include <algorithm>
#include <memory>
#include <random>
#include <cstdlib>
#include <ctime>
#include <thread>
#include <mutex>
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/properties.hpp"
#include "boost/property_map/property_map.hpp"
#include "boost/graph/graphviz.hpp"
#include "boost/format.hpp"

namespace cwd {
	using namespace std;
	const uint KMER_LEN = 31;
	const int thread_i = 24;
	template<typename T1, typename T2, typename T3 = std::hash<T1>>
	using hash = std::unordered_multimap<T1, T2, T3>;
	typedef struct {
		uint readID;
		uint begin;
		// uint end;
	} hashValue_t;
	typedef std::string kmer_t;
	struct kmer_hash {
		uint operator()(const kmer_t& kmer) const
		{
			uint hash = 0;
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
	typedef seqan::StringSet<seqan::Dna5String> seqData_t;
	// typedef std::tuple<uint, uint, uint> hashValue_t; // (readId, begin, end)
	typedef hash<kmer_t, hashValue_t> kmerHashTable_t;
	typedef std::pair<uint, uint> numPair_t;
	typedef numPair_t posPair_t;

	typedef struct {
		uint SP1;
		uint SP2;
		uint EP1;
		uint EP2;
		bool orient = true;
	} overlapInfo_t;

	typedef struct {
		bool orient = true;
		uint SP1;
		uint SP2;
	} alignInfo_t;

	typedef struct assemblyInfo_t {
		uint r1, r2;
		uint SP1;
		uint EP1;
		uint SP2;
		uint EP2;
		bool orient = true;
	} assemblyInfo_t;

	typedef struct AVertex {
		uint r;
		uint SP;
		uint EP;
		bool orient = true;

		friend ostream & operator<<(ostream& out, const AVertex& v)
		{
			out << v.r;
			return out;
		}
	} AVertex;

	typedef struct AEdge {
		enum Adj
		{
			HeadHead, HeadTail, TailHead, TailTail
		} adj;

		assemblyInfo_t ovl;
		int weight = -1;

		friend ostream& operator<<(ostream& out, const AEdge& e)
		{
			out << bitset<2>(e.adj);
			return out;
		}

	} AEdge;

	struct vertex_property_t {
		typedef boost::vertex_property_tag kind;
	};

	struct edge_property_t {
		typedef boost::edge_property_tag kind;
	};

	template <class WeightMap>
	class edge_writer
	{
	public:
		edge_writer(WeightMap w) : wm(w) {}

		template <class Edge>
		void operator()(ostream& out, const Edge& e) const {
			out << "[label=\"" << bitset<2>(wm[e]) << "\"]";
		}
	private:
		WeightMap wm;
	};

	template <class WeightMap>
	inline edge_writer<WeightMap> make_edge_writer(WeightMap w)
	{
		return edge_writer<WeightMap>(w);
	}

	using AGraph = boost::adjacency_list<boost::mapS, boost::vecS, boost::directedS, boost::property<vertex_property_t, AVertex>, AEdge>;

	kmerHashTable_t& createKmerHashTable(const seqData_t& seq);
	vector<shared_ptr<list<alignInfo_t>>> chainFromStart(seqData_t& seq, hash<kmer_t, alignInfo_t>& CKS, int k, int ks, int alpha, int beta, double gamma, int r, int t);
	vector<overlapInfo_t> finalOverlap(vector<shared_ptr<list<alignInfo_t>>>& chain, uint len1, uint len2);
	uint maxKmerFrequency(std::ifstream& kmerFrequency);
	unique_ptr<hash<uint, alignInfo_t>> findSameKmer(kmerHashTable_t& kmerHashTable, seqData_t & seq, uint r);

	template<typename T, typename R>
	unique_ptr<hash<kmer_t, alignInfo_t>> getCommonKmerSet(T range, R read)
	{

		unique_ptr<hash<kmer_t, alignInfo_t>> commonKmerSet( new hash<kmer_t, alignInfo_t>() );
		std::set<kmer_t> keySet;
		while (range.first != range.second)
		{
			uint startPos1 = (range.first->second).SP1;
			uint startPos2 = (range.first->second).SP2;
			bool orient = range.first->second.orient;
			kmer_t kmer = { seqan::begin(read) + startPos1, seqan::begin(read) + startPos1 + KMER_LEN };
			commonKmerSet->insert({ kmer, {orient, startPos1, startPos2} });
			keySet.insert(kmer);
			++range.first;
		}
		for (auto kmer : keySet)
		{
			if (commonKmerSet->count(kmer) > 1)
			{
				commonKmerSet->erase(kmer);
			}
		}
		//int num = count_if(commonKmerSet->begin(), commonKmerSet->end(), [](decltype(*commonKmerSet->begin())& a) {return a.second.orient == true;	});
		//for (auto& x : *commonKmerSet)
		//{
		//	x.second.orient = num > commonKmerSet->size() / 2 ? true : false;
		//}
		return commonKmerSet;
	}

	void loadSeqData(const std::string& seqFileName, seqan::StringSet<seqan::CharString>& ID, seqData_t& seq);
	void filterKmer(kmerHashTable_t& kmerHashTable, const std::string& kfFileName);
	void outputOverlapInfo(uint r, uint i, vector<shared_ptr<list<alignInfo_t>>>& chain_v, seqData_t& seq, seqan::StringSet<seqan::CharString> & ID, ofstream& outFile, int minSize);
	void mainProcess(kmerHashTable_t& kmerHashTable, seqData_t& seq, seqan::StringSet<seqan::CharString> & ID, int block1, int block2, ofstream& outFile);
	void assembler(const seqData_t& seq);
	kmer_t revComp(const kmer_t& kmer);
	bool findSmallerSameKmer(seqData_t& seq, uint r, uint t, uint SKMER_LEN, int s, int s2, int d, bool orient);
	std::string getCurrentDate();
	float jaccard(string& a, string& b);
	float hamming(string& a, string& b);
	void toyAssembly(seqData_t& seq);
}