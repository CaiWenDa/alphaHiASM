#include "proccess.h"

using namespace std;
using namespace seqan;
using namespace cwd;
namespace cwd {
	uint KMER_STEP = 1;
	const int KMER_LIMIT = 51;
	const int CHAIN_LEN = 2;
	const alignInfo_t& F_END = { 0, numeric_limits<uint>::max(), numeric_limits<uint>::max() };
}

kmerHashTable_t& cwd::createKmerHashTable(const seqData_t& seq)
{
	kmerHashTable_t* kmerHashTable = new kmerHashTable_t();
	uint readID = 0;
	std::default_random_engine dre;
	std::uniform_int_distribution<int> di(1, KMER_LIMIT);
	cout << "Creating HashTable...\n";
	for (auto& read : seq)
	{
		uint pos = 0;
		for (auto i = begin(read); i < end(read) - KMER_LEN; i += KMER_STEP, pos += KMER_STEP)
		{
			KMER_STEP = di(dre);
			kmer_t kmer = { i, i + KMER_LEN }; // *kmer = string(left:b, right:s)
			hashValue_t kmerInfo{ readID, pos }; //, pos + KMER_LEN };
			kmerHashTable->insert({ kmer, kmerInfo });
		}
		readID++;
	}
	cout << "HashTable has been created!\n";
	return *kmerHashTable;
}

vector<shared_ptr<list<alignInfo_t>>> cwd::chainFromStart(seqData_t& seq, hash<kmer_t, alignInfo_t>& CKS, int k, int ks, int alpha, int beta, double gamma, int r, int t)
{
	shared_ptr<list<alignInfo_t>> chain( new list<alignInfo_t>() );
	vector<decltype(chain)> chain_v;
	// chain->push_back(CKS.begin()->second);
	vector<alignInfo_t> cks;
	//auto iend = chain->begin();
	// cout << CKS.size();
	for (auto& kmer : CKS)
	{
		cks.push_back(kmer.second);
	}
	sort(cks.begin(), cks.end(), [](alignInfo_t& a, alignInfo_t& b) { return a.SP1 < b.SP1; });
	chain->push_back(*cks.begin());
	for (auto ix = cks.begin(), nextx = next(ix); ix != cks.end() && nextx != cks.end();)
	{
		uint d1 = 0;
		uint d2 = 0;
		if (ix->SP2 < nextx->SP2 && ix->orient && nextx->orient || ix->SP2 > nextx->SP2 && !ix->orient && !nextx->orient)
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
				//chain->push_back(*ix);
				chain->push_back(*nextx);
				//ix = next(nextx);
				//nextx = next(ix);
				ix = nextx;
				nextx++;
			}
			else if (d1 < beta && d2 < beta)
			{
				if ((double)d1 / d2 > 0.9 || (double)d2 / d1 > 0.9)
				{
					int s = min(ix->SP1, nextx->SP1), s2 = min(ix->SP2, nextx->SP2);
					int min_d = min(d1, d2);
					if (findSmallerSameKmer(seq, r, t, ks, s, s2, min_d,ix->orient && nextx->orient))
					{
						chain->push_back(*nextx);
						ix = nextx;
						nextx++;
					}
					else
					{
						if (chain->size() > CHAIN_LEN)
						{
							chain_v.push_back(chain);
							chain = decltype(chain)( new list<alignInfo_t> () );
							//chain->push_back(F_END);
							//iend = prev(chain->end());
							chain->push_back(*nextx);
						}
						else
							chain->erase(chain->begin(), chain->end());
						//chain->erase(next(iend), chain->end());
						ix = nextx;
						nextx++;
					}
				}
			}
			else
			{
				if (chain->size() > CHAIN_LEN)
				{
					chain_v.push_back(chain);
					chain = decltype(chain)(new list<alignInfo_t>());
					//chain->push_back(F_END);
					//iend = prev(chain->end());
					chain->push_back(*nextx);
				}
				else
					//chain->erase(next(iend), chain->end());
					chain->erase(chain->begin(), chain->end());
				ix = nextx;
				nextx++;
				continue;
			}
		}
		else
		{
			if (chain->size() > CHAIN_LEN)
				chain_v.push_back(chain);
			chain = decltype(chain)(new list<alignInfo_t>());
			//chain->push_back(F_END);//
			//iend = prev(chain->end());
			ix = nextx;
			nextx++;
			continue;
		}
	}
	//chain->push_back(F_END);
	if (chain->size() > CHAIN_LEN)
		chain_v.push_back(chain);
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

vector<overlapInfo_t> cwd::finalOverlap(vector<shared_ptr<list<alignInfo_t>>>& chain_v, uint len1, uint len2)
{
	//auto r1 = chain.begin();
	//auto r2 = find_if(r1, chain.end(), [](alignInfo_t& a) { return a.orient == false && a.SP1 == numeric_limits<uint>::max() && a.SP2 == numeric_limits<uint>::max(); });
	vector<overlapInfo_t> res;
	for (auto& ch : chain_v)
	{
		if (ch->size() > 2)
		{
			uint P1 = ch->begin()->SP1; // P1
			uint Q1 = min_element(ch->begin(), ch->end(), [](alignInfo_t& a, alignInfo_t& b) {return a.SP2 < b.SP2;})->SP2;//chain.begin()->SP2; // Q1
			uint Pnk = ch->rbegin()->SP1 + KMER_LEN; // Pn + k
			uint Qnk = max_element(ch->begin(), ch->end(), [](alignInfo_t& a, alignInfo_t& b) {return a.SP2 < b.SP2;})->SP2;//prev(chain.end())->SP2 + KMER_LEN; // Qn + k;

			uint ovl_str1, ovl_str2, ovl_end1, ovl_end2;
			ovl_str1 = P1, ovl_str2 = Q1, ovl_end1 = Pnk, ovl_end2 = Qnk + KMER_LEN;
			res.push_back({ ovl_str1, ovl_str2, ovl_end1, ovl_end2 , ch->begin()->orient });
		}
		//r1 = next(r2);
		//r2 = find_if(r1, chain_v.end(), [](alignInfo_t a) { return a.orient == false && a.SP1 == F_END.SP1 && a.SP2 == F_END.SP2; });
	}
	
	//uint P1 = chain.begin()->SP1; // P1
	//uint Q1 = min_element(chain.begin(), chain.end(), [](alignInfo_t& a, alignInfo_t& b) {return a.SP2 < b.SP2;})->SP2;//chain.begin()->SP2; // Q1
	//uint Pnk = prev(chain.end())->SP1 + KMER_LEN; // Pn + k
	//uint Qnk = max_element(chain.begin(), chain.end(), [](alignInfo_t& a, alignInfo_t& b) {return a.SP2 < b.SP2;})->SP2;//prev(chain.end())->SP2 + KMER_LEN; // Qn + k;

	//uint ovl_str1, ovl_str2, ovl_end1, ovl_end2;
	//ovl_str1 = P1, ovl_str2 = Q1, ovl_end1 = Pnk, ovl_end2 = Qnk + KMER_LEN;
/*	if (P1 > Q1 && len1 - Pnk <= len2 - Qnk)
	{
		ovl_str1 = P1 - Q1;// -1;
		ovl_end1 = len1 - 1;
		ovl_str2 = 0;
		ovl_end2 = Qnk + len1 - Pnk;// -1;
	}
	else if (P1 <= Q1 && len1 - Pnk <= len2 - Qnk)
	{
		ovl_str1 = 0;
		ovl_end1 = len1 - 1;
		ovl_str2 = Q1 - P1;// -1;
		ovl_end2 = Qnk + len1 - Pnk;// -1;
	}
	else if (P1 > Q1 && len1 - Pnk > len2 - Qnk)
	{
		ovl_str1 = P1 - Q1;
		ovl_end1 = Pnk + len2 - Qnk;// -1;
		ovl_str2 = 0;
		ovl_end2 = len2 - 1;
	}
	else if (P1 <= Q1 && len1 - Pnk >= len2 - Qnk)
	{
		ovl_str1 = 0;
		ovl_end1 = Pnk + len2 - Qnk;// - 1;
		ovl_str2 = Q1 - P1;// -1;
		ovl_end2 = len2 - 1;
	}
*/
	return res;
}

unique_ptr<cwd::hash<uint, alignInfo_t>> cwd::findSameKmer(kmerHashTable_t& kmerHashTable, seqData_t & seq, uint r)
{
	//每一个读数一个表，用 ReadID 作为索引，记录 readx 与 readID 之间的相同的 kmer
	unique_ptr<cwd::hash<uint, alignInfo_t>> kmerSet( new cwd::hash<uint, alignInfo_t> ); //[length(seq)];
	auto& read1 = seq[r];
	kmer_t kmer1;
	for (auto i = begin(read1); i < end(read1) - KMER_LEN; i++)
	{
		kmer1 = { i, i + KMER_LEN };
		auto rangeP = kmerHashTable.equal_range(kmer1);
		auto rangeN = kmerHashTable.equal_range(revComp(kmer1));
		kmer_t rkmer1 = revComp(kmer1);
		auto range = rangeP;
		bool orient = true;
		int d = distance(range.first, range.second);
		//{
//			std::default_random_engine dre;
//			std::uniform_int_distribution<int> di(0, 30);
//			int x = di(dre);
//			kmer1[x] = revComp(kmer1.c_str() + x)[0];
//			rangeP = kmerHashTable.equal_range(kmer1);
//			rangeN = kmerHashTable.equal_range(revComp(kmer1));
//			range = rangeP;
		//}
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
			kmerSet->insert({ readID, {orient, startPos1, startPos2 } });
			//else
			//{
			//	kmerSet->insert({ readID, {orient, startPos1, startPos2 } });
			//}
			++range.first; // Increment begin iterator
		}
	}
	return kmerSet;
}

seqData_t* cwd::loadSeqData(const string& seqFileName, StringSet<CharString>& ID, seqData_t& seq)
{
	StringSet<CharString> id;
	//seqData_t* seq = new StringSet<Dna5String>();
	SeqFileIn seqFileIn(seqFileName.c_str());
	cout << "Reading seqFile...\n";
	readRecords(id, seq, seqFileIn);
	cout << "seqFile has been read.\n";
	//ofstream dict("dict3.txt", ios_base::out);
	//int i = 0;
	//for (auto& str : id)
	//{
	//	StringSet<CharString> split;
	//	strSplit(split, str);
	//	//erase(split[1], 0, 3);
	//	//appendValue(ID, split[1]);
	//	dict << split[0] << " " << i++ << endl;
	//}
	return &seq;
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

void cwd::outputOverlapInfo(uint r, uint & i, vector<shared_ptr<list<alignInfo_t>>>& chain_v, seqData_t& seq, StringSet<CharString> & ID, ofstream& outFile, int minSize)
{
	//cout << "seq: " << r << " and: " << i << " ";
	//cout << " has OVERLAP! " << endl;
	auto v_ovl = finalOverlap(chain_v, length(seq[r]), length(seq[i]));
	//cout << ovl.SP1 << " , ";
	//cout << ovl.EP1 << " , ";
	//cout << ovl.SP2 << " , ";
	//cout << ovl.EP2 << endl;
	for (auto& ovl : v_ovl)
	{
		if (ovl.EP1 - ovl.SP1 > minSize && ovl.EP2 - ovl.SP2 > minSize)
		{
			outFile << r << "," << i << ",";
			outFile << ovl.orient << ",";
			outFile << ovl.SP1 << ",";
			outFile << ovl.EP1 << ",";
			outFile << ovl.SP2 << ",";
			outFile << ovl.EP2 << "," << length(seq[r]) << "," << length(seq[i]) << endl;
		}
	}
}

void cwd::mainProcess(cwd::kmerHashTable_t& kmerHashTable, seqData_t& seq, StringSet<CharString> & ID)
{
	// 取出表中的一行 ，放到新的表 commonKmerSet 中，然后再去除重复的 kmer
	ofstream outFile("result_sampled-" + getCurrentDate() + ".csv", ios_base::out);
	for (uint r = 0; r < length(seq); r++)
	{
		//每一个读数一个表，用 ReadID 作为索引，记录 readx 与 readID 之间的相同的 kmer
		auto kmerSet = findSameKmer(kmerHashTable, seq, r);
		for (uint i = r + 1; i < length(seq); i++)
		{
			//if (r == 6423 && i == 6859)
			//	cout << string{ begin(seq[r]), begin(seq[r]) + 5347 } << "\n" << string{ begin(seq[r]), begin(seq[r]) + 5347 } << endl;
			//else continue;
			auto range = kmerSet->equal_range(i);
			if (range.first == range.second) continue;
			auto commonKmerSet = getCommonKmerSet(range, seq[r]);
			//cout << commonKmerSet.size() << endl;
			if (commonKmerSet->size() > 0)
			{
				auto chain_v = chainFromStart(seq, *commonKmerSet, KMER_LEN, 15, 300, 500, 0.2, r, i);
				if (chain_v.size() > 0)
				{
					outputOverlapInfo(r, i, chain_v, seq, ID, outFile, 600);
				}
			}
		}
		seqan::clear(seq[r]);
	}
}

kmer_t cwd::revComp(const kmer_t& kmer)
{
	kmer_t rvcKmer = kmer;
	for (auto& dNtp : rvcKmer)
	{
		switch (dNtp)
		{
		case 'A': dNtp = 'T'; break;
		case 'T': dNtp = 'A'; break;
		case 'C': dNtp = 'G'; break;
		case 'G': dNtp = 'C'; break;
		default:
			break;
		}
	}
	std::reverse(rvcKmer.begin(), rvcKmer.end());
	// cout << rvcKmer;
	return rvcKmer;
}

bool cwd::findSmallerSameKmer(seqData_t& seq, uint r, uint t, uint SKMER_LEN, int s, int s2, int d, bool orient)
{
	//unique_ptr<cwd::hash<uint, alignInfo_t>> kmerSet(new cwd::hash<uint, alignInfo_t>);
	auto& read1 = seq[r];
	auto& read2 = seq[t];
	
	string gap1 = { begin(read1) + s, begin(read1) + s + d };
	string gap2 = { begin(read2) + s2, begin(read2) + s2 + d };
	/*kmer_t kmer1, kmer2;
	int count = 0;
	for (auto i = begin(read1) + s; i < begin(read1) + e - SKMER_LEN; i++)
	{
		for (auto j = begin(read2) + s; j < begin(read2) + e - SKMER_LEN;)
		{
			int a = distance(begin(read2), j), b = distance(begin(read1), i);
			if (a > b)
			{
				t = a - b;
				if (t >= SKMER_LEN)
				{
					break;
				}
			}
			else
			{
				t = b - a;
				if (t >= SKMER_LEN)
				{
					j += t;
				}
			}
			kmer1 = { i, i + SKMER_LEN };
			kmer2 = { j, j + SKMER_LEN };
			if (kmer1 == kmer2)
			{
				count++;
				j += SKMER_LEN;
				i += SKMER_LEN;
			}
			else
			{
				j++;
			}
		}
	}
	if (count > 1)
	{
		float score = float(count) / (float(e - s ) / SKMER_LEN);
		if (score > 0.5)
		{
			return true;
		}
	}*/
	return hamming(gap1, gap2) < 0.5f;
}

std::string cwd::getCurrentDate()
{
	time_t nowtime;
	nowtime = time(NULL); //获取日历时间   
	char tmp[64];
	strftime(tmp, sizeof(tmp), "%Y-%m-%d-%H-%M", localtime(&nowtime));
	return tmp;
}

float cwd::jaccard(string& a, string& b)
{
	if (a == "" && b == "")
	{
		return 1.0f;
	}
	// 都为空相似度为 1
	if (a == "" || b == "")
	{
		return 0.0f;
	}
	std::set<int> aChar(a.begin(), a.end());
	std::set<int> bChar(b.begin(), b.end());
	std::set<int> in;
	std::set<int> un;
	set_intersection(aChar.begin(), aChar.end(), bChar.begin(), bChar.end(), inserter(in, in.begin()));
	set_union(aChar.begin(), aChar.end(), bChar.begin(), bChar.end(), inserter(un, un.begin()));
	int inter = in.size(); // 交集数量
	int uni = un.size(); // 并集数量
	if (inter == 0) return 0;
	return ((float)inter) / (float)uni;
}

float cwd::hamming(string& a, string& b)
{
	if (a == "" || b == "")
	{
		return 0.0f;
	}
	if (a.length() != b.length())
	{
		return 0.0f;
	}

	int disCount = 0;
	for (int i = 0; i < a.length(); i++)
	{
		if (a.at(i) != b.at(i))
		{
			disCount++;
		}
	}
	return (float)disCount / (float)a.length();
}