#include "proccess.h"
#include "main.h"

void mainProcess(cwd::kmerHashTable_t& kmerHashTable, seqan::StringSet<seqan::Dna5String>& seq)
{
	using namespace cwd;
	//每一个读数一个表，用 ReadID 作为索引，记录 readx 与 readID 之间的相同的 kmer
	auto kmerSet = findSameKmer(kmerHashTable, seq);

	// 取出表中的一行 ，放到新的表 commonKmerSet 中，然后再去除重复的 kmer
	for (size_t r = 0; r < length(seq); r++)
	{
		auto& read = kmerSet[r];
		for (size_t i = r + 1; i < length(seq); i++)
		{
			auto range = read.equal_range(i);
			if (range.first != range.second)
			{
				//cout << "seq: " << r << " and: " << i << " ";
				//cin.get();
				//cin.get();
			}
			else continue;
			auto commonKmerSet = getCommonKmerSet(range, seq[r]);
			//cout << commonKmerSet.size() << endl;
			if (commonKmerSet.size() > 0)
			{
				auto chain = chainFromStart(commonKmerSet, 15, 13, 500, 0.5, 0.5);
				if (chain.size() > 2)
				{
					outputOverlapInfo(r, i, chain, seq);
				}
			}
		}
	}
}

int main(int argc, char* argv[])
{
	using namespace std;
	using namespace seqan;
	using namespace cwd;
	//if (argc <= 1)
	//{
	//    return 0;
	//}
	string seqFileName = "/home/caiwenda/software/LROD/test/long_read.fa";
	string kfFileName = "/home/caiwenda/software/LROD/test/kmer_file.txt";
	cout << "seqFile : " << seqFileName << endl;
	cout << "frequencyFile : " << kfFileName << endl;
	auto& seq = *loadSeqData(seqFileName);
	kmerHashTable_t& kmerHashTable = createKmerHashTable(seq);
	filterKmer(kmerHashTable, kfFileName);
	mainProcess(kmerHashTable, seq);
	cout << "done!\n";
	getchar();
	return 0;
}
