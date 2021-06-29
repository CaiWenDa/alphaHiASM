#include "proccess.h"

int main(int argc, char* argv[])
{
	using namespace std;
	using namespace seqan;
	using namespace cwd;
	//if (argc <= 1)
	//{
	//    return 0;
	//}
	//string seqFileName = "/home/caiwenda/software/LROD/test/long_read.fa";
	//string kfFileName = "/home/caiwenda/software/LROD/test/kmer_file.txt";
	string seqFileName = "/publicdata/Reads/HiFi/D.mel/dmel_hifi_40x.fasta";
	// string kfFileName = "/publicdata/Reads/HiFi/D.mel/kmer31.txt";
	cout << "seqFile : " << seqFileName << endl;
	// cout << "frequencyFile : " << kfFileName << endl;
	unique_ptr<seqData_t> p_seq(loadSeqData(seqFileName));
	auto& seq = *p_seq;
	kmerHashTable_t& kmerHashTable = createKmerHashTable(seq);
	// filterKmer(kmerHashTable, kfFileName);
	mainProcess(kmerHashTable, seq);
	cout << "done!\n";
	getchar();
	return 0;
}
