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
	string seqFileName = "/publicdata/Reads/HiFi/D.mel/dmel_hifi_40x_sample_1000.fasta";
	// string kfFileName = "/publicdata/Reads/HiFi/D.mel/kmer31.txt";
	cout << "seqFile : " << seqFileName << endl;
	// cout << "frequencyFile : " << kfFileName << endl;
	seqData_t seq;
	StringSet<CharString> ID;
	loadSeqData(seqFileName, ID, seq);
	kmerHashTable_t& kmerHashTable = createKmerHashTable(seq);
	// filterKmer(kmerHashTable, kfFileName);
	mainProcess(kmerHashTable, seq, ID);
	cout << "done!\n";
	getchar();
	return 0;
}
