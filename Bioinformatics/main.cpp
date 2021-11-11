#include "proccess.h"

void compSeqInRange(cwd::seqData_t seq, uint r1, uint r2, uint start1, uint start2, uint end1, uint end2, uint len, bool orient = true);

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
	string seqFileName = "/publicdata_two/publicdata/publicdata/Reads/HiFi/D.mel/dmel_hifi_40x_sample.fasta";
	// string kfFileName = "/publicdata/Reads/HiFi/D.mel/kmer31.txt";
	cout << "seqFile : " << seqFileName << endl;
	// cout << "frequencyFile : " << kfFileName << endl;
	seqData_t seq;
	StringSet<CharString> ID;
	loadSeqData(seqFileName, ID, seq);
	kmerHashTable_t& kmerHashTable = createKmerHashTable(seq);
	// filterKmer(kmerHashTable, kfFileName);
	//cout << string{ begin(seq[6859]) + 5221, begin(seq[6859]) + 5221 +200 } << endl;
	//cout << string{ begin(seq[1712]) + 21497, begin(seq[1712]) + 21497 + 200 } << endl;
	//compSeqInRange(seq, 25, 9552, 5353, 12059, 8734, 15759, 300);
	mainProcess(kmerHashTable, seq, ID);
	cout << "done!\n";
	getchar();
	return 0;
}

void compSeqInRange(cwd::seqData_t seq, uint r1, uint r2, uint start1, uint start2, uint end1, uint end2, uint len, bool orient)
{
	using namespace std;
	if (!orient)
	{
		cout << string{ begin(seq[r1]) + start1, begin(seq[r1]) + start1 + len } << endl;
		cout << cwd::revComp(string{ begin(seq[r2]) + end2 - len, begin(seq[r2]) + end2 }) << endl;

	}
	else
	{
		cout << string{ begin(seq[r1]) + start1, begin(seq[r1]) + start1 + len } << endl;
		cout << string{ begin(seq[r2]) + start2, begin(seq[r2]) + start2 + len } << endl;
	}
}