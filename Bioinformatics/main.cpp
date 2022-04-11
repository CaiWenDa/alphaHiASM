#include "proccess.h"
#include <malloc.h>

void compSeqInRange(cwd::seqData_t& seq, uint r1, uint r2, uint start1, uint start2, uint end1, uint end2, uint len, bool orient = true);

int main(int argc, char* argv[])
{
	using namespace std;
	using namespace seqan;
	using namespace cwd;

	//if (argc <= 1)
	//{
	//    return 0;
	//}
	string seqFileName = "/home/caiwenda/ecoli_40x.fastq";
	//string kfFileName = "/home/caiwenda/software/LROD/test/kmer_file.txt";
	//string seqFileName = "/public_data/publicdata/publicdata/Reads/HiFi/D.mel/dmel_hifi_40x_sample.fasta";
	// string kfFileName = "/publicdata/Reads/HiFi/D.mel/kmer31.txt";
	cout << "seqFile : " << seqFileName << endl;
	// cout << "frequencyFile : " << kfFileName << endl;
	seqData_t seq;
	StringSet<CharString> ID;
	clock_t start = clock();
	loadSeqData(seqFileName, ID, seq);
	seqan::clear(ID);
	auto kmerHashTable = createKmerHashTable(seq);
	// filterKmer(kmerHashTable, kfFileName);
	int block1 = 0;
	int block2 = 0;
	int b_size = length(seq) / thread_i;
	ofstream outFile("result-" + getCurrentDate() + ".csv", ios_base::out);
	//ofstream outFile;
	vector<thread*> threadPool;
	cout << "Detecting Overlap...\n";
	for (size_t i = 0; i < thread_i; i++)
	{
		block2 += b_size;
		threadPool.push_back(new thread(mainProcess, 
			ref(*kmerHashTable), ref(seq), ref(ID), block1, block2, ref(outFile)));
		block1 = block2;
	}
	if (length(seq) % thread_i != 0)
	{
		threadPool.push_back(new thread(mainProcess, 
			ref(*kmerHashTable), ref(seq), ref(ID), block1, length(seq), ref(outFile)));
	}
	for (auto& th : threadPool)
	{
		th->join();
		delete th;
	}
	delete kmerHashTable;
	threadPool.clear();
	malloc_trim(0);
	//boost::thread_group tg;
	//tg.create_thread(bind(toyAssembly, ref(seq), block1, block2));
	if (1)
	{
		cout << "Creating Overlap Graph...\n";
		extern vector<assemblyInfo_t> overlap;
		toyAssembly(seq, 0, overlap.size());
		cout << "Assembling reads...\n";
		assembler(seq);
	}
	cout << "done!\n";
	cout << "time: " << (clock() - start) / CLOCKS_PER_SEC << " sec(s)\n";
	outFile.close();
	getchar();
	return 0;
}

void compSeqInRange(cwd::seqData_t& seq, uint r1, uint r2, uint start1, uint start2, uint end1, uint end2, uint len, bool orient)
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