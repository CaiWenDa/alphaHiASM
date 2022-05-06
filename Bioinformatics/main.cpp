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
	string seqFileName = "/home/caiwenda/ecoli_20x_2.fastq";
	string outFilePre = getCurrentDate();
	string asmFileName = "result_ecoli_.fasta";
	//string kfFileName = "/home/caiwenda/software/LROD/test/kmer_file.txt";
	//string seqFileName = "/public_data/publicdata/publicdata/Reads/HiFi/D.mel/dmel_hifi_20x.fasta";
	// string kfFileName = "/publicdata/Reads/HiFi/D.mel/kmer31.txt";
	cout << "seqFile : " << seqFileName << endl;
	// cout << "frequencyFile : " << kfFileName << endl;
	seqData_t seq;
	StringSet<CharString> ID;
	clock_t start = clock();
	loadSeqData(seqFileName, ID, seq);
	seqan::clear(ID);
	auto kmerHashTable = createKmerHashTable(seq, true);
	// filterKmer(kmerHashTable, kfFileName);
	int block1 = 0;
	int block2 = 0;
	int b_size = length(seq) / thread_i;
	ofstream outFile;//("result-" + getCurrentDate() + ".csv", ios_base::out);
	ofstream seqOut(asmFileName, ios_base::out);
	vector<thread> threadPool;
	cout << "Detecting Overlap...\n";
	for (size_t i = 0; i < thread_i; i++)
	{
		block2 += b_size;
		threadPool.push_back(thread(mainProcess, 
			ref(*kmerHashTable), ref(seq), ref(ID), block1, block2, ref(outFile), 2, 1000));
		block1 = block2;
	}
	if (length(seq) % thread_i != 0)
	{
		threadPool.push_back(thread(mainProcess, 
			ref(*kmerHashTable), ref(seq), ref(ID), block1, length(seq), ref(outFile), 2, 1000));
	}
	for (auto& th : threadPool)
	{
		th.join();
	}
	kmerHashTable->clear();
	delete kmerHashTable;
	threadPool.clear();
	malloc_trim(0);
	//boost::thread_group tg;
	//tg.create_thread(bind(toyAssembly, ref(seq), block1, block2));
	cout << "Creating Overlap Graph...\n";
	extern vector<assemblyInfo_t> overlap;
	toyAssembly(seq, 0, overlap.size());
	cout << "Assembling reads...\n";
	assembler(seq, seqOut);

	extern seqData_t assemblySeq;
#if 0

	while (length(assemblySeq))
	{
		seq = assemblySeq;
		clear(assemblySeq);
		kmerHashTable = createKmerHashTable(seq, true);
		block1 = 0;
		block2 = 0;
		b_size = length(seq) / thread_i;
		//outFile.open("result-" + getCurrentDate() + ".csv", ios_base::out);
		//ofstream outFile;
		cout << "Detecting Overlap...\n";
		for (size_t i = 0; i < thread_i; i++)
		{
			block2 += b_size;
			threadPool.push_back(thread(mainProcess,
				ref(*kmerHashTable), ref(seq), ref(ID), block1, block2, ref(outFile), 2, 600));
			block1 = block2;
		}
		if (length(seq) % thread_i != 0)
		{
			threadPool.push_back(thread(mainProcess,
				ref(*kmerHashTable), ref(seq), ref(ID), block1, length(seq), ref(outFile), 2, 600));
		}
		for (auto& th : threadPool)
		{
			th.join();
		}
		kmerHashTable->clear();
		threadPool.clear();
		malloc_trim(0);
		//boost::thread_group tg;
		//tg.create_thread(bind(toyAssembly, ref(seq), block1, block2));
		cout << "Creating Overlap Graph...\n";
		extern vector<assemblyInfo_t> overlap;
		//extern vector<assemblyInfo_t> ex_overlap;
		//overlap.insert(overlap.end(), ex_overlap.begin(), ex_overlap.end());
		toyAssembly(seq, 0, overlap.size());
		cout << "Assembling reads...\n";
		assembler(seq, seqOut);
	}
	
	seqOut.close();

	loadSeqData(asmFileName, ID, seq);
	seqan::clear(ID);
	kmerHashTable = createKmerHashTable(seq, false);
	// filterKmer(kmerHashTable, kfFileName);
	block1 = 0;
	block2 = 0;
	b_size = length(seq) / thread_i;
	std::set<size_t> dump;
	threadPool.clear();
	cout << "Detecting Overlap...\n";
	for (size_t i = 0; i < thread_i; i++)
	{
		block2 += b_size;
		threadPool.push_back(thread(mainProcess2,
			ref(*kmerHashTable), ref(seq), ref(ID), block1, block2, ref(outFile), 2, 600, ref(dump)));
		block1 = block2;
	}
	if (length(seq) % thread_i != 0)
	{
		threadPool.push_back(thread(mainProcess2,
			ref(*kmerHashTable), ref(seq), ref(ID), block1, length(seq), ref(outFile), 2, 600, ref(dump)));
	}
	for (auto& th : threadPool)
	{
		th.join();
	}
	kmerHashTable->clear();
	delete kmerHashTable;
	threadPool.clear();
	malloc_trim(0);

	seqOut.open(asmFileName);
	std::set<size_t> all;
	for (size_t i = 0; i < length(seq); i++)
	{
		all.insert(i);
	}
	std::set<size_t> una;
	set_difference(all.begin(), all.end(), dump.begin(), dump.end(), inserter(una, una.end()));
	for (auto& i : una)
	{
		seqOut << boost::format(">contig_%d length=%d; reads=%d; type=dna\n") % i % length(seq[i]) % 1;
		auto a = begin(seq[i]);
		for (; a < prev(end(seq[i]), 1000); std::advance(a, 1000))
		{
			seqOut << string{ a, next(a, 1000) } << endl;
		}
		seqOut << string{ a, end(seq[i]) } << endl;
		seqOut << endl;
	}
#endif // 0
	cout << "done!\n";
	cout << "time: " << (clock() - start) / CLOCKS_PER_SEC << " sec(s)\n";
	outFile.close();
	seqOut.close();
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