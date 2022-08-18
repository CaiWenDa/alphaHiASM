#include "assembly.h"
#include "utility.h"
#include <cstdlib>
#include <malloc.h>
#include <getopt.h>
#include <iostream>

extern int optind, opterr, optopt;
extern char* optargi;

int parseOption(int argc, char* argv[], std::string& file, std::string& outFile, int& genomeSize, int& thread_i)
{
	int index = 0;
	int c = 0; //用于接收选项
	//定义长选项
	static struct option long_options[] =
	{
		{"help", no_argument, NULL, 'h'},
		{"file", required_argument, NULL, 'f'},
		{"output", required_argument, NULL, 'o'},
		{"genomeSize", required_argument, NULL, 'g'},
		{"theads", required_argument, NULL, 't'}
	};

	/*循环处理参数*/
	while (EOF != (c = getopt_long(argc, argv, "hf:o:g:t:", long_options, &index)))
	{
		using std::cerr;
		using std::endl;
		switch (c)
		{
		case 'h':
			cerr << 
			"usage: ToyAssembly \n\
		--file - file_1[file_2 ...]\n\
		--outfile - dir PATH\n\
		--genome - size SIZE\n\
		--threads int\n\
		[--min - overlap SIZE]\n\
		[--help]\n\
		[--read - error float]\n";
			break;
		case 'f':
			file = optarg;
			break;
			//-n选项必须要参数
		case 'o':
			outFile = optarg;
			break;
		case 'g':
			genomeSize = atoi(optarg);
			break;
		case 't':
			thread_i = atoi(optarg);
			break;
			//表示选项不支持
		case '?':
			cerr << "unknown option: " << static_cast<char>(optopt) << endl;
			break;
		default:
			break;
		}
	}
	return 0;
}

void compSeqInRange(cwd::seqData_t& seq, uint r1, uint r2, uint start1, uint start2, uint end1, uint end2, uint len, bool orient = true);

int main(int argc, char* argv[])
{
	using namespace std;
	using namespace seqan;
	using namespace cwd;
	int ovLen = 2000;
	extern int thread_i;
	extern int genomeSize;
	//string seqFileName = "/home/caiwenda/dmel.trimmedReads_20x.fasta";
	string seqFileName;//argv[2];//"/home/caiwenda/SRR11292120_sample.fastq";
	string asmFileName;//argv[3];//"result_chm13_debug713.fasta";
	//string asmFileName = "tmp.fasta";
	string outFilePre = getCurrentDate();
	//string seqFileName = "result_ecoli_debug.fasta";
	string kfFileName = "/home/caiwenda/ecoli_kmerFrequency.txt";
	//string seqFileName = "dmel.trimmedReads_10x.fasta";
	// string kfFileName = "/publicdata/Reads/HiFi/D.mel/kmer31.txt";
	parseOption(argc, argv, seqFileName, asmFileName, genomeSize, thread_i);
	if (argc <= 4)
	{
		return 0;
	}
	cerr << "seqFile : " << seqFileName << endl;
	// cout << "frequencyFile : " << kfFileName << endl;
	seqData_t seq;
	StringSet<CharString> ID;
	clock_t start = clock();

	try
	{
		loadSeqData(seqFileName, ID, seq);
	}
	catch (exception & e)
	{
		cerr << e.what() << endl;
		getchar();
		return 1;
	}
	seqan::clear(ID);
	malloc_trim(0);
	auto kmerHashTable = createKmerHashTable(seq, true);
	//filterKmer(*kmerHashTable, kfFileName);
	uint block1 = 0;
	uint block2 = 0;
	uint b_size = length(seq) / thread_i;
	ofstream outFile("result-" + getCurrentDate() + ".csv", ios_base::out);
	ofstream seqOut(asmFileName, ios_base::out);
	//ofstream seqOut;
	vector<thread> threadPool;
	cerr << "Detecting Overlap...\n";
#if 1
	for (size_t i = 0; i < thread_i; i++)
	{
		block2 += b_size;
		threadPool.emplace_back(mainProcess, 
		                        ref(*kmerHashTable), ref(seq), ref(ID), block1, block2, ref(outFile), 2, ovLen);
		block1 = block2;
	}
	if (length(seq) % thread_i != 0)
	{
		threadPool.emplace_back(mainProcess, 
		                        ref(*kmerHashTable), ref(seq), ref(ID), block1, length(seq), ref(outFile), 2, ovLen);
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
	//tg.create_thread(bind(createOverlapGraph, ref(seq), block1, block2));
	cerr << "Creating Overlap Graph...\n";
	extern vector<assemblyInfo_t> overlap;
	createOverlapGraph(seq, 0, overlap.size());
	cerr << "Assembling reads...\n";
	assembler(seq, seqOut);

	if (false)
	{
		std::set<size_t> all;
		for (size_t i = 0; i < length(seq); i++)
		{
			all.insert(i);
		}
		extern std::set<size_t> delReads;
		std::set<size_t> una;
		set_difference(all.begin(), all.end(), delReads.begin(), delReads.end(), inserter(una, una.end()));
		if (!una.empty())
		{
			for (auto& i : una)
			{
				auto assembly = seq[i];
				seqOut << boost::format(">contig_%d; length=%d; reads=%d; type=dna\n") % (length(seq) + i) % length(assembly) % 1;
				auto a = begin(assembly);
				for (; a < prev(end(assembly), 1000); std::advance(a, 1000))
				{
					seqOut << string{ a, next(a, 1000) } << endl;
				}
				seqOut << string{ a, end(assembly) } << endl;
				seqOut << endl;
			}
		}
	}

	//ifstream overlapFile;
	//assemblyInfo_t ovl;
	//uint len1, len2;
	//while (overlapFile >> ovl.r1 >> ovl.r2 >> ovl.orient >> ovl.SP1 >> ovl.EP1 >> ovl.SP2 >> ovl.EP2 >> len1 >> len2)
	//{
	//	overlap.emplace_back(ovl);
	//}
#endif

#if 0
	seqOut.close();
	seqOut.open("result_ecoli_debug3.fasta");
	extern seqData_t assemblySeq;
	extern int OVL_TIP_LEN;
	OVL_TIP_LEN = 200000;
	ovLen = 1000;
	while (length(assemblySeq) > 1)
	{
		seq = assemblySeq;
		clear(assemblySeq);
		kmerHashTable = createKmerHashTable(seq);
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
				ref(*kmerHashTable), ref(seq), ref(ID), block1, block2, ref(outFile), 2, ovLen));
			block1 = block2;
		}
		if (length(seq) % thread_i != 0)
		{
			threadPool.push_back(thread(mainProcess,
				ref(*kmerHashTable), ref(seq), ref(ID), block1, length(seq), ref(outFile), 2, ovLen));
		}
		for (auto& th : threadPool)
		{
			th.join();
		}
		kmerHashTable->clear();
		threadPool.clear();
		malloc_trim(0);
		//boost::thread_group tg;
		//tg.create_thread(bind(createOverlapGraph, ref(seq), block1, block2));
		cout << "Creating Overlap Graph...\n";
		extern vector<assemblyInfo_t> overlap;
		createOverlapGraph(seq, 0, overlap.size());
		cout << "Assembling reads...\n";
		assembler(seq, seqOut);
	}
#endif

#if 0
	seqOut.close();
	seqan::clear(seq);
	loadSeqData(asmFileName, ID, seq);
	seqan::clear(ID);
	kmerHashTable = createKmerHashTable(seq, false);
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

	seqOut.open("remove_dup.fasta");
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
#endif

	cerr << "done!\n";
	cerr << "time: " << (clock() - start) / CLOCKS_PER_SEC << " sec(s)\n";
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
		cerr << string{ begin(seq[r1]) + start1, begin(seq[r1]) + start1 + len } << endl;
		cerr << cwd::revComp(string{ begin(seq[r2]) + end2 - len, begin(seq[r2]) + end2 }) << endl;

	}
	else
	{
		cerr << string{ begin(seq[r1]) + start1, begin(seq[r1]) + start1 + len } << endl;
		cerr << string{ begin(seq[r2]) + start2, begin(seq[r2]) + start2 + len } << endl;
	}
}