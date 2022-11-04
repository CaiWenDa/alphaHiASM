# alphaHiASM
A developing long reads assmebly

# Before installing
    This project is depended on "seqan" and "boost" libraries. You should download them firstly.
    seqan: http://packages.seqan.de/seqan-library/seqan-library-2.4.0.zip
    boost: https://boostorg.jfrog.io/artifactory/main/release/1.78.0/source/boost_1_78_0.tar.gz

# Install
    Create a main directory (eg:alphaHiASM). 
    Copy all source code to alphaHiASM.
    Copy "seqan-library-2.4.0/*" to alphaHiASM/seqan
    Copy "boost_1_78_0/*" to alphaHiASM/boost.
    cd alphaHiASM
    make all

# Usage
	alphaHiASM
	--file - file
	--output - dir PATH
	--genomeSize SIZE
	--threads int
	[--minOverlap SIZE]
	[--help]
	[--paf inOverlapFile]
	[--overlap outOverlapFile]
