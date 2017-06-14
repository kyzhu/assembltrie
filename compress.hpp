#ifndef COMPRESS_HPP_
#define COMPRESS_HPP_

#include <fstream>
#include <experimental/string_view>
#include <vector>

#include "error.hpp"
#include "read.hpp"

class BitWriter {
public:
	/* for test */
	//int test = 0;
	bool heuristic_flag = false;
	std::ofstream stream;
	int curByte = 0;
	int curBits = 0;
	int L = 100;
	static const int BYTE = 8;
	static const long MAXINT = (1L << 32) - 1;

	std::vector<int> seqSize;

	/* The position of letter N. */
	std::vector<int> Nstream;
	int block_length = 0;
	//int sum_length = 0;
	int nflag = 0;

	BitWriter(const std::string &filename, int);
	~BitWriter() {
	}

	char vote(int*, int);

	void dumpReads(std::vector<Read*> reads, int);

	void writeSequence(int, int, size_t, int, char*, int, std::vector<int>, std::vector<int>,
			std::vector<int>, std::vector<Error>);

	void writeUncompressedReads(std::vector<std::experimental::string_view>, std::vector<int>);

	static int lowBits(int, long);

	static int highBits(int, int, long);

	void writeBits(int, int);

	void writeInt(int);

	void writeInt(int, int);

	void writeTwo(int);

	void writeBase(char);

	void writeErrBase(char, char);

	void writeSuffDistance(int);

	void writeErrNum(int);

	void writeLastSequence();

	void fillNstream(char);

	void sendNstream(BitWriter*);

	void writeNstream(bool);

	void finish();

	static int const suffix_len_codebook[2][254];

	static int const error_num_codebook[2][256];

};

#endif

