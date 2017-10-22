#ifndef DECOMPRESS_HPP_
#define DECOMPRESS_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "util.hpp"

class BitReader {
public:
	//int test = 0;
	int heuristic_flag = 0;
	std::ifstream stream;
	std::string filename;
	int curBits = 0;
	int curByte = 0;
	int L = 100;
	static const int BYTE = 8;

	int stopCode = 2, stopBits = 3;
	int lengthCode = 4, lengthBits = 5;

	std::vector<std::vector<std::string>> sequences;
	std::vector<std::string> htcReads;
	int numhtcReads = 0;
	std::vector<std::vector<int>> rc_flags;

	std::vector<int> Nstream;
	unsigned long long sum_length = 0, scanned_length = 0;
	int nflag = 0;
	size_t ns_index = 0;

	/* Contructors */
	BitReader(int, int, const std::string &cfilename);
	~BitReader() {
	}

	void openCompressedFile();

	std::vector<std::string> readSequence();

	void readNstream();

	void recNstream(std::vector<std::string>&);

	std::vector<std::string> readUncompressedReads();

	void recUncompressedNstream();

	std::string substring(char*, int);

	char readBase();

	char readErrBase(char);

	int readSuffDistance();

	int readErrNum();

	/* Reads a variable-length integer. */
	int readInt();
	int readInt(int);

	/* Reads COUNT bits. */
	int readBits(int);
	//int readBits(int, bool);

	/* Finishes reading from file. */
	void finish();

	/* Write the decompressed reads to a file */
	void wrtReads(std::string);

};

#endif

