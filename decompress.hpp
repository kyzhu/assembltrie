#ifndef DECOMPRESS_HPP_
#define DECOMPRESS_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "util.hpp"

class BitReader {
public:
	std::ifstream stream;
	int curBits = 0;
	int curByte = 0;
	int L = 101;
	static const int BYTE = 8;

	std::vector<std::vector<std::string>> sequences;

	/* Contructors */
	BitReader(const std::string &filename);
	~BitReader() {
	}

	std::vector<std::string> readSequence();

	std::string substring(char*, int);

	char readBase() {
		return Util::intToBase(readBits(3));
	}

	int readDistance();

	/* Reads a variable-length integer. */
	int readInt();

	/* Reads COUNT bits. */
	int readBits(int);

	/* Finishes reading from file. */
	void finish();

	/* Write the decompressed reads to a file */
	void wrtReads();

};

#endif

