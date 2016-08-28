#ifndef COMPRESS_HPP_
#define COMPRESS_HPP_

#include <fstream>
#include <string>
#include <vector>

#include "error.hpp"

class BitWriter {
public:
	std::ofstream stream;
	int curByte = 0;
	int curBits = 0;
	static const int BYTE = 8;
	static const long MAXINT = (1L << 32) - 1;

	BitWriter(const std::string &filename);
	~BitWriter() {
	}

	void print();

	static std::string join();

	void writeSequence(int, int, int, char*, int, std::vector<int>,
			std::vector<int>, std::vector<Error>);

	static int lowBits(int, long);

	void writeBits(int, int);

	void writeInt(int);

	void writeTwo(int);

	void writeBase(char);

	void writeDistance(int);

	void finish();

};

#endif
