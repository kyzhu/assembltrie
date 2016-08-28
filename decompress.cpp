#include <string>
#include <cstring>
#include <vector>

#include "compress.hpp"
#include "decompress.hpp"

BitReader::BitReader(const std::string &filename) {
	try {
		stream.open(filename.c_str(),
				std::ios::in | std::ios::binary | std::ios::ate);
	} catch (std::ifstream::failure e) {
		std::cerr << "File does not exist.\n";
	}
}

std::string BitReader::substring(char *bases, int start) {
	char *out = new char[L + 1];
	out[L] = '\0';
	Util::arraycpy(bases, start, out, 0, L);
	std::string res(out);
	delete[] out;
	return res;
}

int BitReader::readBits(int count) {
	int curShift = 0;
	int value = 0;
	int need = count;

	try {
		while (need > curBits) {
			value |= BitWriter::lowBits(curBits, curByte) << curShift;
			need -= curBits;
			curShift += curBits;
			curBits = BYTE;
			if (!stream.eof()) stream.read((char*) &curByte, 1);
			else curByte = -1;
		}
	} catch (std::ifstream::failure e) {
		std::cerr << "Exception reading file\n";
	}

	value |= (BitWriter::lowBits(curBits, curByte) >> (curBits - need))
			<< curShift;
	curBits -= need;

	return value;
}

int BitReader::readDistance() {
	int total = 0;
	while (true) {
		int dist = readBits(2);
		total += dist;
		if (dist != 3)
			return total;
	}
}

int BitReader::readInt() {
	if (readBits(1) == 0)
		return readBits(4);
	if (readBits(1) == 0)
		return readBits(8);
	if (readBits(1) == 0)
		return readBits(16);
	return readBits(32);
}

std::vector<std::string> BitReader::readSequence() {
	/* To return */
	std::vector<std::string> sReads;
	sReads.clear();

	/* Read lengths */
	bool full = (readBits(2) == 1);
	size_t parentSeq = readInt();
	int parentId = readInt();
	int length = readInt();
	if (length == -1)
		return sReads;
	int count = readInt(); // Number of reads
	int eCount = readInt(); // Number of errors

	/* Read bases */
	char *bases = new char[L + length + 1];
	bases[L + length] = '\0';
	for (int i = 0; i < length; i++)
		bases[L + i] = readBase();

	/* Fill in parent bases */
	if (parentSeq < sequences.size()) {
		const char *pBases = sequences[parentSeq][parentId].c_str();
		Util::arraycpy(pBases, 0, bases, 0, L);
		//delete[] pBases;
	} else if (full) {
		char *newBases = new char[length + 1];
		newBases[length] = '\0';
		Util::arraycpy(bases, L, newBases, 0, length);
		bases = newBases;
		//delete[] newBases;
	}

	/* Read reads */
	std::vector<std::string> reads;
	int curDist = 0;
	for (int i = 0; i < count; i++) {
		curDist += readDistance();
		reads.push_back(substring(bases, curDist));
	}

	/* Read errors */
	int curRead = 0;
	for (int i = 0; i < eCount; i++) {
		curRead += readDistance();
		int location = readBits(7);
		char base = readBase();
		reads[curRead][location] = base;
	}

	for (int i = 0; i < count; i++)
		sReads.push_back(reads[i]);

	/* Free temp allocations */
	delete[] bases;
	return sReads;
}

void BitReader::finish() {
	try {
		stream.close();
	} catch (std::ifstream::failure e) {
		std::cerr << "Exception closing file.\n";
	}
}

void BitReader::wrtReads() {
	/* Need to add sth!! */
	stream.seekg(0, std::ios::beg);

	while (true) {
		std::vector<std::string> reads = readSequence();
		if (reads.empty())
			break;
		sequences.push_back(reads);
	}

	std::ofstream ofile("decompressed.out");
	//std::cout<<sequences.size()<<std::endl;
	for (size_t i = 0; i < sequences.size(); i++) {
		for (size_t j = 0; j < sequences[i].size(); j++)
			ofile << sequences[i][j] << "\n";
		//std::cout<<sequences[i][0]<<std::endl;
	}

	ofile.close();
	finish();
}

