#include <string>
#include <cstring>
#include <vector>
#include <assert.h>

#include "compress.hpp"
#include "decompress.hpp"

BitReader::BitReader(int l, int h, const std::string &cfilename) {
	L = l;
	heuristic_flag = h;
	filename = cfilename;
	try {
		stream.open("part.out",
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
			value = (value << curShift) + BitWriter::lowBits(curBits, curByte);
			need -= curBits;
			curShift = BYTE;
			curBits = BYTE;
			if (!stream.eof())
				stream.read((char*) &curByte, 1);
			else
				curByte = -1;
		}
	} catch (std::ifstream::failure e) {
		std::cerr << "Exception reading file\n";
	}
	value = (value << need)
			+ (BitWriter::lowBits(curBits, curByte) >> (curBits - need));
	curBits -= need;
	return value;
}

int BitReader::readErrNum() {
	int codeword = 0, nbits = 0;
	while (true) {
		nbits++;
		int bit = readBits(1);
		codeword = (codeword << 1) + bit;
		for (int i = 0; i < L + 1; i++)
			if (nbits == BitWriter::error_num_codebook[0][i]
					&& codeword == BitWriter::error_num_codebook[1][i])
				return i;
	}
}

char BitReader::readBase() {
	return Util::intToBase(readBits(2));
}

char BitReader::readErrBase(char u) {
	int bit1 = readBits(1);
	if (bit1 == 0)
		return Util::intToBase(0, u);
	else {
		return Util::intToBase(readBits(1) + 1, u);
	}
}

int BitReader::readSuffDistance() {
	int codeword = 0, nbits = 0;
	while (true) {
		nbits++;
		int bit = readBits(1);
		codeword = (codeword << 1) + bit;
		if (nbits == lengthBits && codeword == lengthCode)
			return L;
		if (nbits == stopBits && codeword == stopCode)
			return -1;
		for (int i = 0; i < L; i++)
			if (nbits == BitWriter::suffix_len_codebook[0][i]
					&& codeword == BitWriter::suffix_len_codebook[1][i])
				return i;
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

int BitReader::readInt(int upper_bound) {
	if (upper_bound < 16)
		return readBits(32 - __builtin_clz(upper_bound));
	if (upper_bound < 256) {
		if (readBits(1) == 0)
			return readBits(4);
		else
			return readBits(32 - __builtin_clz(upper_bound));
	}
	if (upper_bound < 65536) {
		if (readBits(1) == 0)
			return readBits(8);
		else
			return readBits(32 - __builtin_clz(upper_bound));
	}
	if (readBits(1) == 0)
		return readBits(16);
	else
		return readBits(32 - __builtin_clz(upper_bound));
}

std::vector<std::string> BitReader::readUncompressedReads() {
	std::vector<std::string> htcReads;
	numhtcReads = readInt();
	int rc = 0;
	for (int i = 0; i < numhtcReads; i++) {
		if (heuristic_flag > 0)
			rc = readBits(1);
		char *bases = new char[L + 1];
		bases[L] = '\0';
		for (int j = 0; j < L; j++)
			bases[j] = readBase();
		if (rc == 0)
			htcReads.push_back(substring(bases, 0));
		else {
			std::string rc = substring(bases, 0);
			htcReads.push_back(Util::getRCstring(rc));
		}
		delete[] bases;
	}
	return htcReads;
}

/* Recover the positions of symbol N for UNCONNECTED reads. */
void BitReader::recUncompressedNstream() {
	size_t read_index = 0;

	/* Run-length decode. */
	while (true) {
		int next_length = readInt();
		int n_length = readErrNum() + 1;
		sum_length += next_length;

		int j = sum_length % L;
		read_index = sum_length / L;
		if (read_index == htcReads.size())
			return;
		for (int i = 0; i < n_length; i++) {
			htcReads[read_index][j++] = 'N'; /* Recover letter N. */
			if (j >= L) {
				j = 0;
				read_index++;
				if (read_index == htcReads.size())
					return;
			}
		}
		sum_length += n_length;
	}
}

void BitReader::readNstream() {
	int num_blocks = readInt();
	for (int i = 0; i < num_blocks - 1; i++) {
		Nstream.push_back(readInt());
		Nstream.push_back(readErrNum() + 1);
	}
	Nstream.push_back(readInt());
	Nstream.push_back(readErrNum());
}

void BitReader::recNstream(std::vector<std::string> &_reads) {
	int read_index = 0, read_pos = 0;

	int num_reads = _reads.size();
	if (nflag == 0) {
		int remaining = sum_length - scanned_length;
		for (int i = 0; i < remaining; i++) {
			_reads[read_index][read_pos++] = 'N'; /* Recover letter N. */
			if (read_pos >= L) {
				read_pos = 0;
				read_index++;
			}
			if (read_index == num_reads)
				break;
		}
	}
	while (sum_length < scanned_length + num_reads * L) { /* Read in next block. */
		nflag = 1 - nflag;
		if (nflag == 0) {
			int block_length = Nstream[ns_index++];
			read_index = (sum_length - scanned_length) / L;
			read_pos = sum_length % L;
			sum_length += block_length;

			for (int i = 0; i < block_length && read_index < num_reads; i++) {
				_reads[read_index][read_pos++] = 'N'; /* Recover letter N. */
				if (read_pos >= L) {
					read_pos = 0;
					read_index++;
				}
			}
		} else {
			int block_length = Nstream[ns_index++];
			sum_length += block_length;
		}
	}
	scanned_length += num_reads * L;
}

void BitReader::openCompressedFile() {
	finish();
	curBits = 0;
	curByte = 0;
	try {
		stream.open(filename.c_str(),
				std::ios::in | std::ios::binary | std::ios::ate);
	} catch (std::ifstream::failure e) {
		std::cerr << "File does not exist.\n";
	}
}

std::vector<std::string> BitReader::readSequence() {
	/* To return */
	std::vector<std::string> reads;
	reads.clear();

	/* Read parent sequence information. */
	bool full = (readBits(1) == 1);
	bool parentSeqIndicator = (readBits(1) == 1);
	size_t parentSeq = 0;
	if (parentSeqIndicator)
		parentSeq = readInt(sequences.size() + 1);
	else
		parentSeq = sequences.size();

	/* Read suffix lengths of reads. */
	int length = 0;
	int count = 0; /* Number of reads. */
	std::vector<int> distances;
	while (curByte != -1) {
		int curDist = readSuffDistance();
		if (curDist != -1) {
			length += curDist;
			distances.push_back(length);
			count += 1;
		} else
			break;
	}

	/* Decompress completed. */
	if (count == 0) {
		return reads;
	}

	if (full)
		length += L - distances[0];

	/* If heuristic is used to identify the orientation of the reads. */
	std::vector<int> rc_bits;
	if (heuristic_flag > 0)
		for (int i = 0; i < count; i++)
			rc_bits.push_back(readBits(1));

	/* Read parent read information. */
	int parentId = 0;
	if (parentSeq != sequences.size()) {
		parentId = readBits(32 - __builtin_clz(sequences[parentSeq].size()));
	} else if (full) {
		parentId = readBits(32 - __builtin_clz(count));
	}

	/* Read in bases */
	char *bases = new char[L + length + 1];
	bases[L + length] = '\0';
	for (int i = 0; i < length; i++)
		bases[L + i] = readBase();

	/* Fill in parent bases */
	if (parentSeq < sequences.size()) {
		const char *pBases = sequences[parentSeq][parentId].c_str();
		Util::arraycpy(pBases, 0, bases, 0, L);
		for (int i = 0; i < L; i++)
			if (sequences[parentSeq][parentId][i] == 'N')
				bases[i] = 'A';
		//delete[] pBases;
	} else if (full) {
		char *newBases = new char[length + distances[0] + 1];
		newBases[length] = '\0';
		Util::arraycpy(bases, L, newBases, distances[0], length);
		for (int i = 0; i < distances[0]; i++)
			newBases[i] = 'A';
		bases = newBases;
		//delete[] newBases;
	}

	/* Recover reads */
	for (int i = 0; i < count; i++)
		reads.push_back(substring(bases, distances[i]));

	/* Apply read errors */
	recNstream (reads);
	if (full)
		for (int i = 0; i < L; i++)
			if (reads[parentId][i] == 'N')
				bases[i] = 'A';

	for (int i = 0; i < count; i++) {
		int curError = readErrNum();
		int eDist = 0;
		for (int j = 0; j < curError; j++) {
			int location = readBits(32 - __builtin_clz(L - eDist)) + eDist;
			char base = readErrBase(bases[distances[i] + location]);
			reads[i][location] = base;
			eDist = location;
		}
	}

	/* Record the orientation */
	if (heuristic_flag)
		rc_flags.push_back(rc_bits);

	/* Free temp allocations */
	delete[] bases;
	return reads;
}

void BitReader::finish() {
	try {
		stream.close();
	} catch (std::ifstream::failure e) {
		std::cerr << "Exception closing file.\n";
	}
}

void BitReader::wrtReads(std::string outputfile) {
	/* Recover Unconnected Reads. */
	stream.seekg(0, std::ios::beg);
	readNstream();
	htcReads = readUncompressedReads();
	recUncompressedNstream();

	/* Locate the compressed (major) file pointer. */
	openCompressedFile();
	stream.seekg(0, std::ios::beg);
	sum_length = 0;

	/* Reconstruct reads. */
	while (true) {
		std::vector<std::string> reads = readSequence();
		if (reads.empty())
			break;
		sequences.push_back(reads);
	}
	assert(
			(ns_index >= Nstream.size() - 2)
					&& (ns_index <= Nstream.size() - 1));
	assert(sum_length == scanned_length);

	/* Output the reads. */
	std::ofstream ofile(outputfile);
	if (heuristic_flag > 0)
		for (size_t i = 0; i < sequences.size(); i++)
			for (size_t j = 0; j < sequences[i].size(); j++)
				if (rc_flags[i][j] == 0)
					ofile << sequences[i][j] << "\n";
				else
					ofile << Util::getRCstring(sequences[i][j]) << "\n";
	else
		for (size_t i = 0; i < sequences.size(); i++)
			for (size_t j = 0; j < sequences[i].size(); j++)
				ofile << sequences[i][j] << "\n";
	for (size_t i = 0; i < htcReads.size(); i++)
		ofile << htcReads[i] << "\n";

	ofile.close();
	finish();
}

