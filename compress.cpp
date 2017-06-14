#include <iostream>
#include <fstream>
#include <experimental/string_view>
#include <vector>

#include "read.hpp"
#include "compress.hpp"
#include "error.hpp"
#include "util.hpp"

BitWriter::BitWriter(const std::string &filename, int _L) {
	L = _L;
	try {
		stream.open(filename.c_str(), std::ios::out | std::ios::binary);
	} catch (std::ofstream::failure e) {
		std::cerr << "Cannot open file \'" << filename << "\'.\n";
	}
}

char BitWriter::vote(int *votes, int row) {
	int best = 0;
	for (int col = 1; col < 4; col++) {
		if (votes[row * 4 + col] > votes[row * 4 + best])
			best = col;
	}
	return Util::intToBase(best);
}

/* Prepare the sequence to be written. */
void BitWriter::dumpReads(std::vector<Read*> reads, int seq) {
	/* Distances of reads, starting from position 0. */
	std::vector<int> distances;

	/* One bit per read, indicating its orientation. */
	std::vector<int> rc_bits;

	int id = 0;
	int offset = 0;
	for (size_t i = 0; i < reads.size(); i++) {
		reads[i]->seqId = id;
		distances.push_back(reads[i]->suffix);
		rc_bits.insert(rc_bits.end(), reads[i]->rc.begin(), reads[i]->rc.end());
		offset += reads[i]->suffix;
		for (int j = 1; j < reads[i]->count; j++)
			distances.push_back(0);
		id += reads[i]->count;
	}

	/* Total length of added bases */
	int length = offset;

	/* Record #reads. */
	seqSize.push_back(id);

	/* Vote on bases */
	int *votes = new int[length * 4];
	std::memset(votes, 0, length * 4 * sizeof(int));
	offset = 0;
	for (size_t i = 0; i < reads.size(); i++) {
		offset += reads[i]->suffix;
		std::experimental::string_view curBases = reads[i]->bases();
		for (int tmp = 0; tmp < reads[i]->count; tmp++)
			for (int j = 0; j < L; j++) {
				fillNstream(curBases[j]);
				if (offset + j - L >= 0)
					votes[(offset + j - L) * 4 + Util::baseToInt(curBases[j])] +=
							1;
			}
	}

	/* Set voted bases */
	char *newBases = new char[length + 1];
	newBases[length] = '\0';
	for (int i = 0; i < length; i++)
		newBases[i] = vote(votes, i);

	/* Set all bases (with parent bases) */
	char *bases = new char[L + length + 1];
	bases[L + length] = '\0';
	Read *start = reads[0];
	bool fullBases = false;
	int seqIdBits;
	/* Comment: don't touch the first L bases when FULLBASES is false. */
	if (start->parent == NULL) {
		start->parent = start;
		seqIdBits = 0;
	} else {
		int parentSeq = start->parent->seq;
		seqIdBits = 32 - __builtin_clz(seqSize[parentSeq]);
		if (parentSeq == seq)
			fullBases = true;
	}
	Util::arraycpy(start->parent->bases().data(), 0, bases, 0, L);
	Util::arraycpy(newBases, 0, bases, L, length);
	for (int i = 0; i < L; i++)
		if (bases[i] == 'N')
			bases[i] = 'A';

	/* Record errors */
	std::vector<int> errorNums;
	std::vector<Error> errors;
	int eDist = 0, last = 0;
	offset = 0;
	for (size_t i = 0; i < reads.size(); i++) {
		offset += reads[i]->suffix;
		std::experimental::string_view curBases = reads[i]->bases();
		for (int tmp = 0; tmp < reads[i]->count; tmp++) {
			last = eDist;
			for (int j = 0; j < L; j++)
				if (curBases[j] != 'N' && curBases[j] != bases[offset + j]) {
					errors.push_back(Error(j, curBases[j], bases[offset + j]));
					eDist++;
				}
			errorNums.push_back(eDist - last);
		}
	}

	if (fullBases)
		writeSequence(1, seqIdBits, start->parent->seq, start->parent->seqId,
				bases + distances[0], L + length - distances[0], distances,
				rc_bits, errorNums, errors);
	else
		writeSequence(0, seqIdBits, start->parent->seq, start->parent->seqId,
				newBases, length, distances, rc_bits, errorNums, errors);

	delete[] votes;
	delete[] bases;
	delete[] newBases;
}

/* Run-length encoding of the positions of N. */
void BitWriter::fillNstream(char c) {
	if ((c == 'N' && nflag == 0) || (c != 'N' && nflag == 1)) {
		Nstream.push_back(block_length);
		block_length = 1;
		nflag = 1 - nflag;
	} else
		block_length++;
}

/* Writes a fake sequence as the bound of N stream. */
void BitWriter::writeLastSequence() {
	writeBits(5, 2);
}

/* Writes the information of a sequence to file. */
void BitWriter::writeSequence(int full, int seqIdBits, size_t parentSeq,
		int parentId, char *bases, int length, std::vector<int> distances,
		std::vector<int> rc_bits, std::vector<int> errorNums,
		std::vector<Error> errors) {
	/* Write parent sequence information. */
	writeBits(1, full);

	if (parentSeq == seqSize.size() - 1)
		writeBits(1, 0);
	else {
		writeBits(1, 1);
		writeInt(parentSeq, seqSize.size());
	}

	/* Write read distances */
	for (int distance : distances)
		writeSuffDistance(distance);
	/* Stop word: 010 */
	writeBits(3, 2);

	/* Reverse complement indicator */
	if (heuristic_flag)
		for (int rc_bit : rc_bits)
			writeBits(1, rc_bit);

	/* Write parent read information. */
	if (seqIdBits > 0)
		writeBits(seqIdBits, parentId);

	/* Write sequence */
	for (int i = 0; i < length; i++)
		writeBase(bases[i]);

	/* Write errors */
	int eDist = 0;
	for (size_t i = 0; i < errorNums.size(); i++) {
		/* Write error numbers. */
		writeErrNum(errorNums[i]);

		/* Write errors with simple differential code. */
		int readDist = 0;
		for (int j = 0; j < errorNums[i]; j++) {
			Error e = errors[eDist++];
			writeBits(32 - __builtin_clz(L - readDist), e.location - readDist);
			readDist = e.location;
			writeErrBase(e.base, e.underlying);
		}
	}
}

void BitWriter::writeUncompressedReads(
		std::vector<std::experimental::string_view> htcReads,
		std::vector<int> htcBits) {
	/* Calculate and encode the number of unconnected reads. */
	//std::cout << htcReads.size() << " single reads" << std::endl;
	writeInt(htcReads.size());
	Nstream.clear();

	/* Encode the bases. */
	char next_char;
	for (size_t i = 0; i < htcReads.size(); i++) {
		if (heuristic_flag != 0)
			writeBits(1, htcBits[i]);
		for (int j = 0; j < L; j++) {
			next_char = htcReads[i][j];
			fillNstream(next_char);
			writeTwo(Util::baseToInt(next_char));
		}
	}
	Nstream.push_back(block_length);
	if (Nstream.size() % 2 != 0)
		Nstream.push_back(0);
}

/* Send the information of N to another bit stream DEST. */
void BitWriter::sendNstream(BitWriter *dest) {
	Nstream.push_back(block_length);
	if (Nstream.size() % 2 != 0)
		Nstream.push_back(0);
	dest->Nstream = Nstream;
}

/* Encode the run-length code. */
void BitWriter::writeNstream(bool _write_nsize) {
	/* Encode the run-length code. */
	size_t n = Nstream.size() / 2;
	if (_write_nsize)
		writeInt(n);
	for (size_t i = 0; i < n - 1; i++) {
		writeInt (Nstream[2 * i]);
		writeErrNum(Nstream[2 * i + 1] - 1);
	}
	writeInt (Nstream[2 * n - 2]);
	writeErrNum(Nstream[2 * n - 1]);
}

/* Returns the lowest COUNT bits of VALUE. */
int BitWriter::lowBits(int count, long value) {
	return (int) ((int) (value & BitWriter::MAXINT) & ((1L << count) - 1));
}

/* Returns the highest COUNT bits of VALUE. */
int BitWriter::highBits(int begin, int count, long value) {
	return (int) (value >> (begin - count));
}

/* Writes the lowest COUNT bits of VALUE to file. */
void BitWriter::writeBits(int count, int value) {
	/* Write next byte when available */
	try {
		while (curBits + count >= BYTE) {
			int toWrite = BYTE - curBits;
			curByte = (curByte << toWrite) + highBits(count, toWrite, value);
			count -= toWrite;
			value = lowBits(count, value);

			stream.write((char*) &curByte, 1);
			curByte = 0;
			curBits = 0;
		}
	} catch (std::ofstream::failure e) {
		std::cerr << "Exception writing to file\n";
	}

	/* Add to current byte */
	curByte = (curByte << count) + lowBits(count, value);
	curBits += count;
}

/* Writes a variable-length integer to file. */
void BitWriter::writeInt(int value) {
	if (value < 16) {
		writeBits(1, 0);
		writeBits(4, value);
	} else if (value < 256) {
		writeBits(2, 2);
		writeBits(8, value);
	} else if (value < 65536) {
		writeBits(3, 6);
		writeBits(16, value);
	} else {
		writeBits(3, 7);
		writeBits(32, value);
	}
}

void BitWriter::writeInt(int value, int upper_bound) {
	if (upper_bound < 16) {
		writeBits(32 - __builtin_clz(upper_bound), value);
	} else if (upper_bound < 256) {
		if (value < 16) {
			writeBits(1, 0);
			writeBits(4, value);
		} else {
			writeBits(1, 1);
			writeBits(32 - __builtin_clz(upper_bound), value);
		}
	} else if (upper_bound < 65536) {
		if (value < 256) {
			writeBits(1, 0);
			writeBits(8, value);
		} else {
			writeBits(1, 1);
			writeBits(32 - __builtin_clz(upper_bound), value);
		}
	} else {
		if (value < 65536) {
			writeBits(1, 0);
			writeBits(16, value);
		} else {
			writeBits(1, 1);
			writeBits(32 - __builtin_clz(upper_bound), value);
		}
	}
}

/* Writes two bits to file. */
void BitWriter::writeTwo(int value) {
	writeBits(2, value);
}

/* Writes a base C to file. */
void BitWriter::writeBase(char c) {
	writeTwo(Util::baseToInt(c));
}

/* Writes an error base C to file. */
void BitWriter::writeErrBase(char c, char u) {
	int toWrite = Util::baseToIntHuffman(c, u);
	switch (toWrite) {
	case 0:
		writeBits(1, 0);
		break;
	case 2:
		writeBits(2, 2);
		break;
	case 3:
		writeBits(2, 3);
		break;
	}
}

/* Writes suffix distances. */
/* Should throw an error when DISTANCE >= 254 */
void BitWriter::writeSuffDistance(int distance) {
	int nbits = 5, codeword = 4; /* Default codeword: 00100. */
	if (distance < L) {
		nbits = suffix_len_codebook[0][distance];
		codeword = suffix_len_codebook[1][distance];
	}
	writeBits(nbits, codeword);
}

/* Writes error numbers for each read. */
/* Should throw an error when NUM >= 256 */
void BitWriter::writeErrNum(int num) {
	writeBits(error_num_codebook[0][num], error_num_codebook[1][num]);
}

/* Finishes writing to file. */
void BitWriter::finish() {
	/* Write terminator and padding bits */
	writeBits(BYTE, 0);
	/* Close file. */
	try {
		stream.close();
	} catch (std::ofstream::failure e) {
		std::cerr << "Exception closing file.\n";
	}
}

/* Codebooks. */
const int BitWriter::suffix_len_codebook[2][254] = { { 6, 4, 5, 5, 5, 5, 5, 5,
		5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9,
		9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10,
		10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12,
		12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
		12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
		12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
		12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
		12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
		12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
		12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
		12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
		12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13 }, { 46, 8,
		31, 29, 27, 24, 22, 19, 15, 13, 7, 5, 1, 61, 56, 52, 50, 43, 40, 37, 28,
		24, 12, 5, 0, 115, 107, 103, 95, 85, 82, 73, 58, 50, 26, 13, 9, 242,
		229, 212, 205, 189, 169, 167, 145, 119, 103, 55, 30, 28, 486, 483, 481,
		480, 456, 427, 409, 377, 376, 336, 332, 289, 237, 205, 109, 62, 58, 49,
		965, 964, 914, 852, 817, 675, 674, 577, 472, 217, 100, 1949, 1830, 1706,
		946, 817, 432, 252, 206, 141, 3897, 3896, 3663, 3662, 3415, 3414, 3267,
		3266, 3265, 3264, 2671, 2670, 2669, 2668, 2667, 2666, 2665, 2664, 2307,
		2306, 2305, 2304, 1895, 1894, 1639, 1638, 1637, 1636, 1633, 1632, 867,
		866, 511, 510, 509, 508, 507, 506, 479, 478, 477, 476, 475, 474, 473,
		472, 415, 414, 411, 410, 409, 408, 407, 406, 405, 391, 404, 390, 389,
		387, 388, 385, 386, 384, 285, 286, 287, 281, 284, 278, 279, 280, 274,
		275, 276, 277, 269, 270, 271, 272, 273, 261, 262, 263, 264, 265, 266,
		267, 268, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127,
		256, 257, 258, 259, 260, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75,
		76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93,
		94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108,
		109, 110, 111, 112, 113, 114, 115, 7800, 7801, 7802, 7803, 7804, 7805,
		7806, 7807 } };

const int BitWriter::error_num_codebook[2][256] = { { 1, 2, 3, 4, 5, 6, 7, 8, 9,
		10, 11, 12, 13, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 21,
		21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
		21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
		21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
		21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
		21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
		21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
		21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
		21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
		21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
		21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
		21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
		21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
		21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21 }, { 1, 1, 1, 0, 3,
		5, 9, 17, 33, 65, 129, 257, 513, 65536, 65537, 65538, 65539, 65540,
		65541, 65542, 65543, 65544, 65545, 65546, 65547, 65548, 131098, 131099,
		131100, 131101, 131102, 131103, 131104, 131105, 131106, 131107, 131108,
		131109, 131110, 131111, 131112, 131113, 131114, 131115, 131116, 131117,
		131118, 131119, 131120, 131121, 131122, 131123, 131124, 131125, 131126,
		131127, 131128, 131129, 131130, 131131, 131132, 131133, 131134, 131135,
		131136, 131137, 131138, 131139, 131140, 131141, 131142, 131143, 131144,
		131145, 131146, 131147, 131148, 131149, 131150, 131151, 131152, 131153,
		131154, 131155, 131156, 131157, 131158, 131159, 131160, 131161, 131162,
		131163, 131164, 131165, 131166, 131167, 131168, 131169, 131170, 131171,
		131172, 131173, 131174, 131175, 131176, 131177, 131178, 131179, 131180,
		131181, 131182, 131183, 131184, 131185, 131186, 131187, 131188, 131189,
		131190, 131191, 131192, 131193, 131194, 131195, 131196, 131197, 131198,
		131199, 131200, 131201, 131202, 131203, 131204, 131205, 131206, 131207,
		131208, 131209, 131210, 131211, 131212, 131213, 131214, 131215, 131216,
		131217, 131218, 131219, 131220, 131221, 131222, 131223, 131224, 131225,
		131226, 131227, 131228, 131229, 131230, 131231, 131232, 131233, 131234,
		131235, 131236, 131237, 131238, 131239, 131240, 131241, 131242, 131243,
		131244, 131245, 131246, 131247, 131248, 131249, 131250, 131251, 131252,
		131253, 131254, 131255, 131256, 131257, 131258, 131259, 131260, 131261,
		131262, 131263, 131264, 131265, 131266, 131267, 131268, 131269, 131270,
		131271, 131272, 131273, 131274, 131275, 131276, 131277, 131278, 131279,
		131280, 131281, 131282, 131283, 131284, 131285, 131286, 131287, 131288,
		131289, 131290, 131291, 131292, 131293, 131294, 131295, 131296, 131297,
		131298, 131299, 131300, 131301, 131302, 131303, 131304, 131305, 131306,
		131307, 131308, 131309, 131310, 131311, 131312, 131313, 131314, 131315,
		131316, 131317, 131318, 131319, 131320, 131321, 131322, 131323, 131324,
		131325, 131326, 131327 } };

