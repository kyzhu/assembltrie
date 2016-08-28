#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "compress.hpp"
#include "error.hpp"
#include "util.hpp"

BitWriter::BitWriter(const std::string &filename) {
	try {
		stream.open(filename.c_str(), std::ios::out | std::ios::binary);
	} catch (std::ofstream::failure e) {
		std::cerr << "Cannot open file \'" << filename << "\'.\n";
	}
}

void BitWriter::writeSequence(int full, int parentSeq, int parentId,
		char *bases, int length, std::vector<int> distances,
		std::vector<int> errorDists, std::vector<Error> errors) {
	/* Write lengths */
	writeTwo(full);
	writeInt(parentSeq);
	writeInt(parentId);
	writeInt(length);
	writeInt(distances.size());
	writeInt(errorDists.size()); // Same as errors.size()

	/* Write sequence */
	// System.out.printf("Bases: %s%n", new String(bases));
	for (int i = 0; i < length; i++)
		writeBase(bases[i]);

	/* Write read distances */
	for (int distance : distances)
		writeDistance(distance);

	/* Write errors */
	for (size_t i = 0; i < errorDists.size(); i++) {
		/* Write distance */
		writeDistance(errorDists[i]);
		/* Write error */
		Error e = errors[i];
		writeBits(7, e.location);
		writeBase(e.base);
	}
}

/* Returns the lowest COUNT bits of VALUE. */
int BitWriter::lowBits(int count, long value) {
	return (int) ((int) (value & BitWriter::MAXINT) & ((1L << count) - 1));
}

/* Writes the lowest COUNT bits of VALUE to file. */
void BitWriter::writeBits(int count, int value) {
	/* Write next byte when available */
	try {
		while (curBits + count >= BYTE) {
			int toWrite = BYTE - curBits;
			curByte = (curByte << toWrite) + lowBits(toWrite, value);
			value >>= toWrite;
			count -= toWrite;

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
		writeBits(1, 1);
		writeBits(1, 0);
		writeBits(8, value);
	} else if (value < 65536) {
		writeBits(1, 1);
		writeBits(1, 1);
		writeBits(1, 0);
		writeBits(16, value);
	} else {
		writeBits(1, 1);
		writeBits(1, 1);
		writeBits(1, 1);
		writeBits(32, value);
	}
}

/* Writes two bits to file. */
void BitWriter::writeTwo(int value) {
	writeBits(2, value);
}

/* Writes a base C to file. */
void BitWriter::writeBase(char c) {
	//writeTwo(Util::baseToInt(c));
	writeBits(3, Util::baseToInt(c));
}

/* Writes distances */
void BitWriter::writeDistance(int distance) {
	while (distance >= 3) {
		distance -= 3;
		writeTwo(3);
	}
	writeTwo(distance);
}

/* Finishes writing to file. */
void BitWriter::finish() {
	/* Write terminator and padding bits */
	writeInt(-1);
	writeBits(BYTE, 0);

	/* Close file. */
	try {
		stream.close();
	} catch (std::ofstream::failure e) {
		std::cerr << "Exception closing file.\n";
	}
}

/*void BitWriter::print(Object x) {
 std::cout << x << " " << std::endl;
 }

 static string BitWriter::join(List<?> x) {
 StringBuilder s = new StringBuilder();
 for (Object o : x) {
 s.append(o);
 s.append(" ");
 }
 return s.toString();
 }*/

