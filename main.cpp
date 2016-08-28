#include <iostream>

#include "fqreader.hpp"
#include "decompress.hpp"

int main() {
	/* Test Compression */
	FastaReader writer(101, 4, "chr2.fa", "chr2.out");
	writer.run(8);

	/* Test Decompression */
	BitReader reader("chr2.out");
	reader.wrtReads();

	return 0;
}

