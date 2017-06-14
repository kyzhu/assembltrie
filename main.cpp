#include <iostream>
#include <unistd.h>

#include "fqreader.hpp"
#include "decompress.hpp"

int main(int argc, char *argv[]) {
	/* Command line option */
	int opt;

	/* Read length. */
	int read_length = 100;

	int hash_length = 20;

	int satisfactory_suffix_length = read_length / 10;

	/* Maximum allowed errors. */
	int max_errors = 3;

	/* Input and output file names */
	char *inputfile = NULL;
	char *outputfile = NULL;

	/* Number of threads */	
	int num_threads = 8;

	/* Mode: 0 - compression; 1 - decompression */
	int mode = 0;

	/* Whether post-process */
	int postprocess = 1;

	/* Correct the orientation of the reads */
	int rc_heuristic = 0;

	while ((opt = getopt(argc, argv, "cdL:K:h:s:i:o:e:n:p")) != -1)
		switch (opt) {
      			case 'c':
        				mode = 0;
        				break;
      			case 'd':
        				mode = 1;
        				break;
			case 'p':
					postprocess = 1;
					break;
      			case 'L':
        				read_length = atoi(optarg);
					satisfactory_suffix_length = read_length / 10;
        				break;
			case 'K':
					hash_length = atoi(optarg);
					break;
			case 'h':
					rc_heuristic = atoi(optarg);
					break;
			case 's':
					satisfactory_suffix_length = atoi(optarg);
					break;
			case 'i':
        				inputfile = optarg;
        				break;
			case 'o':
        				outputfile = optarg;
        				break;
			case 'e':
        				max_errors = atoi(optarg);
        				break;
			case 'n':
        				num_threads = atoi(optarg);
        				break;
      			default:
        				abort();
      		}
	
	if (mode == 0) {
		/* Test Compression */
		FastaReader encoder(read_length, hash_length, max_errors, satisfactory_suffix_length,
				   inputfile, outputfile, rc_heuristic, postprocess);
		encoder.run(num_threads);
	}
	else {
		/* Test Decompression */
		BitReader decoder(read_length, rc_heuristic, inputfile);
		decoder.wrtReads(outputfile);
	}

	return 0;
}

