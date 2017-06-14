#ifndef FQREADER_HPP_
#define FQREADER_HPP_

#include <string>
#include <experimental/string_view>
#include <vector>
#include <pthread.h>

#include "tbb/concurrent_vector.h"

#include "hash.hpp"
#include "read.hpp"
#include "compress.hpp"

class FastaReader {
private:
	/* Parameters */
	static int L;
	static int K;
	static int MAX_ERRORS;

	/* Heuristic to correct read orientation */
	static int heuristic;

public:
	bool print_info = false;

	/* If the suffix length of a read is less than S, then it will not be considered in addChildren. */
	int S;

	/* Whether post-process. */
	int POST_PROCESS = 0;

	/* Input and output file names */
	std::string INFILE;
	std::string OUTFILE;

	/* Shared variables */
	Hash *starts;
	Hash *ends;
	tbb::concurrent_vector<Read*> reads;
	int sequences = 0;
	int current_bucket = 0;
	std::vector<std::experimental::string_view> htcReads;
	std::vector<int> htcBits;

	int test = 0;

	/* Thread control */
	int num_threads = 1;
	pthread_mutex_t thread_lock;
	pthread_spinlock_t read_access;
	std::vector<std::string> bases_unsplit;
	std::vector<std::vector<std::string>> bases_split;
	long split_size;
	std::vector<int> rc_indicator;

	FastaReader(int, int, int, int, const std::string &ifn,
			const std::string &ofn, int, int);
	~FastaReader() {
		delete starts; 
		delete ends;
		for (size_t i = 0; i < reads.size(); i++)
			if (reads[i])
				delete reads[i];
	}

	/* For encoder. */
	static int getL() {
		return L;
	}

	/* Read in fastq files. */
	void readFastaq(int);
	void readFastaqRC(int);

	/* Returns a new read with bases BASES and correct parent, or null if already exists. */
	Read* addRead(std::experimental::string_view);
	Read* addRead(std::experimental::string_view, int);

	/* Attaches appropriate children to a newly-created read. */
	void addChildren(Read*, std::experimental::string_view);

	/* Postprocess reads such that the width of tries is reduced.*/
	static bool suff_order(Read *r, Read *s) {
		return (r->suffix > s->suffix);
	}
	void postProcess(Read*);

	/* Adds a read to a sequence if not already in one, and dumps sequence to file. */
	void addSequence(Read*, BitWriter*);

	/* Returns base with the most votes. */
	char vote(int*, int);

	/* Returns the Hamming distance between two strings. */
	int diff(std::experimental::string_view, std::experimental::string_view);

	/* Returns the current time (in millisecond) */
	unsigned long currentTimeMillis();

	/* Compress the fasta(q) file */
	void *_thread();
	void *_thread_orcom();
	void *_thread_shortsuf();

	static void *thread(void *obj) {
		switch (heuristic) {
		case 0:
			return ((FastaReader*) obj)->_thread();
		case 1:
			return ((FastaReader*) obj)->_thread_orcom();
		case 2:
			return ((FastaReader*) obj)->_thread_shortsuf();
		default:
			return ((FastaReader*) obj)->_thread();
		}
	}

	void run(int);
};

#endif

