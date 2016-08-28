#include <string>
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
	static int MAX_ERRORS;

public:
	/* Input and output file names */
	std::string INFILE;
	std::string OUTFILE;

	/* Shared variables */
	//Hash starts, ends;
	tbb::concurrent_vector<Read*> reads;
	int sequences = 0;
	int read_id = 0;

	/* Thread control */
	//pthread_mutex_t if_mutex;
	//pthread_spinlock_t read_access;
	int num_threads = 1;
	pthread_mutex_t thread_lock;
	bool *processed;
	std::vector<std::string> bases_unsplit;
	std::vector<std::vector<std::string>> bases_split; 

	FastaReader(int, int, const std::string &ifn, const std::string &ofn);
	~FastaReader() {
		//for (size_t i = 0; i < reads.size(); i++)
		//	delete reads[i];
	}

	static int getL() {
		return L;
	}

	void readFastaq(int);

	/* Returns a new read with bases BASES and correct parent, or null if already exists. */
	Read* addRead(std::string, Hash*);

	/* Attaches appropriate children to a newly-created read. */
	void addChildren(Read*, std::string, Hash*);

	/* Adds a read to a sequence if not already in one, and dumps sequence to file. */
	void addSequence(Read*, BitWriter*);

	/* Returns base with the most votes. */
	char vote(int*, int);

	/* Dumps READS to file. */
	void dumpReads(std::vector<Read*>, BitWriter*);

	/* Returns the Hamming distance between two strings. */
	int diff(std::string, std::string);

	/* Returns the current time (in millisecond) */
	unsigned long currentTimeMillis();

	/* Compress the fasta(q) file */
	void *_thread();
	static void *thread(void *obj) {
		return ((FastaReader*)obj)->_thread();
	}
	void run(int);

};
