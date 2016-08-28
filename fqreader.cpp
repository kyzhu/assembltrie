#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <stdlib.h>
#include <pthread.h>

#include "tbb/concurrent_vector.h"

#include "util.hpp"
#include "hash.hpp"
#include "read.hpp"
#include "compress.hpp"
#include "fqreader.hpp"

int FastaReader::L;
int FastaReader::MAX_ERRORS;

FastaReader::FastaReader(int l, int max_errors, const std::string &ifn,
		const std::string &ofn) {
	pthread_mutex_init(&thread_lock, NULL);
	L = l;
	MAX_ERRORS = max_errors;
	INFILE = ifn;
	OUTFILE = ofn;
}

void FastaReader::readFastaq(int num_threads) {
	std::string bases;
	std::ifstream inputFile;

	for (int i = 0; i < num_threads; i++) {
		std::vector<std::string> newvec;
		bases_split.push_back(newvec);
	}
	/**/
	inputFile.open(INFILE);
	while(std::getline(inputFile, bases)) 
		bases_unsplit.push_back(bases); 
	inputFile.close();

	/**/
	/*std::sort(bases_unsplit.begin(), bases_unsplit.end());
	for (int i = 0; i < 10; i++)
		printf("%s\n", bases_unsplit[i].c_str());
	printf("Sorted\n");*/

	long split_size = bases_unsplit.size() / num_threads;

	/**/
	for (int i = 0; i < num_threads - 1; i++)
		bases_split[i].insert(bases_split[i].end(), bases_unsplit.begin() + i * split_size, bases_unsplit.begin() + (i + 1) * split_size);
	//printf("Sorted\n");

	bases_split[num_threads - 1].insert(bases_split[num_threads - 1].end(), bases_unsplit.begin() + (num_threads - 1) * split_size, bases_unsplit.end());
	//printf("Sorted\n");
}

Read* FastaReader::addRead(std::string bases, Hash *ends) {
	Read *r = new Read(bases, bases);
	for (int i = L - Hash::K; i >= 0; i--) {
		std::string curBases = bases.substr(0, i + Hash::K);
		/* Check matches */
		for (Read* parent : ends->matches(bases.substr(i, Hash::K))) { 
			int diffs = diff(parent->bases(i + Hash::K), curBases);
			if (diffs > MAX_ERRORS)
				continue;

			/* Entire read matches */
			if (i + Hash::K == L && diffs == 0) {
				parent->count += 1;
				delete r;
				return NULL;
			}
			return r->connect(parent, i + Hash::K);
		}
	}
	return r;
}

void FastaReader::addChildren(Read *parent, std::string bases, Hash *starts) {
	for (int i = 1; i <= L - Hash::K; i++) {
		for (Read* child : starts->matches(bases.substr(i, Hash::K))) {
			if (child->suffix.length() < (unsigned int)i)
				continue;
			if (diff(child->bases().substr(0, L - i), bases.substr(i)) > MAX_ERRORS)
				continue;

			child->connect(parent, L - i);
		}
	}
}

void FastaReader::addSequence(Read *end, BitWriter *writer) {
	if (end->seq != -1)
		return;

	Read *read = end;

	/* Compute sequence and length */
	std::vector<Read*> reads;
	while (read != NULL && read->seq == -1) {
		read->seq = sequences;
		reads.push_back(read);
		read = read->parent;
	}

	reverse(reads.begin(), reads.end());

	dumpReads(reads, writer);
	sequences += 1;
}

char FastaReader::vote(int *votes, int row) {
	int best = 0;
	for (int col = 1; col < 5; col++) {
		if (votes[row * 5 + col] > votes[row * 5 + best])
			best = col;
	}
	return Util::intToBase(best);
}

void FastaReader::dumpReads(std::vector<Read*> reads, BitWriter *writer) {
	std::vector<int> distances; // Distances of reads, starting from position 0
	int id = 0;
	int offset = 0;
	for (size_t i = 0; i < reads.size(); i++) {
		reads[i]->seqId = id;
		distances.push_back(reads[i]->suffix.length());
		offset += reads[i]->suffix.length();
		for (int j = 1; j < reads[i]->count; j++)
			distances.push_back(0);
		id += reads[i]->count;
	}
	int length = offset; // Total length of added bases

	/* Vote on bases */
	int *votes = new int[length * 5];
	std::memset(votes, 0, length * 5);
	offset = 0;
	for (size_t i = 0; i < reads.size(); i++) {
		offset += reads[i]->suffix.length();
		std::string curBases = reads[i]->bases();
		for (int j = 0; j < L; j++)
			if (offset + j - L >= 0)
				votes[(offset + j - L) * 5 + Util::baseToInt(curBases[j])] += 1;
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
	if (start->parent == NULL) {
		start->parent = start;
	} else if (start->parent->seq == sequences) {
		fullBases = true;
	}
	Util::arraycpy(start->parent->bases().c_str(), 0, bases, 0, L);
	Util::arraycpy(newBases, 0, bases, L, length);

	/* Record errors */
	std::vector<int> errorDists;
	std::vector<Error> errors;
	int eDist = 0;
	offset = 0;
	for (size_t i = 0; i < reads.size(); i++) {
		offset += reads[i]->suffix.length();
		std::string curBases = reads[i]->bases();

		for (int tmp = 0; tmp < reads[i]->count; tmp++) {
			for (int j = 0; j < L; j++) {
				if (curBases[j] != bases[offset + j]) {
					errorDists.push_back(eDist);
					errors.push_back(Error(j, curBases[j]));
					eDist = 0;
				}
			}
			eDist += 1;
		}
	}

	if (fullBases) {
		writer->writeSequence(1, start->parent->seq, start->parent->seqId,
				bases, L + length, distances, errorDists, errors);
	} else {
		writer->writeSequence(0, start->parent->seq, start->parent->seqId,
				newBases, length, distances, errorDists, errors);
	}

	delete[] votes;
	delete[] bases;
	delete[] newBases;
}

int FastaReader::diff(std::string s1, std::string s2) {
	int count = 0;
	for (size_t i = 0; i < s1.length(); i++)
		if (s1[i] != s2[i]) {
			count += 1;
			if (count > MAX_ERRORS)
				break;
		}
	return count;
}

unsigned long FastaReader::currentTimeMillis() {
	return std::chrono::duration_cast < std::chrono::milliseconds
			> (std::chrono::system_clock::now().time_since_epoch()).count();
}

void* FastaReader::_thread() {
	int i = 0;
	pthread_mutex_lock(&thread_lock);
	for (i = 0; i < num_threads; i++)
		if (!processed[i]) {
			processed[i] = true;
			break;
		}
	pthread_mutex_unlock(&thread_lock);
	printf("Begin to process bucket %d by thread %lu.\n", i, pthread_self());
	Hash *starts = new Hash();
	Hash *ends = new Hash();
	int count = 0;
	for (std::string bases : bases_split[i]) {
		Read *read = addRead(bases, ends);		
		if (read == NULL)
			continue;
		read->ID = read_id++;
		
		/* Find new children */
		addChildren(read, bases, starts);

		/* Add hashes, record read */
		ends->add(bases.substr(L - Hash::K), read);
		starts->add(bases.substr(0, Hash::K), read);		
		reads.push_back(read);

		count++;
		if (count % 100000 == 0) {
			printf("Processed %d reads by thread %lu.\n", count, pthread_self());
		}
	}
	delete starts;
	delete ends;
	return 0;
}

void FastaReader::run(int _num_threads) {
	num_threads = _num_threads;
	readFastaq(num_threads);

	processed = new bool[num_threads];
	for (int i = 0; i < num_threads; i++)
		processed[i] = false;

	/* Initialize hashmaps */
	unsigned long elapsed_time = currentTimeMillis();

	/* Process the reads and construct tries */
	pthread_t threads[num_threads];
	for (int i = 0; i < num_threads; i++)
		pthread_create(&threads[i], NULL, thread, this);
	for(int i = 0; i < num_threads; i++)
		pthread_join(threads[i], NULL);

	std::printf("Processed reads in %lu ms\n", currentTimeMillis() - elapsed_time);
	elapsed_time = currentTimeMillis();

	/* Initialize the encoder */
	BitWriter *writer = new BitWriter(OUTFILE);

	/* Count children */
	std::printf("%lu unique reads\n", reads.size());
	double errors = 0;
	double suffix = 0;
	for (size_t i = 0; i < reads.size(); i++) {
		if (reads[i]->parent != NULL)
			reads[i]->parent->children += 1;
		errors += reads[i]->errors.size();
		suffix += reads[i]->suffix.length();
	}
	std::printf("Counted children\n");
	std::printf("%.2f errors, %.2f suffix\n", errors / reads.size(),
			suffix / reads.size());

	/* First add leaves */
	for (size_t i = 0; i < reads.size(); i++)
		if (reads[i]->children == 0) {
			addSequence(reads[i], writer);
		}
	std::printf("Added %d sequences\n", sequences);

	/* Add cycle if necessary */
	for (size_t i = 0; i < reads.size(); i++)
		addSequence(reads[i], writer);
	writer->finish();

	std::printf("%d total sequences\n", sequences);

	int total = 0;
	for (Read *read : reads)
		total += read->suffix.length();
	std::printf("%d base storage\n", total);

	std::printf("Encode tries in %lu ms\n", currentTimeMillis() - elapsed_time);

	for (size_t i = 0; i < reads.size(); i++)
		if (reads[i]) delete reads[i];
	delete writer;
}

