#include <iostream>
#include <fstream>
#include <string>
#include <experimental/string_view>
#include <vector>
#include <algorithm>
#include <chrono>
#include <stdlib.h>
#include <pthread.h>
#include <unordered_set>

#include "tbb/concurrent_vector.h"
#include "util.hpp"
#include "hash.hpp"
#include "read.hpp"
#include "compress.hpp"
#include "fqreader.hpp"

int FastaReader::L;
int FastaReader::K;
int FastaReader::MAX_ERRORS;
int FastaReader::heuristic;

FastaReader::FastaReader(int l, int k, int max_errors, int s,
		const std::string &ifn, const std::string &ofn, int h, int p) {
	pthread_mutex_init(&thread_lock, NULL);
	L = l;
	K = k;
	MAX_ERRORS = max_errors;
	POST_PROCESS = p;
	S = s;
	INFILE = ifn;
	OUTFILE = ofn;
	heuristic = h;
}

void FastaReader::readFastaq(int num_threads) {
	std::string bases;
	std::ifstream inputFile;

	/* Initialize buffer */
	for (int i = 0; i < num_threads; i++) {
		std::vector<std::string> newbases;
		bases_split.push_back(newbases);
	}

	/* Read in file */
	inputFile.open(INFILE);
	while (std::getline(inputFile, bases)) {
		std::getline(inputFile, bases);
		//std::replace(bases.begin(), bases.end(), 'N', 'A');
		bases_unsplit.push_back(bases);
		std::getline(inputFile, bases);
		std::getline(inputFile, bases);
	}
	inputFile.close();
	std::random_shuffle(bases_unsplit.begin(), bases_unsplit.end());

	/* Randomly partition the reads into different buckets */
	split_size = bases_unsplit.size() / num_threads;
	for (int i = 0; i < num_threads - 1; i++)
		bases_split[i].insert(bases_split[i].end(),
				bases_unsplit.begin() + i * split_size,
				bases_unsplit.begin() + (i + 1) * split_size);
	bases_split[num_threads - 1].insert(bases_split[num_threads - 1].end(),
			bases_unsplit.begin() + (num_threads - 1) * split_size,
			bases_unsplit.end());
}

void FastaReader::readFastaqRC(int num_threads) {
	std::string bases;
	std::ifstream inputFile;

	/* Initialize buffer */
	for (int i = 0; i < num_threads; i++) {
		std::vector<std::string> newbases;
		bases_split.push_back(newbases);
	}

	/* Read in file */
	inputFile.open(INFILE);
	if (heuristic == 1) // Heuristic: ORCOM
		while (std::getline(inputFile, bases)) {
			std::getline(inputFile, bases);
			//std::replace(bases.begin(), bases.end(), 'N', 'A');
			std::string rc_bases = Util::getRCstring(bases);

			if (Util::getMin(bases) <= Util::getMin(rc_bases)) {
				bases_unsplit.push_back(bases);
				rc_indicator.push_back(0);
			} else {
				bases_unsplit.push_back(rc_bases);
				rc_indicator.push_back(1);
				test++;
			}

			std::getline(inputFile, bases);
			std::getline(inputFile, bases);
		}
	else
		// Heuristic: Minimum AC weight
		while (std::getline(inputFile, bases)) {
			std::getline(inputFile, bases);
			std::string rc_bases = Util::getRCstring(bases);

			if (Util::getACweight(bases) <= 0) {
				bases_unsplit.push_back(bases);
				rc_indicator.push_back(0);
			} else {
				bases_unsplit.push_back(rc_bases);
				rc_indicator.push_back(1);
			}

			std::getline(inputFile, bases);
			std::getline(inputFile, bases);
		}
	inputFile.close();
	std::random_shuffle(bases_unsplit.begin(), bases_unsplit.end());	

	/* Randomly partition the reads into different buckets */
	split_size = bases_unsplit.size() / num_threads;
	for (int i = 0; i < num_threads - 1; i++)
		bases_split[i].insert(bases_split[i].end(),
				bases_unsplit.begin() + i * split_size,
				bases_unsplit.begin() + (i + 1) * split_size);
	bases_split[num_threads - 1].insert(bases_split[num_threads - 1].end(),
			bases_unsplit.begin() + (num_threads - 1) * split_size,
			bases_unsplit.end());
}

Read* FastaReader::addRead(std::experimental::string_view bases) {
	Read *r = new Read(bases, L);
	for (int i = L - K; i >= 0; i--) {
		/* Check matches */
		for (Read* parent : ends->matches(bases.substr(i, K))) {
			if (!parent)
				continue;

			int diffs = diff(parent->bases(i + K), bases.substr(0, i + K));
			if (diffs > MAX_ERRORS)
				continue;

			/* Entire read matches */
			if (i + K == L && diffs == 0) {
				parent->count += 1;
				delete r;
				return NULL;
			}
			return r->connect(parent, i + K);
		}
	}
	return r;
}

Read* FastaReader::addRead(std::experimental::string_view bases, int ort) {
	Read *r = new Read(bases, L, ort);
	for (int i = L - K; i >= 0; i--) {
		/* Check matches */
		for (Read* parent : ends->matches(bases.substr(i, K))) {
			if (!parent)
				continue;

			int diffs = diff(parent->bases(i + K), bases.substr(0, i + K));
			if (diffs > MAX_ERRORS)
				continue;

			/* Entire read matches */
			if (i + K == L && diffs == 0) {
				parent->rc.push_back(ort);
				parent->count += 1;
				delete r;
				return NULL;
			}
			return r->connect(parent, i + K);
		}
	}
	return r;
}

void FastaReader::addChildren(Read *parent,
		std::experimental::string_view bases) {
	for (int i = 1; i <= L - K; i++)
		for (Read* child : starts->matches(bases.substr(i, K))) {
			if (!child)
				continue;

			if (child->getSuffixlength() <= i)
				continue;
			if (diff(child->bases().substr(0, L - i), bases.substr(i))
					> MAX_ERRORS)
				continue;

			child->connect(parent, L - i);
		}
}

void FastaReader::postProcess(Read *read) {
	while (read->hasChildren()) {
		if (read->merged)
			return;

		/* Try Merging the children with longest suffices */
		std::vector<Read*> children;
		for (int id : read->children) {
			children.push_back(reads[id]);
		}

		int size = children.size();
		Read *next = NULL;
		for (int i = 0; i < size - 1; i++) {
			/* children with longest suffix */
			std::sort(children.begin(), children.end(), suff_order);
			Read *rightRead = children[i];
			Read *leftRead = children[i + 1];
			next = leftRead;

			int x = leftRead->suffix;
			int y = rightRead->suffix;
			int d = diff(leftRead->bases(y - x, L),
					rightRead->bases(0, L + x - y - 1));
			if (d <= MAX_ERRORS) {
				read->disconnect(rightRead);
				rightRead->connect(leftRead, L + x - y);
				leftRead->children.push_back(rightRead->ID);
			}
		}

		read->merged = true;
		if (!next)
			return;
		read = next;
	}
}

void FastaReader::addSequence(Read *end, BitWriter *writer) {
	if (end->seq != -1)
		return;

	if (end->parent != NULL || !end->children.empty()) {
		Read *read = end;
		/* Compute sequence and length */
		std::vector<Read*> reads;
		while (read != NULL && read->seq == -1) {
			read->seq = sequences;
			reads.push_back(read);
			read = read->parent;
		}
		reverse(reads.begin(), reads.end());

		writer->dumpReads(reads, sequences);
		sequences += 1;
	}
}

int FastaReader::diff(std::experimental::string_view s1,
		std::experimental::string_view s2) {
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
	/* Find the bucket of reads to process */
	int i = 0;
	pthread_mutex_lock (&thread_lock);
	i = current_bucket++;
	pthread_mutex_unlock(&thread_lock);
	if (print_info)
		printf("Begin to process bucket %d by thread u %lu.\n", i + 1,
			pthread_self());
	for (size_t j = 0; j < bases_split[i].size(); j++) {
		std::experimental::string_view sv_bases(bases_split[i][j]);
		Read *read = addRead(sv_bases);
		if (read == NULL)
			continue;

		/* Find new children */
		addChildren(read, sv_bases);

		/* Add hashes, record read */
		ends->add(sv_bases.substr(L - K), read);
		if (read->suffix >= S)
			starts->add(sv_bases.substr(0, K), read);
		reads.push_back(read);
	}
	if (print_info)
		printf("Read Connection completed by thread %lu.\n", pthread_self());
	return 0;
}

void* FastaReader::_thread_orcom() {
	/* Find the bucket of reads to process */
	int i = 0;
	pthread_mutex_lock (&thread_lock);
	i = current_bucket++;
	pthread_mutex_unlock(&thread_lock);
	if (print_info)
		printf("Begin to process bucket %d by thread o %lu.\n", i + 1,
			pthread_self());
	for (size_t j = 0; j < bases_split[i].size(); j++) {
		std::experimental::string_view sv_bases(bases_split[i][j]);
		Read *read = addRead(sv_bases, rc_indicator[i * split_size + j]);
		if (read == NULL)
			continue;

		/* Find new children */
		addChildren(read, sv_bases);

		/* Add hashes, record read */
		ends->add(sv_bases.substr(L - K), read);
		if (read->suffix >= S)
			starts->add(sv_bases.substr(0, K), read);
		reads.push_back(read);
	}
	if (print_info)
		printf("Read Connection completed by thread %lu.\n", pthread_self());
	return 0;
}

void* FastaReader::_thread_shortsuf() {
	/* Find the bucket of reads to process */
	int i = 0;
	pthread_mutex_lock (&thread_lock);
	i = current_bucket++;
	pthread_mutex_unlock(&thread_lock);
	if (print_info)
		printf("Begin to process bucket %d by thread s %lu.\n", i + 1,
			pthread_self());
	for (size_t j = 0; j < bases_split[i].size(); j++) {
		std::string rc_bases = Util::getRCstring(bases_split[i][j]);
		std::experimental::string_view sv_bases(bases_split[i][j]);
		std::experimental::string_view sv_rc_bases(rc_bases);
		Read *read = addRead(sv_bases, 0);
		if (read == NULL)
			continue;
		Read *rc_read = addRead(sv_rc_bases, 1);
		if (rc_read == NULL) {
			delete read;
			continue;
		}
		if (read->suffix <= rc_read->suffix) {
			/* Find new children */
			addChildren(read, sv_bases);

			/* Add hashes, record read */
			ends->add(sv_bases.substr(L - K), read);
			if (read->suffix >= S)
				starts->add(sv_bases.substr(0, K), read);
			reads.push_back(read);
			delete rc_read;
		} else {
			bases_split[i][j] = rc_bases;
			sv_bases = std::experimental::string_view(bases_split[i][j]);
			rc_read->_bases = sv_bases;

			addChildren(rc_read, sv_bases);
			ends->add(sv_bases.substr(L - K), rc_read);
			if (read->suffix >= S)
				starts->add(sv_bases.substr(0, K), rc_read);
			reads.push_back(rc_read);
			delete read;
		}
	}
	if (print_info)
		printf("Read Connection completed by thread %lu.\n", pthread_self());
	return 0;
}

void FastaReader::run(int _num_threads) {
	num_threads = _num_threads;
	if (heuristic == 0 || heuristic == 2)
		readFastaq (num_threads);
	else
		readFastaqRC(num_threads);

	/* Initialize hashmaps */
	starts = new Hash(K);
	ends = new Hash(K);

	/* Process the reads and construct tries */
	pthread_mutex_init(&thread_lock, NULL);
	pthread_spin_init(&read_access, 0);

	pthread_t threads[num_threads];
	for (int i = 0; i < num_threads; i++)
		pthread_create(&threads[i], NULL, thread, this);
	for (int i = 0; i < num_threads; i++)
		pthread_join(threads[i], NULL);

	/* Post processing. */
	for (size_t i = 0; i < reads.size(); i++) {
		reads[i]->ID = i;
		if (reads[i]->parent != NULL)
			reads[i]->parent->children.push_back(i);
	}
	if (POST_PROCESS == 1) {
		if (print_info)
			printf("Post process begins. \n");
		for (size_t i = 0; i < reads.size(); i++)
			postProcess(reads[i]);
		if (print_info)
		printf("Post processing completed. \n");
	}

	/* Initialize the encoder */
	BitWriter *writer = new BitWriter(OUTFILE, L);
	BitWriter *htcWriter = new BitWriter("part.out", L);
	if (heuristic > 0) {
		writer->heuristic_flag = true;
		htcWriter->heuristic_flag = true;
	}
	if (print_info) {
		/* Count children */
		double errors = 0;
		double suffix = 0;
		for (size_t i = 0; i < reads.size(); i++) {
			errors += reads[i]->errors.size();
			suffix += reads[i]->suffix;
		}
		printf("%lu unique reads\n", reads.size());
		printf("Counted children\n");
		printf("%.2f errors, %.2f suffix\n", errors / reads.size(),
			suffix / reads.size());
	}

	/* First add leaves */
	for (size_t i = 0; i < reads.size(); i++)
		if (reads[i]->children.empty()) {
			if (reads[i]->parent != NULL)
				addSequence(reads[i], writer);
			else {
				htcBits.insert(htcBits.end(), reads[i]->rc.begin(),
						reads[i]->rc.end());
				for (int k = 0; k < reads[i]->count; k++)
					htcReads.push_back(reads[i]->bases());
			}
		}
	if (print_info)
		printf("Added %d sequences\n", sequences);

	/* Add cycle if necessary */
	for (size_t i = 0; i < reads.size(); i++)
		addSequence(reads[i], writer);
	writer->writeLastSequence();
	writer->sendNstream(htcWriter);
	writer->finish();
	if (print_info) {
		printf("%d total sequences\n", sequences);
		unsigned long total = 0;
		for (size_t i = 0; i < reads.size(); i++)
			total += reads[i]->suffix;
		printf("%lu base storage\n", total);
	}

	/* Compressing a separating stream. */
	htcWriter->writeNstream(true);
	htcWriter->writeUncompressedReads(htcReads, htcBits);
	htcWriter->writeNstream(false);
	htcWriter->finish();
	delete writer;
	delete htcWriter;
}

