#ifndef READ_HPP_
#define READ_HPP_

#include <vector>
#include <experimental/string_view>
#include <stdlib.h>
#include <pthread.h>

#include "error.hpp"

class Read{
public:
	/* Parent and children */
	Read *parent = NULL;
	std::vector<int> children;

	/* Bases */
	std::experimental::string_view _bases;
	int suffix;

	/* Basic information */
	int ID = -1; // ID
	std::vector<int> rc;
	int count = 1; // Number of times this read appeared in the input
	int trieID = -1; // ID of the trie containing this read
	int seq = -1; // ID of the sequence containing this read, or -1 if none
	int seqId = -1; // ID in sequence, or -1 if none
	bool merged = false; // For post-process

	pthread_spinlock_t access;

	/* Errors */
	std::vector<Error> errors;

	/* Constructors */
	Read(std::experimental::string_view);
	Read(std::experimental::string_view, int);
	Read(std::experimental::string_view, int, int);
	~Read(){
                //if (parent != NULL) delete parent;
        }

	int getSuffixlength();

	/* Sets my parent to PARENT (with overlap LENGTH) and return myself. */
	Read* connect(Read*, int);

	/* Disconnect a child read and return myself. */
	Read* disconnect(Read*);

	/* Test whether I am the ancestor of another read. */
	bool isAncestor(Read*);

	/* Test whether a read has children or not. */
	bool hasChildren() {
		return (!children.empty());
	}

	/* Get all bases. */
	std::experimental::string_view bases();

	/* Get last L bases. */
	std::experimental::string_view bases(size_t l);

	/* Get bases of length L from BEGIN_INDEX (inclusive) */
	std::experimental::string_view bases(size_t begin_index, size_t l) {
		return _bases.substr(begin_index, l);
	}

	std::experimental::string_view toString() {
		return bases();
	}

};

#endif

