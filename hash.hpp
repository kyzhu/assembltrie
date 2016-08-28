#ifndef HASH_HPP_
#define HASH_HPP_

#include <vector>
#include <string>
#include <pthread.h>

#include "tbb/concurrent_unordered_map.h"
#include "tbb/concurrent_vector.h"
#include "read.hpp"

typedef tbb::concurrent_unordered_map<std::string, tbb::concurrent_vector<Read*>> HashMap;

class Hash {
public:
	/* Number of independent hashes */
	static const int D = 4;

	/* Number of indices to hash */
	//static const int M = 20;

	/* Length of strings to hash */
	static const int K = 20;

	/* Indices to hash */
	//int indices[D][M];

	pthread_rwlock_t vector_update;

	/* Hashmaps */
	HashMap maps;

	/* Constructor */
	Hash();
	//Hash(const std::string &filename);
	~Hash();

	/* Computes the H-th hash of string S */
	//int hash(int, const char*);

	/* Computes matches */
	tbb::concurrent_vector<Read*> matches(std::string);

	/* Adds read READ to hashmaps */
	void add(std::string, Read*);
};

#endif

