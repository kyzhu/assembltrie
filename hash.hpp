#ifndef HASH_HPP_
#define HASH_HPP_

#include <experimental/string_view>

#include "tbb/concurrent_unordered_map.h"
#include "tbb/concurrent_vector.h"
#include "read.hpp"

typedef tbb::concurrent_unordered_map<std::experimental::string_view, tbb::concurrent_vector<Read*>,
	std::hash<std::experimental::string_view>, std::equal_to<std::experimental::string_view>> HashMap;

class Hash {
public:
	/* Length of strings to hash. */
	int K = 20;

	tbb::concurrent_vector<Read*> reads;

	/* Hashmaps */
	HashMap maps;

	/* Constructors */
	Hash(int);
	~Hash();

	/* Computes matches. */
	const tbb::concurrent_vector<Read*>& matches(std::experimental::string_view);

	/* Adds read READ to hashmaps. */
	void add(std::experimental::string_view, Read*);
};

#endif

