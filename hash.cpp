#include <experimental/string_view>

#include "tbb/concurrent_unordered_map.h"
#include "tbb/concurrent_vector.h"
#include "hash.hpp"
#include "read.hpp"

/* Constructors */
Hash::Hash(int _K) {
	K = _K;
	reads.clear();
}

Hash::~Hash() {
	//delete[] maps;
}

/* Computes matches. */
const tbb::concurrent_vector<Read*>& Hash::matches(std::experimental::string_view bases) {
	HashMap::iterator iter = maps.find(bases);
	if (iter != maps.end())
		return iter->second;
	return reads;
}

/* Adds read READ to hashmaps. */
void Hash::add(std::experimental::string_view bases, Read *read) {
	HashMap::iterator iter = maps.find(bases);
	if (iter != maps.end())
		iter->second.push_back(read);
	else {
		tbb::concurrent_vector<Read*> newreads;
		newreads.push_back(read);
		std::pair<HashMap::iterator, bool> res = maps.insert(make_pair(bases, newreads));
		if (!res.second)
			res.first->second.push_back(read);
	}
}

