#include <string>
//#include <cstring>
#include <vector>
//#include <algorithm>
//#include <time.h>
//#include <stdlib.h>
#include <iostream>
#include <pthread.h>

#include "tbb/concurrent_unordered_map.h"
#include "tbb/concurrent_vector.h"
//#include "tbb/atomic.h"
#include "hash.hpp"
#include "read.hpp"

/* Constructors */
Hash::Hash() {
	pthread_rwlock_init(&vector_update, NULL);
}

Hash::~Hash() {
	//delete[] maps;
}

/* Computes matches */
tbb::concurrent_vector<Read*> Hash::matches(std::string bases) {
	tbb::concurrent_vector<Read*> reads;
	pthread_rwlock_rdlock(&vector_update);
	HashMap::iterator iter = maps.find(bases);
	if (iter != maps.end())
		reads = iter->second;
		//std::this_thread::yield();
	pthread_rwlock_unlock(&vector_update);
	return reads;
}

/* Adds read READ to hashmaps */
void Hash::add(std::string bases, Read *read) {
	HashMap::iterator iter = maps.find(bases);
	if (iter != maps.end()) {
		pthread_rwlock_wrlock(&vector_update);
		iter->second.push_back(read);
		pthread_rwlock_unlock(&vector_update);
	}
	else {
		tbb::concurrent_vector<Read*> newreads;
		newreads.push_back(read);
		pthread_rwlock_wrlock(&vector_update);
		std::pair<HashMap::iterator, bool> ret = maps.insert(make_pair(bases, newreads));
		if (!ret.second) ret.first->second.push_back(read);
		pthread_rwlock_unlock(&vector_update);
	}
}

