#include <vector>
#include <algorithm>
#include <experimental/string_view>
#include <iostream>
#include <pthread.h>

#include "read.hpp"
#include "error.hpp"
#include "fqreader.hpp"

Read::Read(std::experimental::string_view s) {
	suffix = FastaReader::getL();
	_bases = s;
	pthread_spin_init(&access, 0); 
}

Read::Read(std::experimental::string_view b, int suflength) {
	suffix = suflength;
	_bases = b;
	pthread_spin_init(&access, 0);
}

Read::Read(std::experimental::string_view b, int suflength, int orientation) {
	suffix = suflength;
	_bases = b;
	rc.push_back(orientation);
	pthread_spin_init(&access, 0);
}

int Read::getSuffixlength() {
	pthread_spin_lock(&access);
	int res = suffix;
	pthread_spin_unlock(&access);
	return res;
}

std::experimental::string_view Read::bases() {
	return bases(FastaReader::getL());
}

std::experimental::string_view Read::bases(size_t l) {
	return _bases.substr(FastaReader::getL() - l, l);
}

Read* Read::connect(Read *newParent, int length) {
	std::experimental::string_view pBases = newParent->bases(length);
	pthread_spin_lock(&access);
	parent = newParent;
	suffix = FastaReader::getL() - length;
	
	if (!errors.empty()) errors.clear();
	for (int i = 0; i < length; i++)
		if (pBases.at(i) != _bases.at(i))
			errors.push_back(Error(i, _bases.at(i)));
	pthread_spin_unlock(&access);
	return this;
}

Read* Read::disconnect(Read *c) {
	/* Remove the child from parent */
	children.erase(std::remove(children.begin(), children.end(), c->ID), children.end());

	/* Reset child */
	c->parent = NULL;
	c->suffix = FastaReader::getL();
	c->errors.clear();

	return this;
}

bool Read::isAncestor(Read *dec) {
	while(dec) {
		if (!dec->parent) return false;
		if (dec->parent->ID == ID)
			return true;
		dec = dec->parent;
        }
        return false;
}

