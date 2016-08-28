#include <vector>
#include <string>
#include <iostream>

#include "read.hpp"
#include "error.hpp"
#include "fqreader.hpp"

Read::Read(std::string s) {
	this->suffix = s;
}

Read::Read(std::string s, std::string b) {
	this->suffix = s;
	this->_bases = b;
}

std::string Read::bases() {
	return bases(FastaReader::getL());
}

Read* Read::connect(Read *newParent, int length) {
	std::string pBases = newParent->bases(length);
	std::string bases = Read::bases();
	
	parent = newParent;
	suffix = bases.substr(length, bases.size() - length + 1);
	if (!errors.empty()) errors.clear();
	for (int i = 0; i < length; i++)
		if (pBases.at(i) != bases.at(i))
			errors.push_back(Error(i, bases.at(i)));
	return this;
}

std::string Read::bases(Read *r, size_t l, bool useErrors) {
	if(r->count == 0) printf("Bug");
	return r->_bases.substr(FastaReader::getL() - l, l);
}


