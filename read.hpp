#ifndef READ_HPP_
#define READ_HPP_

#include <vector>
#include <string>
#include <stdlib.h>

#include "error.hpp"

class Read{
public:
	Read *parent = NULL;
	std::string suffix;
	std::string _bases;

	//bool connected = false;
	int ID = -1;
	int count = 1; // Number of times this read appeared in the input
	int seq = -1; // Sequence containing this read, or -1 if none
	int seqId = -1; // ID in sequence, or -1 if none
	int children = 0; // Number of children
	std::vector<Error> errors;

	/* Constructor */
	Read(std::string);
	Read(std::string, std::string);
	~Read(){
                //if (parent != NULL) delete parent;
        }

	/* Sets my parent to PARENT (with overlap LENGTH) and return myself. */
	Read* connect(Read*, int);

	/* Last L bases of read R. */
	static std::string bases(Read *r, size_t l) {
		return bases(r, l, true);
	}

	/* Last L bases of read R. */
	static std::string bases(Read*, size_t, bool);

	/* Last L bases. */
	std::string bases(size_t l) {
		return bases(this, l);
	}

	/* Last L bases. */
	std::string bases(size_t l, bool e) {
		return bases(this, l, e);
	}

	/* All bases. */
	std::string bases();

	std::string toString() {
		return bases();
	}

};

#endif