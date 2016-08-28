#ifndef ERROR_HPP_
#define ERROR_HPP_

class Error {
public:
	int location;
	char base;

	Error(int location, char base) {
		this->location = location;
		this->base = base;
	}

	~Error(){};
};

#endif
