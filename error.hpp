#ifndef ERROR_HPP_
#define ERROR_HPP_

class Error {
public:
	int location;
	char underlying;
	char base;

	Error(int location, char base) {
		this->location = location;
		this->base = base;
		this->underlying = 'X';
	}

	Error(int location, char base, char underlying) {
		this->location = location;
		this->base = base;
		this->underlying = underlying;
	}

	~Error(){};
};

#endif

