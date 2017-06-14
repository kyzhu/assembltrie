#include <string>
#include <cstring>
#include <experimental/string_view>
#include <algorithm>

class Util {

public:
	static int baseToInt(char c) {
		switch (c) {
			case 'A': return 0;
			case 'T': return 1;
			case 'G': return 2;
			case 'C': return 3;
			default: return 0;
		}
	}

	static int baseToWt(char c) {
		switch (c) {
			case 'A': return -1;
			case 'T': return 1;
			case 'G': return 1;
			case 'C': return -1;
			default: return 0;
		}
	}

	static int baseToInt(char c, char u) {
		int _u = baseToInt(u);
		int _c = baseToInt(c);
		return (_c > _u) ? _c : _c - 1;
	}

	/* 4-letter model. */
	/* Underlying symbol is not allowed to be 'N'. */
	static int baseToIntHuffman(char c, char u) {
		int _u = baseToInt(u);
		int _c = baseToInt(c);
		if (_c > _u)
			_c--;
		return ((_c > 0) ? _c + 1 : _c);
	}

	static char intToBase(int val) {
		switch (val) {
			case 0: return 'A';
			case 1: return 'T';
			case 2: return 'G';
			case 3: return 'C';
			default: return 'N';
		}
	}

	static char intToBase(int val, char u) {
		if (val < baseToInt(u))
			return intToBase(val);
		else
			return intToBase(val + 1);
	}

	static char complement(char c) {   
    		switch(c) {   
    			case 'A': return 'T';
    			case 'T': return 'A';
    			case 'G': return 'C';
    			case 'C': return 'G';
			default: return 'N';
    		}
	}

	static void arraycpy(char *src, int srcpos, char *dest, int destpos, int l) {
		strncpy(dest + destpos, src + srcpos, l);
	}

	static void arraycpy(const char *src, int srcpos, char *dest, int destpos, int l) {
		strncpy(dest + destpos, src + srcpos, l);
	}
	
	/* Get the reverse complement string of S. */
	static std::string getRCstring(std::string &s) {
		std::string rc = s;
		std::transform(rc.begin(), rc.end(), rc.begin(), complement);
		std::reverse(rc.begin(), rc.end());
		return rc;
	}

	/* Get the minimum 5-mer of S. */
	static int getMin(std::string &s) {
		int min = 1023;
		int temp = min;
		for (size_t i = 0; i < s.length(); i++) {
			if (s[i] == 'N') {
				temp = 1023;
				continue;
			}
			temp = ((temp & 255) << 2) + baseToInt(s[i]);
			if (temp < min)
				min = temp;
		}
		return min;
	}

	static int getACweight(std::string &s) {
		int sum = 0;
		for (size_t i = 0; i < s.length(); i++) {
			sum += baseToWt(s[i]);
		}
		return sum;
	}
};

