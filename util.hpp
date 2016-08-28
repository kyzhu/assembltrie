#include <cstring>

class Util {
public:
	static int baseToInt(char c) {
		if (c == 'A') return 0;
		if (c == 'T') return 1;
		if (c == 'G') return 2;
		if (c == 'C') return 3;
		return 4;
	}

	static char intToBase(int val) {
		if (val == 0) return 'A';
		if (val == 1) return 'T';
		if (val == 2) return 'G';
		if (val == 3) return 'C';
		return 'N';
	}

	static void arraycpy(char *src, int srcpos, char *dest, int destpos, int l) {
		if (srcpos == 0 && destpos == 0)
			strncpy(dest, src, l);
		else {
			for (int i = 0; i < l; i++)
				dest[destpos + i] = src[srcpos + i];
		}
	}

	static void arraycpy(const char *src, int srcpos, char *dest, int destpos, int l) {
		if (srcpos == 0 && destpos == 0)
			strncpy(dest, src, l);
		else {
			for (int i = 0; i < l; i++)
				dest[destpos + i] = src[srcpos + i];
		}
	}

};