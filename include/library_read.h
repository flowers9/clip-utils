#ifndef _LIBRARY_READ_H
#define _LIBRARY_READ_H

#include "read.h"	// Read

class LibraryRead : public Read {
    public:
	int library;
	bool is_forward;
	bool is_reverse;
	const LibraryRead *pair;
	LibraryRead(void) : library(-1), is_forward(0), is_reverse(0), pair(NULL) { }
	explicit LibraryRead(const std::string &__s) : Read(__s), library(-1), is_forward(0), is_reverse(0), pair(NULL) { }
	~LibraryRead(void) { }
};

#endif // !_LIBRARY_READ_H
