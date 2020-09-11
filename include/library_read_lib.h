#ifndef _LIBRARY_READ_LIB_H
#define _LIBRARY_READ_LIB_H

#include <list>		// list<>
#include <string>	// string

class LibraryRead;

extern int library_read_sequence(const char *, std::list<LibraryRead> &, bool);

#endif // !_LIBRARY_READ_LIB_H
