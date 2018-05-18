#ifndef _LIBRARY_READ_LIB_H
#define _LIBRARY_READ_LIB_H

#include "library_read.h"	/* LibraryRead */
#include <list>		/* list<> */
#include <string>	/* string */

extern int library_read_sequence(const char *, std::list<LibraryRead> &, bool);

#endif /* !_LIBRARY_READ_LIB_H */
