#ifndef _FIND_LIBRARY_H
#define _FIND_LIBRARY_H

#include "library_read.h"	/* LibraryRead */
#include "library_match.h"	/* ProtoLibraryPattern */
#include <string>	/* string */
#include <sys/types.h>	/* size_t */

extern void init_library_patterns(const ProtoLibraryPattern * = NULL, size_t = 0);
extern std::string find_library(const LibraryRead &);

#endif /* !_FIND_LIBRARY_H */
