#ifndef _PARSE_READ_H
#define _PARSE_READ_H

#include "library_read.h"	/* LibraryRead */
#include "read_match.h"	/* ProtoReadPattern */
#include <string>	/* string */
#include <sys/types.h>	/* size_t */

extern void init_read_patterns(const ProtoReadPattern * = NULL, size_t = 0);
extern void parse_read_name(LibraryRead &);
extern std::string make_index_name(const LibraryRead &);
extern std::string make_index_pair_name(const LibraryRead &);

#endif /* !_PARSE_READ_H */
