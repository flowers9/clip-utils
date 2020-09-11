#ifndef _PARSE_READ_H
#define _PARSE_READ_H

#include <string>	// string
#include <sys/types.h>	// size_t

class LibraryRead;
struct ProtoReadPattern;

extern void init_read_patterns(const ProtoReadPattern * = NULL, size_t = 0);
extern void parse_read_name(LibraryRead &);
extern std::string make_index_name(const LibraryRead &);
extern std::string make_index_pair_name(const LibraryRead &);

#endif // !_PARSE_READ_H
