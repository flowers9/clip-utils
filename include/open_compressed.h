#ifndef _OPEN_COMPRESSED_H
#define _OPEN_COMPRESSED_H

#include <string>	// string
#include <sys/types.h>	// size_t, ssize_t

extern void get_suffix(const std::string &, std::string &);
extern int find_suffix(std::string &, std::string &);
extern int open_compressed(const std::string &);
extern void close_compressed(int);
extern ssize_t pfgets(int, std::string &, char delim = '\n');
extern ssize_t skip_next_line(int, char delim = '\n');
extern ssize_t skip_next_chars(int, size_t);
extern ssize_t pfread(int, void *, size_t);
extern ssize_t pfpeek(int, void *, size_t);

#endif // !_OPEN_COMPRESSED_H
