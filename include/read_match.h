#ifndef _READ_MATCH_H
#define _READ_MATCH_H

#include <regex.h>	// regex_t, regmatch_t
#include <string>	// string
#include <sys/types.h>	// size_t

class LibraryRead;

struct ProtoReadPattern {
	int library_hint;
	const char *regexp;
	int subexpressions;	// total number of subexpressions
	int direction;		// subexpression to use for direction
	const char *forward;	// match for forward direction
	const char *reverse;	// match for reverse direction
};

class ReadMatch {
    private:
	regex_t pattern;
	size_t nmatch;
	regmatch_t *pmatch;
	std::string regexp;
	int library_hint;
	size_t direction;
	std::string forward;
	std::string reverse;
    public:
	ReadMatch(const ProtoReadPattern &);
	~ReadMatch() { }
	bool parse_name(LibraryRead &);
	bool index_name(const LibraryRead &, std::string &);
	bool index_pair_name(const LibraryRead &, std::string &);
};

#endif // !_READ_MATCH_H
