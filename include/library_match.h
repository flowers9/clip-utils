#ifndef _LIBRARY_MATCH_H
#define _LIBRARY_MATCH_H

#include "library_read.h"	/* LibraryRead */
#include <regex.h>	/* regex_t, regmatch_t */
#include <string>	/* string */
#include <sys/types.h>	/* size_t */

typedef struct _ProtoLibraryPattern {
	int library_hint;
	const char *regexp;
	const char *name;		/*
					 * name of library (empty if
					 * subexpression should be used
					 */
} ProtoLibraryPattern;

class LibraryMatch {
    private:
	regex_t pattern;
	size_t nmatch;
	regmatch_t *pmatch;
	int library_hint;
	std::string name;	/* get name from subexpression if empty */
    public:
	std::string regexp;
	explicit LibraryMatch(const ProtoLibraryPattern &);
	~LibraryMatch() { }
	bool is_match(const LibraryRead &, std::string &) const;
};

#endif /* !_LIBRARY_MATCH_H */
