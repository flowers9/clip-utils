#include "library_read.h"	/* LibraryRead */
#include "read_match.h"
#include <regex.h>	/* REG_EXTENDED, regcomp(), regerror(), regmatch_t */
#include <stdio.h>	/* fprintf(), stderr */
#include <string>	/* string */
#include <strings.h>	/* strcasecmp() */
#include <sys/types.h>	/* size_t */

ReadMatch::ReadMatch(const ProtoReadPattern &a) {
	int i;
	if ((i = regcomp(&pattern, a.regexp, REG_EXTENDED)) != 0) {
		size_t j = regerror(i, &pattern, NULL, 0);
		if (j == 0) {
			fprintf(stderr, "regcomp: regexp error %d: %s\n", i, a.regexp);
		} else {
			char *err = new char[j];
			regerror(i, &pattern, err, j);
			fprintf(stderr, "regcomp: regexp error %s: %s\n", err, a.regexp);
			delete[] err;
		}
	} else {
		regexp = a.regexp;
		nmatch = a.subexpressions + 1;
		pmatch = new regmatch_t[nmatch];
		library_hint = a.library_hint;
		direction = a.direction;
		forward = a.forward;
		reverse = a.reverse;
	}
}

/*
 * checks to see if the read name matches, and, if so, changes it's library and
 * forward/reverse flags to match the library; returns true if the match was
 * sucessful.
 */

bool ReadMatch::parse_name(LibraryRead &a) {
	if (regexec(&pattern, a.name().c_str(), direction + 1, pmatch, 0) != 0) {
		return 0;
	}
	a.library = library_hint;
	regmatch_t *b = &pmatch[direction];
	if (forward.length() == (size_t)(b->rm_eo - b->rm_so) && strcasecmp(forward.c_str(), a.name().substr(b->rm_so, b->rm_eo - b->rm_so).c_str()) == 0) {
		a.is_forward = 1;
	} else if (reverse.length() == (size_t)(b->rm_eo - b->rm_so) && strcasecmp(reverse.c_str(), a.name().substr(b->rm_so, b->rm_eo - b->rm_so).c_str()) == 0) {
		a.is_reverse = 1;
	}
	return 1;
}

/*
 * checks to see if the given read name matches, and, if so, constructs the
 * index name for it; returns true if the match was sucessful
 */

bool ReadMatch::index_name(const LibraryRead &a, std::string &index) {
	if (regexec(&pattern, a.name().c_str(), nmatch, pmatch, 0) != 0) {
		return 0;
	}
	index.clear();
	size_t i;
	for (i = 1; i < nmatch; ++i) {
		if (i != direction) {
			index += a.name().substr(pmatch[i].rm_so, pmatch[i].rm_eo - pmatch[i].rm_so);
		} else if (a.is_forward) {
			index += 'f';
		} else {
			index += 'r';
		}
	}
	return 1;
}

/*
 * checks to see if the given read name matches, and, if so, constructs the
 * index name for its pair; returns true if the match was sucessful
 */

bool ReadMatch::index_pair_name(const LibraryRead &a, std::string &index) {
	if (regexec(&pattern, a.name().c_str(), nmatch, pmatch, 0) != 0) {
		return 0;
	}
	index.clear();
	size_t i;
	for (i = 1; i < nmatch; ++i) {
		if (i != direction) {
			index += a.name().substr(pmatch[i].rm_so, pmatch[i].rm_eo - pmatch[i].rm_so);
		} else if (a.is_forward) {
			index += 'r';
		} else {
			index += 'f';
		}
	}
	return 1;
}
