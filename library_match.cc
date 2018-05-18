#include "library_match.h"
#include "library_read.h"	/* LibraryRead */
#include <regex.h>	/* REG_EXTENDED, REG_NOSUB, regcomp(), regerror(), regmatch_t */
#include <stdio.h>	/* fprintf(), stderr */
#include <string>	/* string */
#include <sys/types.h>	/* size_t */

LibraryMatch::LibraryMatch(const ProtoLibraryPattern &a) {
	bool has_subexpressions = a.name == NULL;
	int i;
	if (has_subexpressions) {
		i = regcomp(&pattern, a.regexp, REG_EXTENDED);
	} else {
		i = regcomp(&pattern, a.regexp, REG_EXTENDED | REG_NOSUB);
	}
	if (i != 0) {
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
		library_hint = a.library_hint;
		if (has_subexpressions) {
			nmatch = 2;
			pmatch = new regmatch_t[2];
		} else {
			nmatch = 0;
			pmatch = NULL;
			name = a.name;
		}
	}
}

/*
 * see if read name is matched by library pattern; if so, returns library name
 */

bool LibraryMatch::is_match(const LibraryRead &a, std::string &library) const {
	if (library_hint != a.library) {
		return 0;
	}
	if (regexec(&pattern, a.name().c_str(), nmatch, pmatch, 0) != 0) {
		return 0;
	}
	if (name.empty()) {
		library = a.name().substr(pmatch[1].rm_so, pmatch[1].rm_eo - pmatch[1].rm_so);
	} else {
		library = name;
	}
	return 1;
}
