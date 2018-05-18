#include "find_library.h"
#include "library_match.h"	/* ProtoLibraryPattern, LibraryMatch */
#include "library_read.h"	/* LibraryRead */
#include <list>		/* list<> */
#include <stdio.h>	/* fprintf(), stderr */
#include <string>	/* string */
#include <sys/types.h>	/* size_t */

static std::list<LibraryMatch> library_patterns;

static const ProtoLibraryPattern standard_list[] = {
	{ 0, "^([[:alpha:]]{3,4})", NULL },
	{ 0, "^([[:alnum:]]{2,3}[[:alpha:]])", NULL },
	{ 0, "", "Other" },
	{ 1, "^([[:alpha:]]{3,4})", NULL },
	{ 1, "", "JGI" },
	{ 2, "", "Los Alamos" },
	{ 3, "", "stg.pld" },
	{ 4, "^.[[:digit:]]+[^[:digit:]][[:digit:]]{1,2}[^[:digit:]]", "WIBR m13" },
	{ 4, "^.[[:digit:]]+[^[:digit:]]6[[:digit:]]{2,4}[^[:digit:]]", "WIBR 4k" },
	{ 4, "^.[[:digit:]]+[^[:digit:]]9[[:digit:]]{2}[^[:digit:]]", "WIBR shatter" },
	{ 4, "^.[[:digit:]]+[^[:digit:]]5[[:digit:]]{3}[^[:digit:]]", "WIBR 10k" },
	{ 4, "^.[[:digit:]]+[^[:digit:]]8[[:digit:]]{3}[^[:digit:]]", "WIBR fosmid" },
	{ 4, "", "WIBR" },
};

void init_library_patterns(const ProtoLibraryPattern *list, size_t n) {
	if (list == NULL) {
		list = standard_list;
		n = sizeof(standard_list) / sizeof(standard_list[0]);
	}
	library_patterns.clear();
	while (n > 0) {
		--n;
		LibraryMatch a(list[n]);
		library_patterns.push_front(a);
	}
}

/* find what library a read belongs to */

std::string find_library(const LibraryRead &a) {
	std::string library;
	std::list<LibraryMatch>::iterator b = library_patterns.begin();
	std::list<LibraryMatch>::iterator end = library_patterns.end();
	/* find matching library */
	for (; b != end && !b->is_match(a, library); ++b) { }
	if (b == end) {
		fprintf(stderr, "Could not find library for read: %s\n", a.name().c_str());
	}
	return library;
}
