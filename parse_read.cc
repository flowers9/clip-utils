#include "library_read.h"	/* LibraryRead */
#include "parse_read.h"
#include "read_match.h"	/* ProtoReadPattern, ReadMatch */
#include <list>		/* list<> */
#include <string>	/* string */
#include <sys/types.h>	/* size_t */

static std::list<ReadMatch> read_patterns;

static ProtoReadPattern standard_list[] = {
	{ 0, "^(.+\\.TR\\.)([FR])$", 2, 2, "f", "r" },
	{ 0, "^(.+\\.)[xy].(.*-)([FfRr])(.*)$", 4, 3, "f", "r" },
	{ 0, "^(.+\\.)[xy].(.*-)([SsTt])(.*)$", 4, 3, "s", "t" },
	{ 1, "^([[:alnum:]]{2,3}[[:alpha:]]+[[:digit:]]+\\.)([xy])", 2, 2, "x", "y" },
	{ 1, "^([[:alpha:]]{3,}[[:digit:]]+\\.)([xy])", 2, 2, "x", "y" },
	{ 1, "^([[:alnum:]]{2,3}[[:alpha:]]+[[:digit:]]+\\.)([bg])", 2, 2, "b", "g" },
	{ 1, "^([[:alpha:]]{3,}[[:digit:]]+\\.)([bg])", 2, 2, "b", "g" },
	{ 2, "^(.+_.+\\.)([xy])([[:digit:]]{1,2})$", 3, 2, "x", "y" },
	{ 3, "^(.+\\.)([xy])\\.stg\\.pld", 2, 2, "x", "y" },
	{ 3, "^(.+\\.)([pq])", 2, 2, "p", "q" },
	{ 4, "^([LG][[:digit:]]+P[[:digit:]]+)(.*g)([FR])(\\.T[[:digit:]]+\\.scf)$", 4, 3, "f", "r" },
	{ 4, "^([LG][[:digit:]]+P[[:digit:]]+)([FR])(.*\\.T[[:digit:]]+\\.scf)$", 3, 2, "f", "r" },
	{ 0, "^(.+\\.)([fr])$", 2, 2, "f", "r" },
	{ 0, "^(.+\\.)(s)$", 2, 2, "s", "t" },
};

void init_read_patterns(const ProtoReadPattern *list, size_t n) {
	if (list == NULL) {
		list = standard_list;
		n = sizeof(standard_list) / sizeof(standard_list[0]);
	}
	read_patterns.clear();
	while (n > 0) {
		--n;
		ReadMatch a(list[n]);
		read_patterns.push_front(a);
	}
}

/*
 * returns a read's library (although it's really just a hint at this point),
 * and marks it as forward or reverse
 */

void parse_read_name(LibraryRead &a) {
	std::list<ReadMatch>::iterator b = read_patterns.begin();
	std::list<ReadMatch>::iterator end = read_patterns.end();
	for (; b != end && !b->parse_name(a); ++b) { }
}

/* returns the index string for the read */

std::string make_index_name(const LibraryRead &a) {
	std::string c;
	std::list<ReadMatch>::iterator b = read_patterns.begin();
	std::list<ReadMatch>::iterator end = read_patterns.end();
	for (; b != end && !b->index_name(a, c); ++b) { }
	return c;
}

/* returns the index string for a read's theoretical pair */

std::string make_index_pair_name(const LibraryRead &a) {
	std::string c;
	std::list<ReadMatch>::iterator b = read_patterns.begin();
	std::list<ReadMatch>::iterator end = read_patterns.end();
	for (; b != end && !b->index_pair_name(a, c); ++b) { }
	return c;
}
