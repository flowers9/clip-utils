/*
 * this takes a blast output file (well, blat in blast output mode), and
 * compresses it into many fewer lines; =readname starts a section of
 * matches against that readname, >readname are matches from that read,
 * and matches themselves are simply identity, score, start, stop, s_start,
 * s_stop, matched sequence, and (if different) matching sequence
 */

#include "open_compressed.h"	/* close_compressed(), open_compressed(), pfgets() */
#include "strtostr.h"	/* strtostr() */
#include "version.h"	/* VERSION */
#include <algorithm>	/* swap() */
#include <errno.h>	/* errno */
#include <getopt.h>	// getopt(), optarg, optind
#include <stdio.h>	/* EOF, fclose(), fopen(), fprintf(), stderr, stdout */
#include <stdlib.h>	// atoi(), exit()
#include <string.h>	/* strerror() */
#include <string>	/* string */

static std::string opt_output_file;

class Match {
    private:
	FILE *_fp;
    public:
	long length, identity, score, start, stop, s_start, s_stop;
	std::string query_id, id, qs, ss;
	explicit Match(const std::string &__s) : length(0), identity(0), score(0), start(0), stop(0), s_start(0), s_stop(0) {
		if (__s.empty() || __s == "-") {
			_fp = stdout;
		} else {
			_fp = fopen(__s.c_str(), "w");
			if (_fp == NULL) {
				fprintf(stderr, "Error: fopen: %s: %s\n", __s.c_str(), strerror(errno));
				exit(1);
			}
		}
	}
	~Match(void) {
		if (_fp != NULL) {
			fclose(_fp);
		}
	}
	void normalize(void);
	inline void print_header1(void) {
		fprintf(_fp, "=%s\n", query_id.c_str());
	}
	inline void print_header2(void) {
		fprintf(_fp, ">%s %ld\n", id.c_str(), length);
	}
	inline void print_match(void) {
		fprintf(_fp, "%ld %ld %ld %ld %ld %ld\n%s", identity, score, start, stop, s_start, s_stop, qs.c_str());
		if (qs != ss) {
			fprintf(_fp, " %s", ss.c_str());
		}
		fprintf(_fp, "\n");
	}
	
};

class BlatFile {
    private:
	int _fd;
    public:
	explicit BlatFile(const std::string &__s) {
		if ((_fd = open_compressed(__s)) == -1) {
			fprintf(stderr, "Error: open_compressed: %s: %s\n", __s.c_str(), strerror(errno));
			exit(1);
		}
	}
	~BlatFile(void) {
		if (_fd != -1) {
			close_compressed(_fd);
		}
	}
	bool read_next(Match &);
};

// reverse a string in place
static void reverse(std::string &s) {
	if (s.size() < 2) {
		return;
	}
	std::string::size_type i = 0;
	std::string::size_type j = s.size() - 1;
	for (; i < j; ++i, --j) {
		std::swap(s[i], s[j]);
	}
}

// complement the basepairs of a string - acgt -> tgca
static void complement(std::string &s) {
	std::string::iterator a = s.begin();
	std::string::iterator end_a = s.end();
	for (; a != end_a; ++a) {
		switch (*a) {
		    case 'a':
			*a = 't';
			break;
		    case 'c':
			*a = 'g';
			break;
		    case 'g':
			*a = 'c';
			break;
		    case 't':
			*a = 'a';
			break;
		}
	}
}

// complement if needed
void Match::normalize() {
	if (start > stop) {
		std::swap(start, stop);
		std::swap(s_start, s_stop);
		reverse(qs);
		reverse(ss);
		complement(qs);
		complement(ss);
	}
}

// return if a match was found;
// XXX - this could have a bit more input format verification
bool BlatFile::read_next(Match &b) {
	std::string line;
	std::string::size_type i, j;
	while (pfgets(_fd, line) != -1) {
		if (line.empty()) {
		} else if (line[0] == '>') {
			b.id = strtostr(line, &(i = 1));
			if (pfgets(_fd, line) == -1) {	// length
				fprintf(stderr, "Error: unexpected end of file 1: %s\n", b.id.c_str());
				return 0;
			}
			if ((j = line.find("Length = ")) == std::string::npos) {
				fprintf(stderr, "Error: missing length: %s\n", b.id.c_str());
				return 0;
			}
			b.length = atoi(line.substr(j + 9).c_str());
			b.print_header2();
			if (pfgets(_fd, line) == -1) {	// blank
				fprintf(stderr, "Error: unexpected end of file 2: %s\n", b.id.c_str());
				return 0;
			}
			if (pfgets(_fd, line) == -1) {	// Score
				fprintf(stderr, "Error: unexpected end of file 3: %s\n", b.id.c_str());
				return 0;
			}
			if (line.compare(0, 9, " Score = ") != 0) {
				fprintf(stderr, "Error: missing score: %s\n", b.id.c_str());
				return 0;
			}
			break;
		} else if (line.compare(0, 9, " Score = ") == 0) {
			break;
		} else if (line.compare(0, 6, "Query=") == 0) {
			b.query_id = strtostr(line, &(i = 6));
			b.print_header1();
		}
	}
	if (line.empty()) {
		return 0;
	}
	b.score = atoi(line.substr(9).c_str());
	if (pfgets(_fd, line) == -1) {	// identities
		fprintf(stderr, "Error: unexpected end of file 4: %s\n", b.id.c_str());
		return 0;
	}
	if ((j = line.find("Identities = ")) == std::string::npos) {
		fprintf(stderr, "Error: missing identities: %s\n", b.id.c_str());
		return 0;
	}
	b.identity = atoi(line.substr(j + 13).c_str());
	if (pfgets(_fd, line) == -1) {	// strand
		fprintf(stderr, "Error: unexpected end of file 5: %s\n", b.id.c_str());
		return 0;
	}
	if (pfgets(_fd, line) == -1) {	// blank
		fprintf(stderr, "Error: unexpected end of file 6: %s\n", b.id.c_str());
		return 0;
	}
	if (pfgets(_fd, line) == -1) {	// ^(Query:\s+(\d+)\s+)(\S+)\s+(\d+)$
		fprintf(stderr, "Error: unexpected end of file 7: %s\n", b.id.c_str());
		return 0;
	}
	b.start = atoi(strtostr(line, &(i = 6)).c_str());
	std::string::size_type header_length = i = line.find_first_not_of(' ', i + 1);
	if (i == std::string::npos) {
		fprintf(stderr, "Error: malformed query line: %s\n", b.id.c_str());
		return 0;
	}
	b.qs = strtostr(line, &i);
	b.stop = atoi(strtostr(line, &i).c_str());
	if (pfgets(_fd, line) == -1) {	// tag
		fprintf(stderr, "Error: unexpected end of file 8: %s\n", b.id.c_str());
		return 0;
	}
	if (pfgets(_fd, line) == -1) {	// ^Sbjct:\s+(\d+)\s+(\S+)\s+(\d+)$
		fprintf(stderr, "Error: unexpected end of file 9: %s\n", b.id.c_str());
		return 0;
	}
	b.s_start = atoi(strtostr(line, &(i = 6)).c_str());
	b.ss = strtostr(line, &i);
	b.s_stop = atoi(strtostr(line, &i).c_str());
	if (pfgets(_fd, line) == -1) {	// blank
		fprintf(stderr, "Error: unexpected end of file 10: %s\n", b.id.c_str());
		return 0;
	}
	// next line will be blank unless it's a continuation
	while (pfgets(_fd, line) != -1 && line.compare(0, 6, "Query:") == 0) {
		b.qs += strtostr(line, &(i = header_length));
		b.stop = atoi(strtostr(line, &i).c_str());
		if (pfgets(_fd, line) == -1) {	// skip tag line
			fprintf(stderr, "Error: unexpected end of file 11: %s\n", b.id.c_str());
			return 0;
		}
		if (pfgets(_fd, line) == -1) {
			fprintf(stderr, "Error: unexpected end of file 12: %s\n", b.id.c_str());
			return 0;
		}
		b.ss += strtostr(line, &(i = header_length));
		b.s_stop = atoi(strtostr(line, &i).c_str());
		if (pfgets(_fd, line) == -1) {	// blank line
			fprintf(stderr, "Error: unexpected end of file 13: %s\n", b.id.c_str());
			return 0;
		}
	}
	b.normalize();			// complement if needed
	b.print_match();
	return 1;
}

static void print_usage() {
	fprintf(stderr, "usage: compress_blat [-o <output file>] [file]\n");
	exit(1);
}

static void get_opts(const int argc, char * const * argv) {
	int c;
	while ((c = getopt(argc, argv, "o:V")) != EOF) {
		switch (c) {
		    case 'o':
			opt_output_file = optarg;
			break;
		    case 'V':
			fprintf(stderr, "compress_blat version %s\n", VERSION);
			exit(0);
		    default:
			print_usage();
		}
	}
}

int main(const int argc, char * const * argv) {
	get_opts(argc, argv);
	BlatFile blat(optind == argc ? "-" : argv[optind]);
	Match b(opt_output_file);
	while (blat.read_next(b)) { }
	return 0;
}
