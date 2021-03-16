/*
 * This takes a file with a list of readnames, and converts blat output
 * from solexa read pairs into a highly compressed format while also
 * filtering out sufficiently bad matches; filtering is done based on
 * a given minimum identity and offset sets from pairs to the insert;
 * the output is binary, read1, read2, length of longest read (negative
 * if both reads are the same length), length of shortest read minus Ns,
 * match length, and gap offset for the leading and trailing pairs.
 * The reads are 4 byte numbers which are indexes into the list of
 * readnames.  The lengths are 2 bytes, and the offsets 1.
 */

// Note: no longer takes a file of readnames - doesn't appear to be worth
// the effort compressing in that fashion.  Output is also text, not binary.

#include "breakup_line.h"	/* breakup_line() */
#include "open_compressed.h"	/* close_compressed(), open_compressed(), pfgets() */
#include "strtostr.h"	/* strtostr_exact() */
#include "version.h"	/* VERSION */
#include <algorithm>	/* max(), min() */
#include <getopt.h>	// getopt(), optarg, optind
#include <list>		/* list<> */
#include <map>		/* map<> */
#include <stdint.h>	/* uint16_t, uint32_t */
#include <stdio.h>	/* EOF, fprintf(), printf(), stderr */
#include <stdlib.h>	/* atoi(), exit(), strtod() */
#include <string>	/* string */
#include <unistd.h>	// write()
#include <utility>	/* make_pair(), pair<> */
#include <vector>	/* vector<> */

typedef uint32_t read_id_t;
typedef uint16_t read_length_t;

#define INSERT_LENGTH 48
static double opt_read_identity;
static int opt_read_offset;
//static std::map<std::string, size_t> read_name_to_int;

class MatchBin {
    public:
	read_id_t read1, read2;					// 32 bits
	read_length_t read_length, match_length, identity;	// 9 bits
	unsigned char gap, flag;				// 4+1 bits
	MatchBin(void) : read1(0), read2(0), read_length(0), match_length(0), identity(0), gap(0), flag(0) { }
	~MatchBin(void) { }
	void set(read_id_t __i, read_id_t __j, read_length_t __k, char __l, char __m) {
		read1 = __i;
		read2 = __j;
		read_length = __k;
		gap = __l;
		flag = __m;
	}
	void print_match(void) const;
};

void MatchBin::print_match() const {
	if (write(1, &read1, 4) != 4) {
		fprintf(stderr, "Error: write error\n");
		exit(1);
	}
	if (write(1, &read2, 4) != 4) {
		fprintf(stderr, "Error: write error\n");
		exit(1);
	}
	// compress remainder into 4 bytes
	const uint32_t i = (read_length << 23) + (match_length << 14) + (identity << 5) + (gap << 1) + flag;
	if (write(1, &i, 4) != 4) {
		fprintf(stderr, "Error: write error\n");
		exit(1);
	}
}

class MatchText {
    public:
	std::string read1, read2;
	size_t read_length, match_length, identity, gap;
	unsigned char flag;
	MatchText(void) : read_length(0), match_length(0), identity(0), gap(0), flag(0) { }
	~MatchText(void) { }
	void set(const std::string &__i, const std::string &__j, size_t __k, size_t __l, unsigned char __m) {
		read1 = __i;
		read2 = __j;
		read_length = __k;
		gap = __l;
		flag = __m;
	}
	void print_match(void) const;
};

void MatchText::print_match() const {
	printf("%s %s %lu %lu %lu %lu %u\n", read1.c_str(), read2.c_str(), read_length, match_length, identity, gap, flag);
}

static void print_usage() {
	//fprintf(stderr, "usage: filter_blat [opts] <-n read_name_file> <blat_file1> [blat_file2] ...\n");
	fprintf(stderr, "usage: filter_blat [opts] <blat_file1> [blat_file2] ...\n");
	fprintf(stderr, "\t-I\tmatch identity [.98]\n");
	fprintf(stderr, "\t-O\tmatch offset [2]\n");
	exit(0);
}

#if 0
static void read_name_list(const char *filename) {
	int fd = open_compressed(filename);
	size_t i = 0;
	std::string line;
	while (pfgets(fd, line) != -1) {
		read_name_to_int[line] = i++;
	}
	close_compressed(fd);
}
#endif

static void get_opts(int argc, char **argv) {
	/* set option defaults */
	opt_read_identity = .98;
	opt_read_offset = 2;
	/* read in options */
	double x;
	int c;
	//while ((c = getopt(argc, argv, "hI:n:O:V")) != EOF) {
	while ((c = getopt(argc, argv, "hI:O:V")) != EOF) {
		switch (c) {
		    case 'h':
			print_usage();
			return;
		    case 'I':
			x = strtod(optarg, NULL);
			if (x < 0 || 1 < x) {
				fprintf(stderr, "Error: match identity is out of range [0,1]: %lf\n", x);
				print_usage();
			}
			opt_read_identity = x;
			break;
		    //case 'n':
			//read_name_list(optarg);
			//break;
		    case 'O':
			c = atoi(optarg);
			if (c < 0) {
				fprintf(stderr, "Error: match offset is negative: %d\n", c);
				print_usage();
			//} else if (c > 15) {
			//	fprintf(stderr, "Warning: match offset is greater than 15, setting to 15\n");
			//	c = 15;
			}
			opt_read_offset = c;
			break;
		    case 'V':
			fprintf(stderr, "parse_output version %s\n", VERSION);
			exit(0);
		    default:
			print_usage();
		}
	}
}

// convert a comma delimited string of ints into an array
static void read_list(const std::string &list, std::vector<int> &data) {
	std::string::size_type i = 0;
	std::string::size_type j = i;
	std::string s = strtostr_exact(list, ",", &i);
	while (j != i && !s.empty()) {
		j = i;
		data.push_back(atoi(s.c_str()));
		s = strtostr_exact(list, ",", &i);
	}
}

// find overlap between given ranges and the insert
static void get_insert_range(const std::vector<int> &starts, const std::vector<int> &lengths, const int offset, std::list<std::pair<int, int> > &ranges) {
	const int x = offset + INSERT_LENGTH;
	size_t i;
	for (i = 0; i != lengths.size(); ++i) {
		if (starts[i] + lengths[i] <= offset) {
		} else if (starts[i] + lengths[i] < x) {
			ranges.push_back(std::make_pair(std::max(offset, starts[i]), starts[i] + lengths[i]));
		} else {
			if (starts[i] < x) {
				ranges.push_back(std::make_pair(std::max(offset, starts[i]), x));
			}
			return;
		}
	}
}

// subtract the insert ranges for 2 from the match ranges for 1
static void get_match_range(const std::vector<int> &starts1, const std::vector<int> &starts2, const std::vector<int> &lengths, const std::list<std::pair<int, int> > &insert_ranges, std::list<std::pair<int, int> > &ranges) {
	std::list<std::pair<int, int> >::const_iterator a = insert_ranges.begin();
	std::list<std::pair<int, int> >::const_iterator end_a = insert_ranges.end();
	size_t i;
	for (i = 0; i != starts1.size() && a != end_a; ++i) {
		if (starts2[i] + lengths[i] <= a->first) {
			ranges.push_back(std::make_pair(starts1[i], starts1[i] + lengths[i]));
		} else {
			if (starts2[i] != a->first) {
				ranges.push_back(std::make_pair(starts1[i], starts1[i] + a->first - starts2[i]));
			}
			if (starts2[i] + lengths[i] != a->second) {
				ranges.push_back(std::make_pair(starts1[i] + a->second - starts2[i], starts1[i] + lengths[i]));
			}
			++a;
		}
	}
	for (; i != starts1.size(); ++i) {
		ranges.push_back(std::make_pair(starts1[i], starts1[i] + lengths[i]));
	}
}

// find distance from nearest matched basepair to insert, and return
// longest of both sides
static int find_gap(int offset, const std::list<std::pair<int, int> > &ranges) {
	std::list<std::pair<int, int> >::const_iterator a = ranges.begin();
	std::list<std::pair<int, int> >::const_iterator end_a = ranges.end();
	int gap = -1;
	for (; a != end_a && a->first < offset; ++a) {
		if (offset > a->second) {
			gap = offset - a->second;
		} else {
			gap = 0;
			break;
		}
	}
	if (gap == -1) {
		return -1;
	}
	offset += INSERT_LENGTH;
	for (; a != end_a && a->second <= offset; ++a) { }
	if (a == end_a) {
		return -1;
	} else if (a->first <= offset) {
		return gap;
	} else {
		return std::max(gap, a->first - offset);
	}
}

// find longest gap between matched sequence and query/target insert
static int find_overall_gap(const std::string &l, const std::string &qs, const std::string &ts, const int q_offset, const int t_offset) {
	std::vector<int> lengths, q_starts, t_starts;
	read_list(l, lengths);
	read_list(qs, q_starts);
	read_list(ts, t_starts);
	// now subtract target insert from query match, and vice-versa
	std::list<std::pair<int, int> > ti_ranges, q_ranges;
	get_insert_range(t_starts, lengths, t_offset, ti_ranges);
	get_match_range(q_starts, t_starts, lengths, ti_ranges, q_ranges);
	const int i = find_gap(q_offset, q_ranges);
	if (i == -1 || i > opt_read_offset) {
		return -1;
	}
	std::list<std::pair<int, int> > qi_ranges, t_ranges;
	get_insert_range(q_starts, lengths, q_offset, qi_ranges);
	get_match_range(t_starts, q_starts, lengths, qi_ranges, t_ranges);
	const int j = find_gap(t_offset, t_ranges);
	if (j == -1 || j > opt_read_offset) {
		return -1;
	}
	return std::max(i, j);
}

static void parse_output(const char *blat_file) {
	const int fd = open_compressed(blat_file);
	std::string line;
	pfgets(fd, line);	// skip header lines
	pfgets(fd, line);
	pfgets(fd, line);
	pfgets(fd, line);
	pfgets(fd, line);
	while (pfgets(fd, line) != -1) {
		std::vector<std::string> fields;
		breakup_line(line, fields);
		if (fields.size() < 21) {
			fprintf(stderr, "Warning: short line in %s: %lu\n", blat_file, fields.size());
			continue;
		}
		// skip if query name == target name
		if (fields[9] == fields[13]) {
			continue;
		}
		// match needs to include insert
		int n_count = atoi(fields[3].c_str());
		if (n_count < INSERT_LENGTH) {
			continue;
		}
		int query_size = atoi(fields[10].c_str());
#if 0
		if (query_size > 511) {
			fprintf(stderr, "Warning: read length greater than 511, skipping match: %s\n", fields[9].c_str());
			continue;
		}
#endif
		int target_size = atoi(fields[14].c_str());
#if 0
		if (target_size > 511) {
			fprintf(stderr, "Warning: read length greater than 511, skipping match: %s\n", fields[13].c_str());
			continue;
		}
#endif
		//MatchBin m;
		MatchText m;
		m.match_length = std::min(query_size, target_size) - n_count;
		m.identity = atoi(fields[0].c_str()) + atoi(fields[2].c_str());
		// skip if not similar
		if (m.identity < opt_read_identity * m.match_length) {
			continue;
		}
		int target_offset = atoi(fields[13].substr(fields[13].rfind('-') + 1).c_str());
		int query_offset = atoi(fields[9].substr(fields[9].rfind('-') + 1).c_str());
		int gap = find_overall_gap(fields[18], fields[19], fields[20], fields[8] == "+" ? query_offset : query_size - query_offset - INSERT_LENGTH, target_offset);
		if (gap == -1) {
			continue;
		}
		if (query_size >= target_size) {
			//m.set(read_name_to_int[fields[9]], read_name_to_int[fields[13]], query_size, std::max(query_gap, target_gap), query_size == target_size ? 1 : 0);
			m.set(fields[9], fields[13], query_size, gap, query_size == target_size ? 1 : 0);
		} else {
			m.set(fields[13], fields[9], target_size, gap, 0);
		}
		m.print_match();
	}
	close_compressed(fd);
}

int main(int argc, char **argv) {
	get_opts(argc, argv);
	if (optind == argc) {
		parse_output("-");
	} else {
		for (; optind != argc; ++optind) {
			parse_output(argv[optind]);
		}
	}
	return 0;
}
