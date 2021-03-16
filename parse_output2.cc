/*
 * this reads a filtered blat output file (from filter_blat) and groups
 * up reads that match each other, in an attempt to eliminate reads that
 * are sufficiently covered by other reads; it includes a recursive option
 * to redo the set with eliminated reads removed from the set to try and
 * improve accuracy
 */

#include "open_compressed.h"	/* close_compressed(), open_compressed(), pfgets(), pfread() */
#include "strtostr.h"	/* strtostr() */
#include "version.h"	/* VERSION */
#include <algorithm>	/* sort() */
#include <errno.h>	/* errno */
#include <getopt.h>	// getopt(), optarg, optind
#include <list>		/* list<> */
#include <map>		/* map<> */
#include <stdint.h>	/* uint16_t, uint32_t */
#include <stdio.h>	/* EOF, fclose(), fopen(), fprintf(), printf(), stderr */
#include <stdlib.h>	/* atoi(), exit(), strtod() */
#include <string.h>	/* strerror() */
#include <string>	/* string */
#include <vector>	/* vector<> */

typedef uint32_t read_id_t;
typedef uint16_t read_length_t;

#define INSERT_LENGTH 48
static double opt_read_identity;
static unsigned int opt_read_offset;
static std::map<read_id_t, read_id_t> read_to_int;	// convert # to interal
static std::vector<read_id_t> read_to_ext;		// convert # to external
static std::map<read_id_t, read_id_t> duplicates;	// external read #s
static std::map<std::string, read_id_t> read_names;

class Match {
    public:
	read_id_t read1, read2;
	read_length_t read_length, match_length, identity;	// 9 bits
	unsigned char gap, flag;				// 4+1 bits
	Match(void) : read1(0), read2(0), read_length(0), match_length(0), identity(0), gap(0), flag(0) { }
	~Match(void) { }
	bool acceptable(void) const {
		return identity >= opt_read_identity * match_length && gap <= opt_read_offset && duplicates.find(read1) == duplicates.end() && duplicates.find(read2) == duplicates.end();
	}
	int read_match_bin(const int);
	int read_match_text(const int);
	void convert_reads(void);
};

int Match::read_match_bin(const int fd) {
	if (pfread(fd, &read1, 4) != 4) {
		return 0;
	}
	if (pfread(fd, &read2, 4) != 4) {
		return -1;
	}
	// uncompress remainder from 4 bytes
	uint32_t i;
	if (pfread(fd, &i, 4) != 4) {
		return -1;
	}
	flag = i & 1;
	i >>= 1;
	gap = i & 15;
	i >>= 4;
	identity = i & 511;
	i >>= 9;
	match_length = i & 511;
	read_length = i >> 9;
	return 1;
}

int Match::read_match_text(const int fd) {
	std::string line;
	if (pfgets(fd, line) == -1) {
		return 0;
	}
	std::string::size_type i = 0;
	read1 = read_names[strtostr(line, &i)];
	read2 = read_names[strtostr(line, &i)];
	read_length = atoi(strtostr(line, &i).c_str());
	match_length = atoi(strtostr(line, &i).c_str());
	identity = atoi(strtostr(line, &i).c_str());
	gap = atoi(strtostr(line, &i).c_str());
	flag = atoi(strtostr(line, &i).c_str());
	return 1;
}

// convert read numbers from external numbering (which can have gaps in
// the list) to internal (which are contiguous)
void Match::convert_reads() {
	std::map<read_id_t, read_id_t>::const_iterator a = read_to_int.find(read1);
	if (a == read_to_int.end()) {
		read_id_t i = read_to_int.size();
		read_to_ext.push_back(read1);
		read1 = read_to_int[read1] = i;
	} else {
		read1 = a->second;
	}
	a = read_to_int.find(read2);
	if (a == read_to_int.end()) {
		read_id_t i = read_to_int.size();
		read_to_ext.push_back(read2);
		read2 = read_to_int[read2] = i;
	} else {
		read2 = a->second;
	}
}

class ReadScore {
    public:
	read_length_t read_length, match_length, identity;
	// would this be faster as a list<>?
	std::vector<read_id_t> read_list;	// sorted list
	ReadScore(void) : read_length(0), match_length(0), identity(0) { }
	explicit ReadScore(const Match &__m) : read_length(__m.read_length), match_length(__m.match_length), identity(__m.identity) {
		if (!__m.flag) {
			read_list.push_back(__m.read1);
		} else if (__m.read1 < __m.read2) {
			read_list.push_back(__m.read1);
			read_list.push_back(__m.read2);
		} else {
			read_list.push_back(__m.read2);
			read_list.push_back(__m.read1);
		}
	}
	~ReadScore(void) { }
	int cmp(const ReadScore &__a) const {
		if (match_length != __a.match_length) {
			return match_length < __a.match_length ? -1 : 1;
		} else if (identity != __a.identity) {
			return identity < __a.identity ? -1 : 1;
		} else if (read_length != __a.read_length) {
			return read_length < __a.read_length ? -1 : 1;
		} else {
			return 0;
		}
	}
	void add(read_id_t);
};

void ReadScore::add(const read_id_t x) {
	// find position to insert in list
	read_id_t i = 0;
	read_id_t n = read_list.size();
	while (i != n) {
		read_id_t j = (i + n) / 2;
		if (x < read_list[j]) {
			n = j;
		} else if (x > read_list[j]) {
			i = j + 1;
		} else {		// already in list
			return;
		}
	}
	i = read_list.size();
	read_list.push_back(0);		// placeholder
	// shuffle array to open the correct position
	for (; i != n; --i) {
		read_list[i] = read_list[i - 1];
	}
	// insert value
	read_list[i] = x;
}

static std::vector<ReadScore> best_reads;

static void read_blat_files(std::list<std::string> &blat_files, const char *blat_file_list) {
	const int fd = open_compressed(blat_file_list);
	if (fd == -1) {
		fprintf(stderr, "Error: open: %s\n", blat_file_list);
		exit(1);
	}
	std::string line;
	while (pfgets(fd, line) != -1) {
		blat_files.push_back(line);
	}
	close_compressed(fd);
}

static void read_read_names(const std::string &read_name_file) {
	const int fd = open_compressed(read_name_file);
	if (fd == -1) {
		fprintf(stderr, "Error: open: %s\n", read_name_file.c_str());
		exit(1);
	}
	std::string line;
	while (pfgets(fd, line) != -1) {
		read_names[line] = 0;
	}
	close_compressed(fd);
	// number reads lexically (input list is not assumed to be sorted)
	read_id_t i = 0;
	std::map<std::string, read_id_t>::iterator a = read_names.begin();
	std::map<std::string, read_id_t>::iterator end_a = read_names.end();
	for (; a != end_a; ++a, ++i) {
		a->second = i;
	}
}

static void print_usage() {
	fprintf(stderr, "usage: parse_output <-n read_list> <-o output_file> [opts] <blat1> [blat2] ...\n");
	fprintf(stderr, "\t-d\tonly print raw duplicate names\n");
	fprintf(stderr, "\t-I\tmatch identity [.98]\n");
	fprintf(stderr, "\t-l ##\tfile with list of extra blat files\n");
	fprintf(stderr, "\t-m\tprint reads matched against\n");
	fprintf(stderr, "\t-O\tmatch offset [2]\n");
	fprintf(stderr, "\t-r\trecurse matching to weed duplicates\n");
	exit(0);
}

static void get_opts(int argc, char **argv, std::string &output_file, bool &opt_print_matched_read, bool &opt_recurse, bool &opt_raw, std::list<std::string> &blat_files) {
	/* set option defaults */
	opt_print_matched_read = 0;
	opt_raw = 0;
	opt_read_identity = .98;
	opt_read_offset = 2;
	opt_recurse = 0;
	/* read in options */
	double x;
	int c;
	while ((c = getopt(argc, argv, "dhI:l:mn:O:o:rV")) != EOF) {
		switch (c) {
		    case 'd':
			opt_raw = 1;
			break;
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
		    case 'l':
			read_blat_files(blat_files, optarg);
			break;
		    case 'm':
			opt_print_matched_read = 1;
			break;
		    case 'n':
			read_read_names(optarg);
			break;
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
		    case 'o':
			output_file = optarg;
			break;
		    case 'r':
			opt_recurse = 1;
			break;
		    case 'V':
			fprintf(stderr, "parse_output version %s\n", VERSION);
			exit(0);
		    default:
			print_usage();
		}
	}
	if (read_names.empty()) {
		fprintf(stderr, "Error: no read name file given\n");
		print_usage();
	}
	if (output_file.empty()) {
		fprintf(stderr, "Error: no output file given\n");
		print_usage();
	}
	if (optind == argc && blat_files.empty()) {
		fprintf(stderr, "Error: no blat files given\n");
		print_usage();
	}
}

static void update_score(const read_id_t read, const ReadScore &score) {
	if (read >= best_reads.size()) {
		best_reads.push_back(score);
	} else {
		switch (best_reads[read].cmp(score)) {
		    case -1:	// replace score
			best_reads[read] = score;
			break;
		    case 0:	// append reads to list
			if (score.read_list.size() == 1 || score.read_list[0] != read) {
				best_reads[read].add(score.read_list[0]);
			} else {
				best_reads[read].add(score.read_list[1]);
			}
			break;
		}
	}
}

static void parse_output(const std::string &blat_file) {
	const int fd = open_compressed(blat_file);
	if (fd == -1) {
		fprintf(stderr, "Error: open: %s\n", blat_file.c_str());
		exit(1);
	}
	Match m;
	while (m.read_match_text(fd) == 1) {
		if (m.acceptable()) {
			m.convert_reads();
			const ReadScore score(m);
			update_score(m.read1, score);
			update_score(m.read2, score);
		}
	}
	close_compressed(fd);
}

class ScoreHolder {
    public:
	// need to use the external read number for sorting purposes
	read_id_t score, read;
	ScoreHolder(read_id_t __i, read_id_t __j) : score(__i), read(__j) { }
	~ScoreHolder(void) { }
	bool operator<(const ScoreHolder &__a) const {
		if (score != __a.score) {
			return score < __a.score;
		} else {
			return read < __a.read;
		}
	}
};

// reduce list of best reads to single best read for each read
static void reduce_sets() {
	std::vector<read_id_t> score(best_reads.size(), 0);
	std::vector<ReadScore>::const_iterator a = best_reads.begin();
	std::vector<ReadScore>::const_iterator end_a = best_reads.end();
	for (; a != end_a; ++a) {
		std::vector<read_id_t>::const_iterator b = a->read_list.begin();
		std::vector<read_id_t>::const_iterator end_b = a->read_list.end();
		for (; b != end_b; ++b) {
			++score[*b];
		}
	}
	std::vector<ScoreHolder> score_list;
	score_list.reserve(score.size());
	size_t i;
	for (i = 0; i != score.size(); ++i) {
		score_list.push_back(ScoreHolder(score[i], read_to_ext[i]));
	}
	std::sort(score_list.begin(), score_list.end());
	for (i = 0; i != score_list.size(); ++i) {
		score[read_to_int[score_list[i].read]] = i;
	}
	std::vector<ReadScore>::iterator c = best_reads.begin();
	std::vector<ReadScore>::iterator end_c = best_reads.end();
	for (; c != end_c; ++c) {
		std::vector<read_id_t>::const_iterator b = c->read_list.begin();
		std::vector<read_id_t>::const_iterator end_b = c->read_list.end();
		if (b == end_b) {
			continue;
		}
		read_id_t best = *b;
		read_id_t best_score = score[*b];
		for (++b; b != end_b; ++b) {
			if (best_score < score[*b]) {
				best_score = score[*b];
				best = *b;
			}
		}
		c->read_list.clear();
		c->read_list.push_back(best);
	}
}

static void print_output(const std::list<read_id_t> &orphan_count, const std::list<read_id_t> &duplicate_count, const std::string &output_file, bool opt_print_matched_read) {
	std::vector<std::string> read_list(read_names.size());
	std::map<std::string, read_id_t>::const_iterator b = read_names.begin();
	std::map<std::string, read_id_t>::const_iterator end_b = read_names.end();
	for (; b != end_b; ++b) {
		read_list[b->second] = b->first;
	}
	read_names.clear();
	FILE *fp(fopen(std::string(output_file + ".unique").c_str(), "w"));
	if (fp == NULL) {
		fprintf(stderr, "Error: could not write %s.unique: %s\n", output_file.c_str(), strerror(errno));
	} else {
		for (read_id_t i(0); i != read_list.size(); ++i) {
			if (duplicates.find(i) == duplicates.end()) {
				std::string s(read_list[i]);
				const size_t j(s.rfind('-'));
				s.resize(j);
				if (*s.rbegin() == '-') {
					fprintf(fp, "%sR1\n%sR2\n", s.c_str(), s.c_str());
				} else {
					fprintf(fp, "%s/1\n%s/2\n", s.c_str(), s.c_str());
				}
			}
		}
		fclose(fp);
	}
	fp = fopen(std::string(output_file + ".dup").c_str(), "w");
	if (fp == NULL) {
		fprintf(stderr, "Error: could not write %s.dup: %s\n", output_file.c_str(), strerror(errno));
	} else {
		std::map<read_id_t, read_id_t>::const_iterator a(duplicates.begin());
		const std::map<read_id_t, read_id_t>::const_iterator end_a(duplicates.end());
		for (; a != end_a; ++a) {
			std::string s(read_list[a->first]);
			const size_t j(s.rfind('-'));
			s.resize(j);
			if (opt_print_matched_read) {
				if (*s.rbegin() == '-') {
					fprintf(fp, "%sR1\t%s\n%sR2\t%s\n", s.c_str(), read_list[a->second].c_str(), s.c_str(), read_list[a->second].c_str());
				} else {
					fprintf(fp, "%s/1\t%s\n%s/2\t%s\n", s.c_str(), read_list[a->second].c_str(), s.c_str(), read_list[a->second].c_str());
				}
			} else {
				if (*s.rbegin() == '-') {
					fprintf(fp, "%sR1\n%sR2\n", s.c_str(), s.c_str());
				} else {
					fprintf(fp, "%s/1\n%s/2\n", s.c_str(), s.c_str());
				}
			}
		}
		fclose(fp);
	}
	printf("orphans:");
	std::list<read_id_t>::const_iterator a(orphan_count.begin());
	std::list<read_id_t>::const_iterator end_a(orphan_count.end());
	for (; a != end_a; ++a) {
		printf(" %u", *a);
	}
	printf("\nduplicates:");
	a = duplicate_count.begin();
	end_a = duplicate_count.end();
	for (; a != end_a; ++a) {
		printf(" %u", *a);
	}
	printf("\nunique: %lu\n", 2 * (read_list.size() - duplicates.size()));
}

static void print_dups(const std::string &output_file) {
	std::vector<std::string> read_list(read_names.size());
	std::map<std::string, read_id_t>::const_iterator b(read_names.begin());
	const std::map<std::string, read_id_t>::const_iterator end_b(read_names.end());
	for (; b != end_b; ++b) {
		read_list[b->second] = b->first;
	}
	read_names.clear();
	FILE * const fp(output_file == "-" ? stdout : fopen(std::string(output_file).c_str(), "w"));
	if (fp == NULL) {
		fprintf(stderr, "Error: could not write %s: %s\n", output_file.c_str(), strerror(errno));
		return;
	}
	std::map<read_id_t, read_id_t>::const_iterator a(duplicates.begin());
	const std::map<read_id_t, read_id_t>::const_iterator end_a(duplicates.end());
	for (; a != end_a; ++a) {
		fprintf(fp, "%s\n", read_list[a->first].c_str());
	}
	fclose(fp);
}

int main(int argc, char **argv) {
	bool opt_print_matched_read, opt_raw, opt_recurse;
	std::string output_file;
	std::list<std::string> blat_files;
	get_opts(argc, argv, output_file, opt_print_matched_read, opt_recurse, opt_raw, blat_files);
	for (int n = optind; n != argc; ++n) {
		blat_files.push_back(argv[n]);
	}
	std::list<read_id_t> orphan_count, duplicate_count;
	while (orphan_count.empty() || (opt_recurse && orphan_count.size() != 100)) {
		std::list<std::string>::const_iterator a = blat_files.begin();
		const std::list<std::string>::const_iterator end_a = blat_files.end();
		for (; a != end_a; ++a) {
			parse_output(*a);
		}
		reduce_sets();
		read_to_int.clear();
		std::map<read_id_t, read_id_t> new_dups;	// external #s
		size_t orphans = 0;
		read_id_t i;
		for (i = 0; i != best_reads.size(); ++i) {
			read_id_t j = best_reads[i].read_list[0];
			if (i != j) {
				new_dups[read_to_ext[i]] = read_to_ext[j];
				if (j != best_reads[j].read_list[0]) {
					++orphans;
				}
			}
		}
		if (new_dups.empty()) {
			opt_recurse = 0;
		} else {
			orphan_count.push_back(2 * orphans);
			duplicate_count.push_back(2 * (new_dups.size() - orphans));
			duplicates.insert(new_dups.begin(), new_dups.end());
		}
		best_reads.clear();
		read_to_ext.clear();
	}
	if (opt_raw) {
		print_dups(output_file);
	} else {
		print_output(orphan_count, duplicate_count, output_file, opt_print_matched_read);
	}
	return 0;
}
