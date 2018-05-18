#include "breakup_line.h"	/* breakup_line() */
#include "open_compressed.h"	/* close_compressed(), open_compressed(), pfgets() */
#include "strtostr.h"	/* strtostr_exact() */
#include "version.h"	/* VERSION */
#include <algorithm>	/* sort() */
#include <getopt.h>	// getopt(), optarg, optind
#include <map>		/* map<> */
#include <stdio.h>	/* EOF, fprintf(), printf(), stderr */
#include <stdlib.h>	/* atoi(), exit(), strtod() */
#include <string>	/* string */
#include <vector>	/* vector<> */

#define INSERT_LENGTH 48
static double opt_read_identity;
static int opt_read_offset;
static std::vector<std::string> read_names;
static std::map<std::string, size_t> read_name_to_int;

// a is assumed to be sorted; returns position at which to insert value
static ssize_t not_in_list(const std::vector<size_t> &list, size_t x) {
	size_t i = 0;
	size_t n = list.size();
	while (i != n) {
		size_t j = (i + n) / 2;
		if (x < list[j]) {
			n = j;
		} else if (x > list[j]) {
			i = j + 1;
		} else {
			return -1;
		}
	}
	return i;
}

class Read_Score {
    public:
	int match_length;
	int identity;
	int read_length;
	std::vector<size_t> read_list;	// sorted list
	Read_Score(void) : match_length(0), identity(0), read_length(0) { }
	Read_Score(int __i, int __j, int __k, int __l, size_t __s, size_t __t) : match_length(__i), identity(__j) {
		if (__k > __l) {
			read_length = __k;
			read_list.push_back(__s);
		} else {
			read_length = __l;
			if (__k != __l) {
				read_list.push_back(__t);
			} else if (__s < __t) {
				read_list.push_back(__s);
				read_list.push_back(__t);
			} else {
				read_list.push_back(__t);
				read_list.push_back(__s);
			}
		}
	}
	~Read_Score(void) { }
	int cmp(const Read_Score &__a) const {
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
	void add(int __i) {
		ssize_t __n = not_in_list(read_list, __i);
		if (__n != -1) {
			ssize_t __j = read_list.size();
			read_list.push_back(0);		// placeholder
			// shuffle array to open the correct position
			for (; __j != __n; --__j) {
				read_list[__j] = read_list[__j - 1];
			}
			// insert value
			read_list[__j] = __i;
		}
	}
};

static std::vector<Read_Score> best_reads;

static void print_usage() {
	fprintf(stderr, "usage: parse_output [opts] <blat_file1> <blat_file2> ...\n");
	fprintf(stderr, "\t-I\tmatch identity [.98]\n");
	fprintf(stderr, "\t-O\tmatch offset [2]\n");
	exit(0);
}

static void get_opts(int argc, char **argv) {
	/* set option defaults */
	opt_read_identity = .98;
	opt_read_offset = 2;
	/* read in options */
	double x;
	int c;
	while ((c = getopt(argc, argv, "hI:O:V")) != EOF) {
		switch (c) {
		    case 'h':
			print_usage();
		    case 'I':
			x = strtod(optarg, NULL);
			if (x < 0 || 1 < x) {
				fprintf(stderr, "Error: read identity is out of range [0,1]: %lf\n", x);
				print_usage();
			}
			opt_read_identity = x;
			break;
		    case 'O':
			c = atoi(optarg);
			if (c < 0) {
				fprintf(stderr, "Error: read offset is negative: %d\n", c);
				print_usage();
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
	if (optind == argc) {
		fprintf(stderr, "Error: no blat files given\n");
		print_usage();
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

// make sure at least one basepair is within read_offset
// basepairs of each end of the insert is contained in a block
static int not_contains(int offset, const std::string &s, const std::vector<int> &lengths) {
	int x = offset - opt_read_offset - 1;
	int y = offset - 1;
	std::vector<int> starts;
	read_list(s, starts);
	size_t i;
	for (i = 0; i != starts.size(); ++i) {
		if (starts[i] > y) {
			return 1;
		} else if (x < starts[i] + lengths[i]) {
			break;
		}
	}
	x = offset + INSERT_LENGTH;
	y = offset + INSERT_LENGTH + opt_read_offset;
	for (; i != starts.size(); ++i) {
		if (starts[i] > y) {
			return 1;
		} else if (x < starts[i] + lengths[i]) {
			break;
		}
	}
	return i == starts.size();
}

static size_t find_read(const std::string &s) {
	std::map<std::string, size_t>::const_iterator a = read_name_to_int.find(s);
	if (a == read_name_to_int.end()) {
		size_t i = read_name_to_int[s] = read_names.size();
		read_names.push_back(s);
		return i;
	} else {
		return a->second;
	}
}

static void update_score(size_t read, const Read_Score &score) {
	if (read == best_reads.size()) {
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

static void parse_output(char *blat_file) {
	int fd = open_compressed(blat_file);
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
		int target_size = atoi(fields[14].c_str());
		int length = std::min(query_size, target_size) - n_count;
		int identity = atoi(fields[0].c_str()) + atoi(fields[2].c_str());
		// skip if not similar
		if (identity < opt_read_identity * length) {
			continue;
		}
		int query_offset = atoi(fields[9].substr(fields[9].rfind('-') + 1).c_str());
		int target_offset = atoi(fields[13].substr(fields[13].rfind('-') + 1).c_str());
		std::vector<int> blocks;
		read_list(fields[18], blocks);
		if (not_contains(target_offset, fields[20], blocks)) {
			continue;
		}
		if ((fields[8] == "+" && not_contains(query_offset, fields[19], blocks)) || (fields[8] == "-" && not_contains(query_size - query_offset - INSERT_LENGTH, fields[19], blocks))) {
			continue;
		}
		size_t query = find_read(fields[9]);
		size_t target = find_read(fields[13]);
		Read_Score score(length, identity, query_size, target_size, query, target);
		update_score(query, score);
		update_score(target, score);
	}
	close_compressed(fd);
}

class Score_Holder {
    public:
	int score;
	std::string read;
	Score_Holder(int __i, const std::string &__s) : score(__i), read(__s) { }
	~Score_Holder(void) { }
	bool operator<(const Score_Holder &__a) const {
		if (score != __a.score) {
			return score < __a.score;
		} else {
			return read < __a.read;
		}
	}
};

// reduce list of best reads to single best read for each read
static void reduce_sets() {
	std::vector<int> score(best_reads.size(), 0);
	std::vector<Read_Score>::const_iterator a = best_reads.begin();
	std::vector<Read_Score>::const_iterator end_a = best_reads.end();
	for (; a != end_a; ++a) {
		std::vector<size_t>::const_iterator b = a->read_list.begin();
		std::vector<size_t>::const_iterator end_b = a->read_list.end();
		for (; b != end_b; ++b) {
			++score[*b];
		}
	}
	std::vector<Score_Holder> score_list;
	score_list.reserve(score.size());
	size_t i;
	for (i = 0; i != score.size(); ++i) {
		score_list.push_back(Score_Holder(score[i], read_names[i]));
	}
	std::sort(score_list.begin(), score_list.end());
	for (i = 0; i != score_list.size(); ++i) {
		score[read_name_to_int[score_list[i].read]] = i;
	}
	std::vector<Read_Score>::iterator c = best_reads.begin();
	std::vector<Read_Score>::iterator end_c = best_reads.end();
	for (c = best_reads.begin(); c != end_c; ++c) {
		std::vector<size_t>::const_iterator b = c->read_list.begin();
		std::vector<size_t>::const_iterator end_b = c->read_list.end();
		int best = *b;
		int best_score = score[*b];
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

int main(int argc, char **argv) {
	get_opts(argc, argv);
	for (; optind != argc; ++optind) {
		parse_output(argv[optind]);
	}
	reduce_sets();
	size_t orphans = 0;
	size_t i;
	for (i = 0; i != best_reads.size(); ++i) {
		size_t j = best_reads[i].read_list[0];
		if (i != j && j != best_reads[j].read_list[0]) {
			++orphans;
		}
	}
	printf("%lu\n", orphans);
	for (i = 0; i != best_reads.size(); ++i) {
		size_t j = best_reads[i].read_list[0];
		if (i != j) {
			printf("%s %s\n", read_names[i].c_str(), read_names[j].c_str());
		}
	}
	return 0;
}
