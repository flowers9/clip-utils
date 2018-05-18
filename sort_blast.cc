#include "breakup_line.h"	// breakup_line()
#include "open_compressed.h"	// close_compressed(), open_compressed()
#include <exception>	// exception
#include <getopt.h>	// getopt(), optarg, optind
#include <iostream>	// cerr, cout
#include <map>		// map<>, multimap<>
#include <set>		// set<>
#include <sstream>	// istringstream
#include <stdio.h>	// EOF
#include <string>	// string
#include <utility>	// make_pair()
#include <vector>	// vector<>

class LocalException : public std::exception {
    private:
	const std::string error_;
	const int show_usage_;
    public:
	explicit LocalException(const std::string &s) : error_(s), show_usage_(0) { }
	LocalException(const std::string &s, int i) : error_(s), show_usage_(i) { }
	virtual ~LocalException(void) throw() { }
	virtual const char * what() const throw() {
		return error_.c_str();
	}
	const int &show_usage(void) const throw() {
		return show_usage_;
	}
};

static int largest_indel, max_matches, max_mismatches, minimum_score;

static void print_usage() {
	std::cerr <<
		"usage: sort_blast [opts] <blast_file>\n" <<
		"\t-i ##\tlargest index [5]\n" <<
		"\t-M ##\tmaximum matches [6]\n" <<
		"\t-m ##\tmaximum mismatches [3]\n" <<
		"\t-s ##\tminimum score [92]\n";
	}

static void get_opts(int argc, char **argv) {
	largest_indel = 5;
	max_matches = 6;
	max_mismatches = 3;
	minimum_score = 92;
	int c;
	while ((c = getopt(argc, argv, "i:M:m:s:")) != EOF) {
		switch (c) {
		    case 'i':
			std::istringstream(optarg) >> largest_indel;
			break;
		    case 'M':
			std::istringstream(optarg) >> max_matches;
			break;
		    case 'm':
			std::istringstream(optarg) >> max_mismatches;
			break;
		    case 's':
			std::istringstream(optarg) >> minimum_score;
			break;
		    default:
			throw LocalException("bad option: " + static_cast<char>(c), 1);
		}
	}
	if (optind == argc) {
		throw LocalException("no files specified", 1);
	}
}

static int count_good_alignments(const char * const filename, std::map<std::string, size_t> &alignment_count) {
	const int fd(open_compressed(filename));
	if (fd == -1) {
		return -1;
	}
	std::string line;
	while (pfgets(fd, line) != EOF) {
		std::vector<std::string> list;
		breakup_line(line, list);
		// the 17 requirement comes from the second pass
		if (list.size() <= 17) {
			continue;
		}
		int i;
		std::istringstream(list[0]) >> i;
		if (i <= minimum_score) {
			continue;
		}
		std::istringstream(list[1]) >> i;
		if (i >= max_mismatches) {
			continue;
		}
		std::istringstream(list[7]) >> i;
		if (i >= largest_indel) {
			continue;
		}
		++alignment_count[list[9]];
	}
	close_compressed(fd);
	// remove name1 reads with too many matches
	std::map<std::string, size_t>::iterator a(alignment_count.begin());
	const std::map<std::string, size_t>::const_iterator end_a(alignment_count.end());
	while (a != end_a) {
		if (a->second == 1 || a->second > static_cast<size_t>(max_matches)) {
			alignment_count.erase(a++);
		} else {
			++a;
		}
	}
	return 0;
}

class Match {
    public:	// columns 0, 15, and 16 (9 and 13 are stored in the lookup)
	int score_, start_, end_;
	Match(const int score, const int start, const int end) : score_(score), start_(start), end_(end) { }
	~Match(void) { }
};


static int get_good_alignments(const char * const filename, const std::map<std::string, size_t> &alignment_count, std::map<std::string, std::multimap<std::string, Match> > &matches, std::map<std::string, std::multimap<int, std::string> > &match_lookup) {
	const int fd(open_compressed(filename));
	if (fd == -1) {
		return -1;
	}
	std::string line;
	while (pfgets(fd, line) != EOF) {
		std::vector<std::string> list;
		breakup_line(line, list);
		if (list.size() <= 17) {
			continue;
		}
		int score;
		std::istringstream(list[0]) >> score;
		if (score <= minimum_score) {
			continue;
		}
		int i;
		std::istringstream(list[1]) >> i;
		if (i >= max_mismatches) {
			continue;
		}
		std::istringstream(list[7]) >> i;
		if (i >= largest_indel) {
			continue;
		}
		if (alignment_count.find(list[9]) == alignment_count.end()) {
			continue;
		}
		int start, end;
		std::istringstream(list[15]) >> start;
		std::istringstream(list[16]) >> end;
		matches[list[9]].insert(std::make_pair(list[13], Match(score, start, end)));
		match_lookup[list[13]].insert(std::make_pair(start, list[9]));
	}
	close_compressed(fd);
	return 0;
}

static int read_file(const char * const filename, std::map<std::string, std::multimap<std::string, Match> > &matches, std::map<std::string, std::multimap<int, std::string> > &match_lookup) {
	std::map<std::string, size_t> alignment_count;
	if (count_good_alignments(filename, alignment_count) == -1) {
		return -1;
	}
	if (get_good_alignments(filename, alignment_count, matches, match_lookup) == -1) {
		return -1;
	}
	return 0;
}

// in the event a name2 matches multiple times against a name1, choose the
// best match (by score) to be the first one shown

static void winnow_match_lookup(const std::map<std::string, std::multimap<std::string, Match> > &matches, std::map<std::string, std::multimap<int, std::string> > &match_lookup) {
	std::map<std::string, std::multimap<int, std::string> >::iterator a(match_lookup.begin());
	const std::map<std::string, std::multimap<int, std::string> >::const_iterator end_a(match_lookup.end());
	for (; a != end_a; ++a) {
		// find name1s that name2 has multiple alignments to
		std::map<std::string, size_t> count;
		std::multimap<int, std::string>::iterator b(a->second.begin());
		const std::multimap<int, std::string>::const_iterator end_b(a->second.end());
		for (; b != end_b; ++b) {
			++count[b->second];
		}
		std::map<std::string, size_t>::const_iterator c(count.begin());
		const std::map<std::string, size_t>::const_iterator end_c(count.end());
		for (; c != end_c; ++c) {
			if (c->second == 1) {
				continue;
			}
			// find all name1 vs name2 alignments, choose best
			// (highest score, lowest start, highest end)
			const std::multimap<std::string, Match> &list(matches.find(c->first)->second);
			std::multimap<std::string, Match>::const_iterator d(list.lower_bound(a->first));
			const std::multimap<std::string, Match>::const_iterator end_d(list.upper_bound(a->first));
			std::multimap<std::string, Match>::const_iterator best(d);
			for (++d; d != end_d; ++d) {
				if (best->second.score_ < d->second.score_ || (best->second.score_ == d->second.score_ && (best->second.start_ > d->second.start_ || (best->second.start_ == d->second.start_ && best->second.end_ < d->second.end_)))) {
					best = d;
				}
			}
			// now remove all other name2-name1 match lookups
			b = a->second.begin();
			while (b != end_b) {
				if (c->first.compare(b->second) == 0 && best->second.start_ != b->first) {
					a->second.erase(b++);
				} else {
					++b;
				}
			}
		}
	}
}

static void print_pairs(const std::string &name1, const std::string &name2, const std::multimap<std::string, Match> &matches, const int start) {
	std::cout << name1;
	// pull out the match we're looking at - name1 vs name2 at start
	// (if there are more than one, go with the best score, then higher end)
	std::multimap<std::string, Match>::const_iterator a(matches.lower_bound(name2));
	std::multimap<std::string, Match>::const_iterator end_a(matches.upper_bound(name2));
	std::multimap<std::string, Match>::const_iterator best(matches.end());
	for (; a != end_a; ++a) {
		if (a->second.start_ == start && (best == matches.end() || best->second.score_ < a->second.score_ || (best->second.score_ == a->second.score_ && best->second.end_ < a->second.end_))) {
			best = a;
		}
	}
	const Match &first(best->second);
	std::cout << "\t" << first.score_ << "\t" << name2 << "\t" << first.start_ << "\t" << first.end_;
	// now print out all other matches against name1
	a = matches.begin();
	end_a = matches.end();
	for (; a != end_a; ++a) {
		if (a != best) {
			std::cout << "\t" << a->second.score_ << "\t" << a->first << "\t" << a->second.start_ << "\t" << a->second.end_;
		}
	}
	std::cout << "\n";
}

// go through all alignments in order by name2; for each one, pull out
// all alignments against them and sort them by start position, then list
// each one per line followed by all the other alignments against that
// same name1, sorted by name2

static void process_data(const std::map<std::string, std::multimap<std::string, Match> > &matches, const std::map<std::string, std::multimap<int, std::string> > &match_lookup) {
	std::map<std::string, std::multimap<int, std::string> >::const_iterator a(match_lookup.begin());
	const std::map<std::string, std::multimap<int, std::string> >::const_iterator end_a(match_lookup.end());
	for (; a != end_a; ++a) {
		std::cout << a->first << ":\n";
		std::multimap<int, std::string>::const_iterator b(a->second.begin());
		const std::multimap<int, std::string>::const_iterator end_b(a->second.end());
		for (; b != end_b; ++b) {
			print_pairs(b->second, a->first, matches.find(b->second)->second, b->first);
		}
	}
}

int main(int argc, char **argv) {
	int had_error(0);
	try {
		get_opts(argc, argv);
		std::map<std::string, std::multimap<std::string, Match> > matches;
		std::map<std::string, std::multimap<int, std::string> > match_lookup;
		if (read_file(argv[optind], matches, match_lookup) == -1) {
			had_error = 1;
		} else {
			winnow_match_lookup(matches, match_lookup);
			process_data(matches, match_lookup);
		}
	} catch (std::exception &e) {
		std::cerr << "Error: " << e.what() << "\n";
		LocalException *x(dynamic_cast<LocalException *>(&e));
		if (x != NULL && x->show_usage()) {
			print_usage();
		}
		had_error = 1;
	}
	return had_error;
}
