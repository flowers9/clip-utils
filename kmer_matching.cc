#include "breakup_line.h"	// breakup_line()
#include "hash_read_hits.h"	// hash_read_hits::read_type
#include "hist_lib_hash.h"	// init_mer_constants(), opt_mer_length
#include "kmer_lookup_info.h"	// KmerLookupInfo
#include "open_compressed.h"	// close_compressed(), open_compressed()
#include <exception>	// exception
#include <getopt.h>	// getopt(), optarg, optind
#include <iostream>	// cerr, cout
#include <limits.h>	// INT32_MAX
#include <map>		// map<>
#include <readline/history.h>	// might be needed to enable history functions?
#include <readline/readline.h>	// readline()
#include <sstream>	// istringstream
#include <stdio.h>	// EOF
#include <stdlib.h>	// free()
#include <string>	// string
#include <unistd.h>	// STDIN_FILENO
#include <vector>	// vector<>

class LocalException : public std::exception {
    private:
	const std::string error_;
	const int show_usage_;
    public:
	explicit LocalException(const std::string &s) : error_(s), show_usage_(0) { }
	explicit LocalException(const std::string &s, const int i) : error_(s), show_usage_(i) { }
	virtual ~LocalException(void) throw() { }
	virtual const char * what() const throw() {
		return error_.c_str();
	}
	const int &show_usage(void) const throw() {
		return show_usage_;
	}
};

static void print_usage() {
	std::cerr <<
		"usage: kmer_matching [options] kmer_index_file\n"
		"    -h    print this help\n"
		"    -m ## set kmer length (1-32) [24]\n";
}

static void get_opts(int argc, char **argv) {
	opt_mer_length = 24;
	int c;
	while ((c = getopt(argc, argv, "hm:")) != EOF) {
		switch (c) {
		    case 'h':
			throw LocalException("", 1);
			break;
		    case 'm':
			std::istringstream(optarg) >> opt_mer_length;
			if (opt_mer_length < 1 || 32 < opt_mer_length) {
				throw LocalException("bad mer length", 1);
			}
			break;
		    default:
			throw LocalException("bad option: " + static_cast<char>(c), 1);
		}
	}
	if (optind + 1 != argc) {
		throw LocalException("missing kmer index file", 1);
	}
}

class Selection {
    public:
	// cutoffs are normalized to the selection length
	int selection_length;
	// cutoffs don't limit was goes into read_hits, but should limit what comes out
	int lower_cutoff, upper_cutoff;
	std::map<hash_read_hits::read_type, int> read_hits;
	Selection() : selection_length(0), lower_cutoff(0), upper_cutoff(INT32_MAX) { };
	~Selection() { };
	void print_hits() const {
		// XXX
	}
};

static void help_function(const std::vector<std::string> &list, const KmerLookupInfo &kmers, Selection &selection) {
	(void)list;
	(void)kmers;
	(void)selection;
	std::cout <<
		"help              this text\n"
		"quit              quit program\n"
		"search ##         search index for matches against given sequence\n"
		"set bounds[##:##] set upper and/or lower cutoffs for coverage amount\n"
		"show bounds       show current bounds\n"
		"write ##          write current search results to given file\n";
}

// XXX - functions are just placeholder stubs for the moment

static void search_function(const std::vector<std::string> &list, const KmerLookupInfo &kmers, Selection &selection) {
	if (list.size() != 2) {
		std::cout << "Error: search only takes one parameter (the sequence to match against)\n";
		return;
	}
	count_read_hits(list[1], kmers, selection.read_hits);
	if (selection.read_hits.empty()) {
		std::cout << "No matching reads found\n";
	} else {
		selection.print_hits();
	}
}

static void set_function(const std::vector<std::string> &list, const KmerLookupInfo &kmers, Selection &selection) {
	(void)list;
	(void)kmers;
	(void)selection;
	std::cout << "set";
}

static void show_function(const std::vector<std::string> &list, const KmerLookupInfo &kmers, Selection &selection) {
	(void)list;
	(void)kmers;
	(void)selection;
	std::cout << "show";
}

static void write_function(const std::vector<std::string> &list, const KmerLookupInfo &kmers, Selection &selection) {
	(void)kmers;
	(void)selection;
	if (list.size() != 2) {
		std::cout << "Error: write only takes one parameter (the filename to write to)\n";
		return;
	}
	std::cout << "write\n";
}

struct ActionType {
	const char * const action;
	void (*function)(const std::vector<std::string> &, const KmerLookupInfo &, Selection &);
};

static ActionType actions[] = {
	{ "help", help_function },
	{ "quit", 0 },
	{ "search", search_function },
	{ "set", set_function },
	{ "show", show_function },
	{ "write", write_function },
};

static void user_input_loop(const KmerLookupInfo &kmers) {
	const size_t possible_actions(sizeof(actions) / sizeof(ActionType));
	// take list out of loop to avoid extra allocations/deallocations
	std::vector<std::string> list;
	Selection selection;
	char *s;
	while ((s = readline("kmers> ")) != 0) {
		list.clear();
		breakup_line(s, list);
		free(s);
		if (list.empty()) {
			continue;
		}
		size_t i(0);
		for (; i != possible_actions; ++i) {
			if (list[0].compare(actions[i].action) == 0) {
				if (actions[i].function != 0) {
					actions[i].function(list, kmers, selection);
					break;
				}
				return;		// quit action
			}
		}
		if (i == possible_actions) {
			std::cout << "invalid command: " << list[0] << "\n";
		}
	}
}

int main(int argc, char **argv) {
	int had_error(0);
	try {
		get_opts(argc, argv);
		init_mer_constants();
		KmerLookupInfo kmers;
		const int fd(open_compressed(argv[optind]));
		if (fd == -1) {
			throw LocalException("could not open " + std::string(argv[optind]));
		} else if (fd == STDIN_FILENO) {
			throw LocalException("can not read kmer index file from stdin");
		}
		std::cout << "Reading kmer index file\n";
		kmers.restore(fd);
		close_compressed(fd);
		user_input_loop(kmers);
	} catch (std::exception &e) {
		if (e.what()[0] != 0) {
			std::cerr << "Error: " << e.what() << "\n";
		}
		LocalException * const x(dynamic_cast<LocalException *>(&e));
		if (x != NULL && x->show_usage()) {
			print_usage();
		}
		had_error = 1;
	}
	return had_error;
}
