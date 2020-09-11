#include "breakup_line.h"	// breakup_line()
#include "hash_read_hits.h"	// hash_read_hits::read_type
#include "hist_lib_hash.h"	// init_mer_constants(), opt_mer_length
#include "kmer_lookup_info.h"	// KmerLookupInfo
#include "open_compressed.h"	// close_compressed(), open_compressed()
#include "write_fork.h"	// close_fork(), pfputc(), pfputs(), write_fork()
#include <exception>	// exception
#include <iostream>	// cerr, cout
#include <list>		// list<>
#include <map>		// map<>
#include <readline/history.h>	// add_history(), using_history()
#include <readline/readline.h>	// readline(), rl_attempted_completion_function, rl_completion_func_t, rl_completion_matches(), rl_readline_name
#include <sstream>	// istringstream
#include <stdio.h>	// EOF
#include <stdlib.h>	// free()
#include <string.h>	// strdup(), strlen(), strncmp()
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
	std::cerr << "usage: kmer_matching kmer_index_file\n";
}

class Selection {
    public:
	int match_value_min;
	hash_read_hits::value_type kmer_hit_max;
	std::string search_sequence;
	std::map<hash_read_hits::read_type, int> read_hits;
	// kmer_hit_max is unsigned, so -1 should be the max value
	Selection() : match_value_min(0), kmer_hit_max(-1) { };
	~Selection() { };
	void print_hits(const KmerLookupInfo &kmers) const {
		// make reverse map to we can order output by match value
		// (we also filter by match value cutoffs at this point)
		// using list here for quick allocation and low memory usage
		std::map<int, std::list<hash_read_hits::read_type> > list;
		std::map<hash_read_hits::read_type, int>::const_iterator a(read_hits.begin());
		const std::map<hash_read_hits::read_type, int>::const_iterator end_a(read_hits.end());
		for (; a != end_a; ++a) {
			if (a->second >= match_value_min) {
				list[a->second].push_back(a->first);
			}
		}
		if (list.empty()) {
			std::cout << "No matches in selection\n";
			return;
		}
		std::map<int, std::list<hash_read_hits::read_type> >::const_iterator b(list.begin());
		const std::map<int, std::list<hash_read_hits::read_type> >::const_iterator end_b(list.end());
		for (; b != end_b; ++b) {
			std::list<hash_read_hits::read_type>::const_iterator c(b->second.begin());
			const std::list<hash_read_hits::read_type>::const_iterator end_c(b->second.end());
			for (; c != end_c; ++c) {
				std::cout << kmers.read_name(*c) << " " << (double(b->first) / (search_sequence.size() - opt_mer_length)) << "\n";
			}
		}
		std::cout << "\n";
	}
	int write_hits(const KmerLookupInfo &kmers, const std::string &file) const {
		const int fd(write_fork(file.c_str()));
		if (fd == -1) {
			std::cerr << "Error: write failed\n";
			return -1;
		}
		int total_written(0);
		std::map<hash_read_hits::read_type, int>::const_iterator a(read_hits.begin());
		const std::map<hash_read_hits::read_type, int>::const_iterator end_a(read_hits.end());
		for (; a != end_a; ++a) {
			if (a->second >= match_value_min) {
				++total_written;
				pfputs(fd, kmers.read_name(a->first));
				pfputc(fd, '\n');
			}
		}
		close_fork(fd);
		return total_written;
	}
};

static void help_function(const std::vector<std::string> &list, const KmerLookupInfo &kmers, Selection &selection) {
	(void)list;		// avoid warnings
	(void)kmers;
	(void)selection;
	std::cout <<
		"help                    this text\n"
		"quit                    quit program\n"
		"search ##               search index for matches against given sequence\n"
		"set kmer_hit_max ##     set maximum hit count for kmers\n"
		"                        (kmers with more matches than this will be ignored)\n"
		"set match_value_min ##  set minimum match value\n"
		"show                    show current cutoffs\n"
		"write ##                write current search results to given file\n";
}

static void search_function(const std::vector<std::string> &list, const KmerLookupInfo &kmers, Selection &selection) {
	if (list.size() > 2) {
		std::cout << "Error: search takes one parameter (the sequence to match against)\n";
		return;
	// +1 as opt_mer_length is one less than the set mer length
	} else if (list.size() == 2 && list[1].size() < opt_mer_length + 1) {
		std::cout << "Error: search sequence is too short; need to be at least " << (opt_mer_length + 1) << " basepairs long\n";
		return;
	}
	if (list.size() == 2) {
		selection.search_sequence = list[1];
	}
	selection.read_hits.clear();
	count_read_hits(selection.search_sequence, kmers, selection.read_hits, selection.kmer_hit_max);
	if (selection.read_hits.empty()) {
		std::cout << "No matching reads found\n";
	} else {
		selection.print_hits(kmers);
	}
}

static void set_function(const std::vector<std::string> &list, const KmerLookupInfo &kmers, Selection &selection) {
	(void)kmers;
	(void)selection;
	if (list.size() != 3) {
		std::cout << "Error: set takes two parameters (variable and value)\n";
	} else if (list[1].compare("match_value_min") == 0) {
		std::istringstream(list[2]) >> selection.match_value_min;
	} else if (list[1].compare("kmer_hit_max") == 0) {
		std::istringstream(list[2]) >> selection.kmer_hit_max;
	} else {
		std::cout << "Error: you can only set match_value_min and kmer_hit_max\n";
	}
}

static void show_function(const std::vector<std::string> &list, const KmerLookupInfo &kmers, Selection &selection) {
	(void)list;
	(void)kmers;
	std::cout << "match_value_min " << selection.match_value_min << "\n" <<
		"kmer_hit_max " << selection.kmer_hit_max << "\n\n";
}

static void write_function(const std::vector<std::string> &list, const KmerLookupInfo &kmers, Selection &selection) {
	if (list.size() != 2) {
		std::cout << "Error: write only takes one parameter (the filename to write to)\n";
		return;
	}
	selection.write_hits(kmers, list[1]);
}

struct ActionType {
	const char * const action;
	void (*function)(const std::vector<std::string> &, const KmerLookupInfo &, Selection &);
};

// actions ordered by likelihood of coming up, with most common ones first
static ActionType actions[] = {
	{ "search", search_function },
	{ "set", set_function },
	{ "write", write_function },
	{ "show", show_function },
	{ "?", help_function },
	{ "help", help_function },
	{ "quit", 0 },
};

static const size_t actions_size(sizeof(actions) / sizeof(ActionType));

// command completion routines (for use by readline)

// have to malloc (through strdup(), here) returned string
static char *command_generator(const char * const text, const int state) {
	static int list_index, len;
	// if first pass, initialize
	if (state == 0) {
		list_index = 0;
		len = strlen(text);
	}
	// advance to next partial match, if any
	while (list_index != actions_size) {
		const char * const s(actions[list_index++].action);
		if (strncmp(s, text, len) == 0) {
			return strdup(s);
		}
	}
	return 0;	// no more matches
}

// only activate command completion if text is at the start of the line
static char **command_completion(const char * const text, const int start, const int end) {
	(void)end;	// avoid warning
	return start == 0 ? rl_completion_matches(text, command_generator) : 0;
}

static void user_input_loop(const KmerLookupInfo &kmers) {
	// take list out of loop to avoid extra allocations/deallocations
	std::vector<std::string> list;
	Selection selection;
	char *s;
	while ((s = readline("kmers> ")) != 0) {
		list.clear();
		breakup_line(s, list);
		// hold onto s to add to history, if needed
		if (list.empty()) {
			free(s);
			continue;
		}
		size_t i(0);
		for (; i != actions_size; ++i) {
			if (list[0].compare(actions[i].action) == 0) {
				if (actions[i].function != 0) {
					add_history(s);
					actions[i].function(list, kmers, selection);
					break;
				}
				free(s);
				return;		// quit action
			}
		}
		free(s);
		if (i == actions_size) {
			std::cout << "invalid command: " << list[0] << "\n";
		}
	}
}

int main(int argc, char **argv) {
	int had_error(0);
	if (argc != 2) {
		print_usage();
		return 1;
	}
	try {
		// set up readline history and command completion
		using_history();
		rl_readline_name = "kmer_matching";
		rl_attempted_completion_function = (rl_completion_func_t *)command_completion;
		KmerLookupInfo kmers;
		const int fd(open_compressed(argv[1]));
		if (fd == -1) {
			throw LocalException("could not open " + std::string(argv[1]));
		} else if (fd == STDIN_FILENO) {
			throw LocalException("can not read kmer index file from stdin");
		}
		std::cout << "Reading kmer index file\n";
		kmers.restore(fd);
		close_compressed(fd);
		opt_mer_length = kmers.mer_length();
		init_mer_constants();
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
