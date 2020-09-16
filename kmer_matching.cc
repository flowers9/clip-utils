#include "breakup_line.h"	// breakup_line()
#include "get_name.h"	// get_name()
#include "hash_read_hits.h"	// hash_read_hits::read_type
#include "hist_lib_hash.h"	// init_mer_constants(), opt_mer_length
#include "kmer_lookup_info.h"	// KmerLookupInfo
#include "open_compressed.h"	// close_compressed(), open_compressed()
#include "write_fork.h"	// close_fork(), pfputc(), pfputs(), write_fork()
#include <cassert>	// assert()
#include <exception>	// exception
#include <iostream>	// cerr, cout, fixed
#include <list>		// list<>
#include <map>		// map<>
#include <math.h>	// ceil()
#include <readline/history.h>	// add_history(), using_history()
#include <readline/readline.h>	// readline(), rl_attempted_completion_function, rl_completion_func_t, rl_completion_matches(), rl_readline_name
#include <sstream>	// istringstream
#include <stdio.h>	// EOF
#include <stdlib.h>	// free(), malloc()
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

static int in_multiline_search(0);

static void print_usage() {
	std::cerr << "usage: kmer_matching <kmer_index_file> [fastq_file]\n";
}

// names are held in KmerLookupInfo, so don't need to duplicate them here
class RawRead {
    public:
	std::string sequence, quality;
	RawRead() { }
	~RawRead() { }
};

class Selection {
    public:
	double match_value_min;
	hash_read_hits::value_type kmer_hit_max;
	std::string search_sequence;
	std::map<hash_read_hits::read_type, int> read_hits;
	// kmer_hit_max is unsigned, so -1 should be the max value
	Selection() : match_value_min(0), kmer_hit_max(-1) { }
	~Selection() { }
	void print_hits(const KmerLookupInfo &kmers) const {
		// normalize the match value to search length
		const double x(search_sequence.size() - opt_mer_length);
		const int min_match(ceil(match_value_min * x));
		// make reverse map to we can order output by match value
		// (we also filter by match value cutoffs at this point)
		// using list here for quick allocation and low memory usage
		std::map<int, std::list<hash_read_hits::read_type> > list;
		std::map<hash_read_hits::read_type, int>::const_iterator a(read_hits.begin());
		const std::map<hash_read_hits::read_type, int>::const_iterator end_a(read_hits.end());
		for (; a != end_a; ++a) {
			if (a->second >= min_match) {
				list[a->second].push_back(a->first);
			}
		}
		if (list.empty()) {
			std::cout << "No matches in selection\n";
			return;
		}
		// go from highest match to lowest
		std::cout << "\n";
		std::map<int, std::list<hash_read_hits::read_type> >::const_reverse_iterator b(list.rbegin());
		const std::map<int, std::list<hash_read_hits::read_type> >::const_reverse_iterator end_b(list.rend());
		for (; b != end_b; ++b) {
			std::list<hash_read_hits::read_type>::const_iterator c(b->second.begin());
			const std::list<hash_read_hits::read_type>::const_iterator end_c(b->second.end());
			for (; c != end_c; ++c) {
				std::cout << kmers.read_name(*c) << " " << (b->first / x) << "\n";
			}
		}
		std::cout << "\n";
	}
	size_t write_hits(const KmerLookupInfo &kmers, const std::string &file, const std::vector<RawRead> &reads) const {
		const int fd(write_fork(file.c_str()));
		if (fd == -1) {
			std::cerr << "Error: write failed\n";
			return -1;
		}
		size_t total_written(0);
		std::map<hash_read_hits::read_type, int>::const_iterator a(read_hits.begin());
		const std::map<hash_read_hits::read_type, int>::const_iterator end_a(read_hits.end());
		if (reads.empty()) {		// just write the read names
			for (; a != end_a; ++a) {
				if (a->second >= match_value_min) {
					++total_written;
					pfputs(fd, kmers.read_name(a->first));
					pfputc(fd, '\n');
				}
			}
		} else {			// write fastq entries
			for (; a != end_a; ++a) {
				if (a->second >= match_value_min) {
					++total_written;
					pfputc(fd, '@');
					pfputs(fd, kmers.read_name(a->first));
					pfputc(fd, '\n');
					pfputs(fd, reads[a->first].sequence);
					pfputs(fd, "\n+\n");
					pfputs(fd, reads[a->first].quality);
					pfputc(fd, '\n');
				}
			}
		}
		close_fork(fd);
		return total_written;
	}
};

static void read_kmer_index(const char * const file, KmerLookupInfo &kmers) {
	std::cout << "Reading kmer index file\n";
	const int fd(open_compressed(file));
	if (fd == -1) {
		throw LocalException("could not open " + std::string(file));
	} else if (fd == STDIN_FILENO) {
		throw LocalException("can not read kmer index file from stdin");
	}
	kmers.restore(fd);
	close_compressed(fd);
}

static void read_fastq(const char * const file, std::vector<RawRead> &reads, const KmerLookupInfo &kmers) {
	std::cout << "Reading fastq file\n";
	reads.resize(kmers.read_count());
	const int fd(open_compressed(file));
	if (fd == -1) {
		throw LocalException("could not open " + std::string(file));
	} else if (fd == STDIN_FILENO) {
		throw LocalException("can not read fastq file from stdin");
	}
	size_t i(0);
	std::string line;
	for (; pfgets(fd, line) != -1; ++i) {
		if (line[0] != '@') {
			throw LocalException("bad header line in fastq file");
		}
		const std::string read_name(get_name(line));
		if (read_name.compare(kmers.read_name(i)) != 0) {
			throw LocalException("mismatched read names: " + read_name + " != " + kmers.read_name(i));
		}
		if (pfgets(fd, reads[i].sequence) == -1) {
			throw LocalException("truncated read: " + line);
		}
		if (pfgets(fd, line) == -1) {
			throw LocalException("truncated read: missing quality header");
		}
		if (pfgets(fd, reads[i].quality) == -1) {
			throw LocalException("truncated read: missing quality");
		}
		assert(reads[i].sequence.size() == reads[i].quality.size());
	}
	close_compressed(fd);
	assert(i == reads.size());
}

static void help_function(const std::vector<std::string> &list, const KmerLookupInfo &kmers, Selection &selection, const std::vector<RawRead> &reads) {
	(void)list;		// avoid warnings
	(void)kmers;
	(void)selection;
	(void)reads;
	std::cout <<
		"help                    this text\n"
		"quit                    quit program\n"
		"msearch ##              multi-line search index for matches against given sequence\n"
		"                        (following lines continue sequence until blank line)\n"
		"search ##               search index for matches against given sequence\n"
		"set kmer_hit_max ##     set maximum hit count for kmers\n"
		"                        (kmers with more matches than this will be ignored)\n"
		"set match_value_min ##  set minimum match value\n"
		"show                    show current cutoffs\n"
		"write ##                write current search results to given file\n";
}

static void search_function(const std::vector<std::string> &list, const KmerLookupInfo &kmers, Selection &selection, const std::vector<RawRead> &reads) {
	(void)reads;
	if (list.size() > 2) {
		std::cout << "Error: search takes one parameter (the sequence to match against)\n";
		return;
	// +1 as opt_mer_length is one less than the set mer length
	} else if (list.size() == 2 && list[1].size() < opt_mer_length + 1) {
		std::cout << "Error: search sequence is too short; need to be at least " << (opt_mer_length + 1) << " basepairs long\n";
		return;
	}
	if (list.size() == 2 && !list[1].empty()) {
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

// continue to concatenate the search selection until a blank line
static void msearch_function(const std::vector<std::string> &list, const KmerLookupInfo &kmers, Selection &selection, const std::vector<RawRead> &reads) {
	(void)kmers;		// avoid warnings
	(void)reads;
	if (list.size() == 1 && list[0].compare("msearch") != 0) {
		selection.search_sequence += list[0];
		// update history - add extra sequence to last search
		// remove last entry (this msearch)
		const int x(where_history());
		HIST_ENTRY *last_entry(remove_history(x));
		free(last_entry->line);
		free(last_entry);
		// we're removing one entry from the history, so we need
		// to move the current position down one
		history_set_pos(x - 1);
		// get (and replace) search entry
		last_entry = remove_history(x - 1);
		const size_t n(strlen(last_entry->line));
		// +1 to include terminating null
		char * const s(static_cast<char*>(malloc(n + list[1].size() + 1)));
		if (!s) {
			throw LocalException("malloc failed");
		}
		memcpy(s, last_entry->line, n);
		memcpy(s + n, list[1].c_str(), list[1].size() + 1);
		free(last_entry->line);
		free(last_entry);
		add_history(s);
		free(s);
	} else if (list.size() == 2 && !list[1].empty()) {
		in_multiline_search = 1;
		selection.search_sequence = list[1];
		// update history - msearch becomes search
		char * const s(static_cast<char*>(malloc(8 + list[1].size())));
		if (!s) {
			throw LocalException("malloc failed");
		}
		memcpy(s, "search ", 7);
		// +1 to include terminating null
		memcpy(s + 7, list[1].c_str(), list[1].size() + 1);
		// remove last entry (this msearch), replace with plain search
		// we don't change the size of the history, so current position is good
		HIST_ENTRY * const last_entry(remove_history(where_history()));
		free(last_entry->line);
		free(last_entry);
		add_history(s);
		free(s);
	} else {
		std::cout << "Error: msearch takes one parameter (the sequence to match against)\n";
		return;
	}
}

static void set_function(const std::vector<std::string> &list, const KmerLookupInfo &kmers, Selection &selection, const std::vector<RawRead> &reads) {
	(void)kmers;
	(void)reads;
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

static void show_function(const std::vector<std::string> &list, const KmerLookupInfo &kmers, Selection &selection, const std::vector<RawRead> &reads) {
	(void)list;
	(void)kmers;
	(void)reads;
	std::cout << "match_value_min " << selection.match_value_min << "\n" <<
		"kmer_hit_max " << selection.kmer_hit_max << "\n\n";
}

static void write_function(const std::vector<std::string> &list, const KmerLookupInfo &kmers, Selection &selection, const std::vector<RawRead> &reads) {
	if (list.size() != 2) {
		std::cout << "Error: write only takes one parameter (the filename to write to)\n";
		return;
	}
	const size_t hits(selection.write_hits(kmers, list[1], reads));
	std::cout << hits << " reads written to " << list[1] << '\n';
}

struct ActionType {
	const char * const action;
	void (*function)(const std::vector<std::string> &, const KmerLookupInfo &, Selection &, const std::vector<RawRead> &);
};

// actions ordered by likelihood of coming up, with most common ones first
static ActionType actions[] = {
	{ "search", search_function },
	{ "msearch", msearch_function },
	{ "set", set_function },
	{ "write", write_function },
	{ "show", show_function },
	{ "?", help_function },
	{ "help", help_function },
	{ "exit", 0 },
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

static void user_input_loop(const KmerLookupInfo &kmers, const std::vector<RawRead> &reads) {
	// take list out of loop to avoid extra allocations/deallocations
	std::vector<std::string> list;
	Selection selection;
	char *s, *t;	// original line, history expanded line
	while ((s = readline("kmers> ")) != 0) {
		list.clear();
		// do history expansion
		const int result(history_expand(s, &t));
		free(s);
		if (result) {
			std::cout << t << '\n';
			if (result < 0 || result == 2) {	// error, or display only
				free(t);
				continue;
			}
		}
		breakup_line(t, list);
		// don't add blank lines to history
		if (list.empty() || (list.size() == 1 && list[1].empty())) {
			free(t);
			// a blank line is the end delimiter for a multiline search
			if (in_multiline_search) {
				in_multiline_search = 0;
				search_function(list, kmers, selection, reads);
			}
			continue;
		}
		// add history for bad commands, too, as they might just
		// be slightly mispelled
		add_history(t);
		// XXX  for multiline search, history should be one search command
		free(t);
		size_t i(0);
		for (; i != actions_size; ++i) {
			if (list[0].compare(actions[i].action) == 0) {
				if (actions[i].function != 0) {
					actions[i].function(list, kmers, selection, reads);
					break;
				}
				return;		// quit action
			}
		}
		if (i == actions_size) {
			if (in_multiline_search) {
				msearch_function(list, kmers, selection, reads);
			} else {
				std::cout << "invalid command: " << list[0] << "\n";
			}
		}
	}
}

int main(int argc, char **argv) {
	int had_error(0);
	if (argc != 2 && argc != 3) {
		print_usage();
		return 1;
	}
	try {
		// set up readline history and command completion
		using_history();
		rl_readline_name = "kmer_matching";
		rl_attempted_completion_function = (rl_completion_func_t *)command_completion;
		KmerLookupInfo kmers;
		read_kmer_index(argv[1], kmers);
		opt_mer_length = kmers.mer_length();
		init_mer_constants();
		std::vector<RawRead> reads;
		if (argc == 3) {
			read_fastq(argv[2], reads, kmers);
		}
		// prevent switching format based on value of output
		std::cout << std::fixed;
		// probably shouldn't need more precision that this
		std::cout.precision(3);
		user_input_loop(kmers, reads);
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
