#include "hash.h"	// hash
#include "open_compressed.h"	// close_compressed(), open_compressed(), pfread()
#include "version.h"	// VERSION
#include <getopt.h>	// getopt(), optarg, optind
#include <iostream>	// cerr, cout
#include <map>		// map<>
#include <math.h>	// sqrt()
#include <mutex>	// lock_guard<>, mutex
#include <sstream>	// istringstream
#include <stdlib.h>	// exit()
#include <string.h>	// memcmp(), memcpy()
#include <string>	// string
#include <thread>	// thread
#include <vector>	// vector<>

// read in N saved hashes and create shared identity stats for all crosses

static int opt_files_to_cutoff;
static int opt_keep_total_kmer_count;
static int opt_max_kmer_sharing;
static int opt_min_kmer_frequency;
static int opt_screen_shared_kmers;
static int opt_threads;
static std::vector<std::string> opt_reference_list;

// extend hash by adding a few methods

class dhash : public hash {
    public:
	void init_from_file2(int);		// don't load alts, all values 1
	void init_for_intersection(const dhash &);
	void set_intersection(const dhash &);
	void set_subtraction(const dhash &, value_type = 0);
	void set_addition(const dhash &);
	double shared_identity(const dhash &) const;
	double unshared_identity(const dhash &) const;
};

// like init_from_file(), but no alt values, and value is compressed to 1

void dhash::init_from_file2(const int fd) {
	const std::string s(boilerplate());
	char t[s.size()];
	pfread(fd, t, s.size());
	if (memcmp(s.c_str(), t, s.size())) {
		std::cerr << "Error: could not read hash from file: header mismatch\n";
		exit(1);
	}
	pfread(fd, &modulus, sizeof(modulus));
	pfread(fd, &collision_modulus, sizeof(collision_modulus));
	pfread(fd, &used_elements, sizeof(used_elements));
	pfread(fd, &alt_size, sizeof(alt_size));
	alt_size = 0;
	key_list = new key_type[modulus];
	value_list = new small_value_type[modulus];
	alt_list = NULL;
	alt_map = NULL;
	// read in values (they're the smallest size)
	pfread(fd, &value_list[0], sizeof(small_value_type) * modulus);
	// read in keys for non-zero values
	for (offset_type i(0); i < modulus; ++i) {
		if (value_list[i] == 0) {
			key_list[i] = INVALID_KEY;
		} else {
			pfread(fd, &key_list[i], sizeof(key_type));
			if (value_list[i] < opt_min_kmer_frequency) {
				value_list[i] = 0;
			} else {
				value_list[i] = 1;
			}
		}
	}
}

// prepare hash for finding intersections
void dhash::init_for_intersection(const dhash &h) {
	used_elements = h.used_elements;
	modulus = h.modulus;
	collision_modulus = h.collision_modulus;
	alt_size = 0;
	key_list = new key_type[modulus];
	value_list = new small_value_type[modulus];
	alt_list = NULL;
	alt_map = NULL;
	memcpy(key_list, h.key_list, modulus * sizeof(key_type));
	memcpy(value_list, h.value_list, modulus * sizeof(small_value_type));
}

// zero any key not in h
void dhash::set_intersection(const dhash &h) {
	for (offset_type i(0); i < modulus; ++i) {
		if (value_list[i] && h.value(key_list[i]) == 0) {
			// don't set key to INVALID, as that could hork lookups
			value_list[i] = 0;
			// --used_elements;	// technically not true, but effectively true
		}
	}
}

// zero any key where h's value is above cutoff
void dhash::set_subtraction(const dhash &h, const hash::value_type max_value) {
	for (offset_type i(0); i < modulus; ++i) {
		if (value_list[i] && h.value(key_list[i]) > max_value) {
			// don't set key to INVALID, as that could hork lookups
			value_list[i] = 0;
			if (!opt_keep_total_kmer_count) {
				--used_elements;	// technically not true, but effectively true
			}
		}
	}
}

// increment all values from h
void dhash::set_addition(const dhash &h) {
	for (offset_type i(0); i < h.modulus; ++i) {
		if (h.value_list[i]) {
			if (!increment(h.key_list[i])) {
				std::cerr << "Error: ran out of space in hash\n";
				exit(1);
			}
		}
	}
}

// count kmers in common
double dhash::shared_identity(const dhash &h) const {
	uint64_t d(0);
	for (offset_type i(0); i < modulus; ++i) {
		if (value_list[i] && h.value(key_list[i])) {
			++d;
		}
	}
	return d;
}

// count kmers shared with only one reference (h is the shared_kmers set)
double dhash::unshared_identity(const dhash &h) const {
	uint64_t d(0);
	for (offset_type i(0); i < modulus; ++i) {
		if (value_list[i] && h.value(key_list[i]) == 1) {
			++d;
		}
	}
	return d;
}

static void print_usage() {
	std::cerr << "usage: dot_hash saved_hash1 saved_hash2 ...\n"
		"    -h    print this help\n"
		"    -k    when calculating fraction, compare to total unique kmers\n"
		"    -l ## limit min frequency to first ## files [all]\n"
		"    -m ## min kmer frequency (only applies to non-references) [0]\n"
		"    -r    add reference file (may be specified multiple times)\n"
		"    -s    screen shared kmers\n"
		"    -t ## threads [1]\n"
		"    -u ## only count kmers shared with at most ## references [all]\n"
		"    -V    print version\n";
	exit(1);
}

static void get_opts(const int argc, char * const * argv) {
	opt_files_to_cutoff = 0;	// zero == all
	opt_keep_total_kmer_count = 0;
	opt_max_kmer_sharing = hash::max_small_value;
	opt_min_kmer_frequency = 0;
	opt_screen_shared_kmers = 0;
	opt_threads = 1;
	int c, x;
	while ((c = getopt(argc, argv, "hkl:m:st:u:V")) != EOF) {
		switch (c) {
		    case 'h':
			print_usage();
			break;
		    case 'k':
			opt_keep_total_kmer_count = 1;
			break;
		    case 'l':
			std::istringstream(optarg) >> x;
			if (x < 0) {
				std::cerr << "Error: -l requires non-negative value\n";
				exit(1);
			}
			opt_files_to_cutoff = x;
			break;
		    case 'm':
			std::istringstream(optarg) >> x;
			if (x < 1) {
				std::cerr << "Error: -m requires positive value\n";
				exit(1);
			} else if (x > hash::max_small_value) {
				std::cerr << "Error: -m value too large: " << x << " (max " << (unsigned int)(hash::max_small_value) << ")\n";
				exit(1);
			}
			opt_min_kmer_frequency = x;
			break;
		    case 's':
			opt_screen_shared_kmers = 1;
			break;
		    case 't':
			std::istringstream(optarg) >> x;
			if (x < 0) {
				std::cerr << "Error: -t requires non-negative value\n";
				exit(1);
			}
			opt_threads = x;
			break;
		    case 'u':
			std::istringstream(optarg) >> x;
			if (x < 1) {
				std::cerr << "Error: -u requires positive value\n";
				exit(1);
			}
			opt_max_kmer_sharing = x;
			break;
		    case 'V':
			std::cerr << "dot_hash version " << VERSION << "\n";
			exit(0);
		    default:
			std::cerr << "Error: unknown option " << char(c) << "\n";
			print_usage();
		}
	}
}

// declared here so threads can easily use them
static int hash_count;
static std::mutex get_next_mutex;
static std::vector<dhash> mer_list;
static std::vector<std::vector<double> > results;
static dhash mer_screen, shared_kmers;

// return next i; returns 0 if done; make sure to initialize i_out to a
// value other than -1 before starting loop

static int get_next_i(int &i_out) {
	std::lock_guard<std::mutex> lock(get_next_mutex);
	static int i(-1);
	if (i_out == -1) {		// reset count
		i = -1;
		return 0;
	} else if (++i < hash_count) {
		i_out = i;
		return 1;
	} else {
		i = hash_count;		// prevent wrapping
		return 0;
	}
}

static int skip_upper_matrix(1);
static int skip_diagonal(0);

// return next [i][j] pairing; returns 0 if done
static int get_next_pair(int &i_out, int &j_out) {
	std::lock_guard<std::mutex> lock(get_next_mutex);
	static int i(0), j(-1);
	if (++j == i) {
		if (skip_diagonal) {
			++j;
		} else {
			i_out = i;
			j_out = j;
			return 1;
		}
	}
	if (j < (skip_upper_matrix ? i : hash_count)) {
		i_out = i;
		j_out = j;
		return 1;
	} else if (++i < hash_count) {
		i_out = i;
		j_out = j = 0;
		return 1;
	} else {
		i = j = hash_count;	// prevent wrapping
		return 0;
	}
}

static void screen_universal_key(void) {
	int i(0);
	while (get_next_i(i)) {
		mer_list[i].set_subtraction(mer_screen);
	}
}

static void screen_universal_keys(void) {
	std::thread threads[opt_threads];
	int x(-1);
	get_next_i(x);		// reset count
	for (int i(0); i < opt_threads; ++i) {
		threads[i] = std::thread(screen_universal_key);
	}
	for (int i(0); i < opt_threads; ++i) {
		threads[i].join();
	}
}

static void screen_shared_key(void) {
	int i(0);
	while (get_next_i(i)) {
		mer_list[i].set_subtraction(shared_kmers, opt_max_kmer_sharing);
	}
}

static void screen_shared_keys(void) {
	std::thread threads[opt_threads];
	int x(-1);
	get_next_i(x);		// reset count
	for (int i(0); i < opt_threads; ++i) {
		threads[i] = std::thread(screen_shared_key);
	}
	for (int i(0); i < opt_threads; ++i) {
		threads[i].join();
	}
}

static void calculate_shared_identity(void) {
	int i, j;
	while (get_next_pair(i, j)) {
		if (i == j) {
			results[i][j] = mer_list[i].unshared_identity(shared_kmers) / mer_list[i].size();
		} else {
			const double x(mer_list[i].shared_identity(mer_list[j]));
			results[i][j] = x / mer_list[i].size();
			results[j][i] = x / mer_list[j].size();
		}
	}
}

static void calculate_shared_identities(void) {
	std::thread threads[opt_threads];
	for (int i(0); i < opt_threads; ++i) {
		threads[i] = std::thread(calculate_shared_identity);
	}
	for (int i(0); i < opt_threads; ++i) {
		threads[i].join();
	}
}

static void print_results(char * const * argv) {
	for (int i(0); i < hash_count; ++i) {
		for (int j(0); j < hash_count; ++j) {
			double x;
			if (i == j && skip_diagonal) {
				x = 1;
			} else if (j > i && skip_upper_matrix) {
				x = results[j][i];
			} else {
				x = results[i][j];
			}
			std::cout << x << " ";
		}
		std::cout << argv[i + optind] << "\n";
	}
}

static void print_mer_histogram(hash &h) {
	std::map<hash::value_type, unsigned long> counts;
	hash::const_iterator a(h.begin());
	const hash::const_iterator end_a(h.end());
	for (; a != end_a; ++a) {
		++counts[a.value];
	}
	std::map<hash::value_type, unsigned long>::const_iterator c(counts.begin());
	const std::map<hash::value_type, unsigned long>::const_iterator end_c(counts.end());
	for (; c != end_c; ++c) {
		std::cout << c->first << " " << c->second << "\n";
	}
}

int main(const int argc, char * const * argv) {
	get_opts(argc, argv);
	if (argc == optind) {
		print_usage();
	}
	hash_count = argc - optind;
	// load in saved hashes
	mer_list.assign(hash_count, dhash());
	for (int i(0); i < hash_count; ++i) {
		if (opt_files_to_cutoff && i == opt_files_to_cutoff) {
			opt_min_kmer_frequency = 0;
		}
		const int fd(open_compressed(argv[i + optind]));
		if (fd == -1) {
			std::cerr << "Error: could not read saved hash: " << argv[i + optind] << "\n";
			return 1;
		}
		mer_list[i].init_from_file2(fd);
		close_compressed(fd);
	}
	if (opt_screen_shared_kmers) {
		mer_screen.init_for_intersection(mer_list[0]);
		// possibly make this threaded
		for (int i(1); i < hash_count; ++i) {
			mer_screen.set_intersection(mer_list[i]);
		}
		screen_universal_keys();
		// XXX - could delete mer_screen here
	}
	if (opt_max_kmer_sharing) {
		uint64_t x(0);
		for (int i(0); i < hash_count; ++i) {
			x += mer_list[i].size();
		}
		shared_kmers.init(x);
		for (int i(0); i < hash_count; ++i) {
			shared_kmers.set_addition(mer_list[i]);
		}
		screen_shared_keys();
		// XXX - could delete shared_kmers here
	}
	results.assign(hash_count, std::vector<double>(hash_count, 0));
	calculate_shared_identities();
	skip_upper_matrix = 0;
	print_results(argv);
}
