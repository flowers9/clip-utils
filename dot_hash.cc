#include "hash.h"	// hash
#include "open_compressed.h"	// close_compressed(), open_compressed(), pfread()
#include "version.h"	// VERSION
#include <getopt.h>	// getopt(), optarg, optind
#include <iostream>	// cerr, cout
#include <math.h>	// sqrt()
#include <sstream>	// istringstream
#include <stdlib.h>	// exit()
#include <string.h>	// memcmp()
#include <string>	// string

// read in N saved hashes and create dot-products for all crosses
// XXX - consider a chi-squared type scoring function (sum of squares of the difference)

static int opt_files_to_cutoff;
static int opt_min_kmer_frequency;
static int opt_screen_shared_kmers;
static int opt_shared_identity;

// extend hash by adding a few methods

class dhash : public hash {
    private:
	double magnitude_;
	uint64_t total_count_;
    private:
	value_type value_from_offset(offset_type) const;
    public:
	dhash(void) : hash(), magnitude_(0), total_count_(0) { }
	void init_from_file2(int);		// don't load alts, calculate magnitude
	double dot_product(const dhash &) const;
	double shared_identity(const dhash &) const;
	void init_for_intersection(const dhash &);
	void find_intersection(const dhash &);
	void screen_keys(const dhash &);
	uint64_t total_count(void) const { return total_count_; }
};

// like init_from_file(), but no alt values and calculate magnitude

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
				total_count_ += value_list[i];
				magnitude_ += double(value_list[i]) * double(value_list[i]);
			}
		}
	}
	// read in overflow map
	offset_type x;
	pfread(fd, &x, sizeof(x));
	for (; x <= 0; --x) {
		key_type i;
		value_type j;
		pfread(fd, &i, sizeof(i));
		pfread(fd, &j, sizeof(j));
		value_map[i] = j;
		total_count_ += j;
		// have to subtract out 255**2, as we already added it in above
		magnitude_ += double(j + 255) * double(j + 255) - 65025;
	}
	magnitude_ = sqrt(magnitude_);
}

// relies on INVALID keys having zero value (as set during init_from_from2())
hash::value_type dhash::value_from_offset(const offset_type i) const {
	if (value_list[i] != max_small_value) {
		return value_list[i];
	} else {
		// use find() to avoid inserting a value into value_map
		const std::map<key_type, value_type>::const_iterator a(value_map.find(key_list[i]));
		if (a == value_map.end()) {
			return max_small_value;
		} else {
			return a->second + max_small_value;
		}
	}
}

double dhash::dot_product(const dhash &h) const {
	double d(0);
	for (offset_type i(0); i < modulus; ++i) {
		const value_type j(value_from_offset(i));
		if (j) {
			const value_type k(h.value(key_list[i]));
			if (k) {
				d += double(j) * double(k);
			}
		}
	}
	return d / magnitude_ / h.magnitude_;
}

// count kmers in common
double dhash::shared_identity(const dhash &h) const {
	double d(0);
	for (offset_type i(0); i < modulus; ++i) {
		const value_type j(value_from_offset(i));
		if (j) {
			const value_type k(h.value(key_list[i]));
			if (k) {
				const double x1(double(j) / double(total_count_));
				const double x2(double(k) / double(h.total_count_));
				d += x1 < x2 ? x1 : x2;
			}
		}
	}
	return d;
}

// prepare hash for finding intersections (copy all values as 1)
// (relies on INVALID keys having zero value)

void dhash::init_for_intersection(const dhash &h) {
	used_elements = h.used_elements;
	modulus = h.modulus;
	collision_modulus = h.collision_modulus;
	alt_size = 0;
	key_list = new key_type[modulus];
	value_list = new small_value_type[modulus];
	alt_list = NULL;
	alt_map = NULL;
	for (offset_type i(0); i < modulus; ++i) {
		if (h.value_list[i]) {
			key_list[i] = h.key_list[i];
			value_list[i] = 1;
		} else {
			key_list[i] = INVALID_KEY;
			value_list[i] = 0;
		}
	}
}

// zero any key not shared with h
void dhash::find_intersection(const dhash &h) {
	for (offset_type i(0); i < modulus; ++i) {
		if (value_list[i] && h.value(key_list[i]) == 0) {
			// don't set key to INVALID, as that could hork lookups
			value_list[i] = 0;
		}
	}
}

// zero any key in h, and recalculate magnitude and total_count
void dhash::screen_keys(const dhash &h) {
	total_count_ = 0;
	magnitude_ = 0;
	for (offset_type i(0); i < modulus; ++i) {
		if (value_list[i] == 0) {
		} else if (h.value(key_list[i])) {
			// don't set key to INVALID, as that could hork lookups
			// don't need to zero value_map
			value_list[i] = 0;
		} else {
			const value_type j(value_from_offset(i));
			total_count_ += j;
			magnitude_ += double(j) * double(j);
		}
	}
	magnitude_ = sqrt(magnitude_);
}

static void print_usage() {
	std::cerr << "usage: dot_hash saved_hash1 saved_hash2 ...\n"
		"    -h    print this help\n"
		"    -d    dot-product [shared-identity]\n"
		"    -l ## limit cutoffs to first # files [all]\n"
		"    -m ## min kmer frequency [0]\n"
		"    -s    screen shared kmers\n"
		"    -V    print version\n";
	exit(1);
}

static void get_opts(int argc, char **argv) {
	opt_files_to_cutoff = 0;	// zero == all
	opt_min_kmer_frequency = 0;
	opt_screen_shared_kmers = 0;
	opt_shared_identity = 1;
	int c, x;
	while ((c = getopt(argc, argv, "dhl:m:sV")) != EOF) {
		switch (c) {
		    case 'd':
			opt_shared_identity = 0;
			break;
		    case 'h':
			print_usage();
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
			if (x < 0) {
				std::cerr << "Error: -m requires non-negative value\n";
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
		    case 'V':
			std::cerr << "dot_hash version " << VERSION << "\n";
			exit(0);
		    default:
			std::cerr << "Error: unknown option " << char(c) << "\n";
			print_usage();
		}
	}
}

int main(int argc, char ** argv) {
	get_opts(argc, argv);
	if (argc == optind) {
		print_usage();
	}
	const int hash_count(argc - optind);
	// load in saved hashes
	dhash mer_list[hash_count];
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
		dhash mer_screen;
		mer_screen.init_for_intersection(mer_list[0]);
		for (int i(1); i < hash_count; ++i) {
			mer_screen.find_intersection(mer_list[i]);
		}
		for (int i(0); i < hash_count; ++i) {
			mer_list[i].screen_keys(mer_screen);
		}
	}
	if (opt_shared_identity) {
		// calculate shared identity
		for (int i(1); i < hash_count; ++i) {
			for (int j(0); j < i; ++j) {
				const double x(mer_list[i].shared_identity(mer_list[j]));
				std::cerr << x << " ";
			}
			std::cerr << "\n";
		}
	} else {
		// calculate dot products
		for (int i(1); i < hash_count; ++i) {
			for (int j(0); j < i; ++j) {
				const double x(mer_list[i].dot_product(mer_list[j]));
				std::cerr << x << " ";
			}
			std::cerr << "\n";
		}
	}
}
