#include "hashn.h"	// hashn
#include "next_prime.h"	// next_prime()
#include "open_compressed.h"	// close_compressed(), open_compressed(), pfread()
#include "version.h"	// VERSION
#include <getopt.h>	// getopt(), optarg, optind
#include <iomanip>	// ios, setprecision()
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

static int opt_keep_total_kmer_count;
static int opt_max_kmer_sharing;
static int opt_max_kmer_frequency;
static int opt_min_kmer_frequency;
static int opt_threads;
static std::vector<std::string> opt_reference_list;

// extend hash by adding a few methods

class dhash : public hashn {
    public:
	void init_from_file2(int, int);		// don't load alts, all values 1
	void set_subtraction(const dhash &, value_type = 0);
	void set_addition(const dhash &);
	double shared_identity(const dhash &) const;
	value_type value(const key_type_base &) const;
	bool increment(const key_type_base &key);
	void print_kmer_matching(const dhash *, size_t) const;
	void print_key(const dhash::key_type_base &key) const;
};

// like increment(key_type), but doesn't create overflow values
bool dhash::increment(const key_type_base &key) {
	const offset_type i(insert_offset(key));
	if (i == modulus) {	// insert failed
		return 0;
	}
	if (value_list[i] != max_small_value) {
		++value_list[i];
	}
	return 1;
}

// like value(key_type), but doesn't check value_map for overflow values
dhash::value_type dhash::value(const key_type_base &key) const {
	const offset_type i(find_offset(key));
	if (i == modulus) {	// key not found
		return 0;
	}
	return value_list[i];
}

// like init_from_file(), but no alt values, and value is always 1 if key is present
// if file_type is 0, it's a reference file, 1 is a fastq file
void dhash::init_from_file2(const int fd, const int file_type) {
	const std::string s(boilerplate());
	char t[s.size()];
	pfread(fd, t, s.size());
	if (memcmp(t, s.c_str(), s.size()) != 0) {
		std::cerr << "Error: could not read hash from file: header mismatch\n" << s << "\n" << t << "\n";
		exit(1);
	}
	pfread(fd, &modulus, sizeof(modulus));
	pfread(fd, &collision_modulus, sizeof(collision_modulus));
	pfread(fd, &used_elements, sizeof(used_elements));
	pfread(fd, &alt_size, sizeof(alt_size));
	pfread(fd, &bit_width, sizeof(bit_width));
	word_width = (bit_width + 8 * sizeof(base_type) - 1) / (8 * sizeof(base_type));
	const offset_type n((modulus + 1) * word_width);
	alt_size = 0;
	alt_list = NULL;
	alt_map = NULL;
	key_list = new base_type[n];
	value_list = new small_value_type[modulus];
	pfread(fd, &value_list[0], sizeof(small_value_type) * modulus);
	base_type *a(key_list);
	for (offset_type i(0); i != modulus; ++i) {
		if (value_list[i] == 0) {
			base_type * const end_a(a + word_width);
			for (; a != end_a; ++a) {
				*a = INVALID_KEY;
			}
		} else {
			pfread(fd, a, sizeof(base_type) * word_width);
			a += word_width;
			if ((file_type == 0 && opt_max_kmer_frequency < value_list[i]) || (file_type == 1 && value_list[i] < opt_min_kmer_frequency)) {
				// don't change key to INVALID_KEY, or we'd have to rehash
				value_list[i] = 0;
			}
		}
	}
	invalid_key.assign(*this, modulus);
	pfread(fd, a, sizeof(base_type) * word_width);		// invalid_key
}

#if 0
// like init_from_file2(), but resize hash to be 50% full
void dhash::rehash_from_file(const int fd) {
	const std::string s(boilerplate());
	char t[s.size()];
	pfread(fd, t, s.size());
	if (memcmp(s.c_str(), t, s.size())) {
		std::cerr << "Error: could not read hash from file: header mismatch\n";
		exit(1);
	}
	offset_type original_modulus;
	pfread(fd, &original_modulus, sizeof(original_modulus));
	pfread(fd, &collision_modulus, sizeof(collision_modulus));
	pfread(fd, &used_elements, sizeof(used_elements));
	pfread(fd, &alt_size, sizeof(alt_size));
	pfread(fd, &bit_width, sizeof(bit_width));
	word_width = (bit_width + 8 * sizeof(base_type) - 1) / (8 * sizeof(base_type));
	const offset_type original_used_elements(used_elements);
	alt_size = 0;
	alt_list = NULL;
	alt_map = NULL;
	small_value_type *old_value_list;
	used_elements = original_used_elements;
	if (opt_min_kmer_frequency) {	// screen out keys with low values
		old_value_list = new small_value_type[modulus];
		// read in values (they're the smallest size)
		pfread(fd, &old_value_list[0], sizeof(small_value_type) * modulus);
		for (offset_type i(0); i < modulus; ++i) {
			// zero values weren't counted as used elements already
			if (old_value_list[i] && old_value_list[i] < opt_min_kmer_frequency) {
				--used_elements;
			}
		}
	} else {
		skip_next_chars(fd, sizeof(small_value_type) * modulus);
	}
	size_t size_asked(2 * used_elements);
	used_elements = 1;	// to account for minimum of one INVALID_KEYs
	if (size_asked < 3) {	// to avoid collision_modulus == modulus
		size_asked = 3;
	}
	modulus = next_prime(size_asked);
	// collision_modulus just needs to be relatively prime with modulus;
	// since modulus is prime, any value will do - I made it prime for fun
	collision_modulus = next_prime(size_asked / 2);
	const offset_type n((modulus + 1) * word_width);
	key_list = new base_type[n];
	value_list = new small_value_type[modulus];
	// initialize keys; values are initialized as keys are entered
	for (offset_type i(0); i != n; ++i) {
		key_list[i] = INVALID_KEY;
	}
	invalid_key.assign(*this, modulus);
	// subsequent calculations want a zero value for invalid entries
	memset(value_list, 0, sizeof(small_value_type) * modulus);
	// read in keys
	base_type *xbuf(new base_type[word_width]);
	key_type_base x(word_width, xbuf);
	if (opt_min_kmer_frequency) {	// screen out keys with low values
		for (offset_type i(0); i < modulus; ++i) {
			if (old_value_list[i]) {
				pfread(fd, xbuf, sizeof(base_type) * word_width);
				if (old_value_list[i] >= opt_min_kmer_frequency) {
					// don't have to worry about hash filling up
					value_list[insert_offset(x)] = 1;
				}
			}
		}
		delete[] old_value_list;
	} else {
		for (offset_type i(0); i < original_used_elements - 1; ++i) {
			pfread(fd, xbuf, sizeof(base_type) * word_width);
			// don't have to worry about hash filling up
			value_list[insert_offset(x)] = 1;
		}
	}
	delete[] xbuf;
}
#endif

// zero any key where h's value is above cutoff
void dhash::set_subtraction(const dhash &h, const value_type max_value) {
	base_type *j(key_list);
	for (offset_type i(0); i < modulus; ++i, j += word_width) {
		if (value_list[i] && h.value(key_type_internal(*this, j)) > max_value) {
			// don't set key to INVALID_KEY, as that could hork lookups
			value_list[i] = 0;
			if (!opt_keep_total_kmer_count) {
				--used_elements;	// only count non-zero keys
			}
		}
	}
}

// increment all values from h
void dhash::set_addition(const dhash &h) {
	base_type *j(h.key_list);
	for (offset_type i(0); i < h.modulus; ++i, j += h.word_width) {
		if (h.value_list[i] && !increment(key_type_internal(h, j))) {
			std::cerr << "Error: ran out of space in hash - recompile with larger hash size\n";
			exit(1);
		}
	}
}

// count kmers in common
double dhash::shared_identity(const dhash &h) const {
	base_type *j(key_list);
	uint64_t x(0);
	for (offset_type i(0); i < modulus; ++i, j += word_width) {
		if (value_list[i] && h.value(key_type_internal(*this, j))) {
			++x;
		}
	}
	return x;
}

void dhash::print_key(const dhash::key_type_base &key) const {
	const char values[4] = { 'A', 'C', 'G', 'T' };
	for (unsigned long i(bit_width - 2); i != 0; i -= 2) {
		std::cout << values[key.basepair(i)];
	}
	std::cout << values[key.basepair(0)];
}

void dhash::print_kmer_matching(const dhash * const h, const size_t n) const {
	base_type *j(key_list);
	for (offset_type i(0); i < modulus; ++i, j += word_width) {
		if (value_list[i]) {
			const key_type_internal key(*this, j);
			print_key(key);
			std::cout << ' ';
			for (size_t k(0); k < n; ++k) {
				std::cout << h[k].value(key);
			}
			std::cout << "\n";
		}
	}
}

static void print_usage() {
	std::cerr << "usage: dot_hash saved_hash1 saved_hash2 ...\n"
		"    -h    print this help\n"
		"    -k    when calculating fraction, compare to total unique kmers\n"
		"    -m ## min kmer frequency (only applies to non-references)\n"
		"    -M ## max kmer frequency (only applies to references) [1]\n"
		"    -r ## add reference file (may be specified multiple times)\n"
		"    -t ## threads [1]\n"
		"    -u ## only count kmers shared with at most ## references\n"
		"          (negative values mean shared by all but ##)\n"
		"    -V    print version\n";
	exit(1);
}

static void get_opts(const int argc, char * const * argv) {
	opt_keep_total_kmer_count = 0;
	opt_max_kmer_sharing = 0;
	opt_min_kmer_frequency = 0;
	opt_max_kmer_frequency = 1;
	opt_threads = 1;
	int c, x;
	while ((c = getopt(argc, argv, "hkm:M:r:t:u:V")) != EOF) {
		switch (c) {
		    case 'h':
			print_usage();
			break;
		    case 'k':
			opt_keep_total_kmer_count = 1;
			break;
		    case 'm':
			std::istringstream(optarg) >> x;
			if (x < 1) {
				std::cerr << "Error: -m requires positive value\n";
				exit(1);
			} else if (x > dhash::max_small_value) {
				std::cerr << "Error: -m value too large: " << x << " (max " << (unsigned int)(dhash::max_small_value) << ")\n";
				exit(1);
			}
			opt_min_kmer_frequency = x;
			break;
		    case 'M':
			std::istringstream(optarg) >> x;
			if (x < 1) {
				std::cerr << "Error: -m requires positive value\n";
				exit(1);
			} else if (x > dhash::max_small_value) {
				std::cerr << "Error: -m value too large: " << x << " (max " << (unsigned int)(dhash::max_small_value) << ")\n";
				exit(1);
			}
			opt_max_kmer_frequency = x;
			break;
		    case 'r':
			opt_reference_list.push_back(optarg);
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
	if (opt_reference_list.empty() && optind == argc) {
		std::cerr << "Error: no files given\n";
		print_usage();
	} else if (opt_reference_list.size() + argc - optind < 2) {
		std::cerr << "Error: only one file specified\n";
		exit(1);
	}
	if (opt_max_kmer_sharing < 0 && opt_reference_list.size() > 1) {
		if (static_cast<size_t>(-opt_max_kmer_sharing) >= opt_reference_list.size()) {
			std::cerr << "Error: -u out of range\n";
			exit(1);
		}
		opt_max_kmer_sharing = opt_reference_list.size() + opt_max_kmer_sharing;
	}
#if 0
	// for printing shared frequencies
	if (argc - optind != 1 || opt_reference_list.empty()) {
		std::cerr << "Error: only one fastq dump may be specified, and at least one reference must be\n";
		exit(1);
	}
#endif
}

class counter_1d {
    private:
	std::mutex mutex_;
	int i_, end_i_, i_offset_;
    public:
	int get_next(int &i_out) {
		std::lock_guard<std::mutex> lock(mutex_);
		if (++i_ < end_i_) {
			i_out = i_;
			return 1;
		} else {
			i_ = end_i_;	// prevent wrapping
			return 0;
		}
	}
	void set(const int end_x, const int x_offset) {
		i_ = -1;
		end_i_ = end_x;
		i_offset_ = x_offset;
	}
	int i_offset(void) const {
		return i_offset_;
	}
};

class counter_2d {
    private:
	std::mutex mutex_;
	int i_, j_, end_i_, end_j_;
	int i_offset_, j_offset_;
	int skip_upper_half_;		// includes skipping diagonal
    public:
	counter_2d(void) : skip_upper_half_(0) { }
	int get_next(int &i_out, int &j_out) {
		std::lock_guard<std::mutex> lock(mutex_);
		if (++j_ < (skip_upper_half_ ? i_ : end_j_)) {
			i_out = i_;
			j_out = j_;
			return 1;
		} else if (++i_ < end_i_) {
			i_out = i_;
			j_out = j_ = 0;
			return 1;
		} else {
			i_ = end_i_;	// prevent wrapping
			j_ = end_j_;
			return 0;
		}
	}
	void set(const int end_x, const int end_y, const int x_offset, const int y_offset) {
		i_ = 0;
		j_ = -1;
		end_i_ = end_x;
		end_j_ = end_y;
		i_offset_ = x_offset;
		j_offset_ = y_offset;
	}
	void set_skip_upper_half(void) {
		skip_upper_half_ = 1;
	}
	int i_offset(void) const {
		return i_offset_;
	}
	int j_offset(void) const {
		return j_offset_;
	}
	int skip_upper_half(void) const {
		return skip_upper_half_;
	}
};

// declared here so threads can easily use them
static std::vector<dhash> mer_list;
static std::vector<std::vector<double> > results;
static dhash *shared_kmers;
static counter_1d i_counter;
static counter_2d pair_counter;

static void screen_shared_key(void) {
	int i(0);
	while (i_counter.get_next(i)) {
		mer_list[i + i_counter.i_offset()].set_subtraction(*shared_kmers, opt_max_kmer_sharing);
	}
}

static void screen_shared_keys(void) {
	std::thread threads[opt_threads];
	for (int i(0); i < opt_threads; ++i) {
		threads[i] = std::thread(screen_shared_key);
	}
	for (int i(0); i < opt_threads; ++i) {
		threads[i].join();
	}
}

static void calculate_shared_identity(void) {
	int i, j;
	while (pair_counter.get_next(i, j)) {
		const double x(mer_list[i + pair_counter.i_offset()].shared_identity(mer_list[j + pair_counter.j_offset()]));
		results[i][j] = x / mer_list[i + pair_counter.i_offset()].size();
		if (pair_counter.skip_upper_half()) {
			results[j][i] = x / mer_list[j + pair_counter.j_offset()].size();
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

// style: 0 - references only, 1 - fastqs only, 2 - fastqs vs references
static void print_results(const int end_x, const int end_y, const int style, char * const * argv) {
	// scale output so numbers aren't just a bunch of zeroes
	double biggest(0);
	for (int i(0); i < end_x; ++i) {
		for (int j(0); j < end_y; ++j) {
			if (i != j || style == 2) {
				if (biggest < results[i][j]) {
					biggest = results[i][j];
				}
			}
		}
	}
	if (biggest == 0) {
		std::cerr << "Warning: no result is greater than zero\n";
		return;
	}
	double scale(1);
	for (; biggest * scale * 10 < 1; scale *= 10) { }
	if (scale > 1) {
		std::cout << "Results multiplied by " << scale << " for ease of display\n\n";
	}
	std::cout.setf(std::ios::fixed);
	for (int i(0); i < end_x; ++i) {
		for (int j(0); j < end_y; ++j) {
			if (i == j && style != 2) {
				std::cout << " ---  ";
			} else {
				std::cout << std::setprecision(3) << results[i][j] * scale << " ";
			}
		}
		if (style == 0) {
			std::cout << opt_reference_list[i] << "\n";
		} else {
			std::cout << argv[i + optind] << "\n";
		}
	}
	if (style == 2) {
		std::cout << "\nReferences:";
		std::vector<std::string>::const_iterator a(opt_reference_list.begin());
		const std::vector<std::string>::const_iterator end_a(opt_reference_list.end());
		for (; a != end_a; ++a) {
			std::cout << " " << *a;
		}
		std::cout << "\n";
	}
}

static void find_dot_values(const int fastq_count, char * const * argv) {
	// this is the x versus y shared amounts part
	int x_size, y_size, y_offset, style;
	if (fastq_count == 0) {				// reference vs reference matrix
		x_size = y_size = opt_reference_list.size();
		y_offset = 0;
		style = 0;
		pair_counter.set_skip_upper_half();
	} else if (opt_reference_list.empty()) {	// fastq vs fastq matrix
		x_size = y_size = fastq_count;
		y_offset = 0;
		style = 1;
		pair_counter.set_skip_upper_half();
	} else {					// fastq vs reference
		x_size = fastq_count;
		y_size = opt_reference_list.size();
		y_offset = fastq_count;
		style = 2;
	}
	pair_counter.set(x_size, y_size, 0, y_offset);
	results.assign(x_size, std::vector<double>(y_size, 0));
	calculate_shared_identities();
	print_results(x_size, y_size, style, argv);
}

int main(const int argc, char * const * argv) {
	get_opts(argc, argv);
	const int fastq_count(argc - optind);
	mer_list.assign(fastq_count + opt_reference_list.size(), dhash());
	// load in fastq saved hashes
	int i(0);
	for (; i < fastq_count; ++i) {
		std::cerr << "reading " << argv[i + optind] << "\n";
		const int fd(open_compressed(argv[i + optind]));
		if (fd == -1) {
			std::cerr << "Error: could not read saved hash: " << argv[i + optind] << "\n";
			return 1;
		}
		mer_list[i].init_from_file2(fd, 1);
		std::cerr << "size " << mer_list[i].size() << "\n";
		close_compressed(fd);
	}
	// load in reference saved hashes
	std::vector<std::string>::const_iterator a(opt_reference_list.begin());
	const std::vector<std::string>::const_iterator end_a(opt_reference_list.end());
	for (; a != end_a; ++a, ++i) {
		std::cerr << "reading " << *a << "\n";
		const int fd(open_compressed(*a));
		if (fd == -1) {
			std::cerr << "Error: could not read saved hash: " << *a << "\n";
			return 1;
		}
		mer_list[i].init_from_file2(fd, 0);
		std::cerr << "size " << mer_list[i].size() << "\n";
		close_compressed(fd);
	}
	if (opt_reference_list.size() > 1 && opt_max_kmer_sharing != 0) { // && static_cast<size_t>(opt_max_kmer_sharing) < opt_reference_list.size()) {
		uint64_t x(0);
		for (size_t j(0); j < opt_reference_list.size(); ++j) {
			x += mer_list[j + fastq_count].size();
		}
		shared_kmers = new dhash;
		shared_kmers->init(x, mer_list[0].bits());
		std::cerr << "screening shared kmers\n";
		for (size_t j(0); j < opt_reference_list.size(); ++j) {
			shared_kmers->set_addition(mer_list[j + fastq_count]);
		}
		std::cerr << "size " << shared_kmers->size() << "\n";
		i_counter.set(opt_reference_list.size(), fastq_count);
		screen_shared_keys();
		delete shared_kmers;
	}
	std::cerr << "processing kmers\n";
	// find_dot_values(fastq_count, argv);
	// print shared frequencies
	mer_list[0].print_kmer_matching(&mer_list[1], opt_reference_list.size());
}
