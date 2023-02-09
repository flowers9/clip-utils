#include "hashl.h"	// hashl
#include "hashl_metadata.h"	// hashl_metadata
#include "open_compressed.h"	// close_compressed(), get_suffix(), open_compressed(), pfgets()
#include "time_used.h"	// elapsed_time(), start_time()
#include "version.h"	// VERSION
#include "write_fork.h"	// close_fork(), write_fork()
#include <fstream>	// ofstream
#include <getopt.h>	// getopt(), optarg, optind
#include <iomanip>	// fixed, setprecision()
#include <iostream>	// cerr, cout, ostream
#include <list>		// list<>
#include <map>		// map<>
#include <set>		// set<>
#include <sstream>	// istringstream
#include <stdint.h>	// uint64_t
#include <stdlib.h>	// exit()
#include <string.h>	// memcpy()
#include <string>	// string
#include <utility>	// make_pair(), pair<>
#include <vector>	// vector<>

static bool opt_feedback;
static bool opt_print_gc;
static double opt_load_lower_bound;
static double opt_load_upper_bound;
static int opt_histogram_restore;
static size_t opt_frequency_cutoff;
static size_t opt_mer_length;
static size_t opt_nmers;
static std::string opt_save_file;

static void save_memory(const hashl &mer_list) {
	std::string suffix;
	get_suffix(opt_save_file, suffix);
	std::list<std::string> args;
	if (suffix == ".gz") {
		args.push_back("gzip");
		args.push_back("-c");
	} else if (suffix == ".bz2") {
		args.push_back("bzip2");
		args.push_back("-c");
	} else if (suffix == ".xz") {
		args.push_back("xz");
		args.push_back("-c");
	} else if (suffix == ".Z") {
		args.push_back("compress");
		args.push_back("-c");
	}
	const int fd(write_fork(args, opt_save_file));
	if (fd == -1) {
		std::cerr << "Error: could not save memory\n";
		exit(1);
	}
	mer_list.save(fd);
	close_fork(fd);
}

// convert key to sequence

static std::string convert_key(const hashl::key_type &key) {
	const char values[4] = { 'A', 'C', 'G', 'T' };
	std::string sequence;
	sequence.reserve(opt_mer_length);
	// relies on wrap-around for termination
	for (size_t i(opt_mer_length * 2 - 2); i < opt_mer_length * 2; i -= 2) {
		sequence += values[key.basepair(i)];
	}
	return sequence;
}

static void reverse_key(const hashl::key_type &key_in, hashl::key_type &key_out) {
	for (size_t i(0); i < opt_mer_length * 2; i += 2) {
		key_out.push_back(3 - key_in.basepair(i));
	}
}

void print_final_input_feedback(const hashl &mer_list) {
	if (opt_feedback && mer_list.size() != 0) {
		std::cerr << time(0) << ": " << mer_list.size() << " entries used (" << double(100) * mer_list.size() / mer_list.capacity() << ")\n";
	}
}

// print n-mer occurence frequency

static void print_mer_frequency(std::ostream &fp_out, const hashl &mer_list) {
	hashl::key_type key(mer_list), comp_key(mer_list);
	hashl::const_iterator a(mer_list.cbegin());
	const hashl::const_iterator end_a(mer_list.cend());
	for (; a != end_a; ++a) {
		if (a.value() >= opt_frequency_cutoff) {
			a.get_key(key);
			fp_out << convert_key(key) << ' ' << a.value() << "\n";
			reverse_key(key, comp_key);
			if (key != comp_key) {
				fp_out << convert_key(comp_key) << ' ' << a.value() << "\n";
			}
		}
	}
}

// number of gc bases in string

static unsigned long count_gc(const hashl::key_type &key) {
	const std::string s(convert_key(key));
	unsigned long n(0);
	std::string::const_iterator a(s.begin());
	const std::string::const_iterator end_a(s.end());
	for (; a != end_a; ++a) {
		if (*a == 'G' || *a == 'g' || *a == 'C' || *a == 'c') {
			++n;
		}
	}
	return n;
}

// print histogram of n-mer occurrences, potentially with gc percent at given
// frequencies

static void print_mer_histogram(std::ostream &fp_out, const hashl &mer_list) {
	std::map<hashl::value_type, unsigned long> counts;
	std::map<hashl::value_type, unsigned long> gc_counts;
	hashl::key_type key(mer_list), comp_key(mer_list);
	hashl::const_iterator a(mer_list.cbegin());
	const hashl::const_iterator end_a(mer_list.cend());
	for (; a != end_a; ++a) {
		++counts[a.value()];
		if (opt_print_gc) {
			a.get_key(key);
			gc_counts[a.value()] += count_gc(key);
		}
	}
	std::map<hashl::value_type, unsigned long>::const_iterator c(counts.begin());
	const std::map<hashl::value_type, unsigned long>::const_iterator end_c(counts.end());
	// don't include single occurrences in total
	// (hashes don't have non-positive occurrence values)
	if (c != end_c && c->first == 1) {
		++c;
	}
	double total(0);
	for (; c != end_c; ++c) {
		total += static_cast<double>(c->first) * static_cast<double>(c->second);
	}
	c = counts.begin();
	if (c != end_c && c->first == 1) {
		fp_out <<  c->first << ' ' << c->second << "\n";
		++c;
	}
	double i(0);
	for (; c != end_c; ++c) {
		const double x(100 * static_cast<double>(c->first) * static_cast<double>(c->second));
		i += x;
		if (opt_print_gc) {
			fp_out << c->first << ' ' << c->second << ' ' << x / total << ' ' << i / total << ' ' << 100 * static_cast<double>(gc_counts[c->first]) / static_cast<double>(c->second) / static_cast<double>(opt_mer_length) << "\n";
		} else {
			fp_out << c->first << ' ' << c->second << ' ' << x / total << ' ' << i / total << "\n";
		}
	}
}

// return the number represented by s, which may be suffixed by a k, m, or g
// which act as multipliers to the base amount

static size_t get_value(const std::string s) {
	// validate string - digits optionally followed by k, m, or g
	const size_t i(s.find_first_not_of("0123456789"));
	if (i == std::string::npos) {
		size_t x;
		std::istringstream(s) >> x;
		return x;
	} else if (i + 1 == s.length()) {
		size_t x;
		std::istringstream(s) >> x;
		switch (s[i]) {
		    case 'g':
			return x << 30;
		    case 'm':
			return x << 20;
		    case 'k':
			return x << 10;
		    default:
			return 0;
		}
	} else {	// bad value
		return 0;
	}
}

// print usage statement and exit

static void print_usage() {
	std::cerr <<
		"usage: histogram [options] file1 [file2] ...\n"
		"    -g    print percent gc content at each frequency\n"
		"    -h    print this information\n"
		"    -i    turn off status updates\n"
		"    -l ## lower bound for hash fill fraction\n"
		"    -L ## upper bound for hash fill fraction\n"
		"    -m ## set mer length [24]\n"
		"    -o ## print output to file instead of stdout\n"
		"    -s ## save histogram memory structure to file\n"
		"    -S ## load histogram memory dump from given file\n"
		"    -V    print version\n"
		"    -w ## print frequency count instead of histogram, for all n-mers with\n"
		"          a frequency of at least ## [0 (off)]\n"
		"    -z ## number of possible n-mers to allocate memory for (overrides -l/-L)\n"
		"          (k, m, or g may be suffixed)\n";
	exit(1);
}

static std::ostream &get_opts(int argc, char **argv) {
	static std::ofstream fp_out;		// static so we can pass out a reference
	std::set<std::string> used_files;	// make sure there're no overlaps
	std::string opt_output;
	opt_feedback = 1;
	opt_frequency_cutoff = 0;
	opt_histogram_restore = -1;
	opt_mer_length = 24;
	opt_print_gc = 0;
	opt_load_lower_bound = 0;
	opt_load_upper_bound = 1;
	opt_nmers = 0;
	int c;
	while ((c = getopt(argc, argv, "ghil:L:m:o:s:S:Vw:z:")) != EOF) {
		switch (c) {
		    case 'g':
			opt_print_gc = 1;
			break;
		    case 'h':
			print_usage();
			break;
		    case 'i':
			opt_feedback = 0;
			break;
		    case 'l':
			std::istringstream(optarg) >> opt_load_lower_bound;
			break;
		    case 'L':
			std::istringstream(optarg) >> opt_load_upper_bound;
			break;
		    case 'm':
			std::istringstream(optarg) >> opt_mer_length;
			if (opt_mer_length < 1) {
				std::cerr << "Error: bad mer length\n";
				exit(1);
			}
			break;
		    case 'o':
			opt_output = optarg;
			if (!used_files.insert(optarg).second) {
				std::cerr << "Error: duplicate file: " << optarg << "\n";
				exit(1);
			}
			break;
		    case 's':
			opt_save_file = optarg;
			if (!used_files.insert(optarg).second) {
				std::cerr << "Error: duplicate file: " << optarg << "\n";
				exit(1);
			}
			break;
		    case 'S':
			opt_histogram_restore = open_compressed(optarg);
			if (opt_histogram_restore == -1) {
				std::cerr << "Error: could not read histogram dump file\n";
				exit(1);
			} else if (!used_files.insert(optarg).second) {
				std::cerr << "Error: duplicate file: " << optarg << "\n";
				exit(1);
			}
			break;
		    case 'V':
			std::cerr << "histogram_hashl version " << VERSION <<
#ifdef COMPRESS_READS
				" (read compression)"
#else
				""
#endif
				<< "\n";
			exit(0);
		    case 'w':
			std::istringstream(optarg) >> opt_frequency_cutoff;
			break;
		    case 'z':
			opt_nmers = get_value(optarg);
			break;
		    default:
			std::cerr << "Error: unknown option " << static_cast<char>(c) << "\n";
			print_usage();
		}
	}
	if (opt_load_lower_bound > opt_load_upper_bound) {
		std::cerr << "Error: lower hash fill fraction must be less than upper\n";
		print_usage();
	}
	if (optind == argc && opt_histogram_restore == -1) {
		std::cerr << "Error: no files to process\n";
		print_usage();
	}
	for (int i(optind); i < argc; ++i) {
		if (!used_files.insert(argv[i]).second) {
			std::cerr << "Error: duplicate file: " << argv[i] << "\n";
			exit(1);
		}
	}
	if (!opt_output.empty()) {
		fp_out.open(opt_output);
		if (!fp_out.is_open()) {
			std::cerr << "Error: could not write to " << opt_output << "\n";
			exit(1);
		}
	}
	return opt_output.empty() ? std::cout : fp_out;
}

// for non-ACGT basepairs, need to split reads

static void get_subread_sizes(const std::string &seq, hashl_metadata &metadata) {
	size_t i(seq.find_first_of("ACGTacgt", 0));
	while (i != std::string::npos) {
		size_t next(seq.find_first_not_of("ACGTacgt", i));
		if (next == std::string::npos) {
			next = seq.size();
		}
		// reads shorter than the mer length are skipped
		if (next - i >= opt_mer_length) {
			metadata.add_read_range(i, next);
		}
		i = seq.find_first_of("ACGTacgt", next);
	}
}

// read file and get total space needed for data and metadata
static void get_read_sizes(const char * const file, hashl_metadata &metadata) {
	const int fd(open_compressed(file));
	if (fd == -1) {
		std::cerr << "Error: open: " << file << "\n";
		exit(1);
	}
	std::string line, seq;
	if (pfgets(fd, line) == -1) {		// empty file
		std::cerr << "Warning: empty file: " << file << "\n";
	} else if (line[0] == '>') {		// fasta file
		do {
			size_t i(1);
			for (; i < line.size() && !isspace(line[i]); ++i) { }
			metadata.add_read(line.substr(1, i - 1));
			seq.clear();
			while (pfgets(fd, line) != -1 && line[0] != '>') {
				seq += line;
			}
			get_subread_sizes(seq, metadata);
		} while (!line.empty());
	} else if (line[0] == '@') {		// fastq file
		do {
			size_t i(1);
			for (; i < line.size() && !isspace(line[i]); ++i) { }
			metadata.add_read(line.substr(1, i - 1));
			if (pfgets(fd, seq) == -1) {		// sequence
				std::cerr << "Error: truncated fastq file: " << file << "\n";
				exit(1);
			}
			get_subread_sizes(seq, metadata);
			// skip quality header and quality
			// (use seq because it'll be the same length as quality)
			if (pfgets(fd, line) == -1 || pfgets(fd, seq) == -1) {
				std::cerr << "Error: truncated fastq file: " << file << "\n";
				exit(1);
			}
		} while (pfgets(fd, line) != -1);
	} else {
		std::cerr << "Error: unknown file format: " << file << "\n";
		exit(1);
	}
	close_compressed(fd);
}

static void count_nmers(hashl &mer_list, const std::vector<size_t> &read_ends) {
	const std::vector<hashl::base_type> &data(mer_list.get_data());
	hashl::key_type key(mer_list), comp_key(mer_list);
	size_t total_read_ranges(0);
	size_t i(0), j(0), k(sizeof(hashl::base_type) * 8 - 2);
	// iterate over all reads (nmers can't cross read boundaries)
	std::vector<size_t>::const_iterator a(read_ends.begin());
	const std::vector<size_t>::const_iterator end_a(read_ends.end());
	for (; a != end_a; ++a, ++total_read_ranges) {
		// print feedback every 10 minutes
		if (opt_feedback && elapsed_time() >= 600) {
			start_time();
			std::cerr << time(0) << ": " << mer_list.size() << " entries used (" << double(100) * mer_list.size() / mer_list.capacity() << ") (" << total_read_ranges << " read ranges)\n";
		}
		const size_t end_i(i + opt_mer_length - 1);
		// load keys with opt_mer_length - 1 basepairs
		for (; i < end_i; ++i) {
			const hashl::base_type c((data[j] >> k) & 3);
			key.push_back(c);
			comp_key.push_front(3 - c);
			if (k) {
				k -= 2;
			} else {
				k = sizeof(hashl::base_type) * 8 - 2;
				++j;
			}
		}
		// run over all nmers, one basepair at a time
		for (; i < *a; ++i) {
			const hashl::base_type c((data[j] >> k) & 3);
			key.push_back(c);
			comp_key.push_front(3 - c);
			if (k) {
				k -= 2;
			} else {
				k = sizeof(hashl::base_type) * 8 - 2;
				++j;
			}
			// increment with bit offset to start of nmer
			if (!mer_list.increment(key, comp_key, 2 * (i + 1 - opt_mer_length))) {
				std::cerr << "Error: ran out of space in hash\n";
				exit(1);
			}
		}
	}
}

// metadata format: # of files, [ filename, # of reads, read offsets, read names ]
// strings are null delimited, # and read lengths are uint64_t
// hashl boilerplate check will ensure byte order is the same read as written

static void read_in_files(int argc, char **argv, hashl &mer_list) {
	const uint64_t file_count(argc - optind);
	hashl_metadata metadata;
	for (size_t i(0); i < file_count; ++i) {
		if (opt_feedback) {
			std::cerr << time(0) << ": Getting read sizes for " << argv[i + optind] << "\n";
		}
		metadata.add_file(argv[i + optind]);
		get_read_sizes(argv[i + optind], metadata);
	}
	std::vector<hashl::base_type> data;
	metadata.read_data(data, opt_feedback);
	if (opt_feedback) {
		std::cerr << time(0) << ": Initializing n-mer hash\n";
	}
	// put data and metadata into mer_list
	mer_list.init(opt_nmers ? opt_nmers : metadata.max_kmers(opt_mer_length), opt_mer_length * 2, data);
	std::vector<char> packed_metadata;
	metadata.pack(packed_metadata);
	mer_list.set_metadata(packed_metadata);
	if (opt_feedback) {
		std::pair<size_t, size_t> read_count(metadata.total_reads());
		std::cerr << time(0) << ": Counting n-mers for " << read_count.first << " reads (" << read_count.second << " ranges)\n";
		start_time();
	}
	const std::vector<size_t> read_ends(metadata.read_ends());
	count_nmers(mer_list, read_ends);
	if (opt_feedback) {
		print_final_input_feedback(mer_list);
	}
	// make sure the fill rate is acceptable
	const double load(double(mer_list.size()) / mer_list.capacity());
	if (load < opt_load_lower_bound || opt_load_upper_bound < load) {
		if (opt_feedback) {
			std::cerr << time(0) << ": hash fill rate out of range, resizing: " << opt_load_lower_bound << " - " << opt_load_upper_bound << ": " << load << "\n";
		}
		mer_list.resize(mer_list.size() * 2 / (opt_load_lower_bound + opt_load_upper_bound));
		if (opt_feedback) {
			print_final_input_feedback(mer_list);
		}
	}
}

int main(int argc, char **argv) {
	std::cerr << std::fixed << std::setprecision(2);
	std::ostream &fp_out(get_opts(argc, argv));
	fp_out << std::fixed << std::setprecision(2);
	hashl mer_list;
	if (opt_histogram_restore != -1) {
		if (opt_feedback) {
			std::cerr << time(0) << ": Initializing n-mer hash\n";
		}
		mer_list.init_from_file(opt_histogram_restore);
		close_compressed(opt_histogram_restore);
		if (opt_feedback) {
			print_final_input_feedback(mer_list);
		}
	}
	if (optind != argc) {
		read_in_files(argc, argv, mer_list);
	}
	if (opt_feedback) {
		std::cerr << time(0) << ": Printing results\n";
	}
	if (opt_frequency_cutoff == 0) {
		print_mer_histogram(fp_out, mer_list);
	} else {
		print_mer_frequency(fp_out, mer_list);
	}
	fp_out.flush();
	if (!opt_save_file.empty()) {
		save_memory(mer_list);
	}
	return 0;
}
