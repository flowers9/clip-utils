#include "hashl.h"	// hashl
#include "hashl_metadata.h"	// hashl_metadata
#include "open_compressed.h"	// close_compressed(), get_suffix(), open_compressed()
#include "version.h"	// VERSION
#include "write_fork.h"	// close_fork(), write_fork()
#include <getopt.h>	// getopt(), optarg, optind
#include <iomanip>	// fixed, setprecision()
#include <iostream>	// cerr, cout
#include <list>		// list<>
#include <sstream>	// istringstream
#include <stdlib.h>	// exit()
#include <string>	// string
#include <vector>	// vector<>

// combine a set of reference hashes (and save this, if asked) and then
// go through a target's hash seeing which kmers match; multiple target
// files are treated as one big target file

static bool opt_print_histogram;
static int opt_fastq_max_kmer_frequency;
static int opt_fastq_min_kmer_frequency;
static int opt_hash_load;
static int opt_max_kmer_sharing;
static int opt_reference_max_kmer_frequency;
static int opt_reference_min_kmer_frequency;
static size_t opt_nmers;
static std::string opt_hash_save;
static std::string opt_index_save;
static std::string opt_purged_hash_save;
static std::string opt_results_save;
static std::vector<std::string> opt_reference_list;

// return the number represented by s, which may be suffixed by a k, m, or g
// which act as multipliers to the base amount

static size_t get_value(const std::string s) {
	size_t x;
	std::istringstream(s) >> x;
	// validate string - digits optionally followed by k, m, or g
	const size_t i(s.find_first_not_of("0123456789"));
	if (i == std::string::npos) {		// plain number
		return x;
	} else if (i + 1 == s.length()) {	// has multiplier
		switch (s[i]) {
		    case 'g':
			return x << 30;
		    case 'm':
			return x << 20;
		    case 'k':
			return x << 10;
		    default:			// invalid multiplier
			std::cerr << "Error: invalid unit suffix: " << s[i] << '\n';
			exit(1);
		}
	} else {				// bad value
		std::cerr << "Error: invalid number: " << s << '\n';
		exit(1);
	}
}

static void print_usage() {
	std::cerr << "usage: screen_kmers_by_ref [target_hash1 [target_hash2 ...]]\n"
		"    -h    print this help\n"
		"    -H    print histogram of combined reference\n"
		"    -f ## fastq min kmer frequency\n"
		"    -F ## fastq max kmer frequency\n"
		"    -i ## save an index to reference kmers found in target hashes\n"
		"    -m ## reference min kmer frequency\n"
		"    -M ## reference max kmer frequency [1]\n"
		"    -o ## save results to a hash dump for later processing\n"
		"    -p ## save combined hash, but purge invalid values before saving\n"
		"    -r ## add reference file (may be specified multiple times)\n"
		"    -s ## save resulting combined reference hash\n"
		"    -S ## load histogram memory dump from given file\n"
		"    -u ## only count kmers shared with at most ## references\n"
		"          (negative values mean shared by all but ##) [-1]\n"
		"          (only affects target screening, not -o output)\n"
		"    -V    print version\n"
		"    -z ## number of unique kmers to pre-allocate for combined reference hash\n"
		"          (k, m, or g may be suffixed)\n";
	exit(1);
}

static void get_opts(const int argc, char * const * const argv) {
	opt_fastq_max_kmer_frequency = hashl::max_small_value;
	opt_fastq_min_kmer_frequency = 0;
	opt_hash_load = -1;
	opt_max_kmer_sharing = -1;
	opt_nmers = 0;
	opt_print_histogram = 0;
	opt_reference_max_kmer_frequency = 1;
	opt_reference_min_kmer_frequency = 0;
	int c;
	while ((c = getopt(argc, argv, "hHf:F:i:m:M:o:p:r:s:S:u:Vz:")) != EOF) {
		switch (c) {
		    case 'h':
			print_usage();
			break;
		    case 'H':
			opt_print_histogram = 1;
			break;
		    case 'f':
			std::istringstream(optarg) >> opt_fastq_min_kmer_frequency;
			break;
		    case 'F':
			std::istringstream(optarg) >> opt_fastq_max_kmer_frequency;
			break;
		    case 'i':
			opt_index_save = optarg;
			break;
		    case 'm':
			std::istringstream(optarg) >> opt_reference_min_kmer_frequency;
			break;
		    case 'M':
			std::istringstream(optarg) >> opt_reference_max_kmer_frequency;
			break;
		    case 'o':
			opt_results_save = optarg;
			break;
		    case 'r':
			opt_reference_list.push_back(optarg);
			break;
		    case 'p':
			opt_purged_hash_save = optarg;
			break;
		    case 's':
			opt_hash_save = optarg;
			break;
		    case 'S':
			opt_hash_load = open_compressed(optarg);
			if (opt_hash_load == -1) {
				std::cerr << "Error: could not read histogram dump file\n";
				exit(1);
			}
			break;
		    case 'u':
			std::istringstream(optarg) >> opt_max_kmer_sharing;
			break;
		    case 'V':
			std::cerr << "dot_hashl version " << VERSION << '\n';
			exit(0);
		    case 'z':
			opt_nmers = get_value(optarg);
			break;
		    default:
			std::cerr << "Error: unknown option " << char(c) << '\n';
			print_usage();
		}
	}
	if (opt_reference_list.empty() && opt_hash_load == -1) {
		std::cerr << "Error: no reference files given\n";
		print_usage();
	}
	if (opt_hash_load != -1 && !opt_hash_save.empty()) {
		std::cerr << "Warning: ignoring -s option because of -S\n";
	}
}

static void save_hash(const hashl &mer_list, const std::string &filename) {
	std::string suffix;
	get_suffix(filename, suffix);
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
	const int fd(write_fork(args, filename));
	if (fd == -1) {
		std::cerr << "Error: could not save hash " << filename << '\n';
		exit(1);
	}
	mer_list.save(fd);
	close_fork(fd);
}

static void save_index(hashl &mer_list, const std::string &filename) {
	// index file is never compressed, so don't even check the suffix
	const int fd(write_fork(std::list<std::string>(), filename));
	if (fd == -1) {
		std::cerr << "Error: could not save index " << filename << '\n';
		exit(1);
	}
	// note: trashes mer_list
	mer_list.save_index(fd);
	close_fork(fd);
}

// print histogram of n-mer occurrences

static void print_mer_histogram(const hashl &mer_list) {
	std::map<hashl::small_value_type, unsigned long> counts;
	hashl::const_iterator a(mer_list.cbegin());
	const hashl::const_iterator end_a(mer_list.cend());
	for (; a != end_a; ++a) {
		++counts[*a];
	}
	std::cout << std::fixed << std::setprecision(2);
	double i(0);
	double total(mer_list.size());
	std::map<hashl::small_value_type, unsigned long>::const_iterator c(counts.begin());
	const std::map<hashl::small_value_type, unsigned long>::const_iterator end_c(counts.end());
	for (; c != end_c; ++c) {
		const double x(double(100) * c->second);
		i += x;
		std::cout << c->first << ' ' << c->second << ' ' << x / total << ' ' << i / total << "\n";
	}
}

static bool load_and_combine_hashes(hashl &kmer_hash, const std::vector<std::string> &files, const int min_cutoff, const int max_cutoff, const size_t starting_hash_size = 0) {
	hashl tmp_hash;			// declare outside loop so memory can get reused
	// load in reference saved hashes
	for (const auto &file : files) {
		std::cerr << time(0) << ": reading " << file << '\n';
		const int fd(open_compressed(file));
		if (fd == -1) {
			std::cerr << "Error: could not read saved hash: " << file << '\n';
			return 0;
		}
		tmp_hash.init_from_file(fd);
		close_compressed(fd);
		if (&file == &files[0]) {	// set the bit width, possibly preallocate
			std::vector<hashl::base_type> tmp_data;
			kmer_hash.init(starting_hash_size, tmp_hash.bits(), tmp_data);
		}
		if (!kmer_hash.add(tmp_hash, min_cutoff, max_cutoff)) {
			std::cerr << "Error: failed to add hash\n";
			return 0;
		}
		std::cerr << time(0) << ": size " << kmer_hash.size() << ' ' << double(100) * kmer_hash.size() / kmer_hash.capacity() << "% " << kmer_hash.capacity() << '\n';
		close_compressed(fd);
	}
	return 1;
}

static void cross_ref_stdout(const hashl &reference_kmers, const hashl &fastq_kmers) {
	hashl::const_iterator a(fastq_kmers.cbegin());
	const hashl::const_iterator end_a(fastq_kmers.cend());
	hashl::key_type key(fastq_kmers.bits(), fastq_kmers.words()), comp_key(fastq_kmers.bits(), fastq_kmers.words());
	std::string s;
	for (; a != end_a; ++a) {
		if (*a && *a != hashl::invalid_value) {
			a.key(key);
			// .value() returns 0 if key not found
			const hashl::small_value_type x(reference_kmers.value(key));
			if (x && x != hashl::invalid_value && x <= static_cast<unsigned int>(opt_max_kmer_sharing)) {
				key.convert_to_string(s);
				std::cout << s << ' ' << static_cast<unsigned int>(x) << '\n';
				comp_key.make_complement(key);
				if (key != comp_key) {
					comp_key.convert_to_string(s);
					std::cout << s << ' ' << static_cast<unsigned int>(x) << ' ' << '\n';
				}
			}
		}
	}
}

static void cross_ref_save(const hashl &reference_kmers, hashl &fastq_kmers) {
	hashl::iterator a(fastq_kmers.begin());
	const hashl::iterator end_a(fastq_kmers.end());
	hashl::key_type key(fastq_kmers.bits(), fastq_kmers.words());
	for (; a != end_a; ++a) {
		if (*a && *a != hashl::invalid_value) {
			a.key(key);
			// .value() returns 0 if key not found
			const hashl::small_value_type x(reference_kmers.value(key));
			if (!x || x == hashl::invalid_value || x > static_cast<unsigned int>(opt_max_kmer_sharing)) {
				*a = hashl::invalid_value;
			}
		}
	}
	fastq_kmers.purge_invalid_values();
	if (!opt_results_save.empty()) {
		save_hash(fastq_kmers, opt_results_save);
	}
	if (!opt_index_save.empty()) {
		// trashes fastq_kmers
		save_index(fastq_kmers, opt_index_save);
	}
}

void save_purged_hash(hashl &reference_kmers) {
	reference_kmers.purge_invalid_values();
	save_hash(reference_kmers, opt_purged_hash_save);
}

int main(const int argc, char * const * const argv) {
	std::cerr << std::fixed << std::setprecision(2);
	get_opts(argc, argv);
	hashl reference_kmers;
	if (opt_hash_load != -1) {
		reference_kmers.init_from_file(opt_hash_load);
		close_compressed(opt_hash_load);
		hashl_metadata md;
		md.unpack(reference_kmers.get_metadata());
		if (opt_max_kmer_sharing < 0) {
			opt_max_kmer_sharing += md.file_count();
		}
		if (opt_max_kmer_sharing < 1) {
			opt_max_kmer_sharing = 1;
		} else if (static_cast<unsigned int>(opt_max_kmer_sharing) > md.file_count()) {
			opt_max_kmer_sharing = md.file_count();
		}
	} else if (!load_and_combine_hashes(reference_kmers, opt_reference_list, opt_reference_min_kmer_frequency, opt_reference_max_kmer_frequency, opt_nmers)) {
		return 1;
	} else {
		if (opt_max_kmer_sharing < 0) {
			opt_max_kmer_sharing += opt_reference_list.size();
		}
		if (opt_max_kmer_sharing < 1) {
			opt_max_kmer_sharing = 1;
		} else if (static_cast<unsigned int>(opt_max_kmer_sharing) > opt_reference_list.size()) {
			opt_max_kmer_sharing = opt_reference_list.size();
		}
		if (2 * reference_kmers.size() > reference_kmers.capacity() || !opt_hash_save.empty()) {
			// reduce load to 50% load to optimize speed of lookups
			// or increase to 50% to reduce size of save file
			std::cerr << time(0) << ": setting hash to 50% load\n";
			reference_kmers.resize(2 * reference_kmers.size());
			std::cerr << time(0) << ": size " << reference_kmers.size() << ' ' << double(100) * reference_kmers.size() / reference_kmers.capacity() << "% " << reference_kmers.capacity() << '\n';
		}
		if (!opt_hash_save.empty()) {
			save_hash(reference_kmers, opt_hash_save);
		}
	}
	if (opt_print_histogram) {
		print_mer_histogram(reference_kmers);
	}
	// do this after print_mer_histogram so invalid values still show up in the histogram
	if (!opt_purged_hash_save.empty()) {
		save_purged_hash(reference_kmers);
	}
	if (optind == argc) {	// just saving the created hash, presumably
		return 0;
	}
	std::vector<std::string> target_hashes;
	for (int i(optind); i < argc; ++i) {
		target_hashes.push_back(argv[i]);
	}
	hashl fastq_kmers;
	if (!load_and_combine_hashes(fastq_kmers, target_hashes, opt_fastq_min_kmer_frequency, opt_fastq_max_kmer_frequency)) {
		return 1;
	}
	std::cerr << time(0) << ": processing kmers\n";
	if (!opt_results_save.empty() || !opt_index_save.empty()) {
		cross_ref_save(reference_kmers, fastq_kmers);
	} else {
		cross_ref_stdout(reference_kmers, fastq_kmers);
	}
	return 0;
}
