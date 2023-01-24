#include "hashl.h"	// hashl
#include "next_prime.h"	// next_prime()
#include "open_compressed.h"	// close_compressed(), get_suffix(), open_compressed()
#include "write_fork.h"	// close_fork(), write_fork()
#include "version.h"	// VERSION
#include <getopt.h>	// getopt(), optarg, optind
#include <iostream>	// cerr, cout
#include <list>		// list<>
#include <sstream>	// istringstream
#include <stdlib.h>	// exit()
#include <string>	// string
#include <vector>	// vector<>

// create a hash of useful kmers (ones unique within each reference it's
// in and not present in all references) and save it

static int opt_fastq_max_kmer_frequency;
static int opt_fastq_min_kmer_frequency;
static int opt_hash_load;
static int opt_max_kmer_sharing;
static int opt_reference_max_kmer_frequency;
static int opt_reference_min_kmer_frequency;
static std::string opt_hash_save;
static std::vector<std::string> opt_reference_list;

static void print_usage() {
	std::cerr << "usage: dot_hash saved_hash1 saved_hash2 ...\n"
		"    -h    print this help\n"
		"    -f ## fastq min kmer frequency\n"
		"    -F ## fastq max kmer frequency\n"
		"    -m ## reference min kmer frequency\n"
		"    -M ## reference max kmer frequency [1]\n"
		"    -r ## add reference file (may be specified multiple times)\n"
		"    -s ## save resulting combined reference hash\n"
		"    -S ## load histogram memory dump from given file\n"
		"    -u ## only count kmers shared with at most ## references\n"
		"          (negative values mean shared by all but ##) [-1]\n"
		"    -V    print version\n";
	exit(1);
}

static void get_opts(const int argc, char * const * const argv) {
	opt_fastq_max_kmer_frequency = hashl::max_small_value;
	opt_fastq_min_kmer_frequency = 0;
	opt_hash_load = -1;
	opt_max_kmer_sharing = -1;
	opt_reference_max_kmer_frequency = 1;
	opt_reference_min_kmer_frequency = 0;
	int c;
	while ((c = getopt(argc, argv, "hf:F:m:M:r:s:S:u:V")) != EOF) {
		switch (c) {
		    case 'h':
			print_usage();
			break;
		    case 'f':
			std::istringstream(optarg) >> opt_fastq_min_kmer_frequency;
			break;
		    case 'F':
			std::istringstream(optarg) >> opt_fastq_max_kmer_frequency;
			break;
		    case 'm':
			std::istringstream(optarg) >> opt_reference_min_kmer_frequency;
			break;
		    case 'M':
			std::istringstream(optarg) >> opt_reference_max_kmer_frequency;
			break;
		    case 'r':
			opt_reference_list.push_back(optarg);
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
			std::cerr << "dot_hash version " << VERSION << "\n";
			exit(0);
		    default:
			std::cerr << "Error: unknown option " << char(c) << "\n";
			print_usage();
		}
	}
	if (opt_reference_list.empty()) {
		std::cerr << "Error: no reference files given\n";
		print_usage();
	} else if (opt_reference_list.size() + argc - optind < 2) {
		std::cerr << "Error: only one file specified\n";
		exit(1);
	}
	if (opt_max_kmer_sharing < 0) {
		opt_max_kmer_sharing += opt_reference_list.size();
	}
	if (opt_max_kmer_sharing < 1 || opt_reference_list.size() < static_cast<size_t>(opt_max_kmer_sharing)) {
		std::cerr << "Error: -u option out of range: 1-" << opt_reference_list.size() << ": " << opt_max_kmer_sharing << "\n";
		exit(1);
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
		std::cerr << "Error: could not save hash " << filename << "\n";
		exit(1);
	}
	mer_list.save(fd);
	close_fork(fd);
}

static bool load_and_combine_hashes(hashl &kmer_hash, const std::vector<std::string> &files, const int min_cutoff, const int max_cutoff) {
	// load in reference saved hashes
	std::vector<std::string>::const_iterator a(files.begin());
	const std::vector<std::string>::const_iterator end_a(files.end());
	std::cerr << "reading " << *a << "\n";
	int fd(open_compressed(*a));
	if (fd == -1) {
		std::cerr << "Error: could not read saved hash: " << *a << "\n";
		return 0;
	}
	kmer_hash.init_from_file(fd);
kmer_hash.print();
std::cout << "\n";
	kmer_hash.normalize(min_cutoff, max_cutoff);
	hashl tmp_hash;			// declare outside loop so memory can get reused
	for (++a; a != end_a; ++a) {
		std::cerr << "reading " << *a << "\n";
		fd = open_compressed(*a);
		if (fd == -1) {
			std::cerr << "Error: could not read saved hash: " << *a << "\n";
			return 0;
		}
		tmp_hash.init_from_file(fd);
tmp_hash.print();
std::cout << "\n";
		if (!kmer_hash.add(tmp_hash, min_cutoff, max_cutoff)) {
			std::cerr << "Error: failed to add hash\n";
			return 0;
		}
		std::cerr << "size " << kmer_hash.size() << "\n";
		close_compressed(fd);
	}
std::cout << "\n";
kmer_hash.print();
std::cout << "\n";
	return 1;
}

static void cross_ref(const hashl &reference_kmers, const hashl &fastq_kmers) {
	hashl::const_iterator a(fastq_kmers.begin());
	const hashl::const_iterator end_a(fastq_kmers.end());
	hashl::key_type key(fastq_kmers), comp_key(fastq_kmers);
	std::string s;
	for (; a != end_a; ++a) {
		if (a.offset() != hashl::invalid_key && a.value() && a.value() != hashl::invalid_value) {
			a.get_key(key);
			// returns 0 if key not found
			const hashl::small_value_type x(reference_kmers.value(key));
			if (x && x != hashl::invalid_value && x <= opt_max_kmer_sharing) {
				key.convert_to_string(s);
				std::cout << s << ' ' << static_cast<unsigned int>(x) << "\n";
				comp_key.make_complement(key);
				if (key != comp_key) {
					comp_key.convert_to_string(s);
					std::cout << s << ' ' << static_cast<unsigned int>(x) << "\n";
				}
			}
		}
	}
}

int main(const int argc, char * const * const argv) {
	get_opts(argc, argv);
	hashl reference_kmers;
	if (opt_hash_load != -1) {
		reference_kmers.init_from_file(opt_hash_load);
	} else if (!load_and_combine_hashes(reference_kmers, opt_reference_list, opt_reference_min_kmer_frequency, opt_reference_max_kmer_frequency)) {
		return 1;
	} else if (!opt_hash_save.empty()) {
		save_hash(reference_kmers, opt_hash_save);
	}
	if (optind == argc) {	// just saving the created hash, presumably
		return 0;
	}
	std::vector<std::string> fastq_files;
	for (int i(optind); i < argc; ++i) {
		fastq_files.push_back(argv[i]);
	}
	hashl fastq_kmers;
	if (!load_and_combine_hashes(fastq_kmers, fastq_files, opt_fastq_min_kmer_frequency, opt_fastq_max_kmer_frequency)) {
		return 1;
	}
	std::cerr << "processing kmers\n";
	cross_ref(reference_kmers, fastq_kmers);
	return 0;
}
