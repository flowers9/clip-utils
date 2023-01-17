#include "hashl.h"	// hashl
#include "next_prime.h"	// next_prime()
#include "open_compressed.h"	// close_compressed(), get_suffix(), open_compressed(), pfread()
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
		"    -S ## save resulting combined reference hash as ##\n"
		"    -u ## only count kmers shared with at most ## references\n"
		"          (negative values mean shared by all but ##) [-1]\n"
		"    -V    print version\n";
	exit(1);
}

static void get_opts(const int argc, char * const * argv) {
	opt_fastq_max_kmer_frequency = hashl::max_small_value;
	opt_fastq_min_kmer_frequency = 0;
	opt_max_kmer_sharing = -1;
	opt_reference_max_kmer_frequency = 1;
	opt_reference_min_kmer_frequency = 0;
	int c;
	while ((c = getopt(argc, argv, "hf:F:m:M:r:S:u:V")) != EOF) {
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
		    case 'S':
			opt_hash_save = optarg;
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
}

static void save_memory(const hashl &mer_list) {
	std::string suffix;
	get_suffix(opt_hash_save, suffix);
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
	const int fd(write_fork(args, opt_hash_save));
	if (fd == -1) {
		fprintf(stderr, "Error: could not save memory\n");
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
		if (!kmer_hash.add(tmp_hash, min_cutoff, max_cutoff)) {
			std::cerr << "Error: failed to add hash\n";
			return 0;
		}
		std::cerr << "size " << kmer_hash.size() << "\n";
		close_compressed(fd);
	}
	return 1;
}

int main(const int argc, char * const * argv) {
	get_opts(argc, argv);
	hashl reference_kmers;
	if (!load_and_combine_hashes(reference_kmers, opt_reference_list, opt_reference_min_kmer_frequency, opt_reference_max_kmer_frequency)) {
		return 1;
	}
	if (!opt_hash_save.empty()) {
		save_memory(reference_kmers);
	}
	hashl fastq_kmers;
	{	// scope to limit fastq_files
		std::vector<std::string> fastq_files;
		for (int i(optind); i < argc; ++i) {
			fastq_files.push_back(argv[i]);
		}
		if (!load_and_combine_hashes(fastq_kmers, fastq_files, opt_reference_min_kmer_frequency, opt_reference_max_kmer_frequency)) {
			return 1;
		}
	}
	std::cerr << "processing kmers\n";
	mer_list[0].print_kmer_matching(&mer_list[1], opt_reference_list.size());
}
