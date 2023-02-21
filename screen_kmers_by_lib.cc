#include "hashl.h"	// hashl
#include "hashl_metadata.h"	// hashl_metadata
#include "open_compressed.h"	// close_compressed(), get_suffix(), open_compressed()
#include "time_used.h"	// elapsed_time(), start_time()
#include "version.h"	// VERSION
#include "write_fork.h"	// close_fork(), write_fork()
#include <errno.h>	// errno
#include <getopt.h>	// getopt(), optarg, optind
#include <iomanip>	// fixed, setprecision()
#include <iostream>	// cerr, cout
#include <list>		// list<>
#include <sstream>	// istringstream
#include <stdio.h>	// rename()
#include <stdlib.h>	// exit()
#include <string.h>	// strerror()
#include <string>	// string
#include <time.h>	// time()
#include <vector>	// vector<>

// take an existing hash, count hits against a library, and mark
// kmers outside the given ranges as invalid

static bool opt_feedback;
static int opt_max_kmer_frequency;
static int opt_min_kmer_frequency;
static size_t opt_mer_length;
static std::string opt_output_hash;

static void print_usage() {
	std::cerr << "usage: screen_kmers_by_lib reference_hash library.fastx [more_library.fastx [...] ]\n"
		"	   multiple library files are treated as one large file - to screen against\n"
		"          multiple libraries, you have to run this program once per library\n"
		"    -h    print this help\n"
		"    -f ## min kmer frequency [1]\n"
		"    -F ## max kmer frequency [" << static_cast<unsigned int>(hashl::max_small_value) << "]\n"
		"    -o ## output file for resulting hash [overwrite original hash]\n"
		"    -q    don't print status updates\n"
		"    -V    print version\n";
	exit(1);
}

static void get_opts(const int argc, char * const * const argv) {
	opt_feedback = 1;
	opt_max_kmer_frequency = hashl::max_small_value;
	opt_min_kmer_frequency = 1;
	int c;
	while ((c = getopt(argc, argv, "hf:F:o:qV")) != EOF) {
		switch (c) {
		    case 'h':
			print_usage();
			break;
		    case 'f':
			std::istringstream(optarg) >> opt_min_kmer_frequency;
			break;
		    case 'F':
			std::istringstream(optarg) >> opt_max_kmer_frequency;
			break;
		    case 'o':
			opt_output_hash = optarg;
			break;
		    case 'q':
			opt_feedback = 0;
			break;
		    case 'V':
			std::cerr << "dot_hashl version " << VERSION << '\n';
			exit(0);
		    default:
			std::cerr << "Error: unknown option " << char(c) << '\n';
			print_usage();
		}
	}
	if (opt_min_kmer_frequency < 1) {
		std::cerr << "Error: -f less than one\n";
		exit(1);
	} else if (static_cast<unsigned int>(opt_min_kmer_frequency) > hashl::max_small_value) {
		std::cerr << "Error: -f greater than " << static_cast<unsigned int>(hashl::max_small_value) << '\n';
		exit(1);
	} else if (opt_min_kmer_frequency > opt_max_kmer_frequency) {
		std::cerr << "Error: -f greater than -F\n";
		exit(1);
	} else if (static_cast<unsigned int>(opt_max_kmer_frequency) > hashl::max_small_value) {
		std::cerr << "Error: -F greater than " << static_cast<unsigned int>(hashl::max_small_value) << '\n';
		exit(1);
	}
	if (optind + 2 > argc) {
		std::cerr << "Error: incorrect number of parameters\n";
		print_usage();
	}
	// we'll just overwrite the original hash if not given an alternative
	if (opt_output_hash.empty()) {
		opt_output_hash = argv[optind];
	}
}

static void save_hash(const hashl &mer_list, const std::string &filename) {
	const std::string tmp(filename + ".tmp");
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
	const int fd(write_fork(args, tmp));
	if (fd == -1) {
		std::cerr << "Error: could not save hash " << filename << '\n';
		exit(1);
	}
	mer_list.save(fd);
	close_fork_wait(fd);
	if (rename(tmp.c_str(), filename.c_str()) == -1) {
		std::cerr << "Error: rename: " << tmp << ": " << filename << ": " << strerror(errno) << '\n';
	}
}

// return 0-3, exit for bad sequence
static hashl::base_type convert_char(const char c) {
	switch (c) {
	    case 'A':
	    case 'a':
		return 0;
	    case 'C':
	    case 'c':
		return 1;
	    case 'G':
	    case 'g':
		return 2;
	    case 'T':
	    case 't':
		return 3;
	    default:
		std::cerr << "Error: non-ACGT basepair: " << static_cast<char>(c) << '\n';
		exit(1);
	}
}

static void count_sequence_mers(hashl &reference_kmers, const std::string &seq, size_t i, const size_t end) {
	hashl::key_type key(reference_kmers), comp_key(reference_kmers);
	const size_t preload_end(i + opt_mer_length - 1);
	for (; i < preload_end; ++i) {
		const hashl::base_type c(convert_char(seq[i]));
		key.push_back(c);
		comp_key.push_front(3 - c);
	}
	for (; i < end; ++i) {
		const hashl::base_type c(convert_char(seq[i]));
		key.push_back(c);
		comp_key.push_front(3 - c);
		reference_kmers.increment(key, comp_key);
	}
}

// for each range of value basepairs (if at least opt_mer_length in length), count kmers

static void process_sequence(hashl &reference_kmers, const std::string &seq) {
	size_t i(seq.find_first_of("ACGTacgt", 0));
	while (i != std::string::npos) {
		size_t next(seq.find_first_not_of("ACGTacgt", i));
		if (next == std::string::npos) {
			next = seq.size();
		}
		// reads shorter than the mer length are skipped
		if (next - i >= opt_mer_length) {
			count_sequence_mers(reference_kmers, seq, i, next);
		}
		i = seq.find_first_of("ACGTacgt", next);
	}
}

// get total match counts for reference kmers
static void process_library(hashl &reference_kmers, const std::string &library_file) {
	const int fd(open_compressed(library_file));
	if (fd == -1) {
		std::cerr << "Error: open: " << library_file << '\n';
		exit(1);
	}
	if (opt_feedback) {
		std::cerr << time(0) << ": processing " << library_file << '\n';
		start_time();
	}
	size_t read_count(0);
	std::string line, seq;
	if (pfgets(fd, line) == -1) {		// empty file
		std::cerr << "Error: File is empty: " << library_file << '\n';
		exit(1);
	} else if (line[0] == '>') {		// fasta file
		do {
			seq.clear();
			while (pfgets(fd, line) != -1 && line[0] != '>') {
				seq += line;
			}
			process_sequence(reference_kmers, seq);
			++read_count;
			if (opt_feedback && elapsed_time() >= 600) {
				start_time();
				std::cerr << time(0) << ": " << read_count << " reads processed\n";
			}
		} while (!line.empty());
	} else if (line[0] == '@') {		// fastq file
		do {
			if (pfgets(fd, seq) == -1) {		// sequence
				std::cerr << "Error: truncated fastq file: " << library_file << '\n';
				exit(1);
			}
			process_sequence(reference_kmers, seq);
			++read_count;
			if (opt_feedback && elapsed_time() >= 600) {
				start_time();
				std::cerr << time(0) << ": " << read_count << " reads processed\n";
			}
			// skip quality header and quality
			// (use seq because it'll be the same length as quality)
			if (pfgets(fd, line) == -1 || pfgets(fd, seq) == -1) {
				std::cerr << "Error: truncated fastq file: " << library_file << '\n';
				exit(1);
			}
		} while (pfgets(fd, line) != -1);
	} else {
		std::cerr << "Error: unknown file format: " << library_file << '\n';
		exit(1);
	}
	close_compressed(fd);
	if (opt_feedback) {
		std::cerr << time(0) << ": " << read_count << " reads processed\n";
	}
}

int main(const int argc, char * const * const argv) {
	std::cerr << std::fixed << std::setprecision(2);
	get_opts(argc, argv);
	// load reference hash
	hashl reference_kmers;
	const int fd(open_compressed(argv[optind]));
	if (fd == -1) {
		std::cerr << "Error: open: " << argv[optind] << '\n';
		return 1;
	}
	if (opt_feedback) {
		std::cerr << time(0) << ": reading in reference hash\n";
	}
	reference_kmers.init_from_file(fd);
	close_compressed(fd);
	reference_kmers.filtering_prep();
	opt_mer_length = reference_kmers.bits() / 2;
	for (int i(optind + 1); i < argc; ++i) {
		process_library(reference_kmers, argv[i]);
	}
	reference_kmers.filtering_finish(opt_min_kmer_frequency, opt_max_kmer_frequency);
	if (opt_feedback) {
		std::cerr << time(0) << ": saving reference hash\n";
	}
	save_hash(reference_kmers, opt_output_hash);
	if (opt_feedback) {
		std::cerr << time(0) << ": save finished\n";
	}
	return 0;
}
