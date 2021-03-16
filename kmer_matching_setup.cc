// This program is meant to process a set of reads (in fastq/a format) to create
// a list of which kmers are part of which read, in a format that allows easy retrieval
// by kmer; the arrays and hashes necessary are then written out to disk.
//
// This program requires a lot of memory, depending on the size of the fastq/a file.

#include "hash.h"	// hash
#include "hist_lib_hash.h"	// add_sequence_mers(), add_sequence_mers_index(), init_mer_constants(), opt_feedback, opt_include, opt_mer_length, opt_skip_size, print_final_input_feedback()
#include "kmer_lookup_info.h"	// KmerLookupInfo
#include "open_compressed.h"	// close_compressed(), get_suffix(), open_compressed(), pfgets()
#include "read.h"	// Read, opt_clip_quality, opt_clip_vector, opt_quality_cutoff
#include "read_file.h"	// ReadFile, opt_strip_tracename
#include "version.h"	// VERSION
#include "write_fork.h"	// close_fork_wait(), pfputs(), write_fork()
#include <getopt.h>	// getopt(), optarg, optind
#include <new>		// new
#include <regex.h>	// REG_EXTENDED, REG_NOSUB
#include <sstream>	// istringstream
#include <stdio.h>	// EOF, fprintf(), stderr, stdout
#include <stdlib.h>	// exit()
#include <string>	// string
#include <sys/types.h>	// size_t
#include <unistd.h>	// STDOUT_FILENO

static int fd_out(STDOUT_FILENO);
static bool opt_track_dups;
static bool opt_warnings;
static size_t opt_batch_size;
static size_t opt_nmers;

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
	fprintf(stderr,
		"usage: kmer_matching_setup [options] file1 [file2] ...\n"
		"    -B ## process seq & qual file in batches of ## reads\n"
		"    -c    clip low quality\n"
		"    -d    when processing in batches, check for duplicates across whole file\n"
		"    -f ## when clipping quality or vector, use ## as the target quality [20]\n"
		"    -h    print this information\n"
		"    -i    turn off status updates\n"
		"    -k ## skip reads smaller than this\n"
		"    -m ## set mer length (1-32) [24]\n"
		"    -o ## print output to file instead of stdout\n"
		"    -p ## don't touch reads not matching pattern (an extended regex)\n"
		"    -q    turn off all warnings\n"
		"    -t    strip first part of trace id\n"
		"    -v    clip vector\n"
		"    -V    print version\n"
		"    -z ## number of possible n-mers to allocate memory for [200m]\n"
		"          (k, m, or g may be suffixed)\n"
	);
	exit(1);
}

static void get_opts(int argc, char **argv) {
	std::string opt_output;
	opt_batch_size = 0;
	opt_clip_quality = 0;
	opt_clip_vector = 0;
	opt_feedback = 1;
	opt_mer_length = 24;
	opt_nmers = static_cast<size_t>(-1);
	opt_quality_cutoff = 20;
	opt_skip_size = 0;
	opt_strip_tracename = 0;
	opt_track_dups = 0;
	opt_warnings = 1;
	int c;
	while ((c = getopt(argc, argv, "B:cdf:hik:m:o:p:qtvVz:")) != EOF) {
		switch (c) {
		    case 'B':
			std::istringstream(optarg) >> c;
			if (c < 0) {
				print_usage();
			}
			opt_batch_size = c;
			break;
		    case 'c':
			opt_clip_quality = 1;
			break;
		    case 'd':
			opt_track_dups = 1;
			break;
		    case 'f':
			std::istringstream(optarg) >> opt_quality_cutoff;
			if (opt_quality_cutoff < 0) {
				print_usage();
			}
			break;
		    case 'h':
			print_usage();
			break;
		    case 'i':
			opt_feedback = 0;
			break;
		    case 'k':
			std::istringstream(optarg) >> c;
			if (c < 0) {
				fprintf(stderr, "Error: invalid skip size %d\n", c);
				print_usage();
			}
			opt_skip_size = c;
			break;
		    case 'm':
			std::istringstream(optarg) >> opt_mer_length;
			if (opt_mer_length < 1 || 32 < opt_mer_length) {
				fprintf(stderr, "Error: bad mer length\n");
				print_usage();
			}
			break;
		    case 'o':
			opt_output = optarg;
			break;
		    case 'p':
			opt_include.initialize(optarg, 0, REG_NOSUB | REG_EXTENDED);
			break;
		    case 'q':
			opt_warnings = 0;
			break;
		    case 't':
			opt_strip_tracename = 1;
			break;
		    case 'v':
			opt_clip_vector = 1;
			break;
		    case 'V':
			fprintf(stderr, "kmer_matching_setup version %s%s\n", VERSION,
#ifdef COMPRESS_READS
" (read compression)"
#else
""
#endif
);
			exit(0);
		    case 'z':
			opt_nmers = get_value(optarg);
			if (opt_nmers == 0) {
				fprintf(stderr, "Error: bad n-mer count %s\n", optarg);
				print_usage();
			}
			break;
		    default:
			fprintf(stderr, "Error: unknown option %c\n", c);
			print_usage();
		}
	}
	if (opt_nmers == static_cast<size_t>(-1)) {
		opt_nmers = 200 * 1024 * 1024;
	}
	if (optind == argc) {
		fprintf(stderr, "Error: no files to process\n");
		print_usage();
	}
	if (!opt_output.empty()) {
		std::string suffix;
		get_suffix(opt_output, suffix);
		std::list<std::string> args;
		if (suffix == ".gz") {
			args.push_back("gzip");
			args.push_back("-c");
		} else if (suffix == ".bz2") {
			args.push_back("bzip2");
			args.push_back("-c");
		} else if (suffix == ".Z") {
			args.push_back("compress");
			args.push_back("-c");
		}
		fd_out = write_fork(args, opt_output);
		if (fd_out == -1) {
			fprintf(stderr, "Error: could not write to %s\n", opt_output.c_str());
			exit(1);
		}
	}
}

// count number of all kmers in the reads, along with the number of reads and the
// size of all the read names (so we can pre-allocate all those arrays later)

static int count_kmers(char **argv, char ** const end_argv, hash &mer_list, size_t &total_reads, size_t &total_names_size) {
	int err(0);
	mer_list.init(opt_nmers);
	for (; argv != end_argv; ++argv) {
		if (opt_feedback) {
			fprintf(stderr, "Reading in %s\n", *argv);
		}
		ReadFile file(*argv, opt_batch_size, opt_track_dups);
		if (file.seq_file.empty()) {
			++err;
			continue;
		}
		while (file.read_batch(opt_warnings) != -1) {
			if (!add_sequence_mers(file.read_list.begin(), file.read_list.end(), mer_list, total_reads)) {
				fprintf(stderr, "Error: %lu reads processed: n-mer list incomplete - give a larger -z value\n", total_reads);
				return err + 1;
			}
			total_reads += file.read_list.size();
			std::list<Read>::const_iterator a(file.read_list.begin());
			const std::list<Read>::const_iterator end_a(file.read_list.end());
			for (; a != end_a; ++a) {
				total_names_size += a->name().size();
			}
		}
	}
	return err;
}

// go back through files and store read names, read kmer hits
static void index_kmers(char **argv, char ** const end_argv, KmerLookupInfo &kmers, const size_t total_reads) {
	size_t reads_processed(0);
	for (; argv != end_argv; ++argv) {
		if (opt_feedback) {
			fprintf(stderr, "Reading in %s\n", *argv);
		}
		ReadFile file(*argv, opt_batch_size, opt_track_dups);
		if (file.seq_file.empty()) {
			continue;
		}
		while (file.read_batch(opt_warnings) != -1) {
			add_sequence_mers_index(file.read_list.begin(), file.read_list.end(), kmers, reads_processed, total_reads);
			reads_processed += file.read_list.size();
		}
	}
}

int main(int argc, char **argv) {
	get_opts(argc, argv);
	if (opt_feedback) {
		fprintf(stderr, "%lu: Initializing n-mer hash\n", time(NULL));
	}
	init_mer_constants();
	size_t total_reads(0), total_name_size(0);
	// allocate this dynamically so we can deallocate it later and free up the memory
	hash *mer_list(new hash);
	const int err(count_kmers(&argv[optind], &argv[argc], *mer_list, total_reads, total_name_size));
	if (err != 0) {
		return err;
	}
	if (opt_feedback) {
		print_final_input_feedback(*mer_list);
		fprintf(stderr, "Initializing kmer lookups\n");
	}
	KmerLookupInfo kmers(opt_mer_length + 1, total_reads, total_name_size, *mer_list);
	delete mer_list;				// free up memory (hopefully)
	index_kmers(&argv[optind], &argv[argc], kmers, total_reads);
	if (opt_feedback) {
		fprintf(stderr, "%lu: all reads processed\n", time(NULL));
		fprintf(stderr, "Saving kmer lookup info\n");
	}
	kmers.save(fd_out);
	close_fork_wait(fd_out);
	if (opt_feedback) {
		fprintf(stderr, "%lu: kmer lookup info saved\n", time(NULL));
	}
	return 0;
}
