#include "breakup_line.h"	/* break_line_exact() */
#include "hashn.h"	/* hashn */
#include "hist_lib_hashn.h"	/* add_sequence_mers(), init_mer_constants(), opt_exclude, opt_feedback, opt_include, opt_mask_lowercase, opt_phred20_anchor, opt_repeat_coverage, opt_reverse_mask, opt_repeat_threshold, opt_repeat_threshold_upper, opt_skip_size, screen_repeats() */
#include "open_compressed.h"	/* close_compressed(), open_compressed(), pfgets() */
#include "read.h"	/* Read, opt_clip_quality, opt_clip_vector, opt_quality_cutoff */
#include "read_file.h"	/* ReadFile, opt_strip_tracename */
#include "version.h"	/* VERSION */
#include <errno.h>	/* errno */
#include <getopt.h>	// getopt(), optarg, optind
#include <list>		/* list<> */
#include <map>		/* map<> */
#include <regex.h>	/* REG_EXTENDED, REG_NOSUB */
#include <stdio.h>	/* EOF, FILE, fclose(), fopen(), fprintf(), stderr */
#include <stdlib.h>	// atoi(), exit()
#include <string.h>	/* strerror() */
#include <string>	/* string */
#include <sys/types.h>	/* size_t */
#include <utility>	// pair<>
#include <vector>	// vector<>

static bool opt_aggregate;
static bool opt_hash_clean;
static bool opt_limit_printout;
static bool opt_print_percent_masked;
static bool opt_print_range;
static bool opt_split;
static bool opt_track_dups;
static bool opt_warnings;
static int opt_histogram_restore;
static int opt_mer_length;
static size_t opt_batch_size;
static size_t opt_nmers;
static std::list<std::string> hist_files;
static std::string opt_suffix;

static FILE *open_output_file(std::string filename = "") {
	FILE *fp = stdout;
	if (!filename.empty()) {
		filename += opt_suffix;
		fp = fopen(filename.c_str(), "w");
		if (fp == NULL) {
			fprintf(stderr, "Error: could not write to %s\n", filename.c_str());
			return NULL;
		}
	}
	return fp;
}

static void close_output_file(FILE *fp) {
	if (fp != stdout) {
		fclose(fp);
	}
}

/*
 * print full reads, with high repeat regions masked out;
 * output goes to filename + opt_suffix
 */

static void print_unique_sequence(std::list<Read>::iterator a, const std::list<Read>::iterator end_a, const hashn &mer_list, FILE *fp = stdout) {
	for (; a != end_a; ++a) {
		if (opt_limit_printout && opt_exclude.find(a->name()) == opt_exclude.end()) {
			continue;
		}
		a->vector_start = a->quality_start = 0;
		a->vector_stop = a->quality_stop = a->size();
		screen_repeats(*a, mer_list);
		if (opt_print_percent_masked) {
			size_t x = a->count_masked();
			if (x != 0) {
				fprintf(fp, "%s %5.2f%%\n", a->name().c_str(), (double)100 * x / a->size());
			}
		} else if (opt_print_range) {
			std::vector<std::pair<size_t, size_t> > ranges;
			a->make_mask_ranges(ranges);
			if (!ranges.empty()) {
				fprintf(fp, "%s", a->name().c_str());
				std::vector<std::pair<size_t, size_t> >::const_iterator b(ranges.begin());
				const std::vector<std::pair<size_t, size_t> >::const_iterator end_b(ranges.end());
				for (; b != end_b; ++b) {
					fprintf(fp, " %lu-%lu", b->first, b->second);
				}
				fprintf(fp, "\n");
			}
		} else {
			a->print_sequence(fp);
		}
	}
}

/*
 * return the number represented by s, which may be suffixed by a k, m, or g
 * which act as multipliers to the base amount
 */

static size_t get_value(std::string s) {
	/* validate string - digits optionally followed by k, m, or g */
	size_t i = s.find_first_not_of("0123456789");
	if (i == std::string::npos) {
		return atol(s.c_str());
	} else if (i + 1 == s.length()) {
		size_t x = atol(s.c_str());
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
	} else {	/* bad value */
		return 0;
	}
}

static void read_excludes(std::string s) {
	if (s.find(',') != std::string::npos) {	// comma separated list
		std::vector<std::string> list;
		breakup_line_exact(s, ",", list);
		std::vector<std::string>::const_iterator a = list.begin();
		const std::vector<std::string>::const_iterator end_a = list.end();
		for (; a != end_a; ++a) {
			if (!a->empty()) {
				opt_exclude[*a] = 1;
			}
		}
	} else {			// filename
		int fd = open_compressed(s);
		if (fd == -1) {
			fprintf(stderr, "Error: open_compressed %s: %s\n", s.c_str(), strerror(errno));
			exit(1);
		}
		std::string line;
		while (pfgets(fd, line) != -1) {
			opt_exclude[line] = 1;
		}
		close_compressed(fd);
	}
}

/* print usage statement and exit */

static void print_usage() {
	fprintf(stderr,
		"usage: mask_repeats [options] file1 [file2] ...\n"
		"    -a ## number of phred20's on both sides of a repeat that will keep\n"
		"          it from being masked (defaults to off)\n"
		"    -B ## process seq & qual file in batches of ## reads\n"
		"    -c    clip low quality when counting n-mers\n"
		"    -d    when processing in batches, check for duplicates across whole file\n"
		"    -f ## when clipping quality or vector, use ## as the target quality [20]\n"
		"    -F    print percentage of masked bases for each read\n"
		"          (will not print out reads with no masked bases)\n"
		"    -g    aggregate sequence from all files for determining repeat\n"
		"          counts, print output to individual files\n"
		"    -G    create histogram for each read only from the read itself\n"
		"    -h    print this information\n"
		"    -H ## use this sequence file to create histogram data, instead of\n"
		"          the input files (option may be specified multiple times)\n"
		"    -i    turn off status updates\n"
		"    -k ## when counting n-mers, skip reads smaller than this\n"
		"    -l ## a comma separated list of reads to exclude from the histogram\n"
		"          (if no comma is present, a file of read names used for same)\n"
		"    -L    mask by lowercasing instead of X\n"
		"    -m ## set mer length (defaults to 24)\n"
		"    -p ## don't touch reads not matching pattern (an extended regex)\n"
		"    -q    turn off all warnings\n"
		"    -r    print read:masked_range rather than sequence\n"
		"    -R    reverse mask before masking (does not affect phred20)\n"
		"    -s ## suffix for individual files (defaults to .kmermasked)\n"
		"    -S ## load histogram memory dump from given file\n"
		"    -t ## number of repetitions for a n-mer to be highly repetitive\n"
		"          (defaults to 20)\n"
		"    -T    strip first part of trace id\n"
		"    -u ## (upper limit) number of repetitions for a n-mer to\n"
		"          no longer be highly repetitive\n"
		"    -x ## number of highly repetitive n-mers a base pair needs to\n"
		"          be part of to be masked (defaults to 1)\n"
		"    -X    only print reads given in the -l option\n"
		"    -v    clip vector when counting n-mers\n"
		"    -V    print version\n"
		"    -z ## number of possible n-mers to allocate memory for\n"
		"          (defaults to 200m) (k, m, or g may be suffixed)\n"
		"    -Z    clean hash if it fills up\n"
	);
	exit(1);
}

static void get_opts(int argc, char **argv) {
	opt_aggregate = 0;
	opt_batch_size = 0;
	opt_clip_quality = 0;
	opt_clip_vector = 0;
	opt_feedback = 1;
	opt_hash_clean = 0;
	opt_histogram_restore = -1;
	opt_limit_printout = 0;
	opt_mask_lowercase = 0;
	opt_mer_length = 24;
	opt_nmers = 200 * 1024 * 1024;
	opt_phred20_anchor = -1;
	opt_print_percent_masked = 0;
	opt_print_range = 0;
	opt_quality_cutoff = 20;
	opt_repeat_coverage = 1;
	opt_reverse_mask = 0;
	opt_repeat_threshold = 20;
	opt_repeat_threshold_upper = (hashn::value_type)-1;
	opt_skip_size = 0;
	opt_split = 0;
	opt_strip_tracename = 0;
	opt_suffix = ".kmermasked";
	opt_track_dups = 0;
	opt_warnings = 1;
	int c;
	while ((c = getopt(argc, argv, "a:B:cdf:FgGhH:ik:l:Lm:p:qrRs:S:t:Tu:vVx:Xz:Z")) != EOF) {
		switch (c) {
		    case 'a':
			opt_phred20_anchor = atoi(optarg);
			if (opt_phred20_anchor < 0) {
				fprintf(stderr, "Error: invalid anchor length %d\n", opt_phred20_anchor);
				print_usage();
			}
			break;
		    case 'B':
			c = atoi(optarg);
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
			opt_quality_cutoff = atoi(optarg);
			if (opt_quality_cutoff < 0) {
				print_usage();
			}
			break;
		    case 'F':
			opt_print_percent_masked = 1;
			break;
		    case 'g':
			opt_aggregate = 1;
			break;
		    case 'G':
			opt_split = 1;
			break;
		    case 'h':
			print_usage();
			break;
		    case 'H':
			opt_aggregate = 1;
			hist_files.push_back(optarg);
			break;
		    case 'i':
			opt_feedback = 0;
			break;
		    case 'k':
			c = atoi(optarg);
			if (c < 0) {
				fprintf(stderr, "Error: invalid skip size %d\n", c);
				print_usage();
			}
			opt_skip_size = c;
			break;
		    case 'l':
			read_excludes(optarg);
			break;
		    case 'L':
			opt_mask_lowercase = 1;
			break;
		    case 'm':
			opt_mer_length = atoi(optarg);
			if (opt_mer_length < 1) {
				fprintf(stderr, "Error: invalid mer length %d\n", opt_mer_length);
				print_usage();
			}
			break;
		    case 'p':
			opt_include.initialize(optarg, 0, REG_NOSUB | REG_EXTENDED);
			break;
		    case 'q':
			opt_warnings = 0;
			break;
		    case 'r':
			opt_print_range = 1;
			break;
		    case 'R':
			opt_reverse_mask = 1;
			break;
		    case 's':
			opt_suffix = optarg;
			if (opt_suffix.empty()) {
				fprintf(stderr, "Error: empty file suffix\n");
				print_usage();
			}
			break;
		    case 'S':
			opt_histogram_restore = open_compressed(optarg);
			if (opt_histogram_restore == -1) {
				fprintf(stderr, "Error: could not read histogram dump file\n");
				print_usage();
			}
			opt_aggregate = 1;
			break;
		    case 't':
			c = atoi(optarg);
			if (c < 1) {
				fprintf(stderr, "Error: invalid repeat threshold %d\n", c);
				print_usage();
			}
			opt_repeat_threshold = c;
			break;
		    case 'T':
			opt_strip_tracename = 1;
			break;
		    case 'u':
			c = atoi(optarg);
			if (c < 1) {
				fprintf(stderr, "Error: invalid upper repeat threshold %d\n", c);
				print_usage();
			}
			opt_repeat_threshold_upper = c;
			break;
		    case 'v':
			opt_clip_vector = 1;
			break;
		    case 'V':
			fprintf(stderr, "mask_repeats_hashn version %s%s\n", VERSION,
#ifdef COMPRESS_READS
" (read compression)"
#else
""
#endif
);
			exit(0);
		    case 'x':
			c = atoi(optarg);
			if (c < 1) {
				fprintf(stderr, "Error: invalid repeat coverage %d\n", c);
				print_usage();
			}
			opt_repeat_coverage = c;
			break;
		    case 'X':
			opt_limit_printout = 1;
			break;
		    case 'z':
			opt_nmers = get_value(optarg);
			if (opt_nmers == 0) {
				fprintf(stderr, "Error: bad n-mer count %s\n", optarg);
				print_usage();
			}
			break;
		    case 'Z':
			opt_hash_clean = 1;
			break;
		    default:
			fprintf(stderr, "Error: unknown option %c\n", c);
			print_usage();
		}
	}
	if (optind == argc) {
		fprintf(stderr, "Error: no files specified\n");
		print_usage();
	}
	if (opt_histogram_restore != -1) {
		if (opt_split) {
			fprintf(stderr, "Error: -S and -G options cannot both be specified\n");
			exit(1);
		} else if (!hist_files.empty()) {
			fprintf(stderr, "Error: -S and -H options cannot both be specified\n");
			exit(1);
		} else if (opt_nmers != 200 * 1024 * 1024) {
			fprintf(stderr, "Error: -S and -z options cannot both be specified\n");
			exit(1);
		} else if (opt_hash_clean) {
			fprintf(stderr, "Error: -S and -Z options cannot both be specified\n");
			exit(1);
		}
	}
	if (opt_split && opt_aggregate) {
		if (hist_files.empty()) {
			fprintf(stderr, "Error: -G and -g options cannot both be specified\n");
		} else {
			fprintf(stderr, "Error: -G and -H options cannot both be specified\n");
		}
		exit(1);
	}
	if (opt_print_percent_masked && opt_print_range) {
		fprintf(stderr, "Error: -F and -r options cannot both be specified\n");
		exit(1);
	}
	if (opt_limit_printout && opt_exclude.empty()) {
		fprintf(stderr, "Error: printed reads limited to an empty list - nothing would be printed\n");
		fprintf(stderr, "        Perhaps you forgot to include a -l option?\n");
		exit(1);
	}
	if (opt_repeat_coverage > opt_mer_length) {
		opt_repeat_coverage = opt_mer_length;
		if (opt_warnings) {
			fprintf(stderr, "Warning: reducing repeat coverage to mer length\n");
		}
	}
	if (hist_files.empty() && optind + 1 == argc && opt_histogram_restore == -1) {
		opt_aggregate = 0;
	}
}

int main(int argc, char **argv) {
	get_opts(argc, argv);
	if (opt_feedback) {
		fprintf(stderr, "Initializing n-mer hash\n");
	}
	init_mer_constants(opt_mer_length);
	int err = 0;
	hashn mer_list;
	if (opt_hash_clean) {
		mer_list.set_no_space_response(hashn::CLEAN_HASH);
	}
	if (opt_histogram_restore != -1) {
		mer_list.init_from_file(opt_histogram_restore);
		close_compressed(opt_histogram_restore);
	} else {
		mer_list.init(opt_nmers, opt_mer_length * 2);
	}
	std::list<std::string>::const_iterator a = hist_files.begin();
	const std::list<std::string>::const_iterator end_a = hist_files.end();
	for (; a != end_a; ++a) {
		if (opt_feedback) {
			fprintf(stderr, "Reading in %s\n", a->c_str());
		}
		ReadFile file(*a, opt_batch_size, opt_track_dups);
		if (file.seq_file.empty()) {
			++err;
			continue;
		}
		size_t total_reads(0);
		while (file.read_batch(opt_warnings) != -1) {
			if (!add_sequence_mers(file.read_list.begin(), file.read_list.end(), mer_list, total_reads)) {
				fprintf(stderr, "Error: n-mer list incomplete - specify a larger -z value\n");
				return 1;
			}
			total_reads += file.read_list.size();
		}
	}
	const int start = optind;
	if (hist_files.empty() && opt_histogram_restore == -1) {
		for (; optind < argc; ++optind) {
			if (opt_feedback) {
				fprintf(stderr, "Reading in %s\n", argv[optind]);
			}
			ReadFile file(argv[optind], opt_batch_size, opt_track_dups);
			if (file.seq_file.empty()) {
				++err;
				continue;
			}
			size_t total_reads(0);
			while (file.read_batch(opt_warnings) != -1) {
				if (opt_split) {
					std::list<Read>::iterator c = file.read_list.begin();
					const std::list<Read>::const_iterator end_c = file.read_list.end();
					while (c != end_c) {
						const std::list<Read>::iterator b = c++;
						add_sequence_mers(b, c, mer_list, total_reads);
						print_unique_sequence(b, c, mer_list);
						mer_list.clear();
					}
				} else if (!add_sequence_mers(file.read_list.begin(), file.read_list.end(), mer_list, total_reads)) {
					fprintf(stderr, "Error: n-mer list incomplete - give a larger -z value\n");
					return 1;
				}
				total_reads += file.read_list.size();
			}
			if (!opt_aggregate && !opt_split) {
				if (opt_feedback) {
					fprintf(stderr, "Printing masked sequence\n");
				}
				file.reset();
				while (file.read_batch(opt_warnings) != -1) {
					print_unique_sequence(file.read_list.begin(), file.read_list.end(), mer_list);
				}
				mer_list.clear();
			}
		}
	}
	if (opt_aggregate) {
		if (opt_feedback) {
			fprintf(stderr, "Printing masked sequence\n");
		}
		for (optind = start; optind < argc; ++optind) {
			if (opt_feedback) {
				fprintf(stderr, "Reading in %s\n", argv[optind]);
			}
			ReadFile file(argv[optind], opt_batch_size, opt_track_dups);
			if (file.seq_file.empty()) {
				++err;
				continue;
			}
			FILE *fout = open_output_file(argv[optind]);
			while (file.read_batch(opt_warnings) != -1) {
				print_unique_sequence(file.read_list.begin(), file.read_list.end(), mer_list, fout);
			}
			close_output_file(fout);
		}
	}
	return err;
}
