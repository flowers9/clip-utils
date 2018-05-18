#include "breakup_line.h"	/* break_line_exact() */
#include "hashz.h"	/* hashz */
#include "hist_lib_hashz.h"	/* add_sequence_mers(), init_mer_constants(), opt_exclude, opt_feedback, opt_include, opt_mask_lowercase, opt_phred20_anchor, opt_repeat_coverage, opt_repeat_threshold, opt_repeat_threshold_upper, opt_skip_size, screen_repeats() */
#include "open_compressed.h"	/* close_compressed(), open_compressed(), pfgets() */
#include "read.h"	/* Read, opt_clip_quality, opt_clip_vector, opt_quality_cutoff */
#include "read_lib.h"	/* opt_strip_tracename, read_sequence() */
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

static bool opt_aggregate;
static bool opt_limit_printout;
static bool opt_print_percent_masked;
static bool opt_split;
static bool opt_warnings;
static int opt_mer_length;
static size_t opt_nmers;
static std::list<std::string> hist_files;
static std::string opt_suffix;

/*
 * print full reads, with high repeat regions masked out;
 * output goes to filename + opt_suffix
 */

static void print_unique_sequence(std::list<Read>::iterator a, std::list<Read>::iterator end_a, const hashz &mer_list, std::string filename = "") {
	FILE *fp = stdout;
	if (!filename.empty()) {
		filename += opt_suffix;
		fp = fopen(filename.c_str(), "w");
		if (fp == NULL) {
			fprintf(stderr, "Error: could not write to %s\n", filename.c_str());
			return;
		}
	}
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
		} else {
			a->print_sequence(fp);
		}
	}
	if (fp != stdout) {
		fclose(fp);
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
			x *= 1024;
		    case 'm':
			x *= 1024;
		    case 'k':
			x *= 1024;
			return x;
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
		std::vector<std::string>::const_iterator end_a = list.end();
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
	fprintf(stderr, "usage: mask_repeats [options] file1 [file2] ...\n");
	fprintf(stderr, "    -a phred20's  number of phred20's on both sides of a repeat that will keep\n");
	fprintf(stderr, "                  it from being masked (defaults to off)\n");
	fprintf(stderr, "    -c            clip low quality when counting n-mers\n");
	fprintf(stderr, "    -f ## when clipping quality or vector, use ## as the target quality [20]\n");
	fprintf(stderr, "    -F            print percentage of masked bases for each read\n");
	fprintf(stderr, "    -g            aggregate sequence from all files for determining repeat\n");
	fprintf(stderr, "                  counts, print output to individual files\n");
	fprintf(stderr, "    -G            create histogram for each read only from the read itself\n");
	fprintf(stderr, "    -h            print this information\n");
	fprintf(stderr, "    -H ##         use this sequence file to create histogram data, instead of\n");
	fprintf(stderr, "                  the input files (option may be specified multiple times)\n");
	fprintf(stderr, "    -i            turn off status updates\n");
	fprintf(stderr, "    -k ##         when counting n-mers, skip reads smaller than this\n");
	fprintf(stderr, "    -l ##         a comma separated list of reads to exclude from the histogram\n");
	fprintf(stderr, "                  (if no comma is present, a file of read names used for same)\n");
	fprintf(stderr, "    -L            mask by lowercasing instead of X\n");
	fprintf(stderr, "    -m mer        set mer length (defaults to 24)\n");
	fprintf(stderr, "    -p pattern    don't touch reads not matching pattern (an extended regex)\n");
	fprintf(stderr, "    -q            turn off all warnings\n");
	fprintf(stderr, "    -s suffix     suffix for individual files (defaults to .kmermasked)\n");
	fprintf(stderr, "    -t threshold  number of repetitions for a n-mer to be highly repetitive\n");
	fprintf(stderr, "                  (defaults to 20)\n");
	fprintf(stderr, "    -T            strip first part of trace id\n");
	fprintf(stderr, "    -u threshold  (upper limit) number of repetitions for a n-mer to\n");
	fprintf(stderr, "                  no longer be highly repetitive\n");
	fprintf(stderr, "    -x threshold  number of highly repetitive n-mers a base pair needs to\n");
        fprintf(stderr, "                  be part of to be masked (defaults to 1)\n");
	fprintf(stderr, "    -X            only print reads given in the -l option\n");
	fprintf(stderr, "    -v            clip vector when counting n-mers\n");
	fprintf(stderr, "    -V            print version\n");
	fprintf(stderr, "    -z n-mers     number of possible n-mers to allocate memory for\n");
	fprintf(stderr, "                  (defaults to 200m) (k, m, or g may be suffixed)\n");
	exit(1);
}

static void get_opts(int argc, char **argv) {
	opt_aggregate = 0;
	opt_clip_quality = 0;
	opt_clip_vector = 0;
	opt_feedback = 1;
	opt_limit_printout = 0;
	opt_mask_lowercase = 0;
	opt_mer_length = 24;
	opt_nmers = 200 * 1024 * 1024;
	opt_phred20_anchor = -1;
	opt_print_percent_masked = 0;
	opt_quality_cutoff = 20;
	opt_repeat_coverage = 1;
	opt_repeat_threshold = 20;
	opt_repeat_threshold_upper = (hashz::value_type)-1;
	opt_skip_size = 0;
	opt_split = 0;
	opt_strip_tracename = 0;
	opt_suffix = ".kmermasked";
	opt_warnings = 1;
	int c, i;
	while ((c = getopt(argc, argv, "a:cf:FgGhH:ik:l:Lm:p:qs:t:Tu:vVx:Xz:")) != EOF) {
		switch (c) {
		    case 'a':
			opt_phred20_anchor = atoi(optarg);
			if (opt_phred20_anchor < 0) {
				fprintf(stderr, "Error: invalid anchor length %d\n", opt_phred20_anchor);
				print_usage();
			}
			break;
		    case 'c':
			opt_clip_quality = 1;
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
			/* use an int here to catch negative values */
			i = atoi(optarg);
			if (i < 0) {
				fprintf(stderr, "Error: invalid skip size %d\n", i);
				print_usage();
			}
			opt_skip_size = i;
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
		    case 's':
			opt_suffix = optarg;
			if (opt_suffix.empty()) {
				fprintf(stderr, "Error: empty file suffix\n");
				print_usage();
			}
			break;
		    case 't':
			/* use an int here to catch negative values */
			i = atoi(optarg);
			if (i < 1) {
				fprintf(stderr, "Error: invalid repeat threshold %d\n", i);
				print_usage();
			}
			opt_repeat_threshold = i;
			break;
		    case 'T':
			opt_strip_tracename = 1;
			break;
		    case 'u':
			/* use an int here to catch negative values */
			i = atoi(optarg);
			if (i < 1) {
				fprintf(stderr, "Error: invalid upper repeat threshold %d\n", i);
				print_usage();
			}
			opt_repeat_threshold_upper = i;
			break;
		    case 'v':
			opt_clip_vector = 1;
			break;
		    case 'V':
			fprintf(stderr, "mask_repeats_hashz version %s%s\n", VERSION,
#ifdef COMPRESS_READS
" (read compression)"
#else
""
#endif
);
			exit(0);
		    case 'x':
			opt_repeat_coverage = atoi(optarg);
			if (opt_repeat_coverage < 1) {
				fprintf(stderr, "Error: invalid repeat coverage %d\n", opt_repeat_coverage);
				print_usage();
			}
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
		    default:
			fprintf(stderr, "Error: unknown option %c\n", c);
			print_usage();
		}
	}
	if (optind == argc) {
		fprintf(stderr, "Error: no files specified\n");
		print_usage();
	}
	if (opt_split && opt_aggregate) {
		if (hist_files.empty()) {
			fprintf(stderr, "Error: -G and -g options cannot both be specified\n");
		} else {
			fprintf(stderr, "Error: -G and -H options cannot both be specified\n");
		}
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
	if (hist_files.empty() && optind + 1 == argc) {
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
	int start = optind;
		// if first read file is also used for histogram, only read once
	int match_start = 0;
	hashz mer_list(opt_nmers, opt_mer_length * 2);
	std::list<std::string>::const_iterator a = hist_files.begin();
	std::list<std::string>::const_iterator end_a = hist_files.end();
	for (; a != end_a; ++a) {
		if (a->compare(argv[optind]) == 0) {
			match_start = 1;
			continue;
		}
		std::list<Read> read_list;
		if (opt_feedback) {
			fprintf(stderr, "Reading in %s\n", a->c_str());
		}
		if (read_sequence(a->c_str(), read_list, opt_warnings) == -1) {
			++err;
			continue;
		}
		if (opt_feedback) {
			fprintf(stderr, "Adding n-mers\n");
		}
		if (!add_sequence_mers(read_list.begin(), read_list.end(), mer_list)) {
			fprintf(stderr, "Error: n-mer list incomplete - specify a larger -z value\n");
			return 1;
		}
	}
	if (hist_files.empty() || match_start) {
		for (; optind < argc; ++optind) {
			std::list<Read> read_list;
			if (opt_feedback) {
				fprintf(stderr, "Reading in %s\n", argv[optind]);
			}
			if (read_sequence(argv[optind], read_list, opt_warnings) == -1) {
				++err;
				continue;
			}
			if (opt_feedback) {
				fprintf(stderr, "Adding n-mers\n");
			}
			if (opt_split) {
				std::list<Read>::iterator c = read_list.begin();
				std::list<Read>::iterator end_c = read_list.end();
				std::list<Read>::iterator b = c;
				for (++b; c != end_c; ++c, ++b) {
					add_sequence_mers(c, b, mer_list);
					print_unique_sequence(c, b, mer_list);
					mer_list.clear();
				}
				continue;
			}
			if (!add_sequence_mers(read_list.begin(), read_list.end(), mer_list)) {
				fprintf(stderr, "Error: n-mer list incomplete - give a larger -z value\n");
				return 1;
			}
			if (!opt_aggregate || match_start) {
				if (opt_feedback) {
					fprintf(stderr, "Printing masked sequence\n");
				}
				if (match_start) {
					print_unique_sequence(read_list.begin(), read_list.end(), mer_list, argv[optind]);
					break;
				} else {
					print_unique_sequence(read_list.begin(), read_list.end(), mer_list);
				}
				mer_list.clear();
			}
		}
	}
	if (opt_aggregate) {
		for (optind = start; optind < argc; ++optind) {
			if (opt_feedback) {
				fprintf(stderr, "Reading in %s\n", argv[optind]);
			}
			std::list<Read> read_list;
			if (read_sequence(argv[optind], read_list, opt_warnings) != -1) {
				if (opt_feedback) {
					fprintf(stderr, "Printing masked sequence for %s\n", argv[optind]);
				}
				print_unique_sequence(read_list.begin(), read_list.end(), mer_list, argv[optind]);
			}
		}
	}
	return err;
}
