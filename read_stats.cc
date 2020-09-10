#include "hash.h"	/* hash */
#include "hist_lib_hash.h"	/* add_sequence_mers(), clear_mer_list(), count_kmers(), init_mer_constants(), opt_feedback, opt_include, opt_mer_length, print_final_input_feedback() */
#include "read.h"	/* Read, opt_clip_quality, opt_clip_vector, opt_quality_cutoff */
#include "read_lib.h"	/* opt_strip_tracename, read_sequence(), set_default_endpoints() */
#include <getopt.h>	// getopt(), optarg, optind
#include <limits.h>	/* UCHAR_MAX */
#include <list>		/* list<> */
#include <regex.h>	/* REG_EXTENDED, REG_NOSUB */
#include <stdio.h>	/* EOF, fprintf(), printf(), putchar(), stderr */
#include <stdlib.h>	// atoi(), atol(), exit()
#include <string.h>	/* strlen() */
#include <string>	/* string */
#include <sys/types.h>	/* size_t */
#include <time.h>	/* time() */

static bool opt_aggregate;
static bool opt_warnings;
static size_t opt_nmers;

static void print_read_stats(const std::list<Read> &read_list, const hash &mer_list) {
	std::list<Read>::const_iterator a = read_list.begin();
	std::list<Read>::const_iterator end = read_list.end();
	for (; a != end; ++a) {
		size_t kmers, r_kmers, ur_kmers;
		count_kmers(*a, mer_list, kmers, r_kmers, ur_kmers);
		printf("%s %6ld", a->name().c_str(), a->quality_stop - a->quality_start);
		if (kmers == 0) {
			printf("\n");
		} else if (r_kmers == 0) {
			printf("   -0-     -0-\n");
		} else {
			printf(" %6.2f%% %6.2f%%\n", (double)100 * r_kmers / kmers, (double)100 * ur_kmers / r_kmers);
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

/* print usage statement and exit */

static void print_usage() {
	fprintf(stderr, "usage: read_stats [options] file1 [file2] ...\n");
	fprintf(stderr, "    -c     clip low quality\n");
	fprintf(stderr, "    -f ##  when clipping quality or vector, use ## as the target quality [20]\n");
	fprintf(stderr, "    -g     aggregate sequence from all files for determining repeat\n");
	fprintf(stderr, "    -h     print this information\n");
	fprintf(stderr, "    -i     turn off status updates\n");
	fprintf(stderr, "    -m ##  mer length (1-32) [24]\n");
	fprintf(stderr, "    -q     turn off all warnings\n");
	fprintf(stderr, "    -t ##  repeat threshold [20]\n");
	fprintf(stderr, "    -T     don't strip first part of trace id\n");
	fprintf(stderr, "    -v     clip vector\n");
	fprintf(stderr, "    -z ##  number of possible n-mers to allocate memory for [200m]\n");
	fprintf(stderr, "           (k, m, or g may be suffixed)\n");
	exit(1);
}

static void get_opts(int argc, char **argv) {
	opt_aggregate = 0;
	opt_clip_quality = 0;
	opt_clip_vector = 0;
	opt_feedback = 1;
	opt_mer_length = 24;
	opt_nmers = 200 * 1024 * 1024;
	opt_quality_cutoff = 20;
	opt_repeat_threshold = 20;
	opt_strip_tracename = 1;
	opt_warnings = 1;
	int c, i;
	while ((c = getopt(argc, argv, "cf:ghim:p:qt:Tvz:")) != EOF) {
		switch (c) {
		    case 'c':
			opt_clip_quality = 1;
			break;
		    case 'f':
			opt_quality_cutoff = atoi(optarg);
			if (opt_quality_cutoff < 0) {
				print_usage();
			}
			break;
		    case 'g':
			opt_aggregate = 1;
			break;
		    case 'h':
			print_usage();
			break;
		    case 'i':
			opt_feedback = 0;
			break;
		    case 'm':
			opt_mer_length = atoi(optarg);
			if (opt_mer_length < 1 || 32 < opt_mer_length) {
				fprintf(stderr, "Error: bad mer length\n");
				print_usage();
			}
			break;
		    case 'q':
			opt_warnings = 0;
			break;
		    case 't':
			i = atoi(optarg);
			if (i < 1 || i > UCHAR_MAX) {
				print_usage();
			}
			opt_repeat_threshold = i;
			break;
		    case 'T':
			opt_strip_tracename = 0;
			break;
		    case 'v':
			opt_clip_vector = 1;
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
		fprintf(stderr, "Error: no files to process\n");
		print_usage();
	}
	if (optind + 1 == argc) {
		opt_aggregate = 0;
	}
}

int main(int argc, char **argv) {
	get_opts(argc, argv);
	if (opt_feedback) {
		fprintf(stderr, "Initializing n-mer hash\n");
	}
	init_mer_constants();
	int err = 0;
	int start = optind;
	hash mer_list(opt_nmers);
	for (; optind != argc; ++optind) {
		if (opt_feedback) {
			fprintf(stderr, "Reading in %s\n", argv[optind]);
		}
		std::list<Read> read_list;
		if (read_sequence(argv[optind], read_list, opt_warnings) == -1) {
			++err;
			continue;
		}
		if (opt_feedback) {
			fprintf(stderr, "Adding n-mers\n");
		}
		if (!add_sequence_mers(read_list.begin(), read_list.end(), mer_list, 0)) {
			fprintf(stderr, "Error: n-mer list incomplete - give a larger -z value\n");
		}
		if (!opt_aggregate) {
			if (opt_feedback) {
				fprintf(stderr, "Printing read stats\n");
			}
			print_read_stats(read_list, mer_list);
			clear_mer_list(mer_list);
		}
	}
	if (opt_aggregate) {
		print_final_input_feedback(mer_list);
		for (optind = start; optind != argc; ++optind) {
			if (opt_feedback) {
				fprintf(stderr, "Rereading %s\n", argv[optind]);
			}
			std::list<Read> read_list;
			if (read_sequence(argv[optind], read_list, opt_warnings) != -1) {
				if (opt_feedback) {
					fprintf(stderr, "Printing read stats for %s\n", argv[optind]);
				}
				print_read_stats(read_list, mer_list);
			}
		}
	}
	return err;
}
