#include "hash.h"	/* hash */
#include "hist_lib_hash.h"	/* add_sequence_mers(), clear_mer_lst(), count_kmers(), init_mer_constants(), opt_feedback, opt_include, opt_mer_length, opt_repeat_threshold, print_final_input_feedback() */
#include "read.h"	/* Read, opt_clip_quality, opt_clip_vector, opt_quality_cutoff */
#include "read_lib.h"	/* opt_strip_tracename, read_sequence(), set_default_endpoints() */
#include <getopt.h>	// getopt(), optarg, optind
#include <limits.h>	/* UCHAR_MAX */
#include <list>		/* list<> */
#include <map>		/* map<> */
#include <regex.h>	/* REG_EXTENDED, REG_NOSUB */
#include <stdio.h>	/* EOF, fprintf(), printf(), putchar(), stderr */
#include <stdlib.h>	// atoi(), atol(), exit()
#include <string.h>	/* strlen() */
#include <string>	/* string */
#include <sys/types.h>	/* size_t */
#include <time.h>	/* time() */
#include <vector>	/* vector<> */

static bool opt_aggregate;
static bool opt_warnings;
static size_t opt_nmers;
static unsigned int opt_phred_count_cutoff;

/*
static void print_gnuplot_style(size_t total, size_t skipped, const hash::value_type repeat_thresholds[], size_t thresholds, const size_t zero_reads[], const size_t table1[][20], const size_t table2[][20]) {
	printf("## Total reads: %lu\n", total);
	printf("## Skipped reads: %lu\n", skipped);
	size_t i;
	for (i = 0; i != thresholds; ++i) {
		printf("# threshold %lu\n", repeat_thresholds[i]);
		printf("# reads with zero repetitions %lu\n", zero_reads[i]);
		int j;
		for (j = 0; j != 20; ++j) {
			printf("%lu %d %lu\n", i, j, table1[i][j]);
		}
		printf("\n");
	}
	printf("\n");	// to separate the two indexes (table1 and table2)
	for (i = 0; i != thresholds; ++i) {
		printf("# threshold %lu\n", repeat_thresholds[i]);
		printf("# reads with zero repetitions %lu\n", zero_reads[i]);
		int j;
		for (j = 0; j != 20; ++j) {
			printf("%lu %d %lu\n", i, j, table2[i][j]);
		}
		printf("\n");
	}
}
*/

static void print_matrix_style(size_t total, size_t skipped, size_t thresholds, const size_t table1[][20], const size_t table2[][20]) {
	printf("## Total reads: %lu\n", total);
	printf("## Skipped reads: %lu\n", skipped);
	size_t i;
	for (i = 0; i != thresholds; ++i) {
		int j;
		for (j = 0; j != 19; ++j) {
			printf("%lu\t", table1[i][j]);
		}
		printf("%lu\n", table1[i][j]);
	}
	printf("\n");	// to separate the two tables
	for (i = 0; i != thresholds; ++i) {
		int j;
		for (j = 0; j != 19; ++j) {
			printf("%lu\t", table2[i][j]);
		}
		printf("%lu\n", table2[i][j]);
	}
}

static void print_read_hist(const std::list<Read> &read_list, const hash &mer_list) {
	const hash::value_type repeat_thresholds[] = { 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100 };
	const size_t thresholds = sizeof(repeat_thresholds) / sizeof(hash::value_type);
	std::vector<const Read *> reads;
	reads.reserve(read_list.size());
	std::list<Read>::const_iterator a = read_list.begin();
	std::list<Read>::const_iterator end_a = read_list.end();
	for (; a != end_a; ++a) {
		if (a->phred_count >= opt_phred_count_cutoff) {
			reads.push_back(&(*a));
		}
	}
	if (opt_feedback) {
		fprintf(stderr, "Total reads: %lu\n", read_list.size());
		fprintf(stderr, "Skipped reads: %lu\n", read_list.size() - reads.size());
	}
	size_t zero_reads[thresholds];
	size_t table1[thresholds][20];
	size_t table2[thresholds][20];
	size_t i;
	// initialize tables
	for (i = 0; i != thresholds; ++i) {
		zero_reads[i] = 0;
		int j;
		for (j = 0; j != 20; ++j) {
			table1[i][j] = 0;
			table2[i][j] = 0;
		}
	}
	// loop over reads in outer loop to allow more efficient swapping
	std::vector<const Read *>::const_iterator b = reads.begin();
	std::vector<const Read *>::const_iterator end_b = reads.end();
	for (; b != end_b; ++b) {
		const Read &c = **b;
		for (i = 0; i != thresholds; ++i) {
			opt_repeat_threshold = repeat_thresholds[i];
			size_t kmers, r_kmers, ur_kmers;
			count_kmers(c, mer_list, kmers, r_kmers, ur_kmers);
			if (r_kmers == 0) {
				++zero_reads[i];
			} else {
				++table1[i][(20 * r_kmers - 1) / kmers];
				++table2[i][(20 * ur_kmers - 1) / r_kmers];
			}
		}
	}
	//print_gnuplot_style(read_list.size(), read_list.size() - reads.size(), repeat_thresholds, thresholds, zero_reads, table1, table2);
	print_matrix_style(read_list.size(), read_list.size() - reads.size(), thresholds, table1, table2);
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

/* print usage statement and exit */

static void print_usage() {
	fprintf(stderr, "usage: read_histogram [options] file1 [file2] ...\n");
	fprintf(stderr, "    -c     clip low quality\n");
	fprintf(stderr, "    -f ##  when clipping quality or vector, use ## as the target quality [20]\n");
	fprintf(stderr, "    -g     aggregate sequence from all files for determining repeat\n");
	fprintf(stderr, "    -h     print this information\n");
	fprintf(stderr, "    -i     turn off status updates\n");
	fprintf(stderr, "    -m ##  mer length (1-32) [24]\n");
	fprintf(stderr, "    -p ##  ignore reads with less than ## qualities of 20 or more\n");
	fprintf(stderr, "    -q     turn off all warnings\n");
	fprintf(stderr, "    -v     clip vector\n");
	fprintf(stderr, "    -w     don't strip first part of trace id\n");
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
	opt_phred_count_cutoff = 0;
	opt_quality_cutoff = 20;
	opt_strip_tracename = 1;
	opt_warnings = 1;
	int c, i;
	while ((c = getopt(argc, argv, "cf:ghim:p:qvwz:")) != EOF) {
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
			i = atoi(optarg);
			if (i < 1 || 32 < i) {
				fprintf(stderr, "Error: bad mer length: %s\n", optarg);
				print_usage();
			}
			opt_mer_length = i;
			break;
		    case 'p':
			i = atoi(optarg);
			if (i < 0) {
				print_usage();
			}
			opt_phred_count_cutoff = i;
			break;
		    case 'q':
			opt_warnings = 0;
			break;
		    case 'v':
			opt_clip_vector = 1;
			break;
		    case 'w':
			opt_strip_tracename = 0;
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
				fprintf(stderr, "Printing read histogram\n");
			}
			print_read_hist(read_list, mer_list);
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
					fprintf(stderr, "Printing read histogram for %s\n", argv[optind]);
				}
				print_read_hist(read_list, mer_list);
			}
		}
	}
	return err;
}
