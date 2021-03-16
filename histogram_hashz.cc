#include "hashz.h"	/* hashz */
#include "hist_lib_hashz.h"	/* add_sequence_mers(), convert_key(), init_mer_constants(), opt_feedback, opt_include, opt_skip_size, reverse_key() */
#include "open_compressed.h"	/* close_compressed(), open_compressed(), pfgets() */
#include "read.h"	/* Read, opt_clip_quality, opt_clip_vector, opt_quality_cutoff */
#include "read_lib.h"	/* opt_strip_tracename, read_sequence(), set_default_endpoints() */
#include "strtostr.h"	/* strtostr() */
#include "version.h"	/* VERSION */
#include <getopt.h>	// getopt(), optarg, optind
#include <gmp.h>	/* mpz_clear(), mpz_cmp(), mpz_init2() */
#include <list>		/* list<> */
#include <map>		/* map<> */
#include <regex.h>	/* REG_EXTENDED, REG_NOSUB */
#include <stdio.h>	/* EOF, fprintf(), putc(), stderr, stdout */
#include <stdlib.h>	// abs(), atoi(), atol(), exit()
#include <string.h>	/* strlen() */
#include <string>	/* string */
#include <sys/types.h>	/* size_t */
#include <time.h>	/* time() */

static FILE *fp_out = stdout;
static bool opt_aggregate;
static bool opt_warnings;
static int opt_readnames_exclude;
static size_t opt_nmers;
static std::map<std::string, hashz::offset_type> opt_readnames;
static unsigned long opt_frequency_cutoff;
static unsigned long opt_mer_length;

/* print n-mer occurence frequency */

static void print_mer_frequency(const hashz &mer_list) {
	hashz::key_type comp_key;
	mpz_init2(comp_key, opt_mer_length * 2);
	hashz::const_iterator a = mer_list.begin();
	hashz::const_iterator end = mer_list.end();
	for (; a != end; ++a) {
		if (a.value >= opt_frequency_cutoff) {
			fprintf(fp_out, "%s %lu\n", convert_key(a.key).c_str(), a.value);
			reverse_key(a.key, comp_key);
			if (mpz_cmp(a.key, comp_key) != 0) {
				fprintf(fp_out, "%s %lu\n", convert_key(comp_key).c_str(), a.value);
			}
		}
	}
	mpz_clear(comp_key);
}

/* print histogram of n-mer occurrences */

static void print_mer_histogram(const hashz &mer_list) {
	std::map<unsigned long, unsigned long> counts;
	hashz::key_type comp_key;
	mpz_init2(comp_key, opt_mer_length * 2);
	hashz::const_iterator a = mer_list.begin();
	hashz::const_iterator end_a = mer_list.end();
	for (; a != end_a; ++a) {
		reverse_key(a.key, comp_key);
		if (mpz_cmp(a.key, comp_key) == 0) {
			counts[a.value] += 2;
		} else {
			++counts[a.value];
		}
	}
	mpz_clear(comp_key);
	std::map<unsigned long, unsigned long>::const_iterator c = counts.begin();
	std::map<unsigned long, unsigned long>::const_iterator end_c = counts.end();
	double total = 0;
	for (; c != end_c; ++c) {
		if (c->first > 1) {
			total += c->first * c->second;
		}
	}
	double i = 0;
	for (c = counts.begin(); c != end_c; ++c) {
		if (c->first > 1) {
			double x = 100 * c->first * c->second;
			i += x;
			fprintf(fp_out, "%lu %lu %.2f %.2f\n", c->first, c->second, x / total, i / total);
		} else {
			fprintf(fp_out, "%lu %lu\n", c->first, c->second);
		}
	}
}

static void print_mer_histogram_sub(const hashz &mer_list) {
	std::map<hashz::value_type, unsigned long> counts[opt_readnames_exclude];
	hashz::const_iterator a = mer_list.begin();
	hashz::const_iterator end_a = mer_list.end();
	hashz::value_type x[opt_readnames_exclude];
	for (; a != end_a; ++a) {
		a.get_alt_values(x);
		int i;
		hashz::value_type n = a.value;
		for (i = 0; i != opt_readnames_exclude; ++i) {
			n += x[i];
		}
		if (n != x[0]) {
			hashz::value_type m = n;
			for (i = 0; i != opt_readnames_exclude; ++i) {
				m -= x[i];
				counts[i][n] += m;
			}
		}
	}
	int i;
	for (i = 0; i != opt_readnames_exclude; ++i) {
		fprintf(fp_out, "\n");
		std::map<hashz::value_type, unsigned long>::const_iterator c = counts[i].begin();
		std::map<hashz::value_type, unsigned long>::const_iterator end_c = counts[i].end();
		for (; c != end_c; ++c) {
			fprintf(fp_out, "%lu %lu\n", c->first, c->second);
		}
	}
}

static void print_mer_histogram_add(const hashz &mer_list) {
	opt_readnames_exclude *= -1;
	std::map<hashz::value_type, unsigned long> counts[opt_readnames_exclude];
	hashz::const_iterator a = mer_list.begin();
	hashz::const_iterator end_a = mer_list.end();
	hashz::value_type x[opt_readnames_exclude];
	for (; a != end_a; ++a) {
		a.get_alt_values(x);
		int i;
		for (i = 0; i != opt_readnames_exclude; ++i) {
			if (x[i] != 0) {
				counts[i][a.value] += x[i];
			}
		}
	}
	int i;
	for (i = 0; i != opt_readnames_exclude; ++i) {
		fprintf(fp_out, "\n");
		std::map<hashz::value_type, unsigned long>::const_iterator c = counts[i].begin();
		std::map<hashz::value_type, unsigned long>::const_iterator end_c = counts[i].end();
		for (; c != end_c; ++c) {
			fprintf(fp_out, "%lu %lu\n", c->first, c->second);
		}
	}
	opt_readnames_exclude *= -1;
}

/* read in read names from the given file, and add them to list */

static void add_readnames(const char *filename, std::map<std::string, hashz::offset_type> &list) {
	int fd = open_compressed(filename);
	if (fd == -1) {
		fprintf(stderr, "Error: could not read %s\n", filename);
		return;
	}
	hashz::offset_type x = 1 << (abs(opt_readnames_exclude) - 1);
	std::string line;
	while (pfgets(fd, line) != -1) {
		std::string s = strtostr(line);
		if (!s.empty()) {
			if (opt_readnames_exclude < 0) {
				list[s] |= x;
			} else if (list[s] == 0) {
				list[s] = x;
			}
		}
	}
	close_compressed(fd);
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
	fprintf(stderr, "usage: histogram [options] file1 [file2] ...\n");
	fprintf(stderr, "    -a         give combined results for all files\n");
	fprintf(stderr, "    -c         clip low quality\n");
	fprintf(stderr, "    -f ## when clipping quality or vector, use ## as the target quality [20]\n");
	fprintf(stderr, "    -h         print this information\n");
	fprintf(stderr, "    -i         turn off status updates\n");
	fprintf(stderr, "    -k ##      skip reads smaller than this\n");
	fprintf(stderr, "    -l ##      filename containing names of reads to subtract from results\n");
	fprintf(stderr, "               (histogram is given as count*frequency, rather than count)\n");
	fprintf(stderr, "    -L ##      filename containing names of reads to compare with results\n");
	fprintf(stderr, "               (count is by given reads, frequency is by other reads)\n");
	fprintf(stderr, "    -m mer     set mer length (1-32, defaults to 24)\n");
	fprintf(stderr, "    -o outputfile  print output to file instead of stdout\n");
	fprintf(stderr, "    -p pattern don't touch reads not matching pattern (an extended regex)\n");
	fprintf(stderr, "    -q         turn off all warnings\n");
	fprintf(stderr, "    -t         strip first part of trace id\n");
	fprintf(stderr, "    -v         clip vector\n");
	fprintf(stderr, "    -V         print version\n");
	fprintf(stderr, "    -w cutoff  print frequency count instead of histogram, for all n-mers with\n");
	fprintf(stderr, "    -z n-mers  number of possible n-mers to allocate memory for\n");
	fprintf(stderr, "               (k, m, or g may be suffixed; defaults to 200m)\n");
	exit(1);
}

static void get_opts(int argc, char **argv) {
	std::string opt_output;
	opt_aggregate = 0;
	opt_clip_quality = 0;
	opt_clip_vector = 0;
	opt_feedback = 1;
	opt_frequency_cutoff = 0;
	opt_mer_length = 24;
	opt_nmers = 200 * 1024 * 1024;
	opt_quality_cutoff = 20;
	opt_readnames_exclude = 0;
	opt_skip_size = 0;
	opt_strip_tracename = 0;
	opt_warnings = 1;
	int i, c;
	while ((c = getopt(argc, argv, "acf:hik:l:L:m:o:p:qtvVw:z:")) != EOF) {
		switch (c) {
		    case 'a':
			opt_aggregate = 1;
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
		    case 'h':
			print_usage();
			break;
		    case 'i':
			opt_feedback = 0;
			break;
		    case 'k':
			/* use an int here to capture negative values */
			i = atoi(optarg);
			if (i < 0) {
				fprintf(stderr, "Error: invalid skip size %d\n", i);
				print_usage();
			}
			opt_skip_size = i;
			break;
		    case 'l':
			if (opt_readnames_exclude < 0) {
				fprintf(stderr, "Warning: -l and -L options conflict: ignoring -l option\n");
			} else {
				++opt_readnames_exclude;
				add_readnames(optarg, opt_readnames);
			}
			break;
		    case 'L':
			if (opt_readnames_exclude > 0) {
				fprintf(stderr, "Warning: -l and -L options conflict: ignoring -L option\n");
			} else {
				--opt_readnames_exclude;
				add_readnames(optarg, opt_readnames);
			}
			break;
		    case 'm':
			opt_mer_length = atoi(optarg);
			if (opt_mer_length < 1) {
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
			fprintf(stderr, "histogram_hashz version %s%s\n", VERSION,
#ifdef COMPRESS_READS
" (read compression)"
#else
""
#endif
);
			exit(0);
		    case 'w':
			opt_frequency_cutoff = atol(optarg);
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
	if (opt_frequency_cutoff != 0 && opt_readnames_exclude) {
		fprintf(stderr, "Warning: -w and -l/-L options conflict: ignoring -w option\n");
	}
	if (optind == argc) {
		fprintf(stderr, "Error: no files to process\n");
		print_usage();
	}
	if (optind + 1 == argc) {
		opt_aggregate = 1;
	}
	if (!opt_output.empty()) {
		fp_out = fopen(opt_output.c_str(), "w");
		if (fp_out == NULL) {
			fprintf(stderr, "Error: could not write to %s\n", opt_output.c_str());
			fp_out = stdout;
		}
	}
}

int main(int argc, char **argv) {
	get_opts(argc, argv);
	if (opt_feedback) {
		fprintf(stderr, "Initializing n-mer hash\n");
	}
	init_mer_constants(opt_mer_length);
	int err = 0;
	hashz mer_list(opt_nmers, opt_mer_length * 2, abs(opt_readnames_exclude));
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
		if (opt_readnames_exclude) {
			if (!add_sequence_mers(read_list.begin(), read_list.end(), mer_list, opt_readnames)) {
				fprintf(stderr, "Error: n-mer list incomplete - give a larger -z value\n");
			}
		} else if (!add_sequence_mers(read_list.begin(), read_list.end(), mer_list)) {
			fprintf(stderr, "Error: n-mer list incomplete - give a larger -z value\n");
		}
		if (!opt_aggregate) {
			if (opt_feedback) {
				fprintf(stderr, "Printing histogram\n");
			}
			fprintf(fp_out, "%s\n", argv[optind]);
			int i = strlen(argv[optind]);
			for (; i != 0; --i) {
				putc('-', fp_out);
			}
			fprintf(fp_out, "\n");
			if (opt_readnames_exclude > 0) {
				print_mer_histogram_sub(mer_list);
			} else if (opt_readnames_exclude < 0) {
				print_mer_histogram_add(mer_list);
			} else if (opt_frequency_cutoff == 0) {
				print_mer_histogram(mer_list);
			} else {
				print_mer_frequency(mer_list);
			}
			if (optind + 1 != argc) {
				fprintf(fp_out, "\n");
			}
			mer_list.clear();
		}
	}
	if (opt_aggregate) {
		if (opt_feedback) {
			fprintf(stderr, "Printing histogram\n");
		}
		if (opt_readnames_exclude > 0) {
			print_mer_histogram_sub(mer_list);
		} else if (opt_readnames_exclude < 0) {
			print_mer_histogram_add(mer_list);
		} else if (opt_frequency_cutoff == 0) {
			print_mer_histogram(mer_list);
		} else {
			print_mer_frequency(mer_list);
		}
	}
	if (fp_out != stdout) {
		fclose(fp_out);
	}
	return err;
}
