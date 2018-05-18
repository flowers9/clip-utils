#include "hashn.h"	// hashn
#include "hist_lib_hashn.h"	// add_sequence_mers(), convert_key(), init_mer_constants(), opt_feedback, opt_include, opt_skip_size, reverse_key()
#include "open_compressed.h"	// close_compressed(), get_suffix(), open_compressed(), pfgets()
#include "read.h"	// Read, opt_clip_quality, opt_clip_vector, opt_quality_cutoff
#include "read_file.h"	// ReadFile, opt_strip_tracename
#include "strtostr.h"	// strtostr()
#include "version.h"	// VERSION
#include "write_fork.h"	// write_fork()
#include <getopt.h>	// getopt(), optarg, optind
#include <list>		// list<>
#include <map>		// map<>
#include <regex.h>	// REG_EXTENDED, REG_NOSUB
#include <stdio.h>	// EOF, fprintf(), putc(), stderr, stdout
#include <stdlib.h>	// abs(), atoi(), atol(), exit()
#include <string.h>	// strlen()
#include <string>	// string
#include <sys/types.h>	// size_t

static FILE *fp_out(stdout);
static bool opt_aggregate;
static bool opt_hash_clean;
static bool opt_print_gc;
static bool opt_track_dups;
static bool opt_warnings;
static int opt_histogram_restore;
static int opt_mer_length;
static int opt_readnames_exclude;
static size_t opt_batch_size;
static size_t opt_nmers;
static std::map<std::string, hashn::offset_type> opt_readnames;
static std::string opt_save_file;
static std::string opt_tmp_file_prefix;
static unsigned long opt_frequency_cutoff;

static void save_memory(const hashn &mer_list) {
	std::string suffix;
	get_suffix(opt_save_file, suffix);
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
	const int fd(write_fork(args, opt_save_file));
	if (fd == -1) {
		fprintf(stderr, "Error: could not save memory\n");
		exit(1);
	}
	mer_list.save(fd);
	close_fork(fd);
}

// print n-mer occurence frequency

static void print_mer_frequency(hashn &mer_list) {
	hashn::key_type comp_key(mer_list);
	hashn::const_iterator a(mer_list.begin());
	const hashn::const_iterator end_a(mer_list.end());
	for (; a != end_a; ++a) {
		if (a.value >= opt_frequency_cutoff) {
			fprintf(fp_out, "%s %lu\n", convert_key(a.key).c_str(), a.value);
			reverse_key(a.key, comp_key);
			if (a.key != comp_key) {
				fprintf(fp_out, "%s %lu\n", convert_key(comp_key).c_str(), a.value);
			}
		}
	}
}

// number of gc bases in string

static unsigned long count_gc(const hashn::key_type_base &key) {
	const std::string s(convert_key(key));
	unsigned long n(0);
	std::string::const_iterator a(s.begin());
	const std::string::const_iterator end_a(s.end());
	for (; a != end_a; ++a) {
		if (*a == 'G' || *a == 'g' || *a == 'C' || *a == 'c') {
			++n;
		}
	}
	return n;
}

// print histogram of n-mer occurrences, potentially with gc percent at given
// frequencies

static void print_mer_histogram(hashn &mer_list) {
	std::map<hashn::value_type, unsigned long> counts;
	std::map<hashn::value_type, unsigned long> gc_counts;
	hashn::key_type comp_key(mer_list);
	hashn::const_iterator a(mer_list.begin());
	const hashn::const_iterator end_a(mer_list.end());
	for (; a != end_a; ++a) {
		reverse_key(a.key, comp_key);
		if (a.key == comp_key) {
			counts[a.value] += 2;
			if (opt_print_gc) {
				gc_counts[a.value] += 2 * count_gc(a.key);
			}
		} else {
			++counts[a.value];
			if (opt_print_gc) {
				gc_counts[a.value] += count_gc(a.key);
			}
		}
	}
	std::map<hashn::value_type, unsigned long>::const_iterator c(counts.begin());
	const std::map<hashn::value_type, unsigned long>::const_iterator end_c(counts.end());
	// don't include single occurrences in total
	// (hashes don't have non-positive occurrence values)
	if (c != end_c && c->first == 1) {
		++c;
	}
	double total(0);
	for (; c != end_c; ++c) {
		total += static_cast<double>(c->first) * static_cast<double>(c->second);
	}
	c = counts.begin();
	if (c != end_c && c->first == 1) {
		fprintf(fp_out, "%lu %lu\n", c->first, c->second);
		++c;
	}
	double i(0);
	for (; c != end_c; ++c) {
		const double x(100 * static_cast<double>(c->first) * static_cast<double>(c->second));
		i += x;
		if (opt_print_gc) {
			fprintf(fp_out, "%lu %lu %.2f %.2f %.2f\n", c->first, c->second, x / total, i / total, 100 * static_cast<double>(gc_counts[c->first]) / static_cast<double>(c->second) / static_cast<double>(opt_mer_length));
		} else {
			fprintf(fp_out, "%lu %lu %.2f %.2f\n", c->first, c->second, x / total, i / total);
		}
	}
}

static void print_mer_histogram_sub(hashn &mer_list) {
	std::map<hashn::value_type, unsigned long> counts[opt_readnames_exclude];
	hashn::const_iterator a(mer_list.begin());
	const hashn::const_iterator end_a(mer_list.end());
	hashn::value_type x[opt_readnames_exclude];
	for (; a != end_a; ++a) {
		a.get_alt_values(x);
		hashn::value_type n(a.value);
		for (int i(0); i != opt_readnames_exclude; ++i) {
			n += x[i];
		}
		if (n != x[0]) {
			hashn::value_type m(n);
			for (int i(0); i != opt_readnames_exclude; ++i) {
				m -= x[i];
				counts[i][n] += m;
			}
		}
	}
	for (int i(0); i != opt_readnames_exclude; ++i) {
		fprintf(fp_out, "\n");
		std::map<hashn::value_type, unsigned long>::const_iterator c(counts[i].begin());
		const std::map<hashn::value_type, unsigned long>::const_iterator end_c(counts[i].end());
		for (; c != end_c; ++c) {
			fprintf(fp_out, "%lu %lu\n", c->first, c->second);
		}
	}
}

static void print_mer_histogram_add(hashn &mer_list) {
	opt_readnames_exclude = -opt_readnames_exclude;
	std::map<hashn::value_type, unsigned long> counts[opt_readnames_exclude];
	hashn::const_iterator a(mer_list.begin());
	const hashn::const_iterator end_a(mer_list.end());
	hashn::value_type x[opt_readnames_exclude];
	for (; a != end_a; ++a) {
		a.get_alt_values(x);
		for (int i(0); i != opt_readnames_exclude; ++i) {
			if (x[i] != 0) {
				counts[i][a.value] += x[i];
			}
		}
	}
	for (int i(0); i != opt_readnames_exclude; ++i) {
		fprintf(fp_out, "\n");
		std::map<hashn::value_type, unsigned long>::const_iterator c(counts[i].begin());
		const std::map<hashn::value_type, unsigned long>::const_iterator end_c(counts[i].end());
		for (; c != end_c; ++c) {
			fprintf(fp_out, "%lu %lu\n", c->first, c->second);
		}
	}
	opt_readnames_exclude = -opt_readnames_exclude;
}

// read in read names from the given file, and add them to list

static void add_readnames(const char * const filename, std::map<std::string, hashn::offset_type> &list) {
	const int fd(open_compressed(filename));
	if (fd == -1) {
		fprintf(stderr, "Error: could not read %s\n", filename);
		return;
	}
	const hashn::offset_type x(1 << (abs(opt_readnames_exclude) - 1));
	std::string line;
	while (pfgets(fd, line) != -1) {
		const std::string s(strtostr(line));
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

// return the number represented by s, which may be suffixed by a k, m, or g
// which act as multipliers to the base amount

static size_t get_value(const std::string s) {
	// validate string - digits optionally followed by k, m, or g
	const size_t i(s.find_first_not_of("0123456789"));
	if (i == std::string::npos) {
		return atol(s.c_str());
	} else if (i + 1 == s.length()) {
		size_t x(atol(s.c_str()));
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
	} else {	// bad value
		return 0;
	}
}

// print usage statement and exit

static void print_usage() {
	fprintf(stderr,
		"usage: histogram [options] file1 [file2] ...\n"
		"    -a    give combined results for all files\n"
		"    -B ## process seq & qual file in batches of ## reads\n"
		"    -c    clip low quality\n"
		"    -d    when processing in batches, check for duplicates across whole file\n"
		"    -f ## when clipping quality or vector, use ## as the target quality [20]\n"
		"    -g    print percent gc content at each frequency\n"
		"    -h    print this information\n"
		"    -i    turn off status updates\n"
		"    -k ## skip reads smaller than this\n"
		"    -l ## filename containing names of reads to subtract from results\n"
		"          (histogram is given as count * frequency, rather than count)\n"
		"    -L ## filename containing names of reads to compare with results\n"
		"          (count is by given reads, frequency is by other reads)\n"
		"    -m ## set mer length [24]\n"
		"    -o ## print output to file instead of stdout\n"
		"    -p ## don't touch reads not matching pattern (an extended regex)\n"
		"    -q    turn off all warnings\n"
		"    -s ## save histogram memory structure to file\n"
		"    -S ## load histogram memory dump from given file\n"
		"    -t    strip first part of trace id\n"
		"    -T ## if the hash fills, store partial dumps with the given filename prefix\n"
		"    -v    clip vector\n"
		"    -V    print version\n"
		"    -w ## print frequency count instead of histogram, for all n-mers with\n"
		"          a frequency of at least ## [0 (off)]\n"
		"    -z ## number of possible n-mers to allocate memory for [200m]\n"
		"          (k, m, or g may be suffixed)\n"
		"    -Z    clean hash if it fills up\n"
	);
	exit(1);
}

static void get_opts(int argc, char **argv) {
	std::string opt_output;
	opt_aggregate = 0;
	opt_batch_size = 0;
	opt_clip_quality = 0;
	opt_clip_vector = 0;
	opt_feedback = 1;
	opt_frequency_cutoff = 0;
	opt_hash_clean = 0;
	opt_histogram_restore = -1;
	opt_mer_length = 24;
	opt_nmers = 200 * 1024 * 1024;
	opt_print_gc = 0;
	opt_quality_cutoff = 20;
	opt_readnames_exclude = 0;
	opt_skip_size = 0;
	opt_strip_tracename = 0;
	opt_track_dups = 0;
	opt_warnings = 1;
	int c;
	while ((c = getopt(argc, argv, "aB:cdf:ghik:l:L:m:o:p:qs:S:tT:vVw:z:Z")) != EOF) {
		switch (c) {
		    case 'a':
			opt_aggregate = 1;
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
		    case 'g':
			opt_print_gc = 1;
			break;
		    case 'h':
			print_usage();
			break;
		    case 'i':
			opt_feedback = 0;
			break;
		    case 'k':
			// use an int here to capture negative values
			c = atoi(optarg);
			if (c < 0) {
				fprintf(stderr, "Error: invalid skip size %d\n", c);
				print_usage();
			}
			opt_skip_size = c;
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
		    case 's':
			opt_save_file = optarg;
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
			opt_strip_tracename = 1;
			break;
		    case 'T':
			opt_tmp_file_prefix = optarg;
			break;
		    case 'v':
			opt_clip_vector = 1;
			break;
		    case 'V':
			fprintf(stderr, "histogram_hashn version %s%s\n", VERSION,
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
		    case 'Z':
			opt_hash_clean = 1;
			break;
		    default:
			fprintf(stderr, "Error: unknown option %c\n", c);
			print_usage();
		}
	}
	if (opt_histogram_restore != -1) {
		if (opt_nmers != 200 * 1024 * 1024) {
			fprintf(stderr, "Error: -S and -z options cannot both be specified\n");
			exit(1);
		} else if (opt_hash_clean) {
			fprintf(stderr, "Error: -S and -Z options cannot both be specified\n");
			exit(1);
		} else if (optind != argc) {
			fprintf(stderr, "Warning: fasta files being ignored, hash is being read from disk\n");
		}
	} else if (optind == argc) {
		fprintf(stderr, "Error: no files to process\n");
		print_usage();
	}
	if (opt_readnames_exclude != 0 && !opt_tmp_file_prefix.empty()) {
		fprintf(stderr, "Error: cannot use -T option with either -l or -L options\n");
		exit(1);
	}
	if (opt_frequency_cutoff != 0 && opt_readnames_exclude) {
		fprintf(stderr, "Warning: -w and -l/-L options conflict: ignoring -w option\n");
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
	int err(0);
	hashn mer_list;
	if (opt_hash_clean || !opt_tmp_file_prefix.empty()) {
		// have to add one to opt_mer_length,
		// as init_mer_constants() subtracts one
		mer_list.set_no_space_response((opt_hash_clean ? hashn::CLEAN_HASH: 0) | (opt_tmp_file_prefix.empty() ? 0 : hashn::TMP_FILE), opt_tmp_file_prefix);
	}
	if (opt_histogram_restore != -1) {
		mer_list.init_from_file(opt_histogram_restore);
	} else {
		mer_list.init(opt_nmers, opt_mer_length * 2, abs(opt_readnames_exclude));
		for (; optind != argc; ++optind) {
			if (opt_feedback) {
				fprintf(stderr, "Reading in %s\n", argv[optind]);
			}
			ReadFile file(argv[optind], opt_batch_size, opt_track_dups);
			if (file.seq_file.empty()) {
				++err;
				continue;
			}
			while (file.read_batch(opt_warnings) != -1) {
				if (opt_readnames_exclude) {
					if (!add_sequence_mers(file.read_list.begin(), file.read_list.end(), mer_list, opt_readnames)) {
						fprintf(stderr, "Error: n-mer list incomplete - give a larger -z value\n");
					}
				} else if (!add_sequence_mers(file.read_list.begin(), file.read_list.end(), mer_list)) {
					fprintf(stderr, "Error: n-mer list incomplete - give a larger -z value\n");
				}
			}
			if (!opt_aggregate) {
				if (opt_feedback) {
					fprintf(stderr, "Printing histogram\n");
				}
				fprintf(fp_out, "%s\n", argv[optind]);
				int i(strlen(argv[optind]));
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
	if (!opt_save_file.empty()) {
		save_memory(mer_list);
	}
	return err;
}
