#include "breakup_line.h"	// breakup_line()
#include "open_compressed.h"	// close_compressed(), open_compressed(), pfgets()
#include "read.h"	// Read, opt_N_is_vector, opt_add_range, opt_all_p20, opt_base_cutoff, opt_clip_quality, opt_clip_vector, opt_line_length, opt_linker, opt_minimum_clip, opt_quality_cutoff, opt_pacbio, opt_repeat_clip, opt_strict_quality, opt_strip_trailing_zero_qual, read_name_translation
#include "read_file.h"	// ReadFile, opt_readname_match, opt_strip_tracename
#include "strtostr.h"	// strtostr()
#include "version.h"	// VERSION
#include <getopt.h>	// getopt(), optarg, optind
#include <list>		// list<>
#include <map>		// map<>
#include <regex.h>	// REG_EXTENDED, REG_ICASE
#include <sstream>	// istringstream
#include <stdio.h>	// EOF, FILE, fclose(), fopen(), fprintf(), stderr
#include <stdlib.h>	// exit()
#include <string>	// string
#include <strings.h>	// strcasecmp()
#include <sys/types.h>	// size_t
#include <vector>	// vector<>

static bool opt_arachne_output;
static bool opt_print_quality;
static bool opt_print_seq_and_qual;
static bool opt_qual_warning;
static bool opt_track_dups;
static bool opt_warnings;
static size_t opt_batch_size;
static std::string opt_output_file;
static unsigned int opt_phred_count_cutoff;
static unsigned int opt_phred_mask_cutoff;
static unsigned int opt_qual_length_cutoff;

class OutputStreams {
    public:
	FILE *fp_seq, *fp_qual, *fp_seq_raw, *fp_qual_raw;
	OutputStreams(void) : fp_seq(stdout), fp_qual(stdout), fp_seq_raw(NULL), fp_qual_raw(NULL) { }
	explicit OutputStreams(const std::string &__s) : fp_seq(NULL), fp_qual(NULL), fp_seq_raw(NULL), fp_qual_raw(NULL) {
		if (opt_arachne_output || opt_print_seq_and_qual) {
			fp_seq = fopen(__s.c_str(), "w");
			fp_qual = fopen(std::string(__s + ".qual").c_str(), "w");
			if (opt_arachne_output) {
				fp_seq_raw = fopen(std::string(__s + ".raw").c_str(), "w");
				fp_qual_raw = fopen(std::string(__s + ".raw.qual").c_str(), "w");
			}
		} else if (opt_output_file.empty() || opt_output_file == "-") {
			fp_seq = fp_qual = stdout;
		} else {
			fp_seq = fp_qual = fopen(__s.c_str(), "w");
		}
	}
	~OutputStreams(void) {
		if (fp_seq != NULL && fp_seq != stdout) {
			fclose(fp_seq);
		}
		if (fp_qual != NULL && fp_qual != fp_seq) {
			fclose(fp_qual);
		}
		if (fp_seq_raw != NULL) {
			fclose(fp_seq_raw);
		}
		if (fp_qual_raw != NULL) {
			fclose(fp_qual_raw);
		}
	}
};

static std::string make_output_filename(const std::string &input_filename) {
	std::string s(input_filename);
	const size_t i(s.rfind('.'));
	if (i != std::string::npos && (s.substr(i) == ".bz2" || s.substr(i) == ".gz" || s.substr(i) == ".Z" || s.substr(i) == ".xz")) {
		s.resize(i);
	}
	s += ".output";
	return s;
}

// print reads in the same order they appeared in the quality file,
// clipping as specified

static void print_clipped_sequence(const ReadFile &file, OutputStreams &out) {
	std::list<Read>::const_iterator a(file.read_list.begin());
	const std::list<Read>::const_iterator end_a(file.read_list.end());
	for (; a != end_a; ++a) {
		if (a->quality_stop - a->quality_start < opt_qual_length_cutoff) {
			if (opt_qual_warning) {
                        	fprintf(stderr, "Warning: quality sequence too short, skipping %s\n", a->name().c_str());
			}
		} else if (a->phred_count < opt_phred_count_cutoff) {
			if (opt_qual_warning) {
				fprintf(stderr, "Warning: phred20 count too small, skipping %s\n", a->name().c_str());
			}
		} else if (opt_arachne_output || opt_print_seq_and_qual) {
			a->print_sequence(out.fp_seq);
			a->print_quality(out.fp_qual);
		} else if (opt_print_quality) {
			a->print_quality(out.fp_qual);
		} else {
			a->print_sequence(out.fp_seq);
		}
	}
	if (opt_arachne_output) {
		const bool opt_clip_quality_orig(opt_clip_quality);
		const bool opt_clip_vector_orig(opt_clip_vector);
		opt_clip_quality = 0;
		opt_clip_vector = 0;
		for (a = file.read_list.begin(); a != end_a; ++a) {
			a->print_sequence(out.fp_seq_raw);
			a->print_quality(out.fp_qual_raw);
		}
		opt_clip_quality = opt_clip_quality_orig;
		opt_clip_vector = opt_clip_vector_orig;
	}
}

/* read in read names from the given file, and add them to opt_readname_match */

static void add_readnames_match(const char * const filename) {
	const int fd(open_compressed(filename));
	if (fd == -1) {
		fprintf(stderr, "Error: could not read %s\n", filename);
		return;
	}
	std::string line;
	while (pfgets(fd, line) != -1) {
		const std::string s(strtostr(line));
		if (!s.empty()) {
			opt_readname_match[s] = 1;
		}
	}
	close_compressed(fd);
}

static void read_translations(const char * const filename) {
	const int fd(open_compressed(filename));
	if (fd == -1) {
		fprintf(stderr, "Error: could not read %s\n", filename);
		return;
	}
	std::string line;
	while (pfgets(fd, line) != -1) {
		std::vector<std::string> list;
		breakup_line(line, list);
		if (list.size() > 1) {
			read_name_translation[list[0]] = list[1];
		}
	}
	close_compressed(fd);
}

// print usage statement and exit

static void print_usage() {
	fprintf(stderr, "usage: clip [options] file1 [file2] ...\n"
		"    -b    clip vector, and treat N's as X's when finding vector\n"
		"    -B ## process seq & qual file in batches of ## reads\n"
		"    -c ## delete sequences with less than ## basepairs after clipping\n"
		"    -d    when processing in batches, check for duplicates for whole file\n"
		"    -D ## after clipping for quality, clip the end to remove any section with\n"
		"          an average repeat length of this much or greater (fractions okay)\n"
		"    -f ## when clipping quality or vector, use ## as the target quality [20]\n"
		"    -h    print this usage information\n"
		"    -H    use stricter rules when doing quality or vector clipping\n"
		"    -k ## clip linker and any sequence past it\n"
		"    -l ## only process readnames found in the given file\n"
		"    -L ## length to wrap output seq/qual (0 = no wrapping)\n"
		"    -m ## mask printed sequence with quality less than ##\n"
		"    -n    do not clip low quality\n"
		"    -N ## use file to translate names\n"
		"    -o [option]\n"
		"       arachne  create standard and raw versions of seq and qual files\n"
		"       est      equivalent to -c 50 -n -p 50 -r -v\n"
		"       qual     prints quality rather than sequence\n"
		"       seq_and_qual  creates both sequence and quality files\n"
		"       pacbio   modify pacbio-style read name if read is trimmed\n"
		"    -p ## delete sequences with less than ## qualities of 20 or more\n"
		"    -P    when counting phred20s, ignore non-ACGT basepairs\n"
		"    -q    turn off all warnings\n"
		"    -r    add quality clip range\n"
		"    -R ## when clipping vector, don't consider sequence with greater than this\n"
		"           fraction of a single base\n"
		"    -s ## make sure clipping includes at least the first ## basepairs\n"
		"          (this will modify the displayed range, if -r is given, otherwise\n"
		"           this will modify the displayed sequence)\n"
		"    -S ## base name to use when saving the output [stdout, or inputfile.output]\n"
		"    -t    strip first part of trace id\n"
		"    -v    clip vector\n"
		"    -V    print version\n"
		"    -w    turn on short quality sequence warning\n"
		"    -z    strip trailing zero qual\n"
	);
	exit(1);
}

static void get_opts(int argc, char **argv) {
	// set option defaults
	opt_N_is_vector = 0;
	opt_add_range = 0;
	opt_all_p20 = 1;
	opt_arachne_output = 0;
	opt_base_cutoff = 0;
	opt_batch_size = 0;
	opt_clip_quality = 1;
	opt_clip_vector = 0;
	opt_minimum_clip = 0;
	opt_pacbio = 0;
	opt_phred_count_cutoff = 0;
	opt_phred_mask_cutoff = 0;
	opt_print_quality = 0;
	opt_print_seq_and_qual = 0;
	opt_qual_length_cutoff = 0;
	opt_qual_warning = 0;
	opt_quality_cutoff = 20;
	opt_repeat_clip = 0;
	opt_strict_quality = 0;
	opt_strip_tracename = 0;
	opt_strip_trailing_zero_qual = 0;
	opt_track_dups = 0;
	opt_warnings = 1;
	// read in options
	int c;
	while ((c = getopt(argc, argv, "bB:c:dD:f:hHk:l:L:m:nN:o:p:PqrR:s:S:tvVwz")) != EOF) {
		switch (c) {
		    case 'b':
			opt_N_is_vector = 1;
			opt_clip_vector = 1;
			break;
		    case 'B':
			std::istringstream(optarg) >> c;
			if (c < 0) {
				print_usage();
			}
			opt_batch_size = c;
			break;
		    case 'c':
			std::istringstream(optarg) >> c;
			if (c < 0) {
				print_usage();
			}
			opt_qual_length_cutoff = c;
			break;
		    case 'd':
			opt_track_dups = 1;
			break;
		    case 'D':
			std::istringstream(optarg) >> opt_repeat_clip;
			if (opt_repeat_clip < 1) {
				print_usage();
			}
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
		    case 'H':
			opt_strict_quality = 1;
			break;
		    case 'k':
			opt_linker.initialize(optarg, 1, REG_EXTENDED | REG_ICASE);
			break;
		    case 'l':
			add_readnames_match(optarg);
			break;
		    case 'L':
			std::istringstream(optarg) >> c;
			if (c < 0) {
				print_usage();
			}
			opt_line_length = c;
			break;
		    case 'm':
			std::istringstream(optarg) >> c;
			if (c < 0) {
				print_usage();
			}
			opt_phred_mask_cutoff = c;
			break;
		    case 'n':
			opt_clip_quality = 0;
			break;
		    case 'N':
			read_translations(optarg);
			break;
		    case 'o':
			if (strcasecmp(optarg, "est") == 0) {
				opt_qual_length_cutoff = 50;
				opt_clip_quality = 0;
				opt_phred_count_cutoff = 50;
				opt_add_range = 1;
				opt_clip_vector = 1;
			} else if (strcasecmp(optarg, "qual") == 0) {
				opt_print_quality = 1;
			} else if (strcasecmp(optarg, "seq_and_qual") == 0) {
				opt_print_seq_and_qual = 1;
			} else if (strcasecmp(optarg, "arachne") == 0) {
				opt_strip_tracename = 1;
				opt_arachne_output = 1;
			} else if (strcasecmp(optarg, "pacbio") == 0) {
				opt_pacbio = 1;
			} else {
				print_usage();
			}
			break;
		    case 'p':
			std::istringstream(optarg) >> c;
			if (c < 0) {
				print_usage();
			}
			opt_phred_count_cutoff = c;
			break;
		    case 'P':
			opt_all_p20 = 0;
			break;
		    case 'q':
			opt_warnings = 0;
			break;
		    case 'r':
			opt_add_range = 1;
			break;
		    case 'R':
			std::istringstream(optarg) >> opt_base_cutoff;
			if (opt_base_cutoff < 0 || 1 < opt_base_cutoff) {
				print_usage();
			}
			opt_clip_vector = 1;
			break;
		    case 's':
			std::istringstream(optarg) >> c;
			if (c < 0) {
				print_usage();
			}
			opt_minimum_clip = c;
			break;
		    case 'S':
			opt_output_file = optarg;
			break;
		    case 't':
			opt_strip_tracename = 1;
			break;
		    case 'v':
			opt_clip_vector = 1;
			break;
		    case 'V':
			fprintf(stderr, "clip version %s%s\n", VERSION,
#ifdef COMPRESS_READS
				" (read compression)"
#else
				""
#endif
			);
			exit(0);
		    case 'w':
			opt_qual_warning = 1;
			break;
		    case 'z':
			opt_strip_trailing_zero_qual = 1;
			break;
		    default:
			print_usage();
		}
	}
	if (optind == argc) {			// no files specified
		print_usage();
	}
	if (opt_print_quality && (opt_print_seq_and_qual || opt_arachne_output)) {
		fprintf(stderr, "Warning: output is to file, not stdout: ignoring -o qual option\n");
	}
	// check option consistency
	if (opt_add_range && opt_clip_quality) {
		fprintf(stderr, "Warning: quality ranges are enabled, so quality clipping is disabled\n");
	}
	if (opt_strict_quality && !opt_clip_vector && !opt_clip_quality && !opt_add_range) {
		fprintf(stderr, "Warning: strict quality clipping was asked for, but no clipping or ranges were enabled, so disabling\n");
		opt_strict_quality = 0;
	}
	if (opt_pacbio && (!read_name_translation.empty() || opt_add_range)) {
		fprintf(stderr, "Warning: cannot perform pacbio name translations with -N or -r options, so disabling\n");
		opt_pacbio = 0;
	}
}

int main(int argc, char **argv) {
	get_opts(argc, argv);
	if (opt_add_range) {
		opt_clip_quality = 1;
	}
	int err(0);
	for (; optind < argc; ++optind) {
		ReadFile file(argv[optind], opt_batch_size, opt_track_dups);
		if (file.seq_file.empty()) {
			++err;
			continue;
		}
		const std::string output_file(opt_output_file.empty() ? make_output_filename(file.seq_file) : opt_output_file);
		OutputStreams out(output_file);
		while (file.read_batch(opt_warnings) != -1) {
			if (opt_add_range) {
				opt_clip_quality = 0;
			}
			if (opt_phred_mask_cutoff > 0) {
				file.mask_by_phred(opt_phred_mask_cutoff);
			}
			print_clipped_sequence(file, out);
			if (opt_add_range) {
				opt_clip_quality = 1;
			}
		}
	}
	return err;
}
