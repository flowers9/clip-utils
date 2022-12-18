#include "hashl.h"	// hashl
#include "open_compressed.h"	// close_compressed(), get_suffix(), open_compressed(), pfgets()
#include "version.h"	// VERSION
#include "write_fork.h"	// close_fork(), write_fork()
#include <getopt.h>	// getopt(), optarg, optind
#include <list>		// list<>
#include <map>		// map<>
#include <sstream>	// istringstream
#include <stdio.h>	// EOF, fprintf(), stderr, stdout
#include <stdlib.h>	// exit()
#include <string.h>	// strlen()
#include <string>	// string
#include <sys/types.h>	// size_t
#include <vector>	// vector<>
#include <stdint.h>	// uint64_t

static FILE *fp_out(stdout);
static bool opt_print_gc;
static bool opt_warnings;
static double opt_load_lower_bound;
static double opt_load_upper_bound;
static int opt_histogram_restore;
static int opt_mer_length;
static std::string opt_save_file;
static unsigned long opt_frequency_cutoff;

static void save_memory(const hashl &mer_list) {
	std::string suffix;
	get_suffix(opt_save_file, suffix);
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
	const int fd(write_fork(args, opt_save_file));
	if (fd == -1) {
		fprintf(stderr, "Error: could not save memory\n");
		exit(1);
	}
	mer_list.save(fd);
	close_fork(fd);
}

// print n-mer occurence frequency

static void print_mer_frequency(const hashl &mer_list) {
	hashl::key_type key(mer_list), comp_key(mer_list);
	hashl::const_iterator a(mer_list.begin());
	const hashl::const_iterator end_a(mer_list.end());
	for (; a != end_a; ++a) {
		if (a.value >= opt_frequency_cutoff) {
			a.get_key(key);
			fprintf(fp_out, "%s %lu\n", convert_key(key).c_str(), a.value);
			reverse_key(key, comp_key);
			if (key != comp_key) {
				fprintf(fp_out, "%s %lu\n", convert_key(comp_key).c_str(), a.value);
			}
		}
	}
}

// number of gc bases in string

static unsigned long count_gc(const hashl::key_type &key) {
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

static void print_mer_histogram(const hashl &mer_list) {
	std::map<hashl::value_type, unsigned long> counts;
	std::map<hashl::value_type, unsigned long> gc_counts;
	hashl::key_type key(mer_list), comp_key(mer_list);
	hashl::const_iterator a(mer_list.begin());
	const hashl::const_iterator end_a(mer_list.end());
	for (; a != end_a; ++a) {
		a.get_key(mer_list, key);
		reverse_key(key, comp_key);
		if (key == comp_key) {
			counts[a.value] += 2;
			if (opt_print_gc) {
				gc_counts[a.value] += 2 * count_gc(key);
			}
		} else {
			++counts[a.value];
			if (opt_print_gc) {
				gc_counts[a.value] += count_gc(key);
			}
		}
	}
	std::map<hashl::value_type, unsigned long>::const_iterator c(counts.begin());
	const std::map<hashl::value_type, unsigned long>::const_iterator end_c(counts.end());
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
		"usage: histogram [options] file1 [file2] ...\n"
		"    -g    print percent gc content at each frequency\n"
		"    -h    print this information\n"
		"    -i    turn off status updates\n"
		"    -l ## lower bound for load value [.5]\n"
		"    -L ## upper bound for load value [.9]\n"
		"    -m ## set mer length [24]\n"
		"    -o ## print output to file instead of stdout\n"
		"    -q    turn off all warnings\n"
		"    -s ## save histogram memory structure to file\n"
		"    -S ## load histogram memory dump from given file\n"
		"    -V    print version\n"
		"    -w ## print frequency count instead of histogram, for all n-mers with\n"
		"          a frequency of at least ## [0 (off)]\n"
	);
	exit(1);
}

static void get_opts(int argc, char **argv) {
	std::string opt_output;
	opt_feedback = 1;
	opt_frequency_cutoff = 0;
	opt_histogram_restore = -1;
	opt_mer_length = 24;
	opt_print_gc = 0;
	opt_warnings = 1;
	opt_load_lower_bound = .5;
	opt_load_upper_bound = .9;
	int c;
	while ((c = getopt(argc, argv, "ghil:L:o:qs:S:Vw:")) != EOF) {
		switch (c) {
		    case 'g':
			opt_print_gc = 1;
			break;
		    case 'h':
			print_usage();
			break;
		    case 'i':
			opt_feedback = 0;
			break;
		    case 'l':
			std::istringstream(optarg) >> opt_load_lower_bound;
			break;
		    case 'L':
			std::istringstream(optarg) >> opt_load_upper_bound;
			break;
		    case 'm':
			std::istringstream(optarg) >> opt_mer_length;
			if (opt_mer_length < 1) {
				fprintf(stderr, "Error: bad mer length\n");
				exit(1);
			}
			break;
		    case 'o':
			opt_output = optarg;
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
				exit(1);
			}
			break;
		    case 'V':
			fprintf(stderr, "histogram_hashl version %s%s\n", VERSION,
#ifdef COMPRESS_READS
" (read compression)"
#else
""
#endif
);
			exit(0);
		    case 'w':
			std::istringstream(optarg) >> opt_frequency_cutoff;
			break;
		    default:
			fprintf(stderr, "Error: unknown option %c\n", c);
			print_usage();
		}
	}
	if (opt_histogram_restore != -1) {
		if (optind != argc) {
			fprintf(stderr, "Warning: fasta files being ignored, hash is being read from disk\n");
		}
		if (!opt_save_file.empty()) {
			fprintf(stderr, "Warning: ignoring -S option - copy it over yourself!\n");
			opt_save_file.clear();
		}
	} else if (optind == argc) {
		fprintf(stderr, "Error: no files to process\n");
		print_usage();
	} else if (opt_histogram_restore == -1 && opt_load_lower_bound > opt_load_upper_bound) {
		fprintf(stderr, "Error: lower load bound must be less than upper\n");
		exit(1);
	}
	if (!opt_output.empty()) {
		fp_out = fopen(opt_output.c_str(), "w");
		if (fp_out == NULL) {
			fprintf(stderr, "Error: could not write to %s\n", opt_output.c_str());
			exit(1);
		}
	}
}

// read file and get total space needed for data and metadata
static void get_sizes(const char * const file, size_t &metadata_size, size_t &data_size, size_t &max_kmers, std::vector<size_t> &read_ends) {
	const int fd(open_compressed(file));
	if (fd == -1) {
		fprintf(stderr, "Error: open: %s\n", file);
		exit(1);
	}
	std::string line;
	if (pfgets(fd, line) == -1) {		// empty file
		fprintf(stderr, "Warning: empty file: %s\n", file);
	} else if (line[0] == '>') {		// fasta file
		do {
			size_t i(1);
			for (; i < line.size() && !isspace(line[i]); ++i) { }
			size_t seq_size(0);
			while (pfgets(fd, line) != -1 && line[0] != '>') {
				seq_size += line.size();
			}
			// reads shorter than the mer length are skipped
			if (seq_size >= opt_mer_length) {
				metadata_size += i;
				data_size += seq_size;
				max_kmers += seq_size - opt_mer_length + 1;
				read_ends.push_back(data_size);
			}
		} while (!line.empty());
	} else if (line[0] == '@') {		// fastq file
		do {
			size_t i(1);
			for (; i < line.size() && !isspace(line[i]); ++i) { }
			if (pfgets(fd, line) == -1) {		// sequence
				fprintf(stderr, "Error: truncated fastq file: %s\n", file);
				exit(1);
			}
			const size_t seq_size(line.size());
			// reads shorter than the mer length are skipped
			if (seq_size >= opt_mer_length) {
				metadata_size += i;
				data_size += seq_size;
				max_kmers += seq_size - opt_mer_length + 1;
				read_ends.push_back(data_size);
			}
			// skip quality header and quality
			if (pfgets(fd, line) == -1 || pfgets(fd, line) == -1) {
				fprintf(stderr, "Error: truncated fastq file: %s\n", file);
				exit(1);
			}
		} while (pfgets(fd, line) != -1);
	} else {
		fprintf(stderr, "Error: unknown file format: %s\n", file);
		exit(1);
	}
	close_compressed(fd);
	if (!read_ends.empty()) {	// only count if we're not going to skip the file
		metadata_size += strlen(file) + 1 + read_ends.size() * sizeof(uint64_t);
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
		fprintf("Error: non-ACGT basepair: %c\n", c);
		exit(1);
	}
}

// sequence is stored from top bits to bottom bits of each offset
static void copy_data(const std::string &seq, hashl::base_type * const data, size_t data_offset) {
	size_t i((2 * data_offset) / (sizeof(hashl::base_type) * 8));
	size_t j(sizeof(hashl::base_type) * 8 - (2 * data_offset) % (sizeof(hashl::base_type) * 8));
	std::string::const_iterator a(seq.begin());
	const std::string::const_iterator end_a(seq.end());
	for (; a != end_a; ++a) {
		if (j) {
			j -= 2;
			data[i] |= convert_char(*a) << j;
		} else {
			j = sizeof(hashl::base_type) * 8 - 2;
			data[++i] = convert_char(*a) << j;
		}
	}
}

static void read_file(const char * const file, char *metadata, hashl::base_type *data, size_t &metadata_offset, size_t &data_offset) {
	const int fd(open_compressed(file));
	if (fd == -1) {
		fprintf(stderr, "Error: open: %s\n", file);
		exit(1);
	}
	std::string line, seq;
	if (pfgets(fd, line) == -1) {		// empty file
		fprintf(stderr, "Warning: empty file: %s\n", file);
	} else if (line[0] == '>') {		// fasta file
		std::string header;
		do {
			size_t i(1);
			for (; i < line.size() && !isspace(line[i]); ++i) { }
			header = line.substr(1, i - 1);
			seq.clear();
			while (pfgets(fd, line) != -1 && line[0] != '>') {
				seq += line;
			}
			// reads shorter than the mer length are skipped
			if (seq_size >= opt_mer_length) {
				memcpy(metadata + metadata_offset, header.c_str(), header.size() + 1);
				metadata_offset += i;
				copy_data(seq, data, data_offset);
				data_offset += seq.size();
			}
		} while (line[0] == '>');	// safe for empty std::string
	} else if (line[0] == '@') {		// fastq file
		do {
			if (pfgets(fd, seq) == -1) {		// sequence
				fprintf(stderr, "Error: truncated fastq file: %s\n", file);
				exit(1);
			}
			// reads shorter than the mer length are skipped
			if (seq.size() >= opt_mer_length) {
				size_t i(1);
				for (; i < line.size() && !isspace(line[i]); ++i) { }
				memcpy(metadata + metadata_offset, line.data() + 1, i - 1);
				metadata_offset += i;
				metadata[metadata_offset - 1] = 0;
				copy_data(seq, data, data_offset);
				data_offset += seq.size();
			}
			// skip quality header and quality
			// (use seq because it'll be the same length as quality)
			if (pfgets(fd, line) == -1 || pfgets(fd, seq) == -1) {
				fprintf(stderr, "Error: truncated fastq file: %s\n", file);
				exit(1);
			}
		} while (pfgets(fd, line) != -1);
	} else {
		fprintf(stderr, "Error: unknown file format: %s\n", file);
		exit(1);
	}
	close_compressed(fd);
}

static void count_nmers(hashl &mer_list, const hashl::base_type * const data, const std::vector<std::vector<uint64_t> > &read_ends) {
	hashl::key_type key(mer_list), comp_key(mer_list);
	uint64_t i(0);
	// go through all reads (nmers can't cross read boundaries)
	std::vector<std::vector<uint64_t> >::const_iterator a(read_ends.begin());
	const std::vector<std::vector<uint64_t> >::const_iterator end_a(read_ends.end());
	for (; a != end_a; ++a) {
		std::vector<uint64_t>::const_iterator b(a->begin());
		const std::vector<uint64_t>::const_iterator end_b(a->end());
		for (; b != end_b; ++b) {
			uint64_t j((2 * i) / (sizeof(hashl::base_type) * 8));
			uint64_t k(sizeof(hashl::data_type) * 8 - (2 * i) % (sizeof(hashl::base_type) * 8) - 2);
			const uint64_t end_i(i + opt_mer_length - 1);
			// load keys with opt_mer_length - 1 basepairs
			for (; i < end_i; ++i) {
				const hashl::base_type c((data[i] >> k) & 3);
				key.push_back(c);
				comp_key.push_front(3 - c);
				if (k) {
					k -= 2;
				} else {
					++j;
					k = sizeof(hashl::data_type) * 8 - 2;
				}
			}
			// run over all nmers, one basepair at a time
			for (; i < *b; ++i) {
				const hashl::base_type c((data[j] >> k) & 3);
				key.push_back(c);
				comp_key.push_front(3 - c);
				if (k) {
					k -= 2;
				} else {
					++j;
					k = sizeof(hashl::data_type) * 8 - 2;
				}
				if (!mer_list.increment(key < comp_key ? key : comp_key, i)) {
					fprintf(stderr, "Error: ran out of space in hash\n");
					exit(1);
				}
			}
		}
	}
}

// metadata format: # of files, [ filename, # of reads, read offsets, read names ]
// strings are null delimited, # and read lengths are uint64_t
// hashl boilerplate check will ensure byte order is the same read as written

static void read_in_files(int argc, char **argv, hashl &mer_list) {
	const uint64_t file_count(argc - optind);
	// leave room for counts for number of files and reads per file
	size_t metadata_size((file_count + 1) * sizeof(uint64_t)), data_size(0), max_kmers(0);
	std::vector<std::vector<uint64_t> > read_ends(file_count, std::vector<size_t>());
	for (int i(0; i < file_count; ++i) {
		get_sizes(argv[i + optind], metadata_size, data_size, max_kmers, read_ends[j]);
	}
	// allocate arrays for data and metadata
	char *metadata = new char[metadata_size];
	hashl::base_type data = new hashl::base_type[data_size];
	memcpy(metadata, file_count, sizeof(file_count));
	size_t metadata_offset(sizeof(file_count)), data_offset(0);
	for (int i(0); i < file_count; ++i) {
		if (read_ends[i].empty()) {			// skip empty files
			continue;
		}
		const size_t name_length(strlen(argv[i + optind]) + 1);	// include trailing null
		memcpy(metadata + metadata_offset, argv[i + optind], name_length);
		metadata_offset += name_length;
		uint64_t x(read_ends[i].size());
		memcpy(metadata + metadata_offset, &x, sizeof(x));
		metadata_offset += sizeof(x);
		memcpy(metadata + metadata_offset, &read_ends[i][0], x * sizeof(uint64_t));
		metadata_offset += x * sizeof(uint64_t);
		read_file(argv[i + optind], metadata, data, metadata_offset, data_offset);
	}
	// put data and metadata into mer_list
	mer_list.init(max_kmers, opt_mer_length * 2, data, data_size);
	mer_list.set_metadata(metadata, metadata_size);
	count_nmers(mer_list, data, read_ends);
}

int main(int argc, char **argv) {
	get_opts(argc, argv);
	if (opt_feedback) {
		fprintf(stderr, "Initializing n-mer hash\n");
	}
	hashl mer_list;
	if (opt_histogram_restore != -1) {
		mer_list.init_from_file(opt_histogram_restore);
	} else {
		read_in_files(argc, argv, mer_list);
	}
	if (opt_feedback) {
		print_final_input_feedback(mer_list);
		fprintf(stderr, "Printing histogram\n");
	}
	if (opt_frequency_cutoff == 0) {
		print_mer_histogram(mer_list);
	} else {
		print_mer_frequency(mer_list);
	}
	if (fp_out != stdout) {
		fclose(fp_out);
	}
	if (!opt_save_file.empty()) {
		save_memory(mer_list);
	}
	return err;
}
