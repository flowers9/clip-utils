#include "hashl.h"	// hashl
#include "open_compressed.h"	// close_compressed(), get_suffix(), open_compressed(), pfgets()
#include "time_used.h"	// elapsed_time(), start_time()
#include "version.h"	// VERSION
#include "write_fork.h"	// close_fork(), write_fork()
#include <getopt.h>	// getopt(), optarg, optind
#include <list>		// list<>
#include <map>		// map<>
#include <set>		// set<>
#include <sstream>	// istringstream, ostringstream
#include <stdint.h>	// uint64_t
#include <stdio.h>	// EOF, fprintf(), stderr, stdout
#include <stdlib.h>	// exit()
#include <string.h>	// memcpy(), strlen()
#include <string>	// string
#include <utility>	// make_pair(), pair<>
#include <vector>	// vector<>

static FILE *fp_out(stdout);
static bool opt_feedback;
static bool opt_print_gc;
static double opt_load_lower_bound;
static double opt_load_upper_bound;
static int opt_histogram_restore;
static size_t opt_mer_length;
static size_t opt_nmers;
static std::string opt_save_file;
static unsigned long opt_frequency_cutoff;

class hash_metadata {
    private:				// convenience variables for read_data()
	hashl::base_type *data;
	size_t byte_offset, bit_offset;
    private:
	std::vector<std::string> files;
	std::vector<std::vector<std::string> > reads;
	// inclusive start, exclusive end
	std::vector<std::vector<std::vector<std::pair<uint64_t, uint64_t> > > > read_ranges;
    public:
	hash_metadata(void) { }
	~hash_metadata(void) { }
	void add_file(const std::string &);		// add new file
	void add_read(const std::string &);		// add new read at current file
	void add_read_range(uint64_t, uint64_t);	// add new read range at current read
	void finalize(void);				// remove last adds if empty
	const void *pack(size_t &) const;
	void unpack(const char *);
	void print(void) const;
	const hashl::base_type *read_data(size_t &);
	std::pair<size_t, size_t> total_reads(void) const;	// reads & read ranges
	size_t max_kmers(void) const;
	std::vector<size_t> read_ends(void) const;
    private:
	void read_file(size_t);
	void get_subreads(const std::string &, const std::vector<std::pair<uint64_t, uint64_t> > &);
};

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
		fprintf(stderr, "Error: non-ACGT basepair: %c\n", c);
		exit(1);
	}
}

void hash_metadata::get_subreads(const std::string &seq, const std::vector<std::pair<uint64_t, uint64_t> > &ranges) {
	for (size_t i(0); i < ranges.size(); ++i) {
		for (size_t j(ranges[i].first); j < ranges[i].second; ++j) {
			if (bit_offset) {
				bit_offset -= 2;
				data[byte_offset] |= convert_char(seq[j]) << j;
			} else {
				bit_offset = sizeof(hashl::base_type) * 8 - 2;
				data[++byte_offset] = convert_char(seq[j]) << j;
			}
		}
	}
}

void hash_metadata::read_file(const size_t i) {
	const int fd(open_compressed(files[i]));
	if (fd == -1) {
		fprintf(stderr, "Error: open: %s\n", files[i].c_str());
		exit(1);
	}
	size_t j(0);	// which read we're on (but only incremented for reads we use)
	std::string line, seq;
	if (pfgets(fd, line) == -1) {		// empty file
		fprintf(stderr, "Error: File is now empty: %s\n", files[i].c_str());
		exit(1);
	} else if (line[0] == '>') {		// fasta file
		std::string header;
		do {
			const std::string &read_name(reads[i][j]);
			if (read_name.compare(0, read_name.size(), line, 1, read_name.size()) && (line.size() == read_name.size() + 1 || isspace(line[read_name.size() + 1]))) {
				seq.clear();
				while (pfgets(fd, line) != -1 && line[0] != '>') {
					seq += line;
				}
				get_subreads(seq, read_ranges[i][j]);
				++j;
			} else {
				while (pfgets(fd, line) != -1 && line[0] != '>') { }
			}
		} while (j < reads[i].size() && line[0] == '>');
	} else if (line[0] == '@') {		// fastq file
		do {
			if (pfgets(fd, seq) == -1) {		// sequence
				fprintf(stderr, "Error: truncated fastq file: %s\n", files[i].c_str());
				exit(1);
			}
			const std::string &read_name(reads[i][j]);
			if (read_name.compare(0, read_name.size(), line, 1, read_name.size()) && (line.size() == read_name.size() + 1 || isspace(line[read_name.size() + 1]))) {
				get_subreads(seq, read_ranges[i][j]);
				++j;
			}
			// skip quality header and quality
			// (use seq because it'll be the same length as quality)
			if (pfgets(fd, line) == -1 || pfgets(fd, seq) == -1) {
				fprintf(stderr, "Error: truncated fastq file: %s\n", files[i].c_str());
				exit(1);
			}
		} while (j < reads[i].size() && pfgets(fd, line) != -1);
	} else {
		fprintf(stderr, "Error: unknown file format: %s\n", files[i].c_str());
		exit(1);
	}
	if (j < reads[i].size()) {
		fprintf(stderr, "Error: File is shorter than before: %s\n", files[i].c_str());
		exit(1);
	}
	close_compressed(fd);
}

const hashl::base_type *hash_metadata::read_data(size_t &data_size) {
	// length of all stored sequence
	size_t sequence_size(0);
	for (size_t i(0); i < files.size(); ++i) {
		for (size_t j(0); j < reads[i].size(); ++j) {
			for (size_t k(0); k < read_ranges[i][j].size(); ++k) {
				sequence_size += read_ranges[i][j][k].second - read_ranges[i][j][k].first;
			}
		}
	}
	// convert size from basepairs to length of hashl::base_type array
	data_size = (2 * sequence_size + sizeof(hashl::base_type) * 8 - 1) / (sizeof(hashl::base_type) * 8);
	data = new hashl::base_type[data_size];
	data[0] = 0;
	byte_offset = 0;
	bit_offset = sizeof(hashl::base_type) * 8;
	for (size_t i(0); i < files.size(); ++i) {
		if (opt_feedback) {
			fprintf(stderr, "%lu: Reading in %s\n", time(0), files[i].c_str());
		}
		read_file(i);
	}
	return data;
}

std::pair<size_t, size_t> hash_metadata::total_reads() const {
	size_t read_count(0), subread_count(0);
	for (size_t i(0); i < files.size(); ++i) {
		read_count += reads[i].size();
		for (size_t j(0); j < reads[i].size(); ++j) {
			subread_count += read_ranges[i][j].size();
		}
	}
	return std::make_pair(read_count, subread_count);
}

size_t hash_metadata::max_kmers() const {
	size_t x(0);
	for (size_t i(0); i < files.size(); ++i) {
		for (size_t j(0); j < reads[i].size(); ++j) {
			for (size_t k(0); k < read_ranges[i][j].size(); ++k) {
				x += read_ranges[i][j][k].second - read_ranges[i][j][k].first - opt_mer_length + 1;
			}
		}
	}
	return x;
}

std::vector<size_t> hash_metadata::read_ends() const {
	std::vector<size_t> list;
	size_t x(0);
	for (size_t i(0); i < files.size(); ++i) {
		for (size_t j(0); j < reads[i].size(); ++j) {
			for (size_t k(0); k < read_ranges[i][j].size(); ++k) {
				x += read_ranges[i][j][k].second - read_ranges[i][j][k].first;
				list.push_back(x);
			}
		}
	}
	return list;
}

void hash_metadata::add_file(const std::string &file_name) {
	finalize();
	files.push_back(file_name);
	reads.push_back(std::vector<std::string>());
	read_ranges.push_back(std::vector<std::vector<std::pair<uint64_t, uint64_t> > >());
}

void hash_metadata::add_read(const std::string &read_name) {
	reads.back().push_back(read_name);
	read_ranges.back().push_back(std::vector<std::pair<uint64_t, uint64_t> >());
}

void hash_metadata::add_read_range(const uint64_t start, const uint64_t end) {
	read_ranges.back().back().push_back(std::make_pair(start, end));
}

void hash_metadata::finalize() {
	if (!files.empty()) {
		// check if last read of last file had any ranges
		if (read_ranges.back().back().empty()) {
			read_ranges.back().pop_back();
			reads.back().pop_back();
		}
		// check if last file had any reads
		if (read_ranges.back().empty()) {
			read_ranges.pop_back();
			reads.pop_back();
			files.pop_back();
		}
	}
}

const void *hash_metadata::pack(size_t &metadata_size) const {
	// count space needed
	metadata_size = sizeof(uint64_t);				// number of files
	for (size_t i(0); i < files.size(); ++i) {
		metadata_size += files[i].size() + 1;			// file name and null
		metadata_size += sizeof(uint64_t);			// number of reads
		for (size_t j(0); j < reads[i].size(); ++j) {
			metadata_size += reads[i][j].size() + 1;	// read name and null
			metadata_size += sizeof(uint64_t);		// number of ranges
			metadata_size += read_ranges[i][j].size() * sizeof(uint64_t) * 2;
		}
	}
	// allocate space
	char * const d(new char[metadata_size]);
	size_t offset(0);
	uint64_t tmp;
	// fill space with metadata
	memcpy(d + offset, &(tmp = files.size()), sizeof(tmp));
	offset += sizeof(tmp);
	for (size_t i(0); i < files.size(); ++i) {
		memcpy(d + offset, &files[i].c_str()[0], files[i].size() + 1);
		offset += files[i].size() + 1;
		memcpy(d + offset, &(tmp = reads.size()), sizeof(tmp));
		offset += sizeof(tmp);
		for (size_t j(0); j < reads[i].size(); ++j) {
			memcpy(d + offset, &reads[i][j].c_str()[0], reads[i][j].size() + 1);
			offset += reads[i][j].size() + 1;
			memcpy(d + offset, &(tmp = read_ranges[i][j].size()), sizeof(tmp));
			offset += sizeof(tmp);
			for (size_t k(0); k < read_ranges[i][j].size(); ++k) {
				memcpy(d + offset, &read_ranges[i][j][k].first, sizeof(uint64_t));
				offset += sizeof(uint64_t);
				memcpy(d + offset, &read_ranges[i][j][k].second, sizeof(uint64_t));
				offset += sizeof(uint64_t);
			}
		}
	}
	return d;
}

void hash_metadata::unpack(const char *d) {
	uint64_t file_count;
	memcpy(&file_count, d, sizeof(file_count));
	d += sizeof(file_count);
	files.clear();
	files.reserve(file_count);
	reads.assign(file_count, std::vector<std::string>());
	read_ranges.assign(file_count, std::vector<std::vector<std::pair<uint64_t, uint64_t> > >());
	for (size_t i(0); i < file_count; ++i) {
		files.push_back(d);
		d += files.back().size() + 1;
		uint64_t read_count;
		memcpy(&read_count, d, sizeof(read_count));
		d += sizeof(read_count);
		reads[i].reserve(read_count);
		read_ranges[i].assign(read_count, std::vector<std::pair<uint64_t, uint64_t> >());
		for (size_t j(0); j < read_count; ++j) {
			reads[i].push_back(d);
			d += reads[i].back().size() + 1;
			uint64_t read_range_count;
			memcpy(&read_range_count, d, sizeof(read_range_count));
			d += sizeof(read_range_count);
			read_ranges[i][j].reserve(read_range_count);
			for (size_t k(0); k < read_range_count; ++k) {
				uint64_t start, end;
				memcpy(&start, d, sizeof(start));
				d += sizeof(start);
				memcpy(&end, d, sizeof(end));
				d += sizeof(end);
				read_ranges[i][j].push_back(std::make_pair(start, end));
			}
		}
	}
}

void hash_metadata::print() const {
	for (size_t i(0); i < files.size(); ++i) {
		printf("%s\n", files[i].c_str());
		for (size_t j(0); j < reads[i].size(); ++j) {
			printf("\t%s\n", reads[i][j].c_str());
			for (size_t k(0); k < read_ranges[i][j].size(); ++k) {
				printf("\t\t%lu %lu\n", read_ranges[i][j][k].first, read_ranges[i][j][k].second);
			}
		}
	}
}

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

// convert key to sequence

static std::string convert_key(const hashl::key_type &key) {
	const char values[4] = { 'A', 'C', 'G', 'T' };
	std::string sequence;
	sequence.reserve(opt_mer_length);
	// relies on wrap-around for termination
	for (size_t i(opt_mer_length * 2 - 2); i < opt_mer_length * 2; i -= 2) {
		sequence += values[key.basepair(i)];
	}
	return sequence;
}

static void reverse_key(const hashl::key_type &key_in, hashl::key_type &key_out) {
	for (size_t i(0); i < opt_mer_length * 2; i += 2) {
		key_out.push_back(3 - key_in.basepair(i));
	}
}

void print_final_input_feedback(const hashl &mer_list) {
	if (opt_feedback && mer_list.size() != 0) {
		fprintf(stderr, "%lu: %10lu entries used (%5.2f%%), %lu overflow\n", time(0), mer_list.size(), double(100) * mer_list.size() / mer_list.capacity(), mer_list.overflow_size());
	}
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
		++counts[a.value];
		if (opt_print_gc) {
			a.get_key(key);
			gc_counts[a.value] += count_gc(key);
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
		"    -l ## lower bound for hash fill fraction\n"
		"    -L ## upper bound for hash fill fraction\n"
		"    -m ## set mer length [24]\n"
		"    -o ## print output to file instead of stdout\n"
		"    -s ## save histogram memory structure to file\n"
		"    -S ## load histogram memory dump from given file\n"
		"    -V    print version\n"
		"    -w ## print frequency count instead of histogram, for all n-mers with\n"
		"          a frequency of at least ## [0 (off)]\n"
		"    -z ## number of possible n-mers to allocate memory for (overrides -l/-L)\n"
		"          (k, m, or g may be suffixed)\n"
	);
	exit(1);
}

static void get_opts(int argc, char **argv) {
	std::set<std::string> used_files;	// make sure there're no overlaps
	std::string opt_output;
	opt_feedback = 1;
	opt_frequency_cutoff = 0;
	opt_histogram_restore = -1;
	opt_mer_length = 24;
	opt_print_gc = 0;
	opt_load_lower_bound = 0;
	opt_load_upper_bound = 1;
	opt_nmers = 0;
	int c;
	while ((c = getopt(argc, argv, "ghil:L:m:o:s:S:Vw:z:")) != EOF) {
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
			if (!used_files.insert(optarg).second) {
				fprintf(stderr, "Error: duplicate file: %s\n", optarg);
				exit(1);
			}
			break;
		    case 's':
			opt_save_file = optarg;
			if (!used_files.insert(optarg).second) {
				fprintf(stderr, "Error: duplicate file: %s\n", optarg);
				exit(1);
			}
			break;
		    case 'S':
			opt_histogram_restore = open_compressed(optarg);
			if (opt_histogram_restore == -1) {
				fprintf(stderr, "Error: could not read histogram dump file\n");
				exit(1);
			} else if (!used_files.insert(optarg).second) {
				fprintf(stderr, "Error: duplicate file: %s\n", optarg);
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
		    case 'z':
			opt_nmers = get_value(optarg);
			break;
		    default:
			fprintf(stderr, "Error: unknown option %c\n", c);
			print_usage();
		}
	}
	if (opt_load_lower_bound > opt_load_upper_bound) {
		fprintf(stderr, "Error: lower hash fill fraction must be less than upper\n");
		print_usage();
	}
	if (optind == argc && opt_histogram_restore == -1) {
		fprintf(stderr, "Error: no files to process\n");
		print_usage();
	}
	for (int i(optind); i < argc; ++i) {
		if (!used_files.insert(argv[i]).second) {
			fprintf(stderr, "Error: duplicate file: %s\n", argv[i]);
			exit(1);
		}
	}
	if (!opt_output.empty()) {
		fp_out = fopen(opt_output.c_str(), "w");
		if (!fp_out) {
			fprintf(stderr, "Error: could not write to %s\n", opt_output.c_str());
			exit(1);
		}
	}
}

// for non-ACGT basepairs, need to split reads

static void get_subread_sizes(const std::string &seq, hash_metadata &metadata) {
	size_t i(seq.find_first_of("ACGTacgt", 0));
	while (i != std::string::npos) {
		size_t next(seq.find_first_not_of("ACGTacgt", i));
		if (next == std::string::npos) {
			next = seq.size();
		}
		// reads shorter than the mer length are skipped
		if (next - i >= opt_mer_length) {
			metadata.add_read_range(i, next);
		}
		i = seq.find_first_of("ACGTacgt", next);
	}
}

// read file and get total space needed for data and metadata
static void get_read_sizes(const char * const file, hash_metadata &metadata) {
	const int fd(open_compressed(file));
	if (fd == -1) {
		fprintf(stderr, "Error: open: %s\n", file);
		exit(1);
	}
	std::string line, seq;
	if (pfgets(fd, line) == -1) {		// empty file
		fprintf(stderr, "Warning: empty file: %s\n", file);
	} else if (line[0] == '>') {		// fasta file
		do {
			size_t i(1);
			for (; i < line.size() && !isspace(line[i]); ++i) { }
			metadata.add_read(line.substr(1, i - 1));
			seq.clear();
			while (pfgets(fd, line) != -1 && line[0] != '>') {
				seq += line;
			}
			get_subread_sizes(seq, metadata);
		} while (!line.empty());
	} else if (line[0] == '@') {		// fastq file
		do {
			size_t i(1);
			for (; i < line.size() && !isspace(line[i]); ++i) { }
			metadata.add_read(line.substr(1, i - 1));
			if (pfgets(fd, seq) == -1) {		// sequence
				fprintf(stderr, "Error: truncated fastq file: %s\n", file);
				exit(1);
			}
			get_subread_sizes(seq, metadata);
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

static void count_nmers(hashl &mer_list, const hashl::base_type * const data, const std::vector<size_t> &read_ends) {
	hashl::key_type key(mer_list), comp_key(mer_list);
	size_t total_reads(0);
	size_t i(0), j(0), k(sizeof(hashl::base_type) * 8 - 2);
	// iterate over all reads (nmers can't cross read boundaries)
	std::vector<size_t>::const_iterator a(read_ends.begin());
	const std::vector<size_t>::const_iterator end_a(read_ends.end());
	for (; a != end_a; ++a, ++total_reads) {
		// print feedback every 10 minutes
		if (opt_feedback && elapsed_time() >= 600) {
			start_time();
			fprintf(stderr, "%lu: %10lu entries used (%5.2f%%), %lu overflow (%lu reads)\n", time(0), mer_list.size(), double(100) * mer_list.size() / mer_list.capacity(), mer_list.overflow_size(), total_reads);
		}
		const size_t end_i(i + opt_mer_length - 1);
		// load keys with opt_mer_length - 1 basepairs
		for (; i < end_i; ++i) {
			const hashl::base_type c((data[j] >> k) & 3);
			key.push_back(c);
			comp_key.push_front(3 - c);
			if (k) {
				k -= 2;
			} else {
				k = sizeof(hashl::base_type) * 8 - 2;
				++j;
			}
		}
		// run over all nmers, one basepair at a time
		for (; i < *a; ++i) {
			const hashl::base_type c((data[j] >> k) & 3);
			key.push_back(c);
			comp_key.push_front(3 - c);
			if (k) {
				k -= 2;
			} else {
				k = sizeof(hashl::base_type) * 8 - 2;
				++j;
			}
			if (!mer_list.increment(key, comp_key, 2 * (i + 1 - opt_mer_length))) {
				fprintf(stderr, "Error: ran out of space in hash\n");
				exit(1);
			}
		}
	}
}

// metadata format: # of files, [ filename, # of reads, read offsets, read names ]
// strings are null delimited, # and read lengths are uint64_t
// hashl boilerplate check will ensure byte order is the same read as written

static void read_in_files(int argc, char **argv, hashl &mer_list) {
	const uint64_t file_count(argc - optind);
	hash_metadata metadata;
	for (size_t i(0); i < file_count; ++i) {
		if (opt_feedback) {
			fprintf(stderr, "%lu: Getting read sizes for %s\n", time(0), argv[i + optind]);
		}
		metadata.add_file(argv[i + optind]);
		get_read_sizes(argv[i + optind], metadata);
	}
	size_t data_size;
	const hashl::base_type * const data(metadata.read_data(data_size));
	if (opt_feedback) {
		fprintf(stderr, "%lu: Initializing n-mer hash\n", time(0));
	}
	// put data and metadata into mer_list
	mer_list.init(opt_nmers ? opt_nmers : metadata.max_kmers(), opt_mer_length * 2, data, data_size);
	size_t metadata_size;
	const void * const md(metadata.pack(metadata_size));
	mer_list.set_metadata(md, metadata_size);
	if (opt_feedback) {
		std::pair<size_t, size_t> read_count(metadata.total_reads());
		fprintf(stderr, "%lu: Counting n-mers for %lu reads (%lu ranges)\n", time(0), read_count.first, read_count.second);
		start_time();
	}
	const std::vector<size_t> read_ends(metadata.read_ends());
	count_nmers(mer_list, data, read_ends);
	if (opt_feedback) {
		print_final_input_feedback(mer_list);
	}
	// make sure the fill rate is acceptable
	const double load(double(mer_list.size()) / mer_list.capacity());
	if (load < opt_load_lower_bound || opt_load_upper_bound < load) {
		if (opt_feedback) {
			fprintf(stderr, "%lu: hash fill rate out of range, rehashing: %f - %f: %f\n", time(0), opt_load_lower_bound, opt_load_upper_bound, load);
		}
		mer_list.resize(mer_list.size() * 2 / (opt_load_lower_bound + opt_load_upper_bound));
		if (opt_feedback) {
			fprintf(stderr, "%lu: Recounting n-mers\n", time(0));
			start_time();
		}
		count_nmers(mer_list, data, read_ends);
		if (opt_feedback) {
			print_final_input_feedback(mer_list);
		}
	}
}

int main(int argc, char **argv) {
	get_opts(argc, argv);
	hashl mer_list;
	if (opt_histogram_restore != -1) {
		if (opt_feedback) {
			fprintf(stderr, "%lu: Initializing n-mer hash\n", time(0));
		}
		mer_list.init_from_file(opt_histogram_restore);
		if (opt_feedback) {
			print_final_input_feedback(mer_list);
		}
	}
	if (optind != argc) {
		read_in_files(argc, argv, mer_list);
	}
	if (opt_feedback) {
		fprintf(stderr, "%lu: Printing results\n", time(0));
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
	return 0;
}
