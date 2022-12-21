#include "hashl.h"	// hashl
#include "open_compressed.h"	// close_compressed(), get_suffix(), open_compressed(), pfgets()
#include "time_used.h"	// elapsed_time(), start_time()
#include "version.h"	// VERSION
#include "write_fork.h"	// close_fork(), write_fork()
#include <getopt.h>	// getopt(), optarg, optind
#include <list>		// list<>
#include <map>		// map<>
#include <sstream>	// istringstream, ostringstream
#include <stdint.h>	// uint64_t
#include <stdio.h>	// EOF, fprintf(), stderr, stdout
#include <stdlib.h>	// exit()
#include <string.h>	// memcpy(), strlen()
#include <string>	// string
#include <utility>	// make_pair(), pair<>
#include <vector>	// vector<>
#include <deque>	// deque<>

class hash_data {
    private:
	std::vector<std::string> files;
	std::vector<std::vector<std::string> > reads;
	// inclusive start, exclusive end
	std::vector<std::vector<std::vector<std::pair<uint64_t, uint64_t> > > > read_ranges;
	std::deque<hashl::base_type> data;
    public:
	hash_data(void) { }
	~hash_data(void) { }
	void read_file(const std::string &);
	const void *pack_metadata(size_t &) const;
	void unpack_metadata(const char *);
	void print_metadata(void) const;
    private:
	void read_file(size_t, size_t &);
	void get_subreads(const std::vector<std::pair<uint64_t, uint64_t> > &, const std::string &, size_t &);
	void copy_data(const std::string &, size_t, size_t, size_t);
	static hashl::base_type hash_data::convert_char(char) const;
};

// return 0-3, exit for bad sequence
hashl::base_type hash_data::convert_char(const char c) {
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

// sequence is stored from top bits to bottom bits of each offset
void hash_data::copy_data(const std::string &seq, size_t k, const size_t end_k, size_t data_offset) {
	size_t i((2 * data_offset) / (sizeof(hashl::base_type) * 8));
	size_t j(sizeof(hashl::base_type) * 8 - (2 * data_offset) % (sizeof(hashl::base_type) * 8));
	for (; k < end_k; ++k) {
		if (j) {
			j -= 2;
			data[i] |= convert_char(seq[k]) << j;
		} else {
			j = sizeof(hashl::base_type) * 8 - 2;
			data[++i] = convert_char(seq[k]) << j;
		}
	}
}

void hash_data::get_subreads(const std::vector<std::pair<uint64_t, uint64_t> > &ranges, const std::string &seq, size_t &data_offset) {
	for (size_t k(0); k < ranges.size(); ++k) {
		copy_data(seq, ranges[k].first, ranges[k].second, data_offset);
		data_offset += ranges[k].second - ranges[k].first;
	}
}

void hash_data::read_file(const size_t i, size_t &data_offset) {
	const int fd(open_compressed(files[i]));
	if (fd == -1) {
		fprintf(stderr, "Error: open: %s\n", files[i].c_str());
		exit(1);
	}
	size_t j(0);
	std::string line, seq;
	if (pfgets(fd, line) == -1) {		// empty file
		fprintf(stderr, "Warning: empty file: %s\n", files[i].c_str());
	} else if (line[0] == '>') {		// fasta file
		std::string header;
		do {
			const std::string &read_name(reads[i][j]);
			if (read_name.compare(0, read_name.size(), line, 1, read_name.size()) && (line.size() == read_name.size() + 1 || isspace(line[read_name.size() + 1]))) {
				seq.clear();
				while (pfgets(fd, line) != -1 && line[0] != '>') {
					seq += line;
				}
				get_subreads(read_ranges[i][j], seq, data_offset);
				++j;
			} else {
				while (pfgets(fd, line) != -1 && line[0] != '>') { }
			}
		} while (line[0] == '>');	// safe after eof
	} else if (line[0] == '@') {		// fastq file
		do {
			if (pfgets(fd, seq) == -1) {		// sequence
				fprintf(stderr, "Error: truncated fastq file: %s\n", files[i].c_str());
				exit(1);
			}
			const std::string &read_name(reads[i][j]);
			if (read_name.compare(0, read_name.size(), line, 1, read_name.size()) && (line.size() == read_name.size() + 1 || isspace(line[read_name.size() + 1]))) {
				get_subreads(read_ranges[i][j], seq, data_offset);
				++j;
			}
			// skip quality header and quality
			// (use seq because it'll be the same length as quality)
			if (pfgets(fd, line) == -1 || pfgets(fd, seq) == -1) {
				fprintf(stderr, "Error: truncated fastq file: %s\n", files[i].c_str());
				exit(1);
			}
		} while (pfgets(fd, line) != -1);
	} else {
		fprintf(stderr, "Error: unknown file format: %s\n", files[i].c_str());
		exit(1);
	}
	close_compressed(fd);
}

void hash_data::read_data() {
	// count space needed for data
	size_t data_size(0);
	for (size_t i(0); i < files.size(); ++i) {
		for (size_t j(0); j < reads[i].size(); ++j) {
			for (size_t k(0); k < read_ranges[i][j].size(); ++k) {
				data_size += read_ranges[i][j][k].second - read_ranges[i][j][k].first;
			}
		}
	}
	// convert data size from basepairs to length of hashl::base_type array
	data_size = (2 * data_size + sizeof(hashl::base_type) * 8 - 1) / (sizeof(hashl::base_type) * 8);
	data.assign(data_size, 0);
	size_t data_offset(0);
	for (size_t i(0); i < files.size(); ++i) {
		read_file(i, data_offset);
	}
}

void hash_data::add_file(const std::string &file_name) {
	finalize();
	files.push_back(file_name);
	reads.push_back(std::vector<std::string>());
	read_ranges.push_back(std::vector<std::vector<std::pair<uint64_t, uint64_t> > >());
}

void hash_data::add_read(const std::string &read_name) {
	reads.last().push_back(read_name);
	read_ranges.last().push_back(std::vector<std::pair<uint64_t, uint64_t> >());
}

void hash_data::add_read_range(const uint64_t start, const uint64_t end) {
	read_ranges.last().last().push_back(std::make_pair(start, end));
}

void hash_data::finalize() {
	if (!files.empty()) {
		// check if last read of last file had any ranges
		if (read_ranges.last().last().empty()) {
			read_ranges.last().pop_back();
			reads.last().pop_back();
		}
		// check if last file had any reads
		if (read_ranges.last().empty()) {
			read_ranges.pop_back();
			reads.pop_back();
			files.pop_back();
		}
	}
}

void hash_data::print_metadata() const {
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

const void *hash_data::pack_metadata() const {
	// count space needed
	size_t total_size(sizeof(uint64_t));			// number of files
	for (size_t i(0); i < files.size(); ++i) {
		total_size += files[i].size() + 1;		// file name and null
		total_size += sizeof(uint64_t);			// number of reads
		for (size_t j(0); j < reads[i].size(); ++j) {
			total_size += reads[i][j].size() + 1; 	// read name and null
			total_size += sizeof(uint64_t);		// number of ranges
			total_size += read_ranges[i][j].size() * sizeof(uint64_t) * 2;
		}
	}
	// allocate space
	char *d(new char[total_size]);
	size_t offset(0);
	uint64_t tmp;
	// fill space with metadata
	mempcy(d + offset, &(tmp = files.size()), sizeof(tmp));
	offset += sizeof(tmp);
	for (size_t i(0); i < files.size(); ++i) {
		mempcy(d + offset, &files[i].c_str()[0], files[i].size() + 1);
		offset += files[i].size() + 1;
		mempcy(d + offset, &(tmp = reads.size()), sizeof(tmp));
		offset += sizeof(tmp);
		for (size_t j(0); j < reads[i].size(); ++j) {
			mempcy(d + offset, &reads[i][j].c_str()[0], reads[i][j].size() + 1);
			offset += reads[i][j].size() + 1;
			mempcy(d + offset, &(tmp = read_ranges[i][j].size()), sizeof(tmp));
			offset += sizeof(tmp);
			for (size_t k(0); k < read_ranges[i][j].size(); ++k) {
				mempcy(d + offset, &read_ranges[i][j][k].first, sizeof(uint64_t));
				offset += sizeof(uint64_t);
				mempcy(d + offset, &read_ranges[i][j][k].second, sizeof(uint64_t));
				offset += sizeof(uint64_t);
			}
		}
	}
	return d;
}

void hash_data::unpack_metadata(const char * d) {
	uint64_t file_count;
	mempcy(&file_count, d, sizeof(file_count));
	d += sizeof(file_count);
	files.clear();
	files.reserve(file_count);
	reads.assign(file_count, std::vector<std::string>());
	read_ranges.assign(file_count, std::vector<std::vector<std::pair<uint64_t, uint64_t> > >());
	for (size_t i(0); i < file_count; ++i) {
		files.push_back(d);
		d += files.last().size() + 1;
		uint64_t read_count;
		mempcy(&read_count, d, sizeof(read_count));
		d += sizeof(read_count);
		reads[i].reserve(read_count);
		read_ranges[i].assign(read_count, std::vector<std::pair<uint64_t, uint64_t> >());
		for (size_t j(0); j < read_count; ++j) {
			reads[i].push_back(d);
			d += reads[i].last().size() + 1;
			uint64_t read_range_count;
			mempcy(&read_range_count, d, sizeof(read_range_count));
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

void hash_data::print() const {
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
			if (!used_files.insert(optarg)->second) {
				fprintf(stderr, "Error: duplicate file: %s\n", optarg);
				exit(1);
			}
			break;
		    case 's':
			opt_save_file = optarg;
			if (!used_files.insert(optarg)->second) {
				fprintf(stderr, "Error: duplicate file: %s\n", optarg);
				exit(1);
			}
			break;
		    case 'S':
			opt_histogram_restore = open_compressed(optarg);
			if (opt_histogram_restore == -1) {
				fprintf(stderr, "Error: could not read histogram dump file\n");
				exit(1);
			} else if (!used_files.insert(optarg)->second) {
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
		if (!used_files.insert(argv[i])->second) {
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

static void get_subread_sizes(const std::string &seq, size_t &metadata_size, size_t &data_size, size_t &max_kmers, std::vector<size_t> &read_ends, const size_t header_size) {
	size_t i(seq.find_first_of("ACGTacgt", 0));
	while (i != std::string::npos) {
		size_t next(seq.find_first_not_of("ACGTacgt", i));
		if (next == std::string::npos) {
			next = seq.size();
		}
		// reads shorter than the mer length are skipped
		if (next - i >= opt_mer_length) {
			// make room for offset if we're splitting the read
			if (i || next < seq.size()) {
				std::ostringstream x;
				x << i;		// convert offset to text
				// +1 for separator
				metadata_size += header_size + 1 + x.str().size();
			} else {
				// don't need the offset
				metadata_size += header_size;
			}
			data_size += next - i;
			max_kmers += next - i - opt_mer_length + 1;
			read_ends.push_back(data_size);
		}
		i = seq.find_first_of("ACGTacgt", next);
	}
}

// read file and get total space needed for data and metadata
void hash_data::get_read_sizes(hash_data &file_data) {
	const int fd(open_compressed(files.last()));
	if (fd == -1) {
		fprintf(stderr, "Error: open: %s\n", files.last().c_str());
		exit(1);
	}
	std::string line, seq;
	if (pfgets(fd, line) == -1) {		// empty file
		fprintf(stderr, "Warning: empty file: %s\n", file);
	} else if (line[0] == '>') {		// fasta file
		do {
			size_t i(1);
			for (; i < line.size() && !isspace(line[i]); ++i) { }
			seq.clear();
			while (pfgets(fd, line) != -1 && line[0] != '>') {
				seq += line;
			}
			get_subread_sizes(seq);
		} while (!line.empty());
	} else if (line[0] == '@') {		// fastq file
		do {
			size_t i(1);
			for (; i < line.size() && !isspace(line[i]); ++i) { }
			if (pfgets(fd, seq) == -1) {		// sequence
				fprintf(stderr, "Error: truncated fastq file: %s\n", file);
				exit(1);
			}
			get_subread_sizes(seq, metadata_size, data_size, max_kmers, read_ends, i);
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
	if (!read_ends.empty()) {	// only count if we're not going to skip the file
		metadata_size += strlen(file) + 1 + read_ends.size() * sizeof(uint64_t);
	}
}

static void count_nmers(hashl &mer_list, const hashl::base_type * const data, const std::vector<std::vector<uint64_t> > &read_ends) {
	hashl::key_type key(mer_list), comp_key(mer_list);
	uint64_t i(0);		// match type of read_ends
	size_t total_reads(0);
	// iterate over all files
	std::vector<std::vector<uint64_t> >::const_iterator a(read_ends.begin());
	const std::vector<std::vector<uint64_t> >::const_iterator end_a(read_ends.end());
	for (; a != end_a; ++a) {
		// iterate over all reads in file (nmers can't cross read boundaries)
		std::vector<uint64_t>::const_iterator b(a->begin());
		const std::vector<uint64_t>::const_iterator end_b(a->end());
		for (; b != end_b; ++b, ++total_reads) {
			// print feedback every 10 minutes
			if (opt_feedback && elapsed_time() >= 600) {
				start_time();
				fprintf(stderr, "%lu: %10lu entries used (%5.2f%%), %lu overflow (%lu reads)\n", time(0), mer_list.size(), double(100) * mer_list.size() / mer_list.capacity(), mer_list.overflow_size(), total_reads);
			}
			uint64_t j((2 * i) / (sizeof(hashl::base_type) * 8));
			uint64_t k(sizeof(hashl::base_type) * 8 - (2 * i) % (sizeof(hashl::base_type) * 8) - 2);
			const uint64_t end_i(i + opt_mer_length - 1);
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
			for (; i < *b; ++i) {
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
}

// XXX - this is a horrible mess and needs to be at least partially reverted
static void read_files(int argc, char **argv, hashl &mer_list) {
	hash_data file_data;
	const uint64_t file_count(argc - optind);
	for (size_t i(0); i < file_count; ++i) {
		if (opt_feedback) {
			fprintf(stderr, "%lu: Reading %s\n", time(0), argv[i + optind]);
		}
		file_data.read_file(argv[i + optind]);
	}
	if (opt_feedback) {
		fprintf(stderr, "%lu: Initializing n-mer hash\n", time(0));
	}
	// put data and metadata into mer_list
	mer_list.init(opt_nmers ? opt_nmers : file_data.max_kmers(), opt_mer_length * 2, file_data.data());
	size_t metadata_size;
	const void *metadata(file_data.pack(metadata_size));
	mer_list.set_metadata(metadata, metadata_size);
	if (opt_feedback) {
		fprintf(stderr, "%lu: Counting n-mers for %lu reads\n", time(0), total_reads);
		start_time();
	}
	count_nmers(mer_list, file_data);
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
		count_nmers(mer_list, file_data);
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
		read_files(argc, argv, mer_list);
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
