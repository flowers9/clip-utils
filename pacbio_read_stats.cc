#include "open_compressed.h"	// close_compressed(), find_suffix(), open_compressed(), pfgets()
#include "pretty_print.h"	// pretty_print()
#include "version.h"	// VERSION
#include <algorithm>	// max(), sort()
#include <ctype.h>	// isdigit(), isspace()
#include <exception>	// exception
#include <getopt.h>	// getopt(), optind
#include <iomanip>	// setfill()
#include <iostream>	// cerr, cout
#include <list>		// list<>
#include <map>		// map<>
#include <sstream>	// istringstream, ostringstream
#include <stdio.h>	// EOF
#include <stdlib.h>	// exit()
#include <string>	// string
#include <sys/stat.h>	// stat(), struct stat
#include <sys/types.h>	// size_t
#include <vector>	// vector<>

class LocalException : public std::exception {
    private:
	const std::string error_;
	const int show_usage_;
    public:
	explicit LocalException(const std::string &s) : error_(s), show_usage_(0) { }
	LocalException(const std::string &s, int i) : error_(s), show_usage_(i) { }
	virtual ~LocalException(void) throw() { }
	virtual const char * what() const throw() {
		return error_.c_str();
	}
	const int &show_usage(void) const throw() {
		return show_usage_;
	}
};

static int opt_all_read_segments(0);
static int opt_full_histogram(0);
static int opt_strip_trace(0);
static size_t (* get_id_start)(std::string &line);
static size_t opt_genome_size(0);
static size_t opt_max_cutoff(0);

// this class holds and processes the current read data

class CurrentState {
    private:
	std::string seq_id_, last_pacbio_id_;
	size_t seq_length_, best_fragment_length_;
	size_t seq_id_offset_;
	std::vector<size_t> read_lengths_;
    private:
	std::string get_pacbio_read(void) const;
    public:
	CurrentState(void) : seq_length_(0), best_fragment_length_(0), seq_id_offset_(0) { }
	~CurrentState(void) { }
	void save_length(void);
	void flush_seq(void) {
		save_length();
		if (best_fragment_length_) {
			read_lengths_.push_back(best_fragment_length_);
			best_fragment_length_ = 0;
		}
		seq_id_.clear();
		last_pacbio_id_.clear();
	}
	void set_seq(const std::string &id, const size_t &length, const size_t &id_offset) {
		seq_id_ = id;
		seq_length_ = length;
		seq_id_offset_ = id_offset;
	}
	std::vector<size_t> &read_lengths(void) {
		return read_lengths_;
	}
};

static CurrentState current;

// parse out pacbio read name from full readname (part up to second /)
// m150314_203941_00123_c100788912550000001823159708251505_s1_p0/14/7926_10200

std::string CurrentState::get_pacbio_read(void) const {
	size_t i(seq_id_offset_);
	const size_t end_i(seq_id_.size());
	for (; i != end_i && seq_id_[i] != '/'; ++i) { }
	if (i != end_i) {
		for (++i; i != end_i && seq_id_[i] != '/'; ++i) { }
	}
	return seq_id_.substr(seq_id_offset_, i - seq_id_offset_);
}

void CurrentState::save_length(void) {
	if (opt_all_read_segments) {
		read_lengths_.push_back(seq_length_);
	} else if (last_pacbio_id_.empty()) {
		last_pacbio_id_ = get_pacbio_read();
		best_fragment_length_ = seq_length_;
	} else {
		std::string s(get_pacbio_read());
		if (s != last_pacbio_id_) {
			read_lengths_.push_back(best_fragment_length_);
			last_pacbio_id_ = s;
			best_fragment_length_ = seq_length_;
		} else if (best_fragment_length_ < seq_length_) {
			best_fragment_length_ = seq_length_;
		}
	}
}

// line is never modified, but needs to be non-const to match
// get_id_start_convert(), which does modify line

static size_t get_id_start_noconvert(std::string &line) {
	if (opt_strip_trace) {
		size_t i(2);
		const size_t end_i(line.size());
		for (; i != end_i && !isspace(line[i]); ++i) { }
		if (i != end_i) {
			for (++i; i != end_i && isspace(line[i]); ++i) { }
			if (i != end_i) {
				return i;
			}
		}
		return 0;
	}
	return 1;
}

// like get_id_start_noconvert(), but also converts the read name from
// FLHBFN1:224:D0WL7ACXX:1:1101:5802:2244 1:N:0:GCCAAT format to
// FLHBFN1_224_D0WL7ACXX_1_1101_5802_2244-R1; will also remove the
// trace from the header line if opt_strip_trace is given

static size_t get_id_start_convert(std::string &line) {
	size_t i(2);
	if (opt_strip_trace) {				// remove trace
		const size_t end_i(line.size());
		for (; i != end_i && !isspace(line[i]); ++i) { }
		if (i == end_i) {
			return 0;
		}
		for (++i; i != end_i && isspace(line[i]); ++i) { }
		if (i == end_i) {
			return 0;
		}
		line.erase(1, i - 1);
	}
	// rewrite header
	const size_t end_i(line.size());
	for (i = 2; i != end_i; ++i) {
		if (line[i] == ':' || line[i] == '-') {
			line[i] = '_';
		} else if (isspace(line[i])) {
			break;
		}
	}
	if (i == end_i || i + 1 == end_i) {
		return 0;
	}
	const char c(line[i + 1]);
	if (c != '1' && c != '2') {
		return 0;
	}
	line.resize(i);
	line += "-R";	// need to do this in two steps;
	line += c;	// curse you, string conversions!
	return 1;
}

static void print_usage() {
	std::cerr <<
		"usage: pacbio_read_stats [opts] <fasta> [<fasta2> ...]\n" <<
		"    -a    use all read segments (not just best one)\n" <<
		"    -f    print full histogram (not binned histogram)\n" <<
		"    -h    print usage\n" <<
		"    -m ## upper limit for cutoffs, in kb [60]\n" <<
		"    -R    convert readnames from new Illumina form to old\n" <<
		"    -s ## estimated genome size (in MB)\n" <<
		"    -t    strip first part of trace id from read headers\n" <<
		"    -V    print version\n";
}

static int get_opts(int argc, char **argv) {
	int opt_convert_readnames(0);
	opt_max_cutoff = 60;
	// set option defaults
	opt_strip_trace = 0;
	int c;
	while ((c = getopt(argc, argv, "afhm:Rs:tV")) != EOF) {
		switch (c) {
		    case 'a':
			opt_all_read_segments = 1;
			break;
		    case 'f':
			opt_full_histogram = 1;
			break;
		    case 'h':
			print_usage();
			return 1;
		    case 'm':
			std::istringstream(optarg) >> opt_max_cutoff;
			break;
		    case 'R':
			opt_convert_readnames = 1;
			break;
		    case 't':
			opt_strip_trace = 1;
			break;
		    case 's':
			std::istringstream(optarg) >> opt_genome_size;
			break;
		    case 'V':
			std::cerr << "pacbio_read_stats version " << VERSION << "\n";
			exit(0);
			break;
		    default:
			throw LocalException("bad option: " + (char)c, 1);
		}
	}
	if (optind == argc) {
		throw LocalException("no files specified", 1);
	}
	get_id_start = opt_convert_readnames ? get_id_start_convert : get_id_start_noconvert;
	opt_max_cutoff *= 1000;
	return 0;
}

// go through arguments to get input files

static void find_files(int argc, char **argv, std::list<std::string> &file_list) {
	for (; optind < argc; ++optind) {
		std::string fasta(argv[optind]), qual;
		// add suffix, if needed
		std::string dummy;
		find_suffix(fasta, dummy);
		file_list.push_back(fasta);
	}
}

static size_t get_next_header_fasta(int fd, std::string &line) {
	int file_status(0);
	if (line.empty()) {		// set file_status if it's unclear
		file_status = pfgets(fd, line);
	}
	while (file_status != -1) {
		if (line.size() > 1 && line[0] == '>' && !isspace(line[1])) {
			size_t i = get_id_start(line);
			if (i != 0) {
				return i;
			}
		}
		file_status = pfgets(fd, line);
	}
	return 0;
}

static size_t get_next_header_fastq(int fd, std::string &line) {
	int file_status(0);
	if (line.empty()) {		// set file_status if it's unclear
		file_status = pfgets(fd, line);
	}
	while (file_status != -1) {
		if (line.size() < 2 || line[0] != '@' || isspace(line[1])) {
			throw LocalException("bad fastq file: expecting a @ line: " + line);
		}
		const size_t i(get_id_start(line));
		if (i != 0) {
			return i;
		}
		pfgets(fd, line);			// seq line
		pfgets(fd, line);			// + line
		pfgets(fd, line);			// qual line
		file_status = pfgets(fd, line);		// @ line
	}
	return 0;
}

// process header to handle n expansion (n\d+-\d+:\d+\.(\d+)\.(\d+)$)

static int process_n_header(const std::string &line, size_t i, size_t &length) {
	const size_t end_i(line.size());
	if (i == std::string::npos || i >= end_i) {
		return 0;
	}
	if (line[i] != 'n' || ++i == end_i || !isdigit(line[i])) {
		return 0;
	}
	for (++i; i != end_i && isdigit(line[i]); ++i) { }
	if (i == end_i || line[i] != '-' || ++i == end_i || !isdigit(line[i])) {
		return 0;
	}
	for (++i; i != end_i && isdigit(line[i]); ++i) { }
	if (i == end_i || line[i] != ':' || ++i == end_i || !isdigit(line[i])) {
		return 0;
	}
	for (++i; i != end_i && isdigit(line[i]); ++i) { }
	if (i == end_i || line[i] != '.' || ++i == end_i || !isdigit(line[i])) {
		return 0;
	}
	size_t j(i);
	for (++i; i != end_i && isdigit(line[i]); ++i) { }
	size_t k(i);
	if (i == end_i || line[i] != '.' || ++i == end_i || !isdigit(line[i])) {
		return 0;
	}
	for (++i; i != end_i && isdigit(line[i]); ++i) { }
	if (i != end_i) {
		return 0;
	}
	size_t x;
	std::istringstream(line.substr(j, k - j)) >> x;
	length += x;
	return 1;
}

// line contains the current line read, which should be a header line;
// return 0 if we didn't get a proper read (either error or eof)

static int find_next_fastq(int fd, std::string &id, std::string &line, size_t &length, size_t &id_offset) {
	id_offset = get_next_header_fastq(fd, line);
	if (id_offset == 0) {
		return 0;
	}
	length = 0;
	if (process_n_header(line, id_offset, length)) {
		id.clear();	// recombining an n-split: read had a fake id
	} else {
		id = line;
	}
	for (;;) {
		// fastq format is fairly strict - it comes in 4 line sets,
		// with a @ line, the sequence, a + line, and the quality
		if (pfgets(fd, line) == -1) {	// sequence line
			return 0;
		}
		length += line.size();
		if (pfgets(fd, line) == -1) {	// + line
			return 0;
		}
		if (line.empty() || line[0] != '+') {
			throw LocalException("bad fastq file: expecting a + line: " + line);
		}
		if (pfgets(fd, line) == -1) {	// quality line
			return 0;
		}
		if (pfgets(fd, line) == -1) {	// next header line
			return 1;
		}
		if (line.size() < 2 || line[0] != '@' || isspace(line[1])) {
			throw LocalException("bad fastq file: expecting a @ line: " + line);
		}
		// return unless the next header is a fake one - in that case,
		// continue adding to the current read
		if (!process_n_header(line, opt_strip_trace ? line.find(' ') : 1, length)) {
			return 1;
		}
	}
}

static int find_next_fasta(int fd, std::string &id, std::string &line, size_t &length, size_t &id_offset) {
	id_offset = get_next_header_fasta(fd, line);
	if (id_offset == 0) {
		return 0;
	}
	length = 0;
	if (process_n_header(line, id_offset, length)) {
		id.clear();	// recombining an n-split: read had a fake id
	} else {
		id = line;
	}
	while (pfgets(fd, line) != -1) {
		if (line.empty()) {
		} else if (line[0] != '>') {
			length += line.size();
		} else if (!process_n_header(line, opt_strip_trace ? line.find(' ') : 1, length)) {
			return 1;
		}
	}
	return 1;
}

static void read_reads(const std::string &seq_file) {
	const int fd(open_compressed(seq_file));
	if (fd == -1) {
		throw LocalException("could not open " + seq_file);
	}
	std::string last_id;
	while (last_id.empty()) {
		if (pfgets(fd, last_id) == -1) {	// empty file?!
			break;
		}
	}
	if (last_id.empty()) {			// skip empty file
	} else if (last_id[0] == '@') {		// fastq input format
		std::string id;
		size_t length, id_offset;
		while (find_next_fastq(fd, id, last_id, length, id_offset)) {
			if (!id.empty()) {
				current.save_length();
				current.set_seq(id, length, id_offset);
			}
		}
	} else {				// fasta input format
		std::string id;
		size_t length, id_offset;
		while (find_next_fasta(fd, id, last_id, length, id_offset)) {
			if (!id.empty()) {
				current.save_length();
				current.set_seq(id, length, id_offset);
			}
		}
	}
	close_compressed(fd);
}

class ReadStats {
    public:
	size_t reads;
	unsigned long long basepairs;
	size_t median_read_length;
	ReadStats(const size_t r, const unsigned long long b, const size_t m) : reads(r), basepairs(b), median_read_length(m) { }
	~ReadStats(void) { }
	void update(const ReadStats &a) {
		if (reads < a.reads) {
			reads = a.reads;
		}
		if (basepairs < a.basepairs) {
			basepairs = a.basepairs;
		}
		if (median_read_length < a.median_read_length) {
			median_read_length = a.median_read_length;
		}
	}
};

class FormatOutput {
    private:
	const int width_;
	const int precision_;
    public:
	FormatOutput(const int w, const int p = 0) : width_(w), precision_(p) { }
	~FormatOutput(void) { }
	friend std::ostream &operator<<(std::ostream &out, const FormatOutput &format) {
		out.fill(' ');
		out.width(format.width_);
		if (format.precision_) {
			out.setf(std::ios_base::fixed, std::ios_base::floatfield);
			out.precision(format.precision_);
		}
		return out;
	}
};

static int get_width(unsigned long long x, const int min_width) {
	if (x < 10) {
		return std::max(1, min_width);
	}
	int i(1);
	for (x /= 10; x != 0; x /= 10, ++i) { }
	return std::max(i, min_width);
}

// Calculate, for a given cutoff in kb, how many reads remain,
// and what their total basepairs and average length will be

static void print_histogram(void) {
	std::vector<size_t> &read_lengths(current.read_lengths());
	std::sort(read_lengths.begin(), read_lengths.end());
	ReadStats max_values(0, 0, 0);
	std::vector<ReadStats> histogram;
	size_t j(0);
	const size_t end_j(read_lengths.size());
	for (size_t i(0); i != opt_max_cutoff; i += 1000) {
		// advance to cutoff
		for (; j != end_j && read_lengths[j] < i; ++j) { }
		if (j == end_j) {
			break;
		}
		const size_t count(end_j - j);
		unsigned long long bps(0);
		for (size_t k(j); k != end_j; bps += read_lengths[k], ++k) { }
		const size_t median(j + count / 2);
		if ((count & 1)) {
			histogram.push_back(ReadStats(count, bps, read_lengths[median]));
		} else {
			histogram.push_back(ReadStats(count, bps, (read_lengths[median] + read_lengths[median - 1]) / 2));
		}
		max_values.update(histogram.back());
	}
	const FormatOutput col1(std::max(pretty_print(histogram.size() * 1000).size(), static_cast<size_t>(6)));
	const FormatOutput col2(std::max(pretty_print(max_values.reads).size(), static_cast<size_t>(5)));
	const FormatOutput col3(std::max(pretty_print(max_values.basepairs).size(), static_cast<size_t>(9)));
	const FormatOutput col4(std::max(pretty_print(max_values.median_read_length).size(), static_cast<size_t>(16)));
	// add 3 for decimal point and two digits of fixed-format precision
	const int coverage_width(get_width(opt_genome_size ? max_values.basepairs / opt_genome_size / 1000000 : 0, 4) + 3);
	// add 1 to account for trailing 'x'
	const FormatOutput col5a(opt_genome_size ? coverage_width + 1 : 8, 2);
	const FormatOutput col5b(opt_genome_size ? coverage_width : 7, 2);
	std::cout << col1 << "Cutoff" <<
		" " << col2 << "Reads" <<
		" " << col3 << "Basepairs" <<
		" " << col4 << "Median Read Size";
	if (opt_genome_size) {
		std::cout << " " << col5a << "Coverage";
		if (opt_genome_size) {
			std::cout << " (" << opt_genome_size << ")";
		}
	}
	std::cout << '\n';
	std::cout << col1 << std::setfill('-') << "-" <<
		" " << col2 << std::setfill('-') << "-" <<
		" " << col3 << std::setfill('-') << "-" <<
		" " << col4 << std::setfill('-') << "-";
	if (opt_genome_size) {
		std::cout << " " << col5a << std::setfill('-') << "-";
	}
	std::cout << '\n';
	for (size_t i(0); i != histogram.size(); ++i) {
		std::cout << col1 << pretty_print(i * 1000) <<
			" " << col2 << pretty_print(histogram[i].reads) <<
			" " << col3 << pretty_print(histogram[i].basepairs) <<
			" " << col4 << pretty_print(histogram[i].median_read_length);
		if (opt_genome_size) {
			std::cout << " " << col5b << static_cast<double>(histogram[i].basepairs) / opt_genome_size / 1000000 << 'x';
		}
		std::cout << '\n';
	}
}

// print simple read histogram

static void print_full_histogram(void) {
	std::vector<size_t> &read_lengths(current.read_lengths());
	std::map<size_t, size_t> histogram;
	std::vector<size_t>::const_iterator a(read_lengths.begin());
	const std::vector<size_t>::const_iterator end_a(read_lengths.end());
	for (; a != end_a; ++a) {
		++histogram[*a];
	}
	std::map<size_t, size_t>::const_iterator b(histogram.begin());
	const std::map<size_t, size_t>::const_iterator end_b(histogram.end());
	for (; b != end_b; ++b) {
		std::cout << b->first << " " << b->second << "\n";
	}
}

static void process_files(const std::list<std::string> &file_list) {
	std::list<std::string>::const_iterator a(file_list.begin());
	const std::list<std::string>::const_iterator end_a(file_list.end());
	for (; a != end_a; ++a) {
		read_reads(*a);
	}
	current.flush_seq();
	if (opt_full_histogram) {
		print_full_histogram();
	} else {
		print_histogram();
	}
}

int main(int argc, char **argv) {
	int had_error(0);
	try {
		if (get_opts(argc, argv)) {
			return 0;
		}
		std::list<std::string> file_list;
		find_files(argc, argv, file_list);
		process_files(file_list);
	} catch (std::exception &e) {
		std::cerr << "Error: " << e.what() << "\n";
		LocalException *x(dynamic_cast<LocalException *>(&e));
		if (x != 0 && x->show_usage()) {
			print_usage();
		}
		had_error = 1;
	}
	return had_error;
}
