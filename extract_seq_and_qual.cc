#include "breakup_line.h"	// breakup_line()
#include "open_compressed.h"	// close_compressed(), find_suffix(), get_suffix(), open_compressed(), pfgets(), pfpeek()
#include "pattern.h"	// Pattern
#include "version.h"	// VERSION
#include "write_fork.h"	// close_fork(), pfputs(), write_fork()
#include <algorithm>	// min(), sort(), swap()
#include <ctype.h>	// isdigit(), isspace()
#include <errno.h>	// ENOENT, errno
#include <exception>	// exception
#include <getopt.h>	// getopt(), optarg, optind
#include <glob.h>	// GLOB_NOCHECK, glob(), globfree(), glob_t
#include <iostream>	// cerr
#include <list>		// list<>
#include <map>		// map<>
#include <regex.h>	// REG_EXTENDED, REG_NOSUB
#include <set>		// set<>
#include <sstream>	// istringstream, ostringstream
#include <stdio.h>	// EOF
#include <stdlib.h>	// exit()
#include <string.h>	// strchr()
#include <string>	// string
#include <sys/stat.h>	// stat(), struct stat
#include <sys/types.h>	// size_t
#include <unistd.h>	// unlink()
#include <utility>	// make_pair(), pair<>
#include <vector>	// vector<>

// this program lets you to extract specific reads from fasta and qual files,
// with the option of doing a little clean up on them as well
//
// you can also extract reads into multiple output files with the -S option;
// write buffering is used to prevent disk thrashing
//
// when using the -i option, duplicate reads or read ranges are removed
// (although the same read or read range can appear in multiple lists when
// using the -S option, they will still be removed if duplicated in the
// same -i list)
//
// reads are output in the same order as they were in the original file(s)
//
// all input files may be compressed; output files will be compressed if
// the -o or -S option has a compression suffix
//
// in the event both -i and -x options are given, -x options only exclude
// reads from -i options given earlier on the command line, and may be
// overridden by later -i options (put them last if you want them to
// override all -i options)

class LocalException : public std::exception {
    private:
	const std::string error_;
	const int show_usage_;
    public:
	explicit LocalException(const std::string &s) : error_(s), show_usage_(0) { }
	LocalException(const std::string &s, const int i) : error_(s), show_usage_(i) { }
	virtual ~LocalException(void) throw() { }
	virtual const char * what() const throw() {
		return error_.c_str();
	}
	const int &show_usage(void) const throw() {
		return show_usage_;
	}
};

class OutputStream {
    private:
	int fd_, has_args_, qual_with_seq_;
	std::string file_, buffer_;
    public:
	OutputStream(void) : fd_(-1), has_args_(0), qual_with_seq_(0) { }
	explicit OutputStream(const std::string &a, const int i, const size_t j) : fd_(-1), has_args_(0), qual_with_seq_(0) {
		this->open(a, i, j);
	}
	~OutputStream(void) {
		this->close(0);
	}
	void open(const std::string &, int, size_t);
	int close(int);
	void write(const std::string &s) {
		if (!buffer_.capacity()) {
			pfputs(fd_, s);
		} else if (buffer_.size() + s.size() > buffer_.capacity()) {
			pfputs(fd_, buffer_);
			pfputs(fd_, s);
			buffer_.clear();
		} else {
			buffer_ += s;
		}
	}
};

class ReadRange {
    public:
	size_t start, stop;
	int output_file;
	ReadRange(const size_t i, const size_t j, const int k) : start(i), stop(j), output_file(k) { }
	// used for excludes
	explicit ReadRange(const int k) : start(static_cast<size_t>(-1)), stop(0), output_file(k) { }
	~ReadRange(void) { }
	bool is_exclude(void) const {
		return start == static_cast<size_t>(-1);
	}
	bool is_range(void) const {
		return start != 0 || stop != 0;
	}
	bool operator<(const ReadRange &a) const {
		return start < a.start || (start == a.start && (stop < a.stop || (stop == a.stop && output_file < a.output_file)));
	}
	bool operator==(const ReadRange &a) const {
		return start == a.start && stop == a.stop && output_file == a.output_file;
	}
};

static int (* get_id_start)(std::string &line);
static int (* is_desired_read)(const std::string &line);
static int opt_complement, opt_fastq_output, opt_qual_only, opt_regex;
static int opt_read_list, opt_strip_leading_zero, opt_strip_trace;
static int opt_strip_trailing_zero, opt_validate;
static size_t opt_line_size, opt_min_length, opt_max_length;
static std::map<std::string, size_t> read_size;
static std::map<std::string, std::list<ReadRange> > reads;
static std::string opt_output_suffix;
static std::vector<std::pair<Pattern, int> > read_patterns;

static std::string get_header(const std::string &line) {
	size_t i(1);
	const size_t end_i(line.size());
	for (; i != end_i && !isspace(line[i]); ++i) { }
	return line.substr(1, i - 1);
}

// this class holds, processes, and prints the current read data;
// however, as seq data and qual data are handled slightly differently,
// certain methods have a _seq or _qual suffix, so make sure to use the
// appropriate one

class CurrentState {
    private:
	std::vector<OutputStream> seq_output_;
	std::vector<OutputStream *> qual_output_;
	std::string conversion_;		// used for seq complementing
	std::string seq_id_, qual_id_, seq_, qual_;
	size_t seq_length_, qual_length_;
    private:
	std::string get_header_seq(void) const {
		return get_header(seq_id_);
	}
	std::string get_header_qual(void) const {
		return get_header(qual_id_);
	}
	void (* print_qual)(const std::string &, const std::string &, OutputStream &, size_t);
	void (* print_qual_range)(const std::string &, const std::vector<std::string> &, const size_t, const size_t, OutputStream &, size_t);
    public:
	CurrentState(void) : seq_length_(0), qual_length_(0) { }
	~CurrentState(void) { }
	void initialize(void);
	void open_outputs(const std::list<std::pair<std::string, std::string> > &files, int qual_with_seq = 0);
	void close_outputs(int had_error);
	void id_check(const std::string &, const std::string &) const;
	void write_seq(void);
	void write_qual(void);
	void write_fastq(void);
	void flush_seq(void) {
		write_seq();
		seq_id_.clear();
		seq_.clear();
	}
	void flush_qual(void) {
		write_qual();
		qual_id_.clear();
		qual_.clear();
	}
	void flush_fastq(void) {
		write_fastq();
		seq_id_.clear();
		qual_id_.clear();
		seq_.clear();
		qual_.clear();
	}
	void set_seq(const std::string &id, const std::string &data, const size_t &length) {
		seq_ = data;
		seq_id_ = id;
		seq_length_ = length;
	}
	void set_qual(const std::string &id, const std::string &data, const size_t &length) {
		qual_ = data;
		qual_id_ = id;
		qual_length_ = length;
	}
	void set_fastq(const std::string &id, const std::string &seq, const size_t &length, const std::string &qual) {
		seq_ = seq;
		qual_ = qual;
		seq_id_ = qual_id_ = id;
		seq_length_ = qual_length_ = length;
	}
	void add_seq(const std::string &data) {
		seq_ += data;
	}
	void add_qual(const std::string &data) {
		// make sure there is a space between data from different lines
		if (!qual_.empty() && *qual_.rbegin() != ' ' && !data.empty() && data[0] != ' ') {
			qual_ += " ";
		}
		qual_ += data;
	}
	void add_fastq(const std::string &seq, const std::string &qual) {
		seq_ += seq;
		qual_ += qual;
	}
	void complement_seq(std::string &) const;
	void complement_qual(std::string &) const;
};

static CurrentState current;

void OutputStream::open(const std::string &file, const int qual_with_seq, const size_t buffer_size) {
	file_ = file;
	qual_with_seq_ = qual_with_seq;
	if (!file_.empty()) {
		std::string suffix;
		get_suffix(file_, suffix);
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
		has_args_ = !args.empty();
		fd_ = write_fork(args, file_);
		if (fd_ == -1) {
			throw LocalException("could not open " + file_);
		}
		buffer_.reserve(buffer_size);
	}
}

int OutputStream::close(const int had_error) {
	const int i(fd_);
	if (fd_ != -1 && fd_ != -2) {
		if (!buffer_.empty()) {		// flush output buffer
			if (!had_error) {
				pfputs(fd_, buffer_);
			}
			buffer_.clear();
		}
		close_fork(fd_);
		fd_ = -2;		// flag to prevent qual memory leak
		struct stat buf;
		if (!has_args_ && (had_error || qual_with_seq_) && stat(file_.c_str(), &buf) == 0 && buf.st_size == 0) {
			unlink(file_.c_str());
		}
	}
	file_.clear();
	return i;
}

void CurrentState::complement_seq(std::string &seq) const {
	if (seq.empty()) {
		return;
	}
	std::string::iterator a(seq.begin());
	std::string::iterator b(seq.end());
	for (--b; a < b; ++a, --b) {
		const char c(conversion_[*b]);
		*b = conversion_[*a];
		*a = c;
	}
	if (a == b) {
		*a = conversion_[*a];
	}
}

static void print_seq(const std::string &id, const std::string &seq, OutputStream &f, const size_t length) {
	std::string t;
	if (opt_complement) {
		current.complement_seq(t = seq);
	}
	const std::string &x(opt_complement ? t : seq);
	f.write((opt_fastq_output ? "@" : ">") + id.substr(1) + "\n");
	for (size_t i(0); i < x.size(); i += length) {
		// substr() is smart and will not go past end of string
		f.write(x.substr(i, length) + "\n");
	}
}

static void print_qual_fasta(const std::string &id, const std::string &qual, OutputStream &f, const size_t length) {
	std::string t;
	if (opt_complement) {
		current.complement_qual(t = qual);
	}
	const std::string &r(opt_complement ? t : qual);
	f.write(">" + id.substr(1) + "\n");
	size_t i(0);
	const size_t end_i(r.size());
	if (id[0] == '@') {			// fastq format quality
		std::ostringstream x;
		for (; i != end_i; ++i) {
			const size_t stop(std::min(i + length, end_i) - 1);
			for (x.str(""); i != stop; ++i) {
				const int z(r[i] - 33);	// convert to int
				x << z << " ";
			}
			const int z(r[i] - 33);		// convert to int
			x << z;
			f.write(x.str() + "\n");
		}
	} else {
		for (; i != end_i && isspace(r[i]); ++i) { }
		while (i != end_i) {
			size_t j(i + 1);
			size_t last_j(j);
			for (size_t k(0); j != end_i && k != length; ++k) {
				for (; j != end_i && !isspace(r[j]); ++j) { }
				last_j = j;
				for (; j != end_i && isspace(r[j]); ++j) { }
			}
			f.write(r.substr(i, last_j - i) + "\n");
			i = j;
		}
	}
}

static void print_qual_fastq(const std::string &id, const std::string &qual, OutputStream &f, const size_t /*length*/) {
	std::string t;
	if (opt_complement) {
		current.complement_qual(t = qual);
	}
	const std::string &r(opt_complement ? t : qual);
	f.write("+\n");
	if (id[0] == '@') {			// fastq format quality
		f.write(r + "\n");
	} else {
		std::string out;
		size_t i(0);
		const size_t end_i(r.size());
		for (; i != end_i && isspace(r[i]); ++i) { }
		while (i != end_i) {
			const size_t j(i);
			for (++i; i != end_i && !isspace(r[i]); ++i) { }
			unsigned int c;
			std::istringstream(r.substr(j, i - j)) >> c;
			out += static_cast<unsigned char>(c + 33);
			for (; i != end_i && isspace(r[i]); ++i) { }
		}
		f.write(out + "\n");
	}
}

static void print_qual_range_fasta(const std::string &header, const std::vector<std::string> &quals, const size_t start, const size_t stop, OutputStream &f, const size_t length) {
	std::ostringstream x;
	x << (start + 1) << "-" << stop;
	f.write(header + x.str() + "\n");
	if (header[0] == '@') {				// fastq format quality
		std::string t;
		if (opt_complement) {
			current.complement_qual(t = quals[0]);
		}
		const std::string &r(opt_complement ? t : quals[0]);
		for (size_t i(start); i != stop; ++i) {
			const size_t end_i(std::min(i + length, stop) - 1);
			for (x.str(""); i != end_i; ++i) {
				const int z(r[i] - 33);	// convert to int
				x << z << " ";
			}
			const int z(r[i] - 33);		// convert to int
			x << z;
			f.write(x.str() + "\n");
		}
	} else if (opt_complement) {
		for (size_t i(stop - 1); i != start - 1; --i) {
			// be careful to avoid wrapping in max()
			const size_t end_i(std::max(i + 1 > length ? i + 1 - length : 0, start));
			for (; i != end_i; --i) {
				f.write(quals[i] + " ");
			}
			f.write(quals[i] + "\n");
		}
	} else {
		for (size_t i(start); i != stop; ++i) {
			const size_t end_i(std::min(i + length, stop) - 1);
			for (; i != end_i; ++i) {
				f.write(quals[i] + " ");
			}
			f.write(quals[i] + "\n");
		}
	}
}

static void print_qual_range_fastq(const std::string &header, const std::vector<std::string> &quals, const size_t start, const size_t stop, OutputStream &f, const size_t /*length*/) {
	f.write("+\n");
	if (header[0] == '@') {				// fastq format quality
		std::string t;
		if (opt_complement) {
			current.complement_qual(t = quals[0]);
		}
		const std::string &r(opt_complement ? t : quals[0]);
		f.write(r.substr(start, stop - start) + "\n");
	} else if (opt_complement) {
		std::string out;
		for (size_t i(stop - 1); i != start - 1; --i) {
			unsigned int c;
			std::istringstream(quals[i]) >> c;
			out += static_cast<unsigned char>(c + 33);
		}
		f.write(out + "\n");
	} else {
		std::string out;
		for (size_t i(start); i != stop; ++i) {
			unsigned int c;
			std::istringstream(quals[i]) >> c;
			out += static_cast<unsigned char>(c + 33);
		}
		f.write(out + "\n");
	}
}

void CurrentState::initialize() {
	if (opt_complement) {
		conversion_.assign(256, 0);
		for (int i = 1; i != 256; ++i) {
			conversion_[i] = i;
		}
		conversion_['A'] = 'T';
		conversion_['C'] = 'G';
		conversion_['G'] = 'C';
		conversion_['T'] = 'A';
		conversion_['a'] = 't';
		conversion_['c'] = 'g';
		conversion_['g'] = 'c';
		conversion_['t'] = 'a';
	}
	if (opt_fastq_output) {
		print_qual = print_qual_fastq;
		print_qual_range = print_qual_range_fastq;
	} else {
		print_qual = print_qual_fasta;
		print_qual_range = print_qual_range_fasta;
	}
}

// we have to open all outputs (seq and qual) up front to handle
// the possibility of fastq input files

void CurrentState::open_outputs(const std::list<std::pair<std::string, std::string> > &files, const int qual_with_seq) {
	// count total number of each kind of file to be opened
	size_t seq(0);
	size_t qual(0);
	std::list<std::pair<std::string, std::string> >::const_iterator a(files.begin());
	const std::list<std::pair<std::string, std::string> >::const_iterator end_a(files.end());
	for (; a != end_a; ++a) {
		if (!a->first.empty()) {
			++seq;
		}
		if (!a->second.empty()) {
			++qual;
		}
	}
	size_t buffer_size(0);
	if (!opt_output_suffix.empty() && seq + qual != 0) {
		buffer_size = (1 << 30) / (seq + qual);
	}
	// any output files not specified are redirected to /dev/null
	// to prevent crashes from writing to nonexistent streams;
	// allocated as a pointer to allow persistence
	OutputStream *dev_null(new OutputStream("/dev/null", 1, 0));
	int dev_null_used(0);
	seq_output_.resize(files.size());
	qual_output_.resize(files.size());
	a = files.begin();
	for (size_t i(0); a != end_a; ++i, ++a) {
		if (!a->first.empty()) {
			seq_output_[i].open(a->first, 0, buffer_size);
		} else {
			seq_output_[i] = *dev_null;
			dev_null_used = 1;
		}
	}
	if (opt_fastq_output) {
		const size_t end_i(files.size());
		for (size_t i(0); i != end_i; ++i) {
			qual_output_[i] = &seq_output_[i];
		}
	} else {
		a = files.begin();
		for (size_t i(0); a != end_a; ++i, ++a) {
			if (!a->second.empty()) {
				qual_output_[i] = new OutputStream;
				qual_output_[i]->open(a->second, qual_with_seq, buffer_size);
			} else {
				qual_output_[i] = dev_null;
				dev_null_used = 1;
			}
		}
	}
	if (!dev_null_used) {
		delete dev_null;
	}
}

void CurrentState::close_outputs(int had_error) {
	const size_t end_s(seq_output_.size());
	for (size_t i(0); i != end_s; ++i) {
		seq_output_[i].close(had_error);
	}
	seq_output_.clear();
	const size_t end_q(qual_output_.size());
	for (size_t i(0); i != end_q; ++i) {
		if (qual_output_[i]->close(had_error) != -2) {
			delete qual_output_[i];
		}
	}
	qual_output_.clear();
}

void CurrentState::write_seq() {
	if (seq_id_.empty()) {
		return;
	}
	if (opt_validate) {
		read_size[get_header_seq()] = seq_.size();
		return;
	}
	if (seq_.size() < opt_min_length) {
		return;
	}
	if (opt_max_length && opt_max_length < seq_.size()) {
		return;
	}
	if (opt_line_size != 0) {
		seq_length_ = opt_line_size;
	}
	if (opt_read_list == 2) {			// process ranges
		std::string header(get_header_seq());
		const std::list<ReadRange> &a(reads[header]);
		const size_t n(seq_.size());
		header = seq_id_[0] + header + ":";
		std::list<ReadRange>::const_iterator b(a.begin());
		const std::list<ReadRange>::const_iterator end_b(a.end());
		for (; b != end_b; ++b) {
			if (b->start >= n) {
				std::cerr << "Warning: specified range not on read: " << header << " " << (b->start + 1) << "-" << b->stop << "\n";
				continue;
			}
			const size_t stop(b->stop == 0 || b->stop > n ? n : b->stop);
			if (b->start == 0 && stop == n) {	// full read
				print_seq(seq_id_, seq_, seq_output_[b->output_file], seq_length_);
			} else {
				std::ostringstream x;
				x << header << (b->start + 1) << "-" << stop;
				print_seq(x.str(), seq_.substr(b->start, stop - b->start), seq_output_[b->output_file], seq_length_);
			}
		}
	} else if (opt_read_list == 1) {		// process includes
		if (opt_regex) {
			const std::string read_name(get_header_seq());
			std::vector<std::pair<Pattern, int> >::iterator a(read_patterns.begin());
			const std::vector<std::pair<Pattern, int> >::const_iterator end_a(read_patterns.end());
			for (; a != end_a; ++a) {
				if (a->second != -1 && a->first.is_match(read_name)) {
					print_seq(seq_id_, seq_, seq_output_[a->second], seq_length_);
				}
			}
		} else {
			const std::list<ReadRange> &a(reads[get_header_seq()]);
			std::list<ReadRange>::const_iterator b(a.begin());
			const std::list<ReadRange>::const_iterator end_b(a.end());
			for (; b != end_b; ++b) {
				print_seq(seq_id_, seq_, seq_output_[b->output_file], seq_length_);
			}
		}
	} else {
		print_seq(seq_id_, seq_, seq_output_[0], seq_length_);
	}
}

static size_t qual_line_size(const std::string &line) {
	size_t i(0);
	const size_t end_i(line.size());
	for (; i != end_i && isspace(line[i]); ++i) { }
	size_t n(0);
	for (; i != end_i; ++n) {
		for (; i != end_i && !isspace(line[i]); ++i) { }
		for (; i != end_i && isspace(line[i]); ++i) { }
	}
	return n;
}

void CurrentState::complement_qual(std::string &qual) const {
	if (qual.empty()) {
		return;
	}
	// we do this in two passes, so we can do it in place
	// first, reverse string
	std::string::iterator a(qual.begin());
	std::string::iterator b(qual.end());
	for (--b; a < b; ++a, --b) {
		std::swap(*a, *b);
	}
	if (qual_id_[0] != '@') {	// not needed for fastq format quals
		// next, go back and reverse words in string
		size_t i(0);
		const size_t end_i(qual.size());
		while (i != end_i) {
			for (; i != end_i && isspace(qual[i]); ++i) { }
			size_t j(i);
			for (; i != end_i && !isspace(qual[i]); ++i) { }
			size_t k(i);
			for (--k; j < k; ++j, --k) {
				std::swap(qual[j], qual[k]);
			}
		}
	}
}

// truncate string to remove trailing zero (if any)

static void strip_trailing_zero(std::string &qual, const bool fastq_format) {
	size_t i(qual.size() - 1);
	const size_t end_i(static_cast<size_t>(-1));
	if (fastq_format) {
		if (i != end_i && qual[i] == 33) {
			qual.resize(i);
		}
	} else {
		for (; i != end_i && isspace(qual[i]); --i) { }
		if (i != end_i && qual[i] == '0') {
			for (--i; i != end_i && qual[i] == '0'; --i) { }
			if (i == end_i) {
				qual.clear();
			} else if (isspace(qual[i])) {
				for (--i; i != end_i && isspace(qual[i]); --i) { }
				qual.resize(i + 1);
			}
		}
	}
}

// change leading zero to spaces (if any)

static void strip_leading_zero(std::string &qual, const bool fastq_format) {
	if (fastq_format) {
		if (!qual.empty() && qual[0] == 33) {
			qual.erase(0, 1);
		}
	} else {
		size_t i(0);
		const size_t end_i(qual.size());
		for (; i != end_i && isspace(qual[i]); ++i) { }
		if (i != end_i && qual[i] == '0') {
			size_t j(i);
			for (++i; i != end_i && qual[i] == '0'; ++i) { }
			if (i == end_i) {
				qual.clear();
			} else if (isspace(qual[i])) {
				qual.replace(j, i - j, i - j, ' ');
			}
		}
	}
}

void CurrentState::write_qual() {
	if (qual_id_.empty()) {
		return;
	}
	// do processing on the quals, if any
	if (opt_strip_trailing_zero) {
		strip_trailing_zero(qual_, qual_id_[0] == '@');
	} else if (opt_strip_leading_zero) {
		strip_leading_zero(qual_, qual_id_[0] == '@');
	}
	if (opt_validate) {
		const std::string header(get_header_qual());
		std::map<std::string, size_t>::iterator a(read_size.find(header));
		if (a == read_size.end()) {
			std::cerr << header << ": sequence missing\n";
		} else {
			const size_t q(qual_id_[0] == '@' ? qual_.size() : qual_line_size(qual_));
			if (a->second != q) {
				std::cerr << header << ": size mismatch: " << a->second << " != " << q << "\n";
			}
			read_size.erase(a);
		}
		return;
	}
	if (qual_.size() < opt_min_length) {
		return;
	}
	if (opt_max_length && opt_max_length < qual_.size()) {
		return;
	}
	if (opt_line_size != 0) {
		qual_length_ = opt_line_size;
	}
	if (opt_read_list == 2) {		// process ranges
		std::string header(get_header_qual());
		const std::list<ReadRange> &a(reads[header]);
		std::vector<std::string> quals;
		if (qual_id_[0] == '@') {		// fastq format quality
			quals.push_back(qual_);
		} else {
			// overestimate, but quick
			quals.reserve(qual_.size() / 2);
			breakup_line(qual_, quals);
		}
		const size_t n(qual_id_[0] == '@' ? qual_.size() : quals.size());
		header = qual_id_[0] + header + ":";
		std::list<ReadRange>::const_iterator b(a.begin());
		const std::list<ReadRange>::const_iterator end_b(a.end());
		for (; b != end_b; ++b) {
			if (b->start >= n) {
				std::cerr << "Warning: specified range not on read: " << header << " " << (b->start + 1) << "-" << b->stop << "\n";
				continue;
			}
			const size_t stop(b->stop == 0 || b->stop > n ? n : b->stop);
			if (b->start == 0 && stop == n) {	// full read
				print_qual(qual_id_, qual_, *qual_output_[b->output_file], qual_length_);
			} else {
				print_qual_range(header, quals, b->start, stop, *qual_output_[b->output_file], qual_length_);
			}
		}
	} else if (opt_read_list == 1) {	// process includes
		if (opt_regex) {
			const std::string read_name(get_header_qual());
			std::vector<std::pair<Pattern, int> >::iterator a(read_patterns.begin());
			const std::vector<std::pair<Pattern, int> >::const_iterator end_a(read_patterns.end());
			for (; a != end_a; ++a) {
				if (a->second != -1 && a->first.is_match(read_name)) {
					print_qual(qual_id_, qual_, *qual_output_[a->second], qual_length_);
				}
			}
		} else {
			const std::list<ReadRange> &a(reads[get_header_qual()]);
			std::list<ReadRange>::const_iterator b(a.begin());
			const std::list<ReadRange>::const_iterator end_b(a.end());
			for (; b != end_b; ++b) {
				print_qual(qual_id_, qual_, *qual_output_[b->output_file], qual_length_);
			}
		}
	} else {
		print_qual(qual_id_, qual_, *qual_output_[0], qual_length_);
	}
}

void CurrentState::write_fastq() {
	if (seq_id_.empty()) {
		return;
	}
	// do processing on the quals, if any
	if (opt_strip_trailing_zero) {
		strip_trailing_zero(qual_, qual_id_[0] == '@');
	} else if (opt_strip_leading_zero) {
		strip_leading_zero(qual_, qual_id_[0] == '@');
	}
	if (opt_validate) {
		const size_t q(qual_id_[0] == '@' ? qual_.size() : qual_line_size(qual_));
		if (seq_.size() != q) {
			std::cerr << get_header_qual() << ": size mismatch: " << seq_.size() << " != " << q << "\n";
		}
		return;
	}
	if (seq_.size() < opt_min_length) {
		return;
	}
	if (opt_max_length && opt_max_length < seq_.size()) {
		return;
	}
	if (opt_line_size != 0) {
		seq_length_ = opt_line_size;
	}
	if (opt_read_list == 2) {			// process ranges
		std::string header(get_header_seq());
		const std::list<ReadRange> &a(reads[header]);
		const size_t n(seq_.size());
		std::vector<std::string> quals;
		if (seq_id_[0] == '@') {		// fastq format quality
			quals.push_back(qual_);
		} else {
			quals.reserve(n);
			breakup_line(qual_, quals);
		}
		header = seq_id_[0] + header + ":";
		std::list<ReadRange>::const_iterator b(a.begin());
		const std::list<ReadRange>::const_iterator end_b(a.end());
		for (; b != end_b; ++b) {
			if (b->start >= n) {
				std::cerr << "Warning: specified range not on read: " << header << " " << (b->start + 1) << "-" << b->stop << "\n";
				continue;
			}
			const size_t stop(b->stop == 0 || b->stop > n ? n : b->stop);
			if (b->start == 0 && stop == n) {	// full read
				print_seq(seq_id_, seq_, seq_output_[b->output_file], seq_length_);
				print_qual(seq_id_, qual_, *qual_output_[b->output_file], qual_length_);
			} else {
				std::ostringstream x;
				x << header << (b->start + 1) << "-" << stop;
				print_seq(x.str(), seq_.substr(b->start, stop - b->start), seq_output_[b->output_file], seq_length_);
				print_qual_range(header, quals, b->start, stop, *qual_output_[b->output_file], seq_length_);
			}
		}
	} else if (opt_read_list == 1) {		// process includes
		if (opt_regex) {
			const std::string read_name(get_header_seq());
			std::vector<std::pair<Pattern, int> >::iterator a(read_patterns.begin());
			const std::vector<std::pair<Pattern, int> >::const_iterator end_a(read_patterns.end());
			for (; a != end_a; ++a) {
				if (a->second != -1 && a->first.is_match(read_name)) {
					print_seq(seq_id_, seq_, seq_output_[a->second], seq_length_);
					print_qual(seq_id_, qual_, *qual_output_[a->second], seq_length_);
				}
			}
		} else {
			const std::list<ReadRange> &a(reads[get_header_seq()]);
			std::list<ReadRange>::const_iterator b(a.begin());
			const std::list<ReadRange>::const_iterator end_b(a.end());
			for (; b != end_b; ++b) {
				print_seq(seq_id_, seq_, seq_output_[b->output_file], seq_length_);
				print_qual(seq_id_, qual_, *qual_output_[b->output_file], seq_length_);
			}
		}
	} else {
		print_seq(seq_id_, seq_, seq_output_[0], seq_length_);
		print_qual(seq_id_, qual_, *qual_output_[0], seq_length_);
	}
}

// make sure the seq and qual ids match; if not, throw an exception

void CurrentState::id_check(const std::string &id_seq, const std::string &id_qual) const {
	size_t i(1);
	const size_t end_i(id_seq.size());
	for (; i != end_i && !isspace(id_seq[i]); ++i) { }
	size_t j(1);
	const size_t end_j(id_qual.size());
	for (; j != end_j && !isspace(id_qual[j]); ++j) { }
	if (id_seq.substr(1, i - 1) != id_qual.substr(1, j - 1)) {
		throw LocalException("id mismatch between seq and qual: " + id_seq.substr(1, i - 1) + " != " + id_qual.substr(1, j - 1));
	}
}

// look for read names with @##-##,##-##,..., parse as ranges
// i is current offset into read, and is changed to end of next parse
static int parse_range(const std::string &read, size_t &offset, const int is_include, const int file_number) {
	// using regex's for read names precludes ranges;
	// excludes can't be ranges, includes can be (@);
	// otherwise, just look for start of next element (,)
	size_t i(read.find_first_of((!opt_regex && is_include) ? "@," : ",", offset));
	if (!opt_regex && is_include && i != std::string::npos && read[i] == '@') {
		const size_t end_i(read.size());
		if (++i != end_i && isdigit(read[i])) {
			const std::string name(read.substr(offset, i - offset - 1));
			std::list<ReadRange> x;
			do {
				size_t j(i + 1);
				for (; j != end_i && isdigit(read[j]); ++j) { }
				if (j == end_i || read[j] != '-' || j == end_i - 1 || !isdigit(read[j + 1])) {
					break;
				}
				size_t k(j + 2);
				for (; k != end_i && isdigit(read[k]); ++k) { }
				if (k != end_i && read[k] != ',') {
					break;
				}
				offset = k == end_i ? std::string::npos : k + 1;
				size_t m, n;
				std::istringstream(read.substr(i, j - i)) >> m;
				std::istringstream(read.substr(j + 1, k - j - 1)) >> n;
				if (m > n && n != 0) {
					std::cerr << "Warning: improper range: " << m << "-" << n << ", discarding\n";
					if (k == end_i) {
						break;
					}
					continue;
				}
				// change ranges from 1 offset inclusive end to
				// 0 offset exclusive end
				x.push_back(ReadRange(m - 1, n, file_number));
				if (k == end_i) {
					break;
				}
				i = k + 1;
			} while (i != end_i);
			if (!x.empty()) {
				std::map<std::string, std::list<ReadRange> >::iterator a(reads.find(name));
				if (a == reads.end()) {
					reads[name] = x;
				} else if (a->second.front().is_exclude()) {
					// later options override earlier ones
					reads[name] = x;
				} else {
					a->second.splice(a->second.end(), x);
				}
				return offset != std::string::npos;
			}
		}
		i = read.find(',', offset);
	}
	std::string name;
	// part of comma separated list
	if (i != std::string::npos) {
		if (i == offset) {	// empty element - consecutive ,'s
			offset = i + 1;
			return 1;
		}
		name = read.substr(offset, i - offset);
		offset = i + 1;
	} else if (offset == read.size()) {	// empty element - trailing ,
		return 0;
	} else {
		name = read.substr(offset);
		offset = std::string::npos;
	}
	if (opt_regex) {
		const Pattern x(name, 0, REG_EXTENDED | REG_NOSUB);
		if (!x.empty()) {
			read_patterns.push_back(std::make_pair(x, is_include ? file_number : -1));
		}
	} else {
		std::map<std::string, std::list<ReadRange> >::iterator a(reads.find(name));
		if (!is_include) {
			if (a == reads.end()) {
				reads[name].assign(1, ReadRange(file_number));
			} else if (!a->second.front().is_exclude()) {
				// later options override earlier ones
				a->second.assign(1, ReadRange(file_number));
			}
		} else if (a == reads.end()) {
			reads[name].assign(1, ReadRange(0, 0, file_number));
		} else if (a->second.front().is_exclude() || a->second.front().output_file == file_number) {
			// later options override earlier ones;
			// a full range include overrides any previous one, if
			// the results aren't going into separate files
			a->second.assign(1, ReadRange(0, 0, file_number));
		} else {
			a->second.push_back(ReadRange(0, 0, file_number));
		}
	}
	return offset != std::string::npos;
}

static bool read_pattern_cmp(const std::pair<Pattern, int> &a, const std::pair<Pattern, int> &b) {
	if (a.second != b.second) {
		return a.second < b.second;
	} else {
		return a.first < b.first;
	}
};

// turn read_list into a map of read names
static void process_read_lists(const std::list<std::pair<std::string, int> > &read_list) {
	int has_includes(0), has_excludes(0);
	std::list<std::pair<std::string, int> >::const_iterator a(read_list.begin());
	const std::list<std::pair<std::string, int> >::const_iterator end_a(read_list.end());
	for (int file_number(0); a != end_a; ++a) {
		if (a->second) {
			has_includes = 1;
		} else {
			has_excludes = 1;
		}
		if (a->first.find(',') != std::string::npos) {	// list of names
			if (!opt_output_suffix.empty()) {
				throw LocalException("-i options can only be given file names when used with the -S option", 1);
			}
			size_t i(0);
			while (parse_range(a->first, i, a->second, file_number)) { }
		} else {				// file
			const int fd(open_compressed(a->first));
			if (fd != -1) {
				std::string line;
				while (pfgets(fd, line) != -1) {
					size_t i(0);
					while (parse_range(line, i, a->second, file_number)) { }
				}
				close_compressed(fd);
			} else {
				std::cerr << "Warning: failed to open " << a->first << "\n";
			}
		}
		if (a->second && !opt_output_suffix.empty()) {
			++file_number;
		}
	}
	if (has_excludes) {
		if (!has_includes) {
			opt_read_list = 3;
			return;
		}
		// remove remaining excludes (as we have includes, anything
		// not already in the include list won't be matched anyway)
		std::map<std::string, std::list<ReadRange> >::iterator b(reads.begin());
		const std::map<std::string, std::list<ReadRange> >::const_iterator end_b(reads.end());
		while (b != end_b) {
			if (b->second.front().is_exclude()) {
				reads.erase(b++);
			} else {
				++b;
			}
		}
	}
	if (has_includes && ((opt_regex && read_patterns.empty()) || (!opt_regex && reads.empty()))) {
		throw LocalException("empty include list: no reads will be selected", 0);
	}
	if (opt_regex) {
		std::sort(read_patterns.begin(), read_patterns.end(), read_pattern_cmp);
	}
	int has_ranges(0);
	// remove any duplicate ReadRanges, but preserve order of ranges
	std::map<std::string, std::list<ReadRange> >::iterator b(reads.begin());
	const std::map<std::string, std::list<ReadRange> >::const_iterator end_b(reads.end());
	for (; b != end_b; ++b) {
		std::set<ReadRange> x;
		std::list<ReadRange> &c(b->second);
		std::list<ReadRange>::iterator d(c.begin());
		const std::list<ReadRange>::const_iterator end_d(c.end());
		while (d != end_d) {
			if (x.find(*d) == x.end()) {	// first sighting
				if (d->is_range()) {
					has_ranges = 1;
				}
				x.insert(*d);
				++d;
			} else {			// duplicate
				d = c.erase(d);
			}
		}
	}
	opt_read_list = has_ranges ? 2 : (has_includes ? 1 : 0);
}

// make the quality file name from the sequence file name - mainly just
// tacks on .qual, but also handles Z/gz/bz2 endings (since the .qual
// needs to be added before the compression suffix); may change the
// seq_file name if seq_file doesn't exist, but seq_file with a compression
// suffix does

static int find_qual(std::string &seq_file, std::string &qual_file, int new_file) {
	if (seq_file.empty() || seq_file == "-") {
		return 0;
	}
	// first find actual seq name (and suffix)
	std::string suffix;
	if (new_file) {
		get_suffix(seq_file, suffix);
	} else if (find_suffix(seq_file, suffix) == -1) {
		return 0;
	}
	const std::string name(seq_file, 0, seq_file.size() - suffix.size());
	qual_file = name + ".qual";
	std::string qual_suffix;
	if (!new_file && find_suffix(qual_file, qual_suffix) == 0) {
		return 1;
	}
	if (name.size() > 4 && name.compare(name.size() - 4, 4, ".fna") == 0) {
		qual_file = name.substr(0, name.size() - 3) + "qual";
		if (new_file) {
			qual_file += suffix;
			return 1;
		} else if (find_suffix(qual_file, qual_suffix) == 0) {
			return 1;
		}
	}
	if (name.size() > 6 && name.compare(name.size() - 6, 6, ".fasta") == 0) {
		qual_file = name.substr(0, name.size() - 5) + "qual";
		if (new_file) {
			qual_file += suffix;
			return 1;
		} else if (find_suffix(qual_file, qual_suffix) == 0) {
			return 1;
		}
	}
	if (name.size() > 1 && name[0] == 'f') {
		size_t i(1);
		const size_t end_i(name.size());
		for (; i != end_i && isdigit(name[i]); ++i) { }
		if (i == end_i) {
			qual_file = "q" + name.substr(1);
			if (new_file) {
				qual_file += suffix;
				return 1;
			} else if (find_suffix(qual_file, qual_suffix) == 0) {
				return 1;
			}
		}
	}
	if (new_file) {
		qual_file = name + ".qual" + suffix;
		return 1;
	} else {
		qual_file.clear();
		return 0;
	}
}

// returns true if attempt to strip the trace worked
static int strip_trace(std::string &line) {
	size_t i(2);
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
	return 1;
}

// returns true if read is desired
static int is_desired_read_plain(const std::string &line) {
	if (opt_read_list == 0) {
		return 1;
	}
	const std::string header(get_header(line));
	const std::map<std::string, std::list<ReadRange> >::const_iterator a(reads.find(header));
	if ((opt_read_list == 3 && a == reads.end()) || (opt_read_list != 3 && a != reads.end())) {
		return 1;
	}
	return 0;
}

// returns true if read is desired
static int is_desired_read_regex(const std::string &line) {
	if (opt_read_list == 0) {
		return 1;
	}
	const std::string header(get_header(line));
	std::vector<std::pair<Pattern, int> >::iterator a(read_patterns.begin());
	const std::vector<std::pair<Pattern, int> >::const_iterator end_a(read_patterns.end());
	for (; a != end_a; ++a) {
		if (a->first.is_match(header)) {
			return a->second != -1 ? 1 : 0;
		}
	}
	// is list is only excludes, no matches means success
	return opt_read_list == 3 ? 1 : 0;
}

// line may be modified if opt_strip_trace is given;
// returns 1 if read is to be read, 0 otherwise

static int get_id_start_noconvert(std::string &line) {
	if (opt_strip_trace && !strip_trace(line)) {		// remove trace
		return 0;
	}
	return is_desired_read(line);
}

// like get_id_start_noconvert(), but also converts the read name from
// FLHBFN1:224:D0WL7ACXX:1:1101:5802:2244 1:N:0:GCCAAT format to
// FLHBFN1_224_D0WL7ACXX_1_1101_5802_2244-R1; will also remove the
// trace from the header line if opt_strip_trace is given

static int get_id_start_convert(std::string &line) {
	if (opt_strip_trace && !strip_trace(line)) {		// remove trace
		return 0;
	}
	// rewrite header
	size_t i(2);
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
	return is_desired_read(line);
}

static void print_usage() {
	std::cerr <<
		"usage: extract_seq_and_qual [opts] <fasta> [<fasta2> ...]\n" <<
		"\t-b\tproduce fastq output instead of fasta/qual\n" <<
		"\t-c\tcomplement output\n" <<
		"\t-h\tprint usage\n" <<
		"\t-i ##\tlist of read names to include (if list contains a comma, it's\n" <<
		"\t\tinterpreted as a comma separated list, otherwise it's treated\n" <<
		"\t\tas a file name; may be specified multiple times)\n" <<
		"\t-L ##\tminimum read length to include\n" <<
		"\t-M ##\tmaximum read length to include\n" <<
		"\t-o ##\tfile to write output to [fasta to stdout, unless -S specified]\n" <<
		"\t-q\tprocess as qual file\n" <<
		"\t-r\ttreat include/exclude read names as regex patterns\n" <<
		"\t-R\tconvert readnames from new Illumina form to old\n" <<
		"\t-s ##\twhen writing output, basepair count to wrap lines at\n" <<
		"\t-S ##\twhen writing output, write to one file for each -i option, named\n" <<
		"\t\tthe same as the -i file but with this parameter as a suffix\n" <<
		"\t-t\tstrip first part of trace id from read headers\n" <<
		"\t-v\tvalidate seq and qual files against each other\n" <<
		"\t-V\tprint version\n" <<
		"\t-x ##\tlist of read names to exclude (see -i for syntax)\n" <<
		"\t-z\tremove trailing zero quality\n" <<
		"\t-Z\tremove leading zero quality\n";
}

static int get_opts(int argc, char **argv, std::list<std::pair<std::string, std::string> > &outputs) {
	int opt_convert_readnames(0);
	// set option defaults
	opt_complement = 0;
	opt_fastq_output = 0;
	opt_line_size = 0;
	opt_min_length = 0;
	opt_max_length = 0;
	opt_qual_only = 0;
	opt_regex = 0;
	opt_strip_leading_zero = 0;
	opt_strip_trace = 0;
	opt_strip_trailing_zero = 0;
	opt_validate = 0;
	std::list<std::pair<std::string, int> > read_list;
	int c;
	while ((c = getopt(argc, argv, "bchi:l:L:M:o:qrRs:S:tvVx:zZ")) != EOF) {
		switch (c) {
		    case 'b':
			opt_fastq_output = 1;
			break;
		    case 'c':
			opt_complement = 1;
			break;
		    case 'h':
			print_usage();
			return 1;
		    case 'l':	// deprecated
		    case 'i':
			if (strchr(optarg, ',')) {	// list of names
				read_list.push_back(std::make_pair(optarg, 1));
			} else {			// files, so glob
				glob_t pglob;
				glob(optarg, GLOB_NOCHECK, 0, &pglob);
				for (size_t i(0); i != pglob.gl_pathc; ++i) {
					read_list.push_back(std::make_pair(pglob.gl_pathv[i], 1));
				}
				globfree(&pglob);
			}
			break;
		    case 'L':
			std::istringstream(optarg) >> opt_min_length;
			break;
		    case 'M':
			std::istringstream(optarg) >> opt_max_length;
			break;
		    case 'o':
			outputs.assign(1, std::make_pair(optarg, ""));
			break;
		    case 'q':
			opt_qual_only = 1;
			break;
		    case 'R':
			opt_convert_readnames = 1;
			break;
		    case 'r':
			opt_regex = 1;
			break;
		    case 's':
			std::istringstream(optarg) >> opt_line_size;
			break;
		    case 'S':
			opt_output_suffix = optarg;
			break;
		    case 't':
			opt_strip_trace = 1;
			break;
		    case 'v':
			opt_validate = 1;
			break;
		    case 'V':
			std::cerr << "extract_seq_and_qual version " << VERSION << "\n";
			exit(0);
			break;
		    case 'x':
			if (strchr(optarg, ',')) {	// list of names
				read_list.push_back(std::make_pair(optarg, 0));
			} else {			// files, so glob
				glob_t pglob;
				glob(optarg, GLOB_NOCHECK, 0, &pglob);
				for (size_t i(0); i != pglob.gl_pathc; ++i) {
					read_list.push_back(std::make_pair(pglob.gl_pathv[i], 0));
				}
				globfree(&pglob);
			}
			break;
		    case 'z':
			opt_strip_trailing_zero = 1;
			break;
		    case 'Z':
			opt_strip_leading_zero = 1;
			break;
		    default:
			throw LocalException("bad option: " + static_cast<char>(c), 1);
		}
	}
	if (optind == argc) {
		throw LocalException("no files specified", 1);
	} else if (opt_strip_leading_zero && opt_strip_trailing_zero) {
		throw LocalException("-z and -Z are mutually exclusive - choose one or the other", 1);
	}
	if (!opt_output_suffix.empty()) {
		if (!outputs.empty()) {
			throw LocalException("-S and -o options are mutually exclusive", 1);
		}
		// convert include read file names to output file names
		std::list<std::pair<std::string, int> >::const_iterator a(read_list.begin());
		const std::list<std::pair<std::string, int> >::const_iterator end_a(read_list.end());
		for (; a != end_a; ++a) {
			if (a->second) {
				outputs.push_back(std::make_pair(a->first + opt_output_suffix, ""));
			}
		}
		if (outputs.empty()) {
			throw LocalException("must give at least one -i with -S", 1);
		}
	} else if (outputs.empty()) {	// default to stdout
		outputs.push_back(std::make_pair("-", ""));
	}
	if (opt_validate && opt_qual_only) {
		std::cerr << "Warning: ignoring -q option (incompatible with -v option)\n";
		opt_qual_only = 0;
	}
	if (opt_validate && opt_fastq_output) {
		std::cerr << "Warning: ignoring -b option (incompatible with -v option)\n";
		opt_fastq_output = 0;
	}
	if (opt_fastq_output) {
		if (opt_qual_only) {
			std::cerr << "Warning: ignoring -q option (incompatible with -b option)\n";
			opt_qual_only = 0;
		}
		if (opt_line_size != 0) {
			std::cerr << "Warning: ignoring -s option (incompatible with -b option)\n";
		}
		opt_line_size = std::string::npos;
	}
	process_read_lists(read_list);
	if (opt_qual_only) {
		std::list<std::pair<std::string, std::string> >::iterator a(outputs.begin());
		const std::list<std::pair<std::string, std::string> >::const_iterator end_a(outputs.end());
		for (; a != end_a; ++a) {
			std::swap(a->first, a->second);
		}
	} else if (!opt_validate && !opt_fastq_output) {
		std::list<std::pair<std::string, std::string> >::iterator a(outputs.begin());
		const std::list<std::pair<std::string, std::string> >::const_iterator end_a(outputs.end());
		for (; a != end_a; ++a) {
			find_qual(a->first, a->second, 1);
		}
	}
	get_id_start = opt_convert_readnames ? get_id_start_convert : get_id_start_noconvert;
	is_desired_read = opt_regex ? is_desired_read_regex : is_desired_read_plain;
	current.initialize();
	return 0;
}

static int check_fastq(const std::string &file) {
	const int fd(open_compressed(file));
	if (fd == -1) {
		return 0;
	}
	char c;
	const int x(pfpeek(fd, &c, 1) == 1 && c == '@');
	if (fd != 0) {			// don't close stdin
		close_compressed(fd);
	}
	return x;
}

// go through arguments to get input files, and deduce associated qual
// file names if needed

static int find_files(int argc, char **argv, int add_qual_files, std::list<std::pair<std::string, std::string> > &file_list) {
	int has_qual_files(opt_qual_only);
	for (; optind < argc; ++optind) {
		std::string fasta(argv[optind]), qual;
		if (opt_qual_only) {
			std::swap(fasta, qual);
			// add suffix, if needed
			std::string dummy;
			find_suffix(qual, dummy);
		} else if (!add_qual_files) {
			// add suffix, if needed
			std::string dummy;
			find_suffix(fasta, dummy);
		} else if (check_fastq(fasta)) {
			has_qual_files = 1;
		} else if (!find_qual(fasta, qual, 0)) {
			if (opt_validate || opt_fastq_output) {
				throw LocalException("could not find qual file for " + fasta);
			} else if (fasta != "-") {
				std::cerr << "Warning: could not find qual file for " << fasta << "\n";
			}
		} else if (!qual.empty()) {
			has_qual_files = 1;
		}
		file_list.push_back(std::make_pair(fasta, qual));
	}
	return has_qual_files;
}

static int get_next_header_fasta(int fd, std::string &line) {
	int file_status(0);
	if (line.empty()) {		// set file_status if it's unclear
		file_status = pfgets(fd, line);
	}
	while (file_status != -1) {
		if (line.size() > 1 && line[0] == '>' && !isspace(line[1]) && get_id_start(line)) {
			return 1;
		}
		file_status = pfgets(fd, line);
	}
	return 0;
}

static int get_next_header_fastq(int fd, std::string &line) {
	int file_status(0);
	if (line.empty()) {		// set file_status if it's unclear
		file_status = pfgets(fd, line);
	}
	while (file_status != -1) {
		if (line.size() < 2 || line[0] != '@' || isspace(line[1])) {
			throw LocalException("bad fastq file: expecting a @ line: " + line);
		}
		if (get_id_start(line)) {
			return 1;
		}
		pfgets(fd, line);			// seq line
		pfgets(fd, line);			// + line
		pfgets(fd, line);			// qual line
		file_status = pfgets(fd, line);		// @ line
	}
	return 0;
}

// process header to handle n expansion (n\d+-\d+:\d+\.(\d+)\.(\d+)$)

static int process_n_header(const std::string &line, size_t i, std::string &seq, std::string &qual, int for_seq) {
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
	if (for_seq == 0) {
		// fasta quality is separated by spaces, so make sure it has one
		if (!qual.empty() && *qual.rbegin() != ' ') {
			qual += " ";
		}
		const std::string s(line.substr(k + 1, end_i - k - 1) + " ");
		qual.reserve(qual.size() + x * s.size());
		for (size_t y(0); y != x; ++y) {
			qual += s;
		}
	} else {
		seq.append(x, 'N');
		if (for_seq == 2) {		// add fastq format quality, too
			unsigned int c;	// using a char here causes problems
			std::istringstream(line.substr(k + 1, end_i - k - 1)) >> c;
			c += 33;
			qual.append(x, c);
		}
	}
	return 1;
}

static int find_next_fastq(int fd, std::string &id, std::string &seq, std::string &qual, std::string &line, size_t &length) {
	if (!get_next_header_fastq(fd, line)) {
		return 0;
	}
	seq.clear();
	qual.clear();
	if (process_n_header(line, 1, seq, qual, 2)) {
		id.clear();	// recombining an n-split: read had a fake id
	} else {
		id = line;
	}
	length = 0;
	for (;;) {
		// fastq format is fairly strict - it comes in 4 line sets,
		// with a @ line, the sequence, a + line, and the quality
		if (pfgets(fd, line) == -1) {
			return 0;
		}
		seq += line;
		if (length == 0) {
			length = line.size();
		}
		if (pfgets(fd, line) == -1) {
			return 0;
		}
		if (line.empty() || line[0] != '+') {
			throw LocalException("bad fastq file: expecting a + line: " + line);
		}
		if (pfgets(fd, line) == -1) {
			return 0;
		}
		qual += line;
		if (pfgets(fd, line) == -1) {
			return 1;
		}
		if (line.size() < 2 || line[0] != '@' || isspace(line[1])) {
			throw LocalException("bad fastq file: expecting a @ line: " + line);
		}
		if (!process_n_header(line, opt_strip_trace ? line.find(' ') : 1, seq, qual, 2)) {
			return 1;
		}
	}
}

static int find_next_seq(int fd, std::string &id, std::string &data, std::string &line, size_t &length) {
	if (!get_next_header_fasta(fd, line)) {
		return 0;
	}
	data.clear();
	std::string dummy;
	if (process_n_header(line, 1, data, dummy, 1)) {
		id.clear();	// recombining an n-split: read had a fake id
	} else {
		id = line;
	}
	length = 0;
	while (pfgets(fd, line) != -1) {
		if (line.empty()) {
		} else if (line[0] != '>') {
			data += line;
			if (length == 0) {
				length = line.size();
			}
		} else if (!process_n_header(line, opt_strip_trace ? line.find(' ') : 1, data, dummy, 1)) {
			return 1;
		}
	}
	return 1;
}

static int find_next_qual(int fd, std::string &id, std::string &data, std::string &line, size_t &length) {
	if (!get_next_header_fasta(fd, line)) {
		return 0;
	}
	data.clear();
	std::string dummy;
	if (process_n_header(line, 1, dummy, data, 0)) {
		id.clear();	// recombining an n-split: read had a fake id
	} else {
		id = line;
	}
	length = 0;
	while (pfgets(fd, line) != -1) {
		if (line.empty()) {
		} else if (line[0] != '>') {
			if (!data.empty() && *data.rbegin() != ' ' && line[0] != ' ') {
				data += " ";
			}
			data += line;
			if (length == 0) {
				length = qual_line_size(line);
			}
		} else if (!process_n_header(line, opt_strip_trace ? line.find(' ') : 1, dummy, data, 0)) {
			return 1;
		}
	}
	return 1;
}

static void output_fastq(const std::string &seq_file, const std::string &qual_file) {
	const int fd_seq(open_compressed(seq_file));
	if (fd_seq == -1) {
		throw LocalException("could not open " + seq_file);
	}
	std::string last_id_seq;
	while (last_id_seq.empty()) {
		if (pfgets(fd_seq, last_id_seq) == -1) {	// empty file?!
			return;
		}
	}
	if (last_id_seq[0] == '@') {		// fastq input format
		std::string id, seq, qual;
		size_t length;
		while (find_next_fastq(fd_seq, id, seq, qual, last_id_seq, length)) {
			if (id.empty()) {
				current.add_fastq(seq, qual);
			} else {
				current.write_fastq();
				current.set_fastq(id, seq, length, qual);
			}
		}
	} else {				// fasta input format
		const int fd_qual(open_compressed(qual_file));
		if (fd_qual == -1) {
			throw LocalException("could not open " + qual_file);
		}
		std::string last_id_qual;
		for (;;) {
			std::string id_seq, id_qual, seq, qual;
			size_t seq_length, qual_length;
			const int found_seq(find_next_seq(fd_seq, id_seq, seq, last_id_seq, seq_length));
			const int found_qual(find_next_qual(fd_qual, id_qual, qual, last_id_qual, qual_length));
			if (found_seq && found_qual) {
			} else if (found_seq) {
				throw LocalException("missing qual: " + id_seq);
			} else if (found_qual) {
				throw LocalException("missing seq: " + id_qual);
			} else {
				break;
			}
			if (id_seq.empty() && id_qual.empty()) {
				current.add_seq(seq);
				current.add_qual(qual);
			} else {
				current.write_fastq();
				current.id_check(id_seq, id_qual);
				current.set_fastq(id_seq, seq, seq_length, qual);
			}
		}
		close_compressed(fd_qual);
	}
	close_compressed(fd_seq);
}

static void output_fasta(const std::string &seq_file, const std::string &qual_file) {
	if (!seq_file.empty()) {
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
			std::string id, seq, qual;
			size_t length;
			while (find_next_fastq(fd, id, seq, qual, last_id, length)) {
				if (id.empty()) {
					current.add_fastq(seq, qual);
				} else {
					current.write_seq();
					current.write_qual();
					current.set_fastq(id, seq, length, qual);
				}
			}
		} else {				// fasta input format
			std::string id, data;
			size_t length;
			while (find_next_seq(fd, id, data, last_id, length)) {
				if (id.empty()) {
					current.add_seq(data);
				} else {
					current.write_seq();
					current.set_seq(id, data, length);
				}
			}
		}
		close_compressed(fd);
	}
	if (!qual_file.empty()) {
		const int fd(open_compressed(qual_file));
		if (fd == -1) {
			throw LocalException("could not open " + qual_file);
		}
		std::string last_id;
		while (last_id.empty()) {
			if (pfgets(fd, last_id) == -1) {	// empty file?!
				break;
			}
		}
		if (last_id.empty()) {			// skip empty file
		} else if (last_id[0] == '@') {		// fastq input format
			std::string id, seq, qual;
			size_t length;
			while (find_next_fastq(fd, id, seq, qual, last_id, length)) {
				// only need quals from fastq
				if (id.empty()) {
					current.add_fastq("", qual);
				} else {
					current.write_qual();
					current.set_fastq(id, "", length, qual);
				}
			}
		} else {				// fasta input format
			std::string id, data;
			size_t length;
			while (find_next_qual(fd, id, data, last_id, length)) {
				if (id.empty()) {
					current.add_seq(data);
				} else {
					current.write_qual();
					current.set_qual(id, data, length);
				}
			}
		}
		close_compressed(fd);
	}
}

static void process_files(const std::list<std::pair<std::string, std::string> > &file_list, const std::list<std::pair<std::string, std::string> > &outputs) {
	std::list<std::pair<std::string, std::string> >::const_iterator a(file_list.begin());
	const std::list<std::pair<std::string, std::string> >::const_iterator end_a(file_list.end());
	if (opt_fastq_output) {
		current.open_outputs(outputs);
		for (; a != end_a; ++a) {
			output_fastq(a->first, a->second);
		}
		// the flush is only done after all files are processed to
		// allow ids to be split across file boundaries (usually
		// because of n-splits)
		current.flush_fastq();
		current.close_outputs(0);
	} else {
		current.open_outputs(outputs, !outputs.front().first.empty());
		for (; a != end_a; ++a) {
			output_fasta(a->first, a->second);
		}
		current.flush_seq();
		current.flush_qual();
		current.close_outputs(0);
	}
}

int main(int argc, char **argv) {
	int had_error(0);
	try {
			// fasta, qual output file pairs
		std::list<std::pair<std::string, std::string> > outputs;
		if (get_opts(argc, argv, outputs)) {
			return 0;
		}
		std::list<std::pair<std::string, std::string> > file_list;
		if (!find_files(argc, argv, opt_validate || opt_fastq_output || !outputs.front().second.empty(), file_list)) {
			// no qual files to read
			std::list<std::pair<std::string, std::string> >::iterator a(outputs.begin());
			const std::list<std::pair<std::string, std::string> >::const_iterator end_a(outputs.end());
			for (; a != end_a; ++a) {
				a->second.clear();
			}
		}
		process_files(file_list, outputs);
		if (opt_validate) {
			std::map<std::string, size_t>::const_iterator a(read_size.begin());
			const std::map<std::string, size_t>::const_iterator end_a(read_size.end());
			for (; a != end_a; ++a) {
				std::cerr << a->first << ": qual missing\n";
			}
		}
	} catch (std::exception &e) {
		std::cerr << "Error: " << e.what() << "\n";
		LocalException *x(dynamic_cast<LocalException *>(&e));
		if (x != NULL && x->show_usage()) {
			print_usage();
		}
		had_error = 1;
	}
	current.close_outputs(had_error);
	return had_error;
}
