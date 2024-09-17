// using a barcode lookup table, separates reads in paired (r1/r2)
// fastq files into separate paired fastq files by sequence barcodes

#include "breakup_line.h"	// breakup_line()
#include "open_compressed.h"	// close_compressed(), open_compressed(), pfgets()
#include "write_fork.h"		// write_fork()
#include <algorithm>	// transform()
#include <ctype.h>	// toupper()
#include <exception>	// exception
#include <iostream>	// cerr
#include <list>		// list<>
#include <map>		// map<>
#include <regex>	// regex, regex_search(), gex_constants::optimize, smatch
#include <string>	// string
#include <utility>	// make_pair(), pair<>
#include <vector>	// vector<>

// check for more than one barcode pair matching each sequence, and sort them separately
// #define CHECK_MULTI
// check beyond start of sequence for matches
// #define FULL_SEQ

const std::list<std::string> gzip_args = {"gzip", "-c"};	// for write_fork()
const std::string r1_suffix = ".R1.fastq.gz";			// for output filenames
const std::string r2_suffix = ".R2.fastq.gz";
// store output fds by barcode name in case multiple barcode pairs have the same name
std::map<std::string, size_t> barcode_lookup;			// [barcode_name] = fd_list offset
std::vector<std::pair<int, int>> fd_list;

class LocalException : public std::exception {
    private:
	const std::string error_;
	const int show_usage_;
    public:
	explicit LocalException(const std::string &s) : error_(s), show_usage_(0) { }
	LocalException(const std::string &s, int i) : error_(s), show_usage_(i) { }
	virtual ~LocalException() throw() { }
	virtual const char * what() const throw() {
		return error_.c_str();
	}
	const int &show_usage() const throw() {
		return show_usage_;
	}
};

class BarcodeSubmap {
    public:
	void add(const std::string &name, const std::string &bc2) {
		// multiple barcode pairs can have the same name - open each *once*
		const std::map<std::string, size_t>::const_iterator a = barcode_lookup.find(name);
		if (a == barcode_lookup.end()) {
			const int fd1 = write_fork(gzip_args, name + r1_suffix);
			if (fd1 == -1) {
				throw LocalException("write_fork: " + name + r1_suffix);
			}
			const int fd2 = write_fork(gzip_args, name + r2_suffix);
			if (fd2 == -1) {
				throw LocalException("write_fork: " + name + r2_suffix);
			}
			fd_list.push_back(std::make_pair(fd1, fd2));
			bc2_lookup_[bc2] = barcode_lookup[name] = fd_list.size() - 1;

		} else {
			bc2_lookup_[bc2] = a->second;
		}
	}
	size_t output_offset(const std::string &bc2) const {
		return bc2_lookup_.at(bc2);
	}
	void finalize() {
		std::string bc2_list;
		for (const auto &a : bc2_lookup_) {
			bc2_list += a.first;
			bc2_list += '|';
		}
		bc2_list.pop_back();
		bc2_re_.assign(bc2_list, std::regex_constants::optimize);
	}
	const std::regex &bc2_re() const {
		return bc2_re_;
	}
    private:
	std::map<std::string, size_t> bc2_lookup_;	// [bc2] = fd_list offset
	std::regex bc2_re_;				// regex for all bc2's we can match

};

class FastqEntry {
    public:
	// return -1 if EOF, otherwise returns entry
	bool read(const int fd) {		// fd needs to be from open_compressed()
		if (pfgets(fd, header_) == -1) {
			return 0;
		}
		if (pfgets(fd, seq_) == -1) {
			throw LocalException("read missing sequence: " + header_);
		}
		if (pfgets(fd, qual_header_) == -1) {
			throw LocalException("read missing quality header: " + header_);
		}
		if (pfgets(fd, qual_) == -1) {
			throw LocalException("read missing quality: " + header_);
		}
#ifndef FULL_SEQ
		target_ = seq_.substr(0, 10);
#endif
		return 1;
	}
	void write(const int fd) const {	// fd needs to be from write_fork()
		pfputs(fd, header_);
		pfputc(fd, '\n');
		pfputs(fd, seq_);
		pfputc(fd, '\n');
		pfputs(fd, qual_header_);
		pfputc(fd, '\n');
		pfputs(fd, qual_);
		pfputc(fd, '\n');
	}
#ifdef FULL_SEQ
	bool search(const std::regex &re, const bool first, std::smatch &match) const {
		if (first) {	// start at the beginning
			return std::regex_search(seq_, match, re);
		} else {	// shift one space from start of last match and look again
			return std::regex_search(match[0].first + 1, seq_.cend(), match, re);
		}
	}
#else
	bool search(const std::regex &re, std::smatch &match) const {
		return std::regex_match(target_, match, re);
	}
#endif
	const std::string header() const {
		return header_;
	}
    private:
	std::string header_, seq_, qual_header_, qual_;
#ifndef FULL_SEQ
	std::string target_;
#endif
};

static void print_usage() {
	std::cerr << "usage: barcode_separation <fastq_r1> <fastq_r2> <barcode_file>\n";
}

// read in the barcode list and make the lookups for it, plus open lots of files
static void prepare_barcodes(const std::string &barcode_file, std::map<std::string, BarcodeSubmap> &barcode_dict, std::regex &bc1_re) {
	const int fd(open_compressed(barcode_file));
	if (fd == -1) {
		throw LocalException("could not open " + barcode_file);
	}
	// line format is name, bc1, bc2
	std::vector<std::string> list;
	// parse the file and open all output files
	std::string line;
	while (pfgets(fd, line) != -1) {
		list.clear();
		breakup_line(line, list);
		if (list.size() != 3) {
			throw LocalException("could not parse line: " + barcode_file + ": " + line);
		}
		// uppercase, just in case some lowercase sequence snuck in
		std::transform(list[1].begin(), list[1].end(), list[1].begin(), toupper);
		std::transform(list[2].begin(), list[2].end(), list[2].begin(), toupper);
		barcode_dict[list[1]].add(list[0], list[2]);
	}
	close_compressed(fd);
	if (barcode_dict.empty()) {
		throw LocalException("barcode file contains no barcodes");
	}
	// make regex's for the submaps that match all included bc2's
	// (plus one for all the bc1's)
        std::string bc1_list;
        for (auto &a : barcode_dict) {
                bc1_list += a.first;
                bc1_list += '|';
                a.second.finalize();
        }
        bc1_list.pop_back();
        bc1_re.assign(bc1_list, std::regex_constants::optimize);
}

static void process_sequence(const std::string &reads_1, const std::string &reads_2, const std::map<std::string, BarcodeSubmap> &barcode_dict, const std::regex &bc1_re) {
	const int r1_fd(open_compressed(reads_1));
	if (r1_fd == -1) {
		throw LocalException("could not open " + reads_1);
	}
	const int r2_fd(open_compressed(reads_2));
	if (r2_fd == -1) {
		throw LocalException("could not open " + reads_2);
	}
	const int nm1_fd(write_fork(gzip_args, "no_match" + r1_suffix));
	if (nm1_fd == -1) {
		throw LocalException("could not open no_match" + r1_suffix);
	}
	const int nm2_fd(write_fork(gzip_args, "no_match" + r2_suffix));
	if (nm2_fd == -1) {
		throw LocalException("could not open no_match" + r2_suffix);
	}
	FastqEntry r1_entry, r2_entry;		// read buffers
	std::smatch r1_match, r2_match;		// match buffers
#ifdef CHECK_MULTI
	const int mm1_fd = write_fork(gzip_args, "multi_match" + r1_suffix);
	if (mm1_fd == -1) {
		throw LocalException("could not open multi_match" + r1_suffix);
	}
	const int mm2_fd = write_fork(gzip_args, "multi_match" + r2_suffix);
	if (mm2_fd == -1) {
		throw LocalException("could not open multi_match" + r2_suffix);
	}
	size_t match_offset;
	while (r1_entry.read(r1_fd) && r2_entry.read(r2_fd)) {
		int matches = 0;
		bool first = 1;
		while (matches < 2 && r1_entry.search(bc1_re, first, r1_match)) {
			const BarcodeSubmap &bc1 = barcode_dict.at(r1_match.str());
			if (r2_entry.search(bc1.bc2_re(), 1, r2_match)) {
				const size_t new_match1 = bc1.output_offset(r2_match.str());
				if (matches == 0) {
					++matches;
					match_offset = new_match1;
				} else if (match_offset != new_match1) {
					++matches;
					break;
				}
				while (r2_entry.search(bc1.bc2_re(), 0, r2_match)) {
					const size_t new_match2 = bc1.output_offset(r2_match.str());
					if (match_offset != new_match2) {
						++matches;
						break;
					}
				}
			}
			first = 0;
		}
		switch (matches) {
		    case 1:
			r1_entry.write(fd_list[match_offset].first);
			r2_entry.write(fd_list[match_offset].second);
			break;
		    case 0:			// no matches
			r1_entry.write(nm1_fd);
			r2_entry.write(nm2_fd);
			break;
		    default:			// multiple matches
			r1_entry.write(mm1_fd);
			r2_entry.write(mm2_fd);
		}
	}
#else
#ifdef FULL_SEQ
	NEXT: while (r1_entry.read(r1_fd) && r2_entry.read(r2_fd)) {
		bool first = 1;
		while (r1_entry.search(bc1_re, first, r1_match)) {
			const BarcodeSubmap &bc1 = barcode_dict.at(r1_match.str());
			if (r2_entry.search(bc1.bc2_re(), 1, r2_match)) {
				const std::pair<int, int> &output_fds = fd_list[bc1.output_offset(r2_match.str())];
				r1_entry.write(output_fds.first);
				r2_entry.write(output_fds.second);
				goto NEXT;
			}
			first = 0;
		}
		r1_entry.write(nm1_fd);
		r2_entry.write(nm2_fd);
	}
#else
	while (r1_entry.read(r1_fd) && r2_entry.read(r2_fd)) {
		if (r1_entry.search(bc1_re, r1_match)) {
			const BarcodeSubmap &bc1 = barcode_dict.at(r1_match.str());
			if (r2_entry.search(bc1.bc2_re(), r2_match)) {
				const std::pair<int, int> &output_fds = fd_list[bc1.output_offset(r2_match.str())];
				r1_entry.write(output_fds.first);
				r2_entry.write(output_fds.second);
				continue;
			}
		}
		r1_entry.write(nm1_fd);
		r2_entry.write(nm2_fd);
	}
#endif
#endif
	// close input/output files
	close_compressed(r1_fd);
	close_compressed(r2_fd);
	close_fork(nm1_fd);
	close_fork(nm2_fd);
#ifdef CHECK_MULTI
	close_fork(mm1_fd);
	close_fork(mm2_fd);
#endif
	for (const auto &a : fd_list) {
		close_fork(a.first);
		close_fork(a.second);
	}
}

int main(int argc, char **argv) {
	int had_error(0);
	try {
		if (argc != 4) {
			throw LocalException("incorrect number of parameters", 1);
		}
		// have to convert these to std::string at some point, anyway
		const std::string reads_1(argv[1]);
		const std::string reads_2(argv[2]);
		const std::string barcode_file(argv[3]);
		std::map<std::string, BarcodeSubmap> barcode_dict;
		std::regex bc1_re;
		prepare_barcodes(barcode_file, barcode_dict, bc1_re);
		process_sequence(reads_1, reads_2, barcode_dict, bc1_re);
	} catch (std::exception &e) {
		std::cerr << "Error: " << e.what() << "\n";
		LocalException *x(dynamic_cast<LocalException *>(&e));
		if (x != NULL && x->show_usage()) {
			print_usage();
		}
		had_error = 1;
	}
	return had_error;
}
