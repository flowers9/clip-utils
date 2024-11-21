// using a barcode lookup table, separates reads in a fastq file into
// separate fastq files by sequence barcodes

// XXX - does not search for revcomp matches
// XXX - could add option to toggle checking for multiple matches

#include "breakup_line.h"	// breakup_line()
#include "open_compressed.h"	// close_compressed(), open_compressed(), pfgets()
#include "write_fork.h"		// close_fork(), close_fork_wait(), pfputc(), pfputs(), write_fork()
#include <exception>	// exception
#include <iostream>	// cerr
#include <list>		// list<>
#include <map>		// map<>
#include <regex>	// regex, regex_search(), regex_constants::optimize, smatch
#include <string>	// string
#include <vector>	// vector<>

const std::list<std::string> gzip_args({"gzip", "-c"});	// for write_fork()
// store output fds by barcode name in case multiple barcode pairs have the same name
std::map<std::string, int> output_fds;			// [barcode_name] = fd_out

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

// all the 5 prime ends (and associated output files) for a given 3 prime barcode
class BarcodeSubmap {
    private:
	static unsigned char comp(unsigned char c) {
		switch (c) {
		    case 'A':
			return 'T';
		    case 'C':
			return 'G';
		    case 'G':
			return 'C';
		    case 'T':
			return 'A';
		    default:
			return c;
		}
	}
    public:
	void add(const std::string &name, std::string &p5_bc) {
		// revcomp 5 prime barcodes
		std::string::iterator b = p5_bc.begin();
		std::string::iterator c = p5_bc.end() - 1;
		for (; b < c; ++b, --c) {
			const unsigned char d = comp(*b);
			*b = comp(*c);
			*c = d;
		}
		if (b == c) {
			*b = comp(*b);
		}
		if (p5_fds_.find(p5_bc) != p5_fds_.end()) {
			throw LocalException("duplicate 5' barcode (" + name + "): " + p5_bc);
		}
		// multiple barcode pairs can have the same name - open each *once*
		const std::map<std::string, int>::const_iterator a = output_fds.find(name);
		if (a == output_fds.end()) {
			const int fd = write_fork(gzip_args, name + ".fastq.gz");
			if (fd == -1) {
				throw LocalException("write_fork: " + name + ".fastq.gz");
			}
			output_fds[name] = fd;
			p5_fds_[p5_bc] = fd;
		} else {
			p5_fds_[p5_bc] = a->second;
		}
	}
	int output_fd(const std::string &p5_bc) const {
		return p5_fds_.at(p5_bc);
	}
	// make re for all 5 prime barcodes that this can pair with
	void finalize() {
		std::string p5_list;
		for (const auto &a : p5_fds_) {
			p5_list += a.first;
			p5_list += '|';
		}
		p5_list.pop_back();
		p5_re_.assign(p5_list, std::regex_constants::optimize);
	}
	const std::regex &p5_re() const {
		return p5_re_;
	}
    private:
	std::map<std::string, int> p5_fds_;	// [name] = output_fd
	std::regex p5_re_;			// regex for all 5 prime barcodes we can match
};

class FastqEntry {
    public:
	bool read(const int fd) {		// fd needs to be from open_compressed()
		if (pfgets(fd, header_) == -1) {	// check for EOF
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
	bool search_3p(const std::regex &re, const bool first, std::smatch &p3_match) const {
		if (first) {	// start at the beginning
			return std::regex_search(seq_, p3_match, re);
		} else {	// shift one space from start of last match and look again
			return std::regex_search(p3_match[0].first + 1, seq_.cend(), p3_match, re);
		}
	}
	// look for a 5 prime barcode following the 3 prime barcode
	bool search_5p(const std::regex &re, std::smatch &p5_match, const std::smatch &p3_match) const {
		return std::regex_search(p3_match[0].second, seq_.cend(), p5_match, re);
	}
	// continue looking for other 5 prime matches
	bool search_5p(const std::regex &re, std::smatch &p5_match) const {
		return std::regex_search(p5_match[0].first + 1, seq_.cend(), p5_match, re);
	}
    private:
	std::string header_, seq_, qual_header_, qual_;
};

static void print_usage(void) {
	std::cerr << "usage: barcode_separation <fastq> <barcode_file>\n";
}

// read in the barcode list and make the lookups for it, plus open all output files
// format for barcode file: name 3prime_barcode 5prime_barcode
static void prepare_barcodes(const std::string &barcode_file, std::map<std::string, BarcodeSubmap> &barcode_dict, std::regex &p3_re) {
	const int fd = open_compressed(barcode_file);
	if (fd == -1) {
		throw LocalException("could not open " + barcode_file);
	}
	// parse the file and open all output files
	std::vector<std::string> list;
	std::string line;
	while (pfgets(fd, line) != -1) {
		list.clear();
		breakup_line(line, list);
		if (list.size() != 3) {
			throw LocalException("could not parse line: " + barcode_file + ": " + line);
		}
		barcode_dict[list[1]].add(list[0], list[2]);
	}
	close_compressed(fd);
	if (barcode_dict.empty()) {
		throw LocalException("barcode file contains no barcodes");
	}
	// make regex's for the submaps that match all included 5 prime barcodes
	// (plus one for all the 3 prime barcodes)
	std::string p3_list;
	for (auto &a : barcode_dict) {
		p3_list += a.first;
		p3_list += '|';
		a.second.finalize();
	}
	p3_list.pop_back();
	p3_re.assign(p3_list, std::regex_constants::optimize);
}

static void process_sequence(const std::string &reads, const std::map<std::string, BarcodeSubmap> &barcode_dict, const std::regex &p3_re) {
	const int reads_fd = open_compressed(reads);
	if (reads_fd == -1) {
		throw LocalException("could not open " + reads);
	}
	const int nomatch_fd = write_fork(gzip_args, "no_match.fastq.gz");
	if (nomatch_fd == -1) {
		throw LocalException("could not open no_match.fastq.gz");
	}
	const int multimatch_fd = write_fork(gzip_args, "multi_match.fastq.gz");
	if (multimatch_fd == -1) {
		throw LocalException("could not open multi_match.fastq.gz");
	}
	FastqEntry entry;			// read buffer
	std::smatch p3_match, p5_match;		// match buffers
	std::vector<int> matches;
	while (entry.read(reads_fd)) {
		matches.clear();
		bool first = 1;
		while (entry.search_3p(p3_re, first, p3_match)) {
			const BarcodeSubmap &p3_entry = barcode_dict.at(p3_match.str());
			if (entry.search_5p(p3_entry.p5_re(), p5_match, p3_match)) {
				matches.push_back(p3_entry.output_fd(p5_match.str()));
				while (entry.search_5p(p3_entry.p5_re(), p5_match)) {
					matches.push_back(p3_entry.output_fd(p5_match.str()));
				}
			}
			first = 0;
		}
		// could not find a matched pair of barcodes
		if (matches.size() == 1) {
			entry.write(matches[0]);
		} else if (matches.size() == 0) {
			entry.write(nomatch_fd);
		} else {
			entry.write(multimatch_fd);
		}
	}
	// close input/output files
	close_compressed(reads_fd);
	close_fork(multimatch_fd);
	for (const auto &a : output_fds) {
		close_fork(a.second);
	}
	close_fork_wait(nomatch_fd);
}

int main(int argc, char **argv) {
	int had_error = 0;
	try {
		if (argc != 3) {
			throw LocalException("incorrect number of parameters", 1);
		}
		const std::string reads = argv[1];
		const std::string barcode_file = argv[2];
		std::map<std::string, BarcodeSubmap> barcode_dict;
		std::regex p3_re;
		prepare_barcodes(barcode_file, barcode_dict, p3_re);
		process_sequence(reads, barcode_dict, p3_re);
	} catch (std::exception &e) {
		std::cerr << "Error: " << e.what() << "\n";
		LocalException *x = dynamic_cast<LocalException *>(&e);
		if (x != NULL && x->show_usage()) {
			print_usage();
		}
		had_error = 1;
	}
	return had_error;
}
