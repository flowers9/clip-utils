// using a barcode lookup table, separates reads in paired (r1/r2)
// fastq files into separate paired fastq files by sequence barcodes;
// r1 barcodes are only matched against the 10 bp start of the sequence

#include "breakup_line.h"	// breakup_line()
#include "open_compressed.h"	// close_compressed(), open_compressed(), pfgets()
#include "write_fork.h"		// write_fork()
#include <algorithm>	// transform()
#include <ctype.h>	// toupper()
#include <exception>	// exception
#include <iostream>	// cerr
#include <list>		// list<>
#include <map>		// map<>
#include <regex>	// regex, regex_search(), smatch
#include <string>	// string
#include <vector>	// vector<>

const std::list<std::string> gzip_args({"gzip", "-c"});	// for write_fork()
const std::string r1_suffix(".R1.fastq.gz");		// for output filenames
const std::string r2_suffix(".R2.fastq.gz");

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

class BarcodeSubmap {
    public:
	std::map<std::string, std::vector<int> > bc2;	// [bc2] = fd_r[12]
	std::regex r2bc_re;			// regex for all bc2's we can match
    public:
	void open(const std::string &name, const std::string &r2_bc) {
		std::vector<int> &fd_list(bc2[r2_bc]);
		fd_list.push_back(write_fork(gzip_args, name + r1_suffix));
		if (fd_list.back() == -1) {
			throw LocalException("write_fork: " + name + r1_suffix);
		}
		fd_list.push_back(write_fork(gzip_args, name + r2_suffix));
		if (fd_list.back() == -1) {
			throw LocalException("write_fork: " + name + r2_suffix);
		}
	}
	void make_re() {
		std::map<std::string, std::vector<int> >::const_iterator b(bc2.begin());
		const std::map<std::string, std::vector<int> >::const_iterator end_b(bc2.end());
		std::string r2bc_list(b->first);
		for (++b; b != end_b; ++b) {
			r2bc_list += '|';
			r2bc_list += b->first;
		}
		r2bc_re.assign(r2bc_list);
	}
	void close() const {
		std::map<std::string, std::vector<int> >::const_iterator a(bc2.begin());
		const std::map<std::string, std::vector<int> >::const_iterator end_a(bc2.end());
		for (; a != end_a; ++a) {
			close_fork(a->second[0]);
			close_fork(a->second[1]);
		}
	}
};

class FastqEntry {
    private:
	std::string target_;
    public:
	std::string header, seq, qual_header, qual;
    public:
	// return -1 if EOF, otherwise returns entry
	int read(const int fd) {		// fd needs to be from open_compressed()
		if (pfgets(fd, header) == -1) {
			return -1;
		}
		if (pfgets(fd, seq) == -1) {
			return -1;
		}
		if (pfgets(fd, qual_header) == -1) {
			return -1;
		}
		if (pfgets(fd, qual) == -1) {
			return -1;
		}
		// we need a permanent string for regex_search(), so make it here
		target_ = seq.substr(0, 10);
		return 0;
	}
	void write(const int fd) const {	// fd needs to be from write_fork()
		pfputs(fd, header);
		pfputc(fd, '\n');
		pfputs(fd, seq);
		pfputc(fd, '\n');
		pfputs(fd, qual_header);
		pfputc(fd, '\n');
		pfputs(fd, qual);
		pfputc(fd, '\n');
	}
	bool search(const std::regex &re, std::smatch &match) const {
		return std::regex_search(target_, match, re);
	}
};

static void print_usage(void) {
	std::cerr << "usage: barcode_separation <fastq_r1> <fastq_r2> <barcode_file>\n";
}

// read in the barcode list and make the lookups for it, plus open lots of files
static void prepare_barcodes(const std::string &barcode_file, std::map<std::string, BarcodeSubmap> &barcode_dict, std::regex &r1bc_re) {
	const int fd(open_compressed(barcode_file));
	if (fd == -1) {
		throw LocalException("could not open " + barcode_file);
	}
	// parse the file and open all output files
	std::string line;
	while (pfgets(fd, line) != -1) {
		// line format is name, bc1, bc2
		std::vector<std::string> list;
		breakup_line(line, list);
		if (list.size() != 3) {
			throw LocalException("could not parse line: " + barcode_file + ": " + line);
		}
		// uppercase, just in case some lowercase sequence snuck in
		std::transform(list[1].begin(), list[1].end(), list[1].begin(), toupper);
		std::transform(list[2].begin(), list[2].end(), list[2].begin(), toupper);
		barcode_dict[list[1]].open(list[0], list[2]);
	}
	close_compressed(fd);
	// make regex's for the submaps that match all included bc2's
	// (plus one for all the bc1's)
	std::map<std::string, BarcodeSubmap>::iterator a(barcode_dict.begin());
	const std::map<std::string, BarcodeSubmap>::const_iterator end_a(barcode_dict.end());
	std::string r1bc_list;
	for (; a != end_a; ++a) {
		a->second.make_re();
		if (a != barcode_dict.begin()) {
			r1bc_list += '|';
		}
		r1bc_list += a->first;
	}
	r1bc_re.assign(r1bc_list);
}

static void process_sequence(const std::string &reads_1, const std::string &reads_2, const std::map<std::string, BarcodeSubmap> &barcode_dict, const std::regex &r1bc_re) {
	const int r1_fd(open_compressed(reads_1));
	if (r1_fd == -1) {
		throw LocalException("could not open " + reads_1);
	}
	const int r2_fd(open_compressed(reads_2));
	if (r2_fd == -1) {
		throw LocalException("could not open " + reads_2);
	}
	const int nu1_fd(write_fork(gzip_args, "newUndetermined.R1.fastq.gz"));
	if (nu1_fd == -1) {
		throw LocalException("could not open newUndetermined.R1.fastq.gz");
	}
	const int nu2_fd(write_fork(gzip_args, "newUndetermined.R2.fastq.gz"));
	if (nu2_fd == -1) {
		throw LocalException("could not open newUndetermined.R2.fastq.gz");
	}
	FastqEntry r1_entry, r2_entry;		// read buffers
	std::smatch match;			// match buffer
	for (;;) {
		if (r1_entry.read(r1_fd) == -1 || r2_entry.read(r2_fd) == -1) {
			break;
		}
		if (r1_entry.search(r1bc_re, match)) {
			const BarcodeSubmap &bc1(barcode_dict.at(*match.begin()));
			if (r2_entry.search(bc1.r2bc_re, match)) {
				const std::vector<int> &fds(bc1.bc2.at(*match.begin()));
				r1_entry.write(fds[0]);
				r2_entry.write(fds[1]);
				continue;
			}
		}
		// could not find bc1 or bc2
		r1_entry.write(nu1_fd);
		r2_entry.write(nu2_fd);
	}
	// close input/output files
	close_compressed(r1_fd);
	close_compressed(r2_fd);
	std::map<std::string, BarcodeSubmap>::const_iterator a(barcode_dict.begin());
	const std::map<std::string, BarcodeSubmap>::const_iterator end_a(barcode_dict.end());
	for (; a != end_a; ++a) {
		a->second.close();
	}
	close_fork(nu1_fd);
	close_fork_wait(nu2_fd);
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
		std::regex r1bc_re;
		prepare_barcodes(barcode_file, barcode_dict, r1bc_re);
		process_sequence(reads_1, reads_2, barcode_dict, r1bc_re);
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
