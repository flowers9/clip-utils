// combinatorial inline barcode demultiplexer

// using a barcode lookup table, separates reads in paired (r1/r2)
// fastq files into separate paired fastq files by sequence barcodes

// barcodes must be exactly 10 basepairs long, as this looks only at the
// leading 10 basepairs of each read, but is pretty fast

// barcode file format is
//	barcode_name r1_barcode_sequence r2_barcode_sequence
// barcode names can appear multiple times, and all matching barcode
// pairs will go into the same file

#include "breakup_line.h"	// breakup_line()
#include "open_compressed.h"	// close_compressed(), open_compressed(), pfgets()
#include "write_fork.h"		// clsoe_fork, close_fork_wait(), write_fork()
#include <algorithm>	// transform()
#include <ctype.h>	// toupper()
#include <exception>	// exception
#include <iostream>	// cerr
#include <list>		// list<>
#include <map>		// map<>
#include <string>	// string
#include <utility>	// make_pair(), pair<>
#include <vector>	// vector<>

const std::list<std::string> gzip_args = {"gzip", "-c"};	// for write_fork()
const std::string r1_suffix = ".R1.fastq.gz";			// for output filenames
const std::string r2_suffix = ".R2.fastq.gz";
std::map<std::string, std::pair<int, int>> barcode_name_fds;	// [barcode_name] = output_fds

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

class FastqEntry {
    public:
	// returns true if entry is read
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
	std::string lead_seq() const {
		return seq_.substr(0, 10);
	}
    private:
	std::string header_, seq_, qual_header_, qual_;
};

static void print_usage() {
	std::cerr << "usage: barcode_separation <fastq_r1> <fastq_r2> <barcode_file>\n";
}

// read in the barcode list and make the lookups for it, plus open lots of files
static void prepare_barcodes(const std::string &barcode_file, std::map<std::string, std::pair<int, int>> &barcode_fds) {
	const int fd(open_compressed(barcode_file));
	if (fd == -1) {
		throw LocalException("could not open " + barcode_file);
	}
	// line format is "barcode_name r1_barcode r2_barcode"
	std::vector<std::string> list;
	// parse the file and open all output files
	std::string line;
	while (pfgets(fd, line) != -1) {
		list.clear();
		breakup_line(line, list);
		if (list.size() != 3) {
			throw LocalException("could not parse line: " + barcode_file + ": " + line);
		}
		if (list[1].size() != 10 || list[2].size() != 10) {
			throw LocalException("barcode length != 10: " + barcode_file + ": " + line);
		}
		// uppercase, just in case some lowercase sequence snuck in
		std::transform(list[1].begin(), list[1].end(), list[1].begin(), toupper);
		std::transform(list[2].begin(), list[2].end(), list[2].begin(), toupper);
		const std::string bc1bc2 = list[1] + list[2];
		if (barcode_fds.find(bc1bc2) != barcode_fds.end()) {
			throw LocalException("duplicate barcode pair: " + barcode_file + ": " + line);
		}
		// multiple barcode pairs can have the same name - open each *once*
		const auto a = barcode_name_fds.find(list[0]);
		if (a == barcode_name_fds.end()) {
			const int fd1 = write_fork(gzip_args, list[0] + r1_suffix);
			if (fd1 == -1) {
				throw LocalException("write_fork: " + list[0] + r1_suffix);
			}
			const int fd2 = write_fork(gzip_args, list[0] + r2_suffix);
			if (fd2 == -1) {
				throw LocalException("write_fork: " + list[0] + r2_suffix);
			}
			barcode_fds[bc1bc2] = barcode_name_fds[list[0]] = std::make_pair(fd1, fd2);
		} else {
			barcode_fds[bc1bc2] = a->second;
		}
	}
	close_compressed(fd);
	if (barcode_name_fds.empty()) {
		throw LocalException("barcode file contains no barcodes");
	}
}

static void process_sequence(const std::string &reads_1, const std::string &reads_2, const std::map<std::string, std::pair<int, int>> &barcode_fds) {
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
	while (r1_entry.read(r1_fd) && r2_entry.read(r2_fd)) {
		std::map<std::string, std::pair<int, int>>::const_iterator a = barcode_fds.find(r1_entry.lead_seq() + r2_entry.lead_seq());
		if (a != barcode_fds.end()) {
			r1_entry.write(a->second.first);
			r2_entry.write(a->second.second);
		} else {
			r1_entry.write(nm1_fd);
			r2_entry.write(nm2_fd);
		}
	}
	// close input/output files
	close_compressed(r1_fd);
	close_compressed(r2_fd);
	for (const auto &a : barcode_name_fds) {
		close_fork(a.second.first);
		close_fork(a.second.second);
	}
	close_fork(nm1_fd);
	close_fork_wait(nm2_fd);
}

int main(int argc, char **argv) {
	int had_error = 0;
	try {
		if (argc != 4) {
			throw LocalException("incorrect number of parameters", 1);
		}
		// have to convert these to std::string at some point, anyway
		const std::string reads_1 = argv[1];
		const std::string reads_2 = argv[2];
		const std::string barcode_file = argv[3];
		std::map<std::string, std::pair<int, int>> barcode_fds;	// [bc1bc2] = output_fds
		prepare_barcodes(barcode_file, barcode_fds);
		process_sequence(reads_1, reads_2, barcode_fds);
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
