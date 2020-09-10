#include "open_compressed.h"	// close_compressed(), open_compressed(), pfgets()
#include <errno.h>	// errno
#include <exception>	// exception
#include <getopt.h>	// getopt(), optarg, optind
#include <iostream>	// cerr
#include <sstream>	// istringstream
#include <stdio.h>	// EOF
#include <string.h>	// strerror()
#include <string>	// string
#include <sys/types.h>	// ssize_t
#include <unistd.h>	// write()
#include <vector>	// vector<>

// this program reads a pair of files and writes N lines of each in succession

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

static int opt_lines;		// number of lines per pass

static void print_usage() {
	std::cerr <<
		"usage: interleave [-l #] <file1a> <file1b> [<file2a> <file2b [...]]\n" <<
		"\t-h\tprint usage\n" <<
		"\t-l ##\tnumber of lines per pass [1]\n";
}

static int get_opts(int argc, char **argv) {
	opt_lines = 1;
	int c;
	while ((c = getopt(argc, argv, "hl:")) != EOF) {
		switch (c) {
		    case 'h':
			print_usage();
			return 1;
		    case 'l':
			std::istringstream(optarg) >> opt_lines;
			break;
		    default:
			throw LocalException("bad option: " + static_cast<char>(c), 1);
		}
	}
	if (optind == argc) {
		throw LocalException("no files specified", 1);
	} else if ((argc - optind) % 2 == 1) {
		throw LocalException("odd number of files specified", 1);
	}
	return 0;
}

// make sure entire line is written; use write() for speed, checking error conditions
static void write_stdout(const std::string &line) {
	const char *buf(line.data());
	ssize_t i(line.size());
	while (i != 0) {
		const ssize_t j(write(1, buf, i));
		if (j == -1) {
			throw LocalException("write failed: " + std::string(strerror(errno)));
		}
		i -= j;
		buf += j;
	}
}

static void interleave(const char* const file1, const char* const file2, std::vector<std::string> &lines1, std::vector<std::string> &lines2) {
	const int f1(open_compressed(file1));
	if (f1 == -1) {
		throw LocalException("could not open " + std::string(file1));
	}
	const int f2(open_compressed(file2));
	if (f2 == -1) {
		throw LocalException("could not open " + std::string(file2));
	}
	int i, j;
	for (;;) {
		for (i = 0; i != opt_lines && pfgets(f1, lines1[i]) != -1; ++i) { }
		for (j = 0; j != opt_lines && pfgets(f2, lines2[j]) != -1; ++j) { }
		if (i != opt_lines || j != opt_lines) {
			break;
		}
		for (i = 0; i != opt_lines; ++i) {
			lines1[i] += "\n";	// pfgets() strings the newline
			write_stdout(lines1[i]);
		}
		for (j = 0; j != opt_lines; ++j) {
			lines2[j] += "\n";	// pfgets() strings the newline
			write_stdout(lines2[j]);
		}
	}
	close_compressed(f1);
	close_compressed(f2);
	if (i == 0 && j == 0) {		// no problems
	} else if (i == 0) {		// various problems
		if (j == opt_lines) {
			throw LocalException("different length files: " + std::string(file2) + " > " + std::string(file1));
		} else {
			throw LocalException("truncated record: " + std::string(file2));
		}
	} else if (j == 0) {
		if (i == opt_lines) {
			throw LocalException("different length files: " + std::string(file1) + " > " + std::string(file2));
		} else {
			throw LocalException("truncated record: " + std::string(file1));
		}
	} else {
		throw LocalException("truncated records: " + std::string(file1) + ", " + std::string(file2));
	}
}

int main(int argc, char **argv) {
	int had_error(0);
	try {
		if (get_opts(argc, argv)) {
			return 0;
		}
		std::vector<std::string> lines1(opt_lines);	// input buffers
		std::vector<std::string> lines2(opt_lines);
		for (; optind != argc; optind += 2) {
			interleave(argv[optind], argv[optind + 1], lines1, lines2);
		}
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
