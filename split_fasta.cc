#include "open_compressed.h"	// close_compressed(), find_suffix(), open_compressed(), pfgets()
#include "write_fork.h"	// close_fork(), pfputs(), write_fork()
#include <ctype.h>	// isspace()
#include <errno.h>	// errno
#include <exception>	// exception
#include <getopt.h>	// getopt(), optarg, optind
#include <glob.h>	// GLOB_NOCHECK, glob(), glob_t
#include <iostream>	// cerr
#include <list>		// list<>
#include <map>		// map<>
#include <stdio.h>	// EOF
#include <string.h>	// strerror()
#include <string>	// string
#include <unistd.h>	// unlink()
#include <vector>	// vector<>

// this program takes a given fasta or qual file and splits it into an
// arbitrary number of files; which reads go into which files are specified
// by files with lists of the desired reads

class LocalException : public std::exception {
    private:
	const std::string error;
	const int show_usage_;
    public:
	explicit LocalException(const std::string &s) : error(s), show_usage_(0) { }
	LocalException(const std::string &s, int i) : error(s), show_usage_(i) { }
	virtual ~LocalException(void) throw() { }
	virtual const char * what() const throw() {
		return error.c_str();
	}
	const int &show_usage(void) const throw() {
		return show_usage_;
	}
};

static int opt_strip_trace(0);
static std::map<std::string, int> read_list;	// read name -> output stream
static std::vector<std::string> output_buffer;

static void print_usage() {
	std::cerr <<
		"usage: split_fasta [opts] -o output_file [input_file1] [input_file2] ...\n" <<
		"\t-i ##\tfile of read names to extract (may be specified multiple times;\n" <<
		"\t\tglobs will be expanded and treated as multiple -i arguments)\n" <<

		"\t-o ##\toutput file name suffix - this is appended to the name of the -i\n" <<
		"\t\tfiles to get the corresponding output file; if a compression\n" <<
		"\t\tsuffix is given, the output will be appropriately compressed\n" <<
		"\t-t\tstrip first part of trace id from query reads\n" <<
		"\n" <<
		"\tif no input files are specified, stdin is read\n";
}

static void get_read_lists(const std::list<std::string> &read_files) {
	std::list<std::string>::const_iterator a(read_files.begin());
	const std::list<std::string>::const_iterator end_a(read_files.end());
	for (int i = 0; a != end_a; ++a, ++i) {
		int fd = open_compressed(*a);
		if (fd == -1) {
			throw LocalException("Error: could not open " + *a + ": " + strerror(errno), 0);
		}
		std::string line;
		while (pfgets(fd, line) != -1) {
			if (read_list.find(line) != read_list.end()) {
				std::cerr << "Warning: read specified in multiple lists: " << line << "\n";
			} else {
				read_list[line] = i;
			}
		}
		close_compressed(fd);
	}
}

static void get_opts(int argc, char **argv, std::list<std::string> &input_files, std::list<std::string> &output_files) {
	std::string output_suffix;
	int c;
	while ((c = getopt(argc, argv, "i:l:o:t")) != EOF) {
		switch(c) {
		    case 'l':	// deprecated; on loan from extract_seq_and_qual
		    case 'i':
			{
				glob_t pglob;
				glob(optarg, GLOB_NOCHECK, 0, &pglob);
				for (size_t i = 0; i != pglob.gl_pathc; ++i) {
					output_files.push_back(pglob.gl_pathv[i]);
				}
			}
			break;
		    case 'o':
			output_suffix = optarg;
			break;
		    case 't':
			opt_strip_trace = 1;
			break;
		    default:
			throw LocalException("bad option: " + (char)c, 1);
		}
	}
	if (output_files.empty()) {
		throw LocalException("no read list files specified", 1);
	} else if (output_suffix.empty()) {
		throw LocalException("-o option not given", 1);
	}
	for (; optind != argc; ++optind) {
		input_files.push_back(argv[optind]);
	}
	if (input_files.empty()) {
		input_files.push_back("");	// default to stdin
	}
	get_read_lists(output_files);
	// convert read file names to output file names
	std::list<std::string>::iterator a(output_files.begin());
	const std::list<std::string>::const_iterator end_a(output_files.end());
	for (; a != end_a; ++a) {
		*a += output_suffix;
	}
}

static int open_input(const std::string &input_file) {
	int fd = open_compressed(input_file);
	if (fd == -1) {
		throw LocalException("Error: could not open " + input_file + ": " + strerror(errno), 0);
	}
	return fd;
}

static void close_outputs(const std::vector<int> &fd_out, int flush = 1) {
	std::vector<int>::const_iterator a(fd_out.begin());
	const std::vector<int>::const_iterator end_a(fd_out.end());
	for (; a != end_a; ++a) {
		if (flush) {
			std::string &buf = output_buffer[*a];
			pfputs(*a, buf);
			buf.clear();
		}
		close_fork(*a);
	}
}

// in addition to opening output descriptors, it also remaps the numbering
// scheme in read_list so the number is the output file descriptor for
// each read

static void open_outputs(std::list<std::string> &output_files, std::vector<int> &fd_out) {
	fd_out.reserve(output_files.size());
	int max_fd = -1;
	std::list<std::string>::iterator a(output_files.begin());
	const std::list<std::string>::const_iterator end_a(output_files.end());
	for (; a != end_a; ++a) {
		std::string suffix;
		find_suffix(*a, suffix);
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
		int fd = write_fork(args, *a);
		if (fd == -1) {
			// clean up any output files we already opened
			close_outputs(fd_out, 0);
			std::list<std::string>::const_iterator b(output_files.begin());
			for (; b != a; ++b) {
				unlink(b->c_str());
			}
			throw LocalException("Error: could not open " + *a + ": " + strerror(errno), 0);
		}
		fd_out.push_back(fd);
		if (max_fd < fd) {
			max_fd = fd;
		}
	}
	output_buffer.resize(max_fd + 1);
	const size_t buf_size((1 << 30) / fd_out.size());
	// now remap read list to use file descriptors
	std::map<std::string, int>::iterator b(read_list.begin());
	const std::map<std::string, int>::const_iterator end_b(read_list.end());
	for (; b != end_b; ++b) {
		b->second = fd_out[b->second];
		output_buffer[b->second].reserve(buf_size);
	}
}

static std::string get_read_name(const std::string &line) {
	if (line.size() < 2) {
		return "";
	}
	std::string::size_type i(2);
	const std::string::size_type end_i(line.size());
	for (; i != end_i && !isspace(line[i]); ++i) { }
	if (!opt_strip_trace) {
		return line.substr(1, i - 1);
	}
	if (i == end_i) {
		return "";
	}
	for (++i; i != end_i && isspace(line[i]); ++i) { }
	if (i == end_i) {
		return "";
	}
	std::string::size_type j(i);
	for (++i; i != end_i && !isspace(line[i]); ++i) { }
	return line.substr(j, i - j);
}

static void process_file(int fd_in) {
	int fd_out(-1);
	std::string line;
	while (pfgets(fd_in, line) != -1) {
		if (line.empty()) {
			continue;
		}
		if (line[0] == '>') {
			std::string name(get_read_name(line));
			const std::map<std::string, int>::iterator a(read_list.find(name));
			if (a == read_list.end()) {
				fd_out = -1;
			} else {
				fd_out = a->second;
				read_list.erase(a);
			}
		}
		if (fd_out != -1) {
			std::string &buf = output_buffer[fd_out];
			if (buf.size() + line.size() + 1 > buf.capacity()) {
				pfputs(fd_out, buf);
				pfputs(fd_out, line + "\n");
				buf.clear();
			} else {
				buf += line + "\n";
			}
		}
	}
}

int main(int argc, char **argv) {
	int had_error(0);
	try {
		std::list<std::string> input_files, output_files;
		get_opts(argc, argv, input_files, output_files);
		std::vector<int> fd_out;
		open_outputs(output_files, fd_out);
		std::list<std::string>::const_iterator a(input_files.begin());
		const std::list<std::string>::const_iterator end_a(input_files.end());
		for (; a != end_a; ++a) {
			const int fd_in(open_input(*a));
			process_file(fd_in);
			close_compressed(fd_in);
		}
		close_outputs(fd_out);
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
