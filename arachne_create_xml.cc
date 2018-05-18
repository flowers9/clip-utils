#include "open_compressed.h"	// close_compressed(), find_suffix(), get_suffix(), open_compressed(), pfgets()
#include "parse_readnames.h"	// ReadNameParser, pick_readname_parser
#include "write_fork.h"	// close_fork(), pfputs(), write_fork()
#include <ctype.h>	// isalnum(), isdigit(), islower(), isupper()
#include <exception>	// exception
#include <getopt.h>	// getopt(), optarg, optind
#include <iostream>	// cerr
#include <sstream>	// istringstream, ostringstream
#include <string>	// string
#include <sys/types.h>	// size_t

class LocalException : public std::exception {
    private:
	const std::string error;
    public:
	LocalException(const std::string &s) : error(s) { }
	virtual ~LocalException(void) throw() { }
	virtual const char * what() const throw() {
		return error.c_str();
	}
};

static int opt_verbose(0), opt_strip_trailing_zero(0), opt_454_3well(0);
static size_t opt_limit(0), opt_split(0);

static std::string get_first_readname(const std::string &fasta) {
	const int fd(open_compressed(fasta));
	if (fd == -1) {
		throw LocalException("could not open " + fasta);
	}
	std::string line;
	if (pfgets(fd, line) == -1) {
		throw LocalException("could not read " + fasta);
	}
	close_compressed(fd);
	return line;
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
	if (name.size() > 1 && name[0] == 'f' && name.find_first_not_of("0123456789", 1) == std::string::npos) {
		qual_file = "q" + name.substr(1);
		if (new_file) {
			qual_file += suffix;
			return 1;
		} else if (find_suffix(qual_file, qual_suffix) == 0) {
			return 1;
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

// sets opt_strip_trailing_zero if needed
static void check_qual_for_trailing_zero(const std::string &fasta, const std::string &qual) {
	// first get length of first read with sequence
	int fd(open_compressed(fasta));
	if (fd == -1) {
		throw LocalException("could not open " + fasta);
	}
	std::string line, data, seq_header;
	while (pfgets(fd, line) != -1) {
		if (line.empty()) {
		} else if (line[0] == '>') {
			if (!data.empty()) {
				break;
			}
			seq_header = line;
		} else {
			data += line;
		}
	}
	close_compressed(fd);
	size_t seq_length(data.size());
	data.clear();
	// next, get quality of first read with quality
	fd = open_compressed(qual);
	if (fd == -1) {
		throw LocalException("could not open " + qual);
	}
	std::string qual_header;
	while (pfgets(fd, line) != -1) {
		if (line.empty()) {
		} else if (line[0] == '>') {
			if (!data.empty()) {
				break;
			}
			qual_header = line;
		} else {
			if (!data.empty() && data[data.size() - 1] != ' ' && line[0] != ' ') {
				data += " ";
			}
			data += line;
		}
	}
	close_compressed(fd);
	if (seq_header != qual_header) {
		// silently fail check if they're different reads - there's
		// no reasonable way to check in this case, so assume the
		// flag isn't needed
		return;
	}
	// now, check to see if the lengths mismatch, and, if they do so,
	// if the quality is exactly one longer and ending in a zero
	size_t qual_length(0);
	std::string::size_type i(0);
	for (; i != data.size() && !isdigit(data[i]); ++i) { }
	std::string::size_type j(i);
	for (; i != data.size(); ++qual_length) {
		j = i;
		for (; i != data.size() && isdigit(data[i]); ++i) { }
		for (; i != data.size() && !isdigit(data[i]); ++i) { }
	}
	if (qual_length == seq_length + 1 && data[j] == '0' && (j + 1 == data.size() || !isdigit(data[j + 1])) && (j == 0 || !isdigit(data[j - 1]))) {
		opt_strip_trailing_zero = 1;
	}
}

static int print_seq(const int fd_out, const std::string &data) {
	std::string::size_type i(0);
	for (; i + 60 < data.size(); i += 60) {
		pfputs(fd_out, data.substr(i, 60) + "\n");
	}
	pfputs(fd_out, data.substr(i) + "\n");
	static size_t count_seq(0);
	++count_seq;
	if (opt_verbose && (count_seq & 0x1FFFF) == 0) {
		std::cerr << "Fasta " << count_seq << "\n";
	}
	if (count_seq == opt_limit) {
		return 2;
	} else if (opt_split != 0 && count_seq % opt_split == 0) {
		return 1;
	} else {
		return 0;
	}
}

// increment the numeric part of .##

static void increment_counter(std::string &s) {
	if (!s.empty()) {
		int x;
		std::istringstream(s.substr(1)) >> x;
		++x;
		std::ostringstream t;
		t << x;
		s.replace(1, s.size() - 1, t.str());
	}
}

static void open_files(int &fd_out, int &fd_xml, const std::string &file) {
	std::list<std::string> args;
	args.push_back("gzip");
	args.push_back("-c");
	fd_out = write_fork(args, file + ".fasta.gz");
	if (fd_out == -1) {
		throw LocalException("could not open " + file + ".fasta.gz");
	}
	fd_xml = write_fork(args, file + ".xml.gz");
	if (fd_xml == -1) {
		throw LocalException("could not open " + file + ".xml.gz");
	}
	pfputs(fd_xml, "<?xml version=\"1.0\"?>\n<trace_volume>\n");
}

static void close_files(const int fd_out, const int fd_xml) {
	pfputs(fd_xml, "</trace_volume>\n");
	close_fork(fd_xml);
	close_fork(fd_out);
}

static void print_fasta(const std::string &fasta, ReadNameParser &parser) {
	const int fd(open_compressed(fasta));
	if (fd == -1) {
		throw LocalException("could not open " + fasta);
	}
	std::string suffix(opt_split == 0 ? "" : ".0");
	int fd_out, fd_xml;
	open_files(fd_out, fd_xml, fasta + suffix);
	std::string line, data;
	while (pfgets(fd, line) != -1) {
		if (line.empty()) {
		} else if (line[0] == '>' && parser.parse(line.substr(1))) {
			if (!data.empty()) {
				const int x(print_seq(fd_out, data));
				if (x == 0) {
				} else if (x == 1) {
					close_files(fd_out, fd_xml);
					increment_counter(suffix);
					open_files(fd_out, fd_xml, fasta + suffix);
				} else if (x == 2) {
					data.clear();
					break;
				}
				data.clear();
			}
			pfputs(fd_out, ">gnl|ti|15447 " + parser.trace() + "\n");
			pfputs(fd_xml, "<trace>\n<CENTER_NAME>SHGC</CENTER_NAME>\n<CHEMISTRY_TYPE>T</CHEMISTRY_TYPE>\n<PLATE_ID>unknown</PLATE_ID>\n<PROGRAM_ID>PHRED-0.961028.I</PROGRAM_ID>\n<RUN_LANE>5</RUN_LANE>\n<SOURCE_TYPE>G</SOURCE_TYPE>\n<SPECIES_CODE>HOMO SAPIENS</SPECIES_CODE>\n<SUBMISSION_TYPE>UPDATE</SUBMISSION_TYPE>\n <SUBSPECIES_ID>JULIO</SUBSPECIES_ID>\n<SVECTOR_CODE>POT</SVECTOR_CODE>\n  <TEMPLATE_ID>" + parser.id() + "</TEMPLATE_ID>\n<TI>11394</TI>\n<TRACE_DIRECTION>" + parser.direction() + "</TRACE_DIRECTION>\n <TRACE_END>" + parser.direction() + "</TRACE_END>\n<TRACE_FORMAT>SCF</TRACE_FORMAT>\n<TRACE_NAME>" + parser.trace() + "</TRACE_NAME>\n <TRACE_TYPE_CODE>WGS</TRACE_TYPE_CODE>\n<WELL_ID>" + parser.well() + "</WELL_ID>\n</trace>\n");
		} else if (line[0] == '>') {
			throw LocalException("unable to parse trace name: " + line);
		} else {
			data += line;
		}
	}
	if (!data.empty()) {
		print_seq(fd_out, data);
	}
	close_files(fd_out, fd_xml);
	close_compressed(fd);
}

// data is passed non-const so it can strip trailing zero if needed
static int print_qual(int fd_out, std::string &data) {
	std::string::size_type i;
	if (opt_strip_trailing_zero) {
		for (i = data.size() - 1; i != std::string::npos && !isdigit(data[i]); --i) { }
		for (; i != std::string::npos && data[i] == '0'; --i) { }
		if (i == std::string::npos || !isdigit(data[i])) {
			data.resize(i + 1);
		}
	}
	for (i = 0; i != data.size() && !isdigit(data[i]); ++i) { }
	while (i != data.size()) {
		std::string::size_type j(i), l(i);
		int k(0);
		for (; i != data.size() && k != 60; ++k) {
			for (; i != data.size() && isdigit(data[i]); ++i) { }
			l = i;
			for (; i != data.size() && !isdigit(data[i]); ++i) { }
		}
		if (k != 0) {
			pfputs(fd_out, data.substr(j, l - j) + "\n");
		}
	}
	static size_t count_qual(0);
	++count_qual;
	if (opt_verbose && (count_qual & 0x1FFFF) == 0) {
		std::cerr << "Qual " << count_qual << "\n";
	}
	if (count_qual == opt_limit) {
		return 2;
	} else if (opt_split != 0 && count_qual % opt_split == 0) {
		return 1;
	} else {
		return 0;
	}
}

static void print_qualfile(const std::string &fasta, ReadNameParser &parser, const std::string &qual) {
	const int fd(open_compressed(qual));
	if (fd == -1) {
		throw LocalException("could not open " + qual);
	}
	std::string suffix(opt_split == 0 ? "" : ".0");
	std::list<std::string> args;
	args.push_back("gzip");
	args.push_back("-c");
	int fd_out = write_fork(args, fasta + suffix + ".fasta.qual.gz");
	if (fd_out == -1) {
		throw LocalException("could not open " + fasta + suffix + ".fasta.qual.gz");
	}
	std::string line, data;
	while (pfgets(fd, line) != -1) {
		if (line.empty()) {
		} else if (parser.parse(line)) {
			if (!data.empty()) {
				const int x(print_qual(fd_out, data));
				if (x == 0) {
				} else if (x == 1) {
					close_fork(fd_out);
					increment_counter(suffix);
					fd_out = write_fork(args, fasta + suffix + ".fasta.qual.gz");
				} else if (x == 2) {
					data.clear();
					break;
				}
				data.clear();
			}
			pfputs(fd_out, ">gnl|ti|15447 " + parser.trace() + "\n");
		} else if (line[0] == '>') {
			throw LocalException("unable to parse trace name: " + line);
		} else {
			if (!data.empty() && data[data.size() - 1] != ' ' && line[0] != ' ') {
				data += " ";
			}
			data += line;
		}
	}
	if (!data.empty()) {
		print_qual(fd_out, data);
	}
	close_fork(fd_out);
	close_compressed(fd);
}

static void print_usage() {
	std::cerr <<
		"usage: arachne_create_xml [opts] <fasta> <newlib>\n"
		"\t-3\tconvert 454 wells from alphanumeric to number\n"
		"\t-l ##\tonly print first ## reads\n"
		"\t-s ##\tsplit output files into ## reads\n"
		"\t-v\tgive user feedback\n";
}

static int get_opts(int argc, char **argv) {
	opt_454_3well = 0;
	opt_limit = 0;
	opt_split = 0;
	opt_verbose = 0;
	int c;
	while ((c = getopt(argc, argv, "3l:s:v")) != -1) {
		switch(c) {
		    case '3':
			opt_454_3well = 1;
			break;
		    case 'l':
			std::istringstream(optarg) >> opt_limit;
			break;
		    case 's':
			std::istringstream(optarg) >> opt_split;
			break;
		    case 'v':
			opt_verbose = 1;
			break;
		   default:
			print_usage();
			return 1;
		}
	}
	if (argc != optind + 2) {	// not enough parameters
		print_usage();
		return 1;
	}
	return 0;
}

int main(int argc, char **argv) {
	if (get_opts(argc, argv)) {
		return 1;
	}
	std::string fasta(argv[optind]);
	const std::string lib(argv[++optind]);
	std::string qual;
	if (!find_qual(fasta, qual, 0)) {
		std::cerr << "Error: could not find qual file for " << fasta << "\n";
		return 1;
	}
	if (opt_verbose) {
		std::cerr << "Using " << qual << "\n";
	}
	try {
		check_qual_for_trailing_zero(fasta, qual);
		const std::string read(get_first_readname(fasta));
		ReadNameParser * const parser(pick_readname_parser(read, opt_454_3well));
		parser->reset_filename(lib);
		print_fasta(fasta, *parser);
		print_qualfile(fasta, *parser, qual);
		delete parser;
	} catch (std::exception &e) {
		std::cerr << "Error: " << e.what() << "\n";
		return 1;
	}
	return 0;
}
