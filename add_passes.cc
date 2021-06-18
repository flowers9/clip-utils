#include "open_compressed.h"	// close_compressed(), open_compressed(), pfgets()
#include "pbbam/BamReader.h"	// BamReader
#include "pbbam/BamRecord.h"	// BamRecord
#include "write_fork.h"		// close_fork(), pfputs(), write_fork()
#include <getopt.h>	// getopt(), optarg, optind
#include <iostream>	// cerr
#include <list>		// list<>
#include <map>		// map<>
#include <sstream>	// ostringstream
#include <stdlib.h>	// exit()
#include <string>	// string
#include <sys/types.h>	// size_t

std::string opt_output_file;

static void print_usage(void) {
	std::cerr <<
		"usage: add_passes [-o output] <ccs_bam_file> <fasta/fastq>\n" <<
		"    -o ## file to store output in [stdout]\n";
	exit(0);
}

static void get_passes(const std::string ccs_bam, std::map<std::string, int> &read_passes) {
	// BamReader throws on error
	PacBio::BAM::BamReader f(ccs_bam);
	PacBio::BAM::BamRecord r;
	while (f.GetNext(r)) {
		std::string s(r.FullName());
		// strip /ccs from end
		if (s.length() > 3 && s.substr(s.length() - 4).compare("/ccs") == 0) {
			s = s.substr(0, s.length() - 4);
		} else {
			std::cerr << "Warning: bad read name: " << s << "\n";
		}
		read_passes[s] = r.NumPasses();
	}
}

static void add_pass(std::string &line, const std::map<std::string, int> &read_passes) {
	size_t i(line.find('/', 1));
	if (i == std::string::npos) {
		std::cerr << "Warning: non-pacbio read name1: " << line << "\n";
		return;
	}
	i = line.find('/', i + 1);
	if (i == std::string::npos) {
		std::cerr << "Warning: non-pacbio read name2: " << line << "\n";
		return;
	}
	const std::map<std::string, int>::const_iterator a(read_passes.find(line.substr(1, i - 1)));
	if (a == read_passes.end()) {
		std::cerr << "Warning: read not found: " << line << "\n";
		return;
	}
	if (line.find(' ', 1) == std::string::npos) {
		// got a raw read name, so add space
		line += " passes=";
	} else {	// otherwise, tack onto comments
		line += ";passes=";
	}
	std::ostringstream x;
	x << a->second;
	line += x.str();
}

// handles fasta & fastq files
static void process_fastx(const std::string fastx, const std::map<std::string, int> &read_passes) {
	int fd_in(open_compressed(fastx));
	if (fd_in == -1) {
		std::cerr << "Error: open: " << fastx << "\n";
		exit(1);
	}
	std::list<std::string> fork_args(1, "gzip");
	int fd_out(write_fork(fork_args, opt_output_file));
	if (fd_out == -1) {
		std::cerr << "Error: could not write output file: " << opt_output_file << "\n";
		close_compressed(fd_in);
		exit(1);
	}
	std::string line;
	while (pfgets(fd_in, line) != -1) {
		if (line.length() > 0 && (line[0] == '>' || line[0] == '@')) {
			add_pass(line, read_passes);
			if (line[0] == '@') {	// fastq - have to skip past quality scores
				// 0 - print header, get seq
				// 1 - print seq, get quality header
				// 2 - print quality header, get quality
				for (int i(0); i != 3; ++i) {
					line += "\n";
					pfputs(fd_out, line);
					if (pfgets(fd_in, line) == -1) {
						std::cerr << "Error: reached eof while on pass " << i << " of fastq entry\n";
						close_compressed(fd_in);
						exit(1);
					}
				}
			}
		}
		line += "\n";
		pfputs(fd_out, line);
	}
	close_compressed(fd_in);
	close_fork(fd_out);
}

static void get_opts(int argc, char **argv) {
	opt_output_file = "-";
	int c;
	while ((c = getopt(argc, argv, "o:")) != EOF) {
		switch (c) {
		    case 'o':
			opt_output_file = optarg;
			break;
		    default:
			std::cerr << "Error: bad option: " << c << "\n";
			print_usage();
		}
	}
	if (optind + 2 != argc) {
		std::cerr << "Error: incorrect number of arguments\n";
		print_usage();
	}
}

int main(int argc, char **argv) {
	get_opts(argc, argv);
	std::map<std::string, int> read_passes;
	get_passes(argv[optind], read_passes);
	process_fastx(argv[optind + 1], read_passes);
	return 0;
}
