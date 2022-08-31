// add quality scores from an associated bam file to a convert a fasta
// file to a fastq file; also adds pass data if not already present

// (also verifies sequence, quality, and pass information against bam
// file, if already present)

#include "open_compressed.h"	// close_compressed(), open_compressed(), pfgets()
#include "pbbam/BamReader.h"	// BamReader
#include "pbbam/BamRecord.h"	// BamRecord
#include "write_fork.h"		// close_fork(), pfputs(), write_fork()
#include <algorithm>	// equal()
#include <getopt.h>	// getopt(), optarg, optind
#include <iostream>	// cerr
#include <list>		// list<>
#include <map>		// map<>
#include <sstream>	// istringstream, ostringstream
#include <stdlib.h>	// exit()
#include <string>	// string
#include <sys/types.h>	// int32_t, size_t

static std::string opt_output_file;

static std::string comp_lookup(256, 0);

// lookup table is to speed up comping reads; here's where we set it up
static void init_comp_lookup() {
	comp_lookup[static_cast<int>('a')] = 'T';
	comp_lookup[static_cast<int>('A')] = 'T';
	comp_lookup[static_cast<int>('c')] = 'G';
	comp_lookup[static_cast<int>('C')] = 'G';
	comp_lookup[static_cast<int>('g')] = 'C';
	comp_lookup[static_cast<int>('G')] = 'C';
	comp_lookup[static_cast<int>('t')] = 'A';
	comp_lookup[static_cast<int>('T')] = 'A';
}

static int comp_not_equal(std::string::const_reverse_iterator a, const std::string::const_reverse_iterator end_a, std::string::const_iterator b) {
	for (; a != end_a; ++a, ++b) {
		if (*a != comp_lookup[static_cast<int>(*b)]) {
			return 1;
		}
	}
	return 0;
}

class Read {
    public:
	int32_t passes;
	std::string seq, qual;
	explicit Read() : passes(-1) { }
	~Read() { }
	void clear() {
		passes = -1;
		seq.clear();
		qual.clear();
	}
	void set_and_comp_seq(const std::string& a) {
		seq.reserve(a.size());
		std::string::const_reverse_iterator b(a.rbegin());
		const std::string::const_reverse_iterator end_b(a.rend());
		for (; b != end_b; ++b) {
			seq.push_back(comp_lookup[static_cast<int>(*b)]);
		}
	}
	void set_and_comp_qual(const std::string& a) {
		qual.reserve(a.size());
		std::string::const_reverse_iterator b(a.rbegin());
		const std::string::const_reverse_iterator end_b(a.rend());
		for (; b != end_b; ++b) {
			qual.push_back(*b);
		}
	}
};

class FullRead : Read {
    private:
	void parse_passes() {
		const size_t i(header_extras.find("passes="));
		if (i == std::string::npos) {
			return;
		}
		std::istringstream(header_extras.substr(i + 7)) >> passes;
	}
	void parse_name(const std::string& line) {
		// header pattern: ^[>@][^ ]+(?: (.+))$
		// name pattern: ^[^/]+/\d+/(\d+)_(\d+)_CCS$, with ($start, $stop)
		size_t i(line.find(' ', 1));
		if (i == std::string::npos) {
			name = line;
		} else {
			name = line.substr(1, i - 1);
			header_extras = line.substr(i + 1);
			parse_passes();
		}
		i = name.find('/');
		if (i == std::string::npos) {
			std::cerr << "Warning: non-pacbio read name1: " << line << "\n";
			return;
		}
		i = name.find('/', i + 1);
		if (i == std::string::npos) {
			std::cerr << "Warning: non-pacbio read name2: " << line << "\n";
			return;
		}
		const size_t j(name.find('_', ++i));
		if (j == std::string::npos) {
			std::cerr << "Warning: non-pbtranscript read name: " << line << "\n";
			return;
		}
		std::istringstream(name.substr(i, j - i)) >> start_;
		std::istringstream(name.substr(j + 1)) >> stop_;
		name.erase(i);					// keep second / as delimiter
	}
    public:
	std::string name, header_extras;
	explicit FullRead() : start_(-1), stop_(-1) { }
	~FullRead() { }
	// with line containing the header line, read in rest of fastx and store in self;
	// returns true if eof is reached
	int read_in(const int fd, std::string& line) {
		const int is_fasta(line[0] == '>');
		if (!is_fasta && line[0] != '@') {
			std::cerr << "Error: could not parse header line: " << line << "\n";
			exit(1);
		}
		parse_name(line);
		if (pfgets(fd, line) == -1) {
			std::cerr << "Error: unexpected eof on read file: " << name << "\n";
			exit(1);
		} else if (is_fasta) {		// just read (possibly multi-line) sequence
			while (line[0] != '>') {
				seq += line;
				if (pfgets(fd, line) == -1) {
					return 1;
				}
			}
			return 0;
		}
		seq = line;			// fastq's never have multi-line sequence
		if (pfgets(fd, line) == -1) {
			std::cerr << "Error: unexpected eof on read file: " << name << "\n";
			exit(1);
		} else if (line[0] != '+') {
			std::cerr << "Error: missing qual header line: " << name << "\n";
			exit(1);
		} else if (pfgets(fd, line) == -1) {
			std::cerr << "Error: unexpected eof on read file: " << name << "\n";
			exit(1);
		}
		qual = line;
		return pfgets(fd, line) == -1 ? 1 : 0;

	}
	// add pass, quality, and/or sequence data from "a" to self
	void update(const Read& a) {
		if (passes == -1) {
			passes = a.passes;
			if (!header_extras.empty()) {
				header_extras += ";";
			}
			std::ostringstream x;
			x << passes;
			header_extras += "passes=" + x.str();
		} else if (a.passes != -1 && a.passes != passes) {
			std::cerr << "Warning: non-matching pass counts: " << name << ": " << a.passes << " != " << passes << "\n";
		}
		if (qual.empty() && !a.qual.empty()) {
			if (start_ == -1) {
				std::cerr << "Error: cannot add quality without read start and stop position: " << name << "\n";
				exit(1);
			} else if (start_ < stop_) {
				qual = a.qual.substr(start_, stop_ - start_);
			} else {
				set_and_comp_qual(a.qual.substr(stop_, start_ - stop_));
			}
		} else if (!a.qual.empty() && start_ != -1) {
			if (start_ < stop_) {
				if (a.qual.substr(start_, stop_ - start_).compare(qual) != 0) {
					std::cerr << "Warning: non-equal quals: " << name << "\n";
				}
			} else {
				// use std::equal to avoid having to copy and reverse a.qual
				if (static_cast<size_t>(start_ - stop_) != qual.size() || !std::equal(qual.rbegin(), qual.rend(), a.qual.substr(stop_, start_ - stop_).begin())) {
					std::cerr << "Warning: non-equal quals: " << name << "\n";
				}
			}
		}
		if (seq.empty() && !a.seq.empty()) {
			if (start_ == -1) {
				std::cerr << "Error: cannot add sequence without read start and stop position: " << name << "\n";
				exit(1);
			} else if (start_ < stop_) {
				seq = a.seq.substr(start_, stop_ - start_);
			} else {
				set_and_comp_seq(a.seq.substr(stop_, start_ - stop_));
			}
		} else if (!a.seq.empty() && start_ != -1) {
			if (start_ < stop_) {
				if (a.seq.substr(start_, stop_ - start_).compare(seq) != 0) {
					std::cerr << "Warning: non-equal seqs: " << name << "\n";
				}
			} else if (static_cast<size_t>(start_ - stop_) != seq.size() || comp_not_equal(seq.rbegin(), seq.rend(), a.seq.substr(stop_, start_ - stop_).begin())) {
				std::cerr << "Warning: non-equal seqs: " << name << "\n";
			}
		}
	}
	// print fastq entry
	void print(const int fd) const {
		pfwrite(fd, "@", 1);
		pfputs(fd, name);
		if (start_ != -1) {
			std::ostringstream x;
			x << start_ << "_" << stop_;
			pfputs(fd, x.str());
			pfwrite(fd, "_CCS", 4);
		}
		if (!header_extras.empty()) {
			pfwrite(fd, " ", 1);
			pfputs(fd, header_extras);
		}
		pfwrite(fd, "\n", 1);
		pfputs(fd, seq);
		pfwrite(fd, "\n+\n", 3);
		pfputs(fd, qual);
		pfwrite(fd, "\n", 1);
	}
	void clear() {
		Read::clear();
		header_extras.clear();
		start_ = stop_ = -1;
	}
    private:
	int start_, stop_;	// if start_ > stop_, bam entry needs complementing
};

static void print_usage(void) {
	std::cerr <<
		"usage: add_quality [-o output] <ccs_bam_file> <fasta/fastq>\n" <<
		"    -o ## file to store output in [stdout]\n";
	exit(0);
}

static void read_bam(const std::string ccs_bam, std::map<std::string, Read>& reads) {
	// BamReader throws on error
	PacBio::BAM::BamReader f(ccs_bam);
	PacBio::BAM::BamRecord r;
	while (f.GetNext(r)) {
		std::string s(r.FullName());
		// strip ccs from end
		if (s.length() > 3 && s.substr(s.length() - 4).compare("/ccs") == 0) {
			s = s.erase(s.length() - 3);
		} else {
			std::cerr << "Warning: bad read name: " << s << "\n";
		}
		Read& x(reads[s]);
		x.passes = r.NumPasses();
		x.seq = r.Sequence();
		x.qual = r.Qualities().Fastq();
	}
}

static void process_fastx(const std::string& read_file, const std::map<std::string, Read>& reads) {
	int fd_in(open_compressed(read_file));
	if (fd_in == -1) {
		std::cerr << "Error: open: " << read_file << "\n";
		exit(1);
	}
	int fd_out(write_fork(opt_output_file));
	if (fd_out == -1) {
		std::cerr << "Error: could not write output file: " << opt_output_file << "\n";
		close_compressed(fd_in);
		exit(1);
	}
	FullRead read;
	std::string line;
	if (pfgets(fd_in, line) == -1) {
		std::cerr << "Warning: empty file: " << read_file << "\n";
	} else {
		for (;;) {
			const int last_read(read.read_in(fd_in, line));
			std::map<std::string, Read>::const_iterator a(reads.find(read.name));
			if (a == reads.end()) {
				std::cerr << "Error: bam is missing read: " << read.name << "\n";
				exit(1);
			}
			read.update(a->second);
			read.print(fd_out);
			if (last_read) {
				break;
			}
			read.clear();
		}
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
	init_comp_lookup();
	get_opts(argc, argv);
	std::map<std::string, Read> reads;
	read_bam(argv[optind], reads);
	process_fastx(argv[optind + 1], reads);
	return 0;
}
