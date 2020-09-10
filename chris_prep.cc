// a streaming pipeline tool to convert raw data into data useful for
// various pipelines

#include "breakup_line.h"	// breakup_line()
#include "open_compressed.h"	// close_compressed(), open_compressed(), pfgets()
#include "write_fork.h"	// close_fork(), close_fork_wait(), pfputc(), pfputs(), pfwrite(), write_fork()
#include <algorithm>	// sort()
#include <ctype.h>	// isspace()
#include <errno.h>	// errno
#include <exception>	// exception
#include <fstream>	// ofstream
#include <getopt.h>	// getopt(), optarg, optind
#include <iostream>	// cerr
#include <list>		// list<>
#include <locale>	// locale, numpunct
#include <map>		// map<>
#include <set>		// set<>
#include <sstream>	// istringstream, ostringstream
#include <stdio.h>	// EOF, rename()
#include <stdlib.h>	// system()
#include <string>	// string
#include <unistd.h>	// STDOUT_FILENO, unlink()
#include <utility>	// make_pair(), move(), pair<>
#include <vector>	// vector<>

// this program reads a pair of files and writes N lines of each in succession
// Note to self: be *very* careful with the size_t's - do not subtract
//	when doing comparisons if it could wrap (and check other subtractions)!

static std::string comp_lookup(256, 0);

class LocalException : public std::exception {
    private:
	const std::string error_;
	const int show_usage_;
    public:
	explicit LocalException(const std::string &s) : error_(s), show_usage_(0) { }
	LocalException(const std::string &s, const int i) : error_(s), show_usage_(i) { }
	virtual ~LocalException() throw() { }
	virtual const char * what() const throw() {
		return error_.c_str();
	}
	const int &show_usage() const throw() {
		return show_usage_;
	}
};

class Options {
    public:
	std::string contaminant_fasta, linker_file, het_rate;
	std::string project_path, library;
	size_t minimum_read_length, max_reads;
	// boolean options, other than mer_size
	int mer_size, paired_reads;
	int diversity, no_simple_filter, output_fasta, print_to_stdout;
	Options() : minimum_read_length(-1), max_reads(-1), mer_size(-1), paired_reads(1), diversity(0), no_simple_filter(0), output_fasta(0), print_to_stdout(0) { }
	~Options() { }
};

class Library {
    public:
	enum class Types { NEXP, JGIFRAG, RNASEQ, CLRS, LFPE, NEXT, NGEN, smRNA, UNSET };
	std::vector<std::string> input_files;
	size_t minimum_read_length;
	Types library_type;
	int is_paired;
	Library() : minimum_read_length(-1), library_type(Types::UNSET), is_paired(0) { }
	~Library() { }
	int is_smrna() const {
		return library_type == Types::smRNA;
	}
	int is_rnaseq() const {
		return library_type == Types::RNASEQ;
	}
	void set_type(const std::string &s) {
		if (s == "NEXP") {
			library_type = Types::NEXP;
		} else if (s == "JGIFRAG") {
			library_type = Types::JGIFRAG;
		} else if (s == "RNASEQ") {
			library_type = Types::RNASEQ;
		} else if (s == "CLRS") {
			library_type = Types::CLRS;
		} else if (s == "LFPE") {
			library_type = Types::LFPE;
		} else if (s == "NEXT") {
			library_type = Types::NEXT;
		} else if (s == "NGEN") {
			library_type = Types::NGEN;
		} else if (s == "smRNA") {
			library_type = Types::smRNA;
		} else {
			throw LocalException("unknown library type");
		}
	}
};

struct OrderContaminants {
    private:
	const std::map<std::string, size_t> &contaminant_count_;
    public:
	explicit OrderContaminants(const std::map<std::string, size_t> &contaminant_count) : contaminant_count_(contaminant_count) { }
	bool operator()(const std::string &s, const std::string &t) const {
		std::map<std::string, size_t>::const_iterator a(contaminant_count_.find(s));
		std::map<std::string, size_t>::const_iterator b(contaminant_count_.find(t));
		if (b->second != a->second) {
			return b->second < a->second;
		} else {
			return s < t;
		}
	}
};

// facet to add commas every 3 number positions
class comma_numpunct : public std::numpunct<char> {
    protected:
	virtual std::string do_grouping() const {
		return "\03";
	}
};

class Counts {
    public:
	std::map<std::string, size_t> contaminant_count;	// [read_name] = count
	size_t reads_extracted, seq_extracted, reads_prepped, seq_prepped, reads_flipped;
	size_t reads_lost_to_ns, reads_lost_to_lq, reads_lost_to_vector, reads_lost_to_polya;
	size_t reads_lost_to_simple, seq_lost_to_simple, reads_lost_to_contaminant;
	size_t seq_lost_to_contaminant, reads_singleton;
	Counts() : reads_extracted(0), seq_extracted(0), reads_prepped(0), seq_prepped(0), reads_flipped(0), reads_lost_to_ns(0), reads_lost_to_lq(0), reads_lost_to_vector(0), reads_lost_to_polya(0), reads_lost_to_simple(0), seq_lost_to_simple(0), reads_lost_to_contaminant(0), seq_lost_to_contaminant(0), reads_singleton(0) { }
	~Counts() { }
	void print_summary(const std::string &library_name, const Library &library) const {
		std::string file(library_name + ".extractionStats");
		std::ofstream out(file.c_str());
		if (!out.is_open()) {
			std::cerr << "Warning: could not write stats file\n";
			return;
		}
		std::string read_file;
		const size_t i(library.input_files[0].rfind('/'));
		read_file = library.input_files[0].substr(i == std::string::npos ? 0 : i + 1);
		if (library.input_files.size() == 2) {
			read_file += ", ";
			const size_t j(library.input_files[1].rfind('/'));
			read_file += library.input_files[1].substr(j == std::string::npos ? 0 : j + 1);
		}
		// create a new local with comma separated numbers
		const std::locale comma_locale(std::locale(), new comma_numpunct());
		// and have out use it
		out.imbue(comma_locale);
		out << "\n\n"
			"=============================\n"
			"CUMULATIVE EXTRACTION SUMMARY\n"
			"=============================\n"
			"Reads File:        " << read_file << "\n"
			"Reads Extracted:   " << reads_extracted << "\n"
			"Seq   Extracted:   " << seq_extracted << " bp\n"
			"---\n"
			"N Lost Reads:      " << reads_lost_to_ns << "\n"
			"---\n"
			"Low Quality Reads: " << reads_lost_to_lq << "\n"
			"---\n"
			"Vector trimmed:    " << reads_lost_to_vector << "\n";
		if (library.is_rnaseq()) {
			out << "Reads Flipped:     " << reads_flipped << "\n";
		}
		out << "---\n"
			"Simple Reads:      " << reads_lost_to_simple << "\n"
			"Simple Seq:        " << seq_lost_to_simple << " bp\n"
			"---\n"
			"Singleton Reads:   " << reads_singleton << "\n"
			"Prepped Reads:     " << reads_prepped << "\n"
			"Prepped Seq:       " << seq_prepped << " bp\n";
		if (!contaminant_count.empty()) {
			out << "---\n"
				"Contaminant Reads: " << reads_lost_to_contaminant << "\n"
				"Contaminant Seq:   " << seq_lost_to_contaminant << " bp\n"
				"Contaminant Summary:\n";
			// print contaminant reads, most hits to least
			std::vector<std::string> list;
			list.reserve(contaminant_count.size());
			std::map<std::string, size_t>::const_iterator a(contaminant_count.begin());
			const std::map<std::string, size_t>::const_iterator end_a(contaminant_count.end());
			for (; a != end_a; ++a) {
				list.push_back(a->first);
			}
			std::sort(list.begin(), list.end(), OrderContaminants(contaminant_count));
			std::vector<std::string>::const_iterator b(list.begin());
			const std::vector<std::string>::const_iterator end_b(list.end());
			for (; b != end_b; ++b) {
				out << "\t\t" << *b << ":\t" << contaminant_count.find(*b)->second << '\n';
			}
		}
		out << "=============================\n\n\n";
		out.close();
	}
};

class Read {
    public:
	std::string name, seq, qual;
	size_t hq_start, hq_end, lfpe_start, lfpe_end;
	Read() : hq_start(-1), hq_end(-1), lfpe_start(-1), lfpe_end(-1) { }
	explicit Read(std::string&& name_in, std::string&& seq_in, std::string&& qual_in) : name(name_in), seq(seq_in), qual(qual_in), hq_start(-1), hq_end(-1), lfpe_start(-1), lfpe_end(-1) { }
	~Read() { }
	void set_limits() {
		hq_start = lfpe_start = 0;
		hq_end = lfpe_end = seq.size();
	}
	size_t hq_length() const {
		// watch that wrapping!
		return hq_start < hq_end ? hq_end - hq_start : 0;
	}
	void create_from_pair(const Read &a) {
		name = a.name;
		// swap -R1 <-> -R2
		name[name.size() - 1] = name[name.size() - 1] == '1' ? '2' : '1';
		// reverse compliment from pair
		seq.clear();
		seq.reserve(a.seq.size());
		std::string::const_reverse_iterator b(a.seq.rbegin());
		const std::string::const_reverse_iterator end_b(a.seq.rend());
		for (; b != end_b; ++b) {
			seq += comp_lookup[static_cast<int>(*b)];
		}
		qual.assign(a.qual.rbegin(), a.qual.rend());	// just reverse the quality
		hq_start = seq.size() - a.hq_end;		// add hq clipping
		hq_end = seq.size() - a.hq_start;
		lfpe_start = seq.size() - a.lfpe_end;		// add lfpe clipping
		lfpe_end = seq.size() - a.lfpe_start;
	}
};

static std::ostringstream output_buffer;	// one buffer, to reduce memory allocations
class Outputs {
    public:
	void (*write_output)(const Outputs &outputs, const std::vector<Read> &reads, Counts &counts);
	std::vector<std::string> output_files;
	std::string contaminant_file, simple_file, singleton_file;
	int fd1, fd2, fd_contaminant, fd_simple, fd_singleton;
	Outputs() : write_output(0), fd1(-1), fd2(-1), fd_contaminant(-1), fd_simple(-1), fd_singleton(-1) { }
	~Outputs() { }
	void close_all() {
		if (fd1 != -1) {
			close_fork_wait(fd1);
			fd1 = -1;
		}
		if (fd2 != -1) {
			close_fork_wait(fd2);
			fd2 = -1;
		}
		if (fd_contaminant != -1) {
			close_fork_wait(fd_contaminant);
			fd_contaminant = -1;
		}
		if (fd_simple != -1) {
			close_fork_wait(fd_simple);
			fd_simple = -1;
		}
		if (fd_singleton != -1) {
			close_fork_wait(fd_singleton);
			fd_singleton = -1;
		}
	}
	void rename_all(const Counts &counts) {
		for (size_t i(0); i != output_files.size(); ++i) {
			const std::string s(output_files[i] + ".tmp");
			rename(s.c_str(), output_files[i].c_str());
		}
		output_files.clear();
		// remove files if they're empty, otherwise rename from tmp
		if (!contaminant_file.empty()) {
			const std::string s(contaminant_file + ".tmp");
			if (counts.reads_lost_to_contaminant == 0) {
				unlink(s.c_str());
			} else {
				rename(s.c_str(), contaminant_file.c_str());
			}
		}
		if (!simple_file.empty()) {
			const std::string s(simple_file + ".tmp");
			if (counts.reads_lost_to_simple == 0) {
				unlink(s.c_str());
			} else {
				rename(s.c_str(), simple_file.c_str());
			}
		}
		if (!singleton_file.empty()) {
			const std::string s(singleton_file + ".tmp");
			if (counts.reads_singleton == 0) {
				unlink(s.c_str());
			} else {
				rename(s.c_str(), singleton_file.c_str());
			}
		}
	}
};

// lookup table is to speed up comping reads; here's where we set it up
static void init_comp_lookup() {
	// unknown stuff just stays the same
	for (unsigned int i(0); i != comp_lookup.size(); ++i) {
		comp_lookup[i] = i;
	}
	comp_lookup[static_cast<int>('A')] = 'T';
	comp_lookup[static_cast<int>('a')] = 't';
	comp_lookup[static_cast<int>('C')] = 'G';
	comp_lookup[static_cast<int>('c')] = 'g';
	comp_lookup[static_cast<int>('G')] = 'C';
	comp_lookup[static_cast<int>('g')] = 'c';
	comp_lookup[static_cast<int>('T')] = 'A';
	comp_lookup[static_cast<int>('t')] = 'a';
}

static void print_usage() {
	std::cerr <<
		"usage: chris_prep [options] <project_path> <library_base_name>\n" <<
		"    -C     print output to stdout\n" <<
		"    -c ##  contaminant fasta file [none]\n" <<
		"    -d     Diversity run\n" <<
		"    -f     output fasta & qual files (instead of fastq)\n" <<
		"    -h     print this help\n" <<
		"    -m ##  set mer size [8/10/14, depends on library]\n" <<
		"    -n ##  number of reads to extract [all]\n" <<
		"    -p ##  minimum read length after clip & trim [50/75 for R<250/R>=250]\n" <<
		"    -s     don't filter simple sequence\n" <<
		"    -u     allow unpaired reads\n" <<
		"    -v ##  fasta file with linker\n";
}

static int get_opts(int argc, char **argv, Options &opts) {
	int c;
	while ((c = getopt(argc, argv, "Cc:dfhm:n:p:suv:")) != EOF) {
		switch (c) {
		    case 'C':
			opts.print_to_stdout = 1;
			break;
		    case 'c':
			opts.contaminant_fasta = optarg;
			break;
		    case 'd':
			opts.diversity = 1;
			break;
		    case 'f':
			opts.output_fasta = 1;
			break;
		    case 'h':
			print_usage();
			return 1;
		    case 'm':
			std::istringstream(optarg) >> opts.mer_size;
			break;
		    case 'n':
			std::istringstream(optarg) >> opts.max_reads;
			break;
		    case 'p':
			std::istringstream(optarg) >> opts.minimum_read_length;
			break;
		    case 'r':
			opts.het_rate = optarg;
			break;
		    case 's':
			opts.no_simple_filter = 1;
			break;
		    case 'u':
			opts.paired_reads = 0;
			break;
		    case 'v':
			opts.linker_file = optarg;
			break;
		    default:
			throw LocalException("bad option: " + static_cast<char>(c), 1);
		}
	}
	if (optind + 2 != argc) {
		throw LocalException("incorrect number of arguments", 1);
	}
	opts.project_path = argv[optind];
	opts.library = argv[optind + 1];
	return 0;
}

static void read_config_file(const Options &opts, Library &library) {
	std::string config_file(opts.project_path + "/unProcessed/lib.config");
	int fd(open_compressed(config_file));
	if (fd == -1) {
		config_file = opts.project_path + "/CONFIG/lib.config";
		fd = open_compressed(config_file);
		if (fd == -1) {
			throw LocalException("could not read configuration file");
		}
	}
	std::string line;
	while (pfgets(fd, line) != -1) {
		const size_t i(line.find(opts.library));
		if (i != static_cast<size_t>(-1) && isspace(line[i + opts.library.size()])) {
			close_compressed(fd);
			std::vector<std::string> list;
			breakup_line(line, list);
			if (list[4] != "UNPROCESSED") {
				throw LocalException("library status is not UNPROCESSED");
			}
			library.set_type(list[1]);
			if (library.library_type == Library::Types::NEXP || library.library_type == Library::Types::CLRS || library.library_type == Library::Types::LFPE) {
				library.is_paired = 1;
			}
			breakup_line(list[5], library.input_files, ',');
			std::vector<std::string> list2;
			breakup_line(list[3], list2, 'x');
			if (list2.size() != 2) {
				throw LocalException("read length format incorrect in config file");
			}
			size_t library_read_length;
			std::istringstream(list2[1]) >> library_read_length;
			library.minimum_read_length = library.is_paired || library_read_length < 250 ? 50 : 75;
			return;
		}
	}
	close_compressed(fd);
	throw LocalException("could not find library in configuration file");
}

static void update_config_file(const Options &opts) {
	std::string config_file(opts.project_path + "/unProcessed/lib.config");
	int fd(open_compressed(config_file));
	if (fd == -1) {
		config_file = opts.project_path + "/CONFIG/lib.config";
		fd = open_compressed(config_file);
		if (fd == -1) {
			std::cerr << "Warning: could not update configuration file: could not read file\n";
		}
	}
	std::string line;
	std::vector<std::string> lines;
	while (pfgets(fd, line) != -1) {
		lines.push_back(line);
	}
	close_compressed(fd);
	std::vector<std::string>::iterator a(lines.begin());
	const std::vector<std::string>::const_iterator end_a(lines.end());
	for (; a != end_a; ++a) {
		const size_t i(a->find(opts.library));
		if (i != std::string::npos && isspace((*a)[i + 1])) {
			const size_t j(a->find("UNPROCESSED"));
			if (j == std::string::npos) {
				std::cerr << "Warning: library is no longer unprocessed\n";
				// looks like we're not changing it after all
				return;
			}
			a->erase(j, 2);
			break;
		}
	}
	fd = write_fork(std::list<std::string>(), config_file.c_str());
	if (fd == -1) {
		std::cerr << "Warning: could not rewrite configuration file\n";
	}
	for (a = lines.begin(); a != end_a; ++a) {
		if (pfputs(fd, *a) == -1 || pfputc(fd, '\n') == -1) {
			std::cerr << "Error writing configuration file: " << *a;
		}
	}
	close_fork(fd);
}

static void apply_library_defaults(Options &opts, const Library &library) {
	if (opts.minimum_read_length == static_cast<size_t>(-1)) {
		opts.minimum_read_length = library.minimum_read_length;
	}
	if (opts.mer_size == -1) {
		switch (library.library_type) {
		    case Library::Types::NEXP:
		    case Library::Types::JGIFRAG:
		    case Library::Types::CLRS:
		    case Library::Types::LFPE:
		    case Library::Types::NEXT:
		    case Library::Types::NGEN:
			opts.mer_size = 10;
			break;
		    case Library::Types::RNASEQ:
			opts.mer_size = 14;
			break;
		    case Library::Types::smRNA:
			opts.mer_size = 11;
			break;
		    case Library::Types::UNSET:		// can't happen
			break;
		}
	}
	if (opts.linker_file.empty()) {
		switch (library.library_type) {
		    case Library::Types::NEXP:
			opts.linker_file = "/home/raid2/SEQ/sharedPythonLibrary/prep_scripts_cbp/linkerSeq/nexteraAdapter.fasta";
			break;
		    case Library::Types::JGIFRAG:
			opts.linker_file = "/home/raid2/SEQ/sharedPythonLibrary/prep_scripts_cbp/linkerSeq/illuminaLinker.fasta";
			break;
		    case Library::Types::RNASEQ:
			opts.linker_file = "/home/raid2/SEQ/sharedPythonLibrary/prep_scripts_cbp/linkerSeq/rnaSeqLinker.fasta";
			break;
		    case Library::Types::CLRS:
			opts.linker_file = "/home/raid2/SEQ/sharedPythonLibrary/prep_scripts_cbp/linkerSeq/CRELOX_linker.fasta";
			break;
		    case Library::Types::LFPE:
			opts.linker_file = "/home/raid2/SEQ/sharedPythonLibrary/prep_scripts_cbp/linkerSeq/LFPE_linker.fasta";
			break;
		    case Library::Types::NEXT:
			opts.linker_file = "/home/raid2/SEQ/sharedPythonLibrary/prep_scripts_cbp/linkerSeq/nexteraAdapter.fasta";
			break;
		    case Library::Types::NGEN:
			opts.linker_file = "/home/raid2/SEQ/sharedPythonLibrary/prep_scripts_cbp/linkerSeq/nugenAdapter.fasta";
			break;
		    case Library::Types::smRNA:
			opts.linker_file = "/global/dna/projectdirs/plant/geneAtlas/HAGSC_TOOLS/PREP_TESTING/adapters.fa";
			break;
		    case Library::Types::UNSET:		// can't happen
			break;
		}
	}
	if (!opts.no_simple_filter && library.is_rnaseq()) {
		opts.no_simple_filter = 1;
	}
}

// reverse compliment in place
static void reverse_compliment(std::string &s) {
	size_t i(0), j(s.size() - 1);
	for (; i < j; ++i, --j) {
		const char c(comp_lookup[static_cast<int>(s[i])]);
		s[i] = comp_lookup[static_cast<int>(s[j])];
		s[j] = c;
	}
	if (i == j) {
		s[i] = comp_lookup[static_cast<int>(s[j])];
	}
}

static void get_linker_kmers(const std::string &linker_file, const int mer_size, const int is_rnaseq, std::set<std::string> &linker_mers, std::set<std::string> &linker_7mers) {
	const int fd(open_compressed(linker_file));
	if (fd == -1) {
		throw LocalException("could not open linker file");
	}
	std::string line;
	while (pfgets(fd, line) != -1) {
		if (line.empty() || line[0] != '>') {
			throw LocalException("incorrect header line in linker file");
		}
		if (pfgets(fd, line) == -1) {
			throw LocalException("truncated linker file");
		}
		for (size_t i(0); i + mer_size <= line.size(); ++i) {
			std::string mer(line.substr(i, mer_size));
			linker_mers.insert(mer);
			reverse_compliment(mer);
			linker_mers.insert(mer);
		}
		if (is_rnaseq) {
			for (size_t i(0); i + 7 <= line.size(); ++i) {
				std::string mer(line.substr(i, 7));
				linker_7mers.insert(mer);
				reverse_compliment(mer);
				linker_7mers.insert(mer);
			}
		}
	}
	close_compressed(fd);
}

// use first 1k reads to figure out quality ranges
static int get_qual_offset(const std::string &file) {
	const int fd(open_compressed(file));
	if (fd == -1) {
		throw LocalException("could not open input file(s)");
	}
	std::string line;
	if (pfgets(fd, line) == -1) {
		throw LocalException("input file is empty");
	}
	if (line.empty() || line[0] != '@') {
		throw LocalException("bad read header format in input file");
	}
	pfgets(fd, line);	// skip seq
	pfgets(fd, line);	// skip qual header
	if (pfgets(fd, line) == -1) {
		throw LocalException("input file truncated in first read");
	}
	int min_qual(line[0]), max_qual(line[0]);
	int count(0);
	for (;;) {
		std::string::const_iterator a(line.begin());
		const std::string::const_iterator end_a(line.end());
		for (; a != end_a; ++a) {
			if (min_qual > *a) {
				min_qual = *a;
			} else if (max_qual < *a) {
				max_qual = *a;
			}
		}
		pfgets(fd, line);	// skip header
		pfgets(fd, line);	// skip seq
		pfgets(fd, line);	// skip qual header
		if (++count == 1000 || pfgets(fd, line) == -1) {
			break;
		}
	}
	close_compressed(fd);
	if (63 < min_qual && max_qual < 77) {
		throw LocalException("could not distinguish quality score encoding");
	} else if (32 < min_qual && max_qual < 77) {
		return 33;
	} else if (63 < min_qual && max_qual < 115) {
		return 64;
	} else {
		throw LocalException("unknown quality score encoding");
	}
}

// dunno if it's faster to write output to a buffer and then print,
// or to print each part separately

static void write_fastq_no_count(const int fd, const std::vector<Read> &reads) {
	std::vector<Read>::const_iterator a(reads.begin());
	const std::vector<Read>::const_iterator end_a(reads.end());
	for (; a != end_a; ++a) {
		output_buffer.str("");
		output_buffer << '@' << a->name << '\n' << a->seq << "\n+\n" << a->qual << '\n';
		pfputs(fd, output_buffer.str());
	}
}

// this is "unclipped", which ignores any hq clipping,
// but still uses lfpe clipping

static void write_fastq(const Outputs &outputs, const std::vector<Read> &reads, Counts &counts) {
	std::vector<Read>::const_iterator a(reads.begin());
	const std::vector<Read>::const_iterator end_a(reads.end());
	for (; a != end_a; ++a) {
		const size_t n(a->lfpe_end - a->lfpe_start);
		output_buffer.str("");
		output_buffer << '@' << a->name << '\n' << a->seq.substr(a->lfpe_start, n) << "\n+\n" << a->qual.substr(a->lfpe_start, n) << '\n';
		pfputs(outputs.fd1, output_buffer.str());
		counts.seq_prepped += n;
	}
	counts.reads_prepped += reads.size();
}

static void write_fastq_clipped(const Outputs &outputs, const std::vector<Read> &reads, Counts &counts) {
	std::vector<Read>::const_iterator a(reads.begin());
	const std::vector<Read>::const_iterator end_a(reads.end());
	for (; a != end_a; ++a) {
		const size_t n(a->hq_end - a->hq_start);
		output_buffer.str("");
		output_buffer << '@' << a->name << '\n' << a->seq.substr(a->hq_start, n) << "\n+\n" << a->qual.substr(a->hq_start, n) << '\n';
		pfputs(outputs.fd1, output_buffer.str());
		counts.seq_prepped += n;
	}
	counts.reads_prepped += reads.size();
}

static void write_fastq_split_clipped(const Outputs &outputs, const std::vector<Read> &reads, Counts &counts) {
	std::vector<Read>::const_iterator a(reads.begin());
	const std::vector<Read>::const_iterator end_a(reads.end());
	for (; a != end_a; ++a) {
		const size_t n(a->hq_end - a->hq_start);
		const char c(a->name[a->name.size() - 1]);
		output_buffer.str("");				// clear buffer
		// convert -R# notation to /#
		output_buffer << '@' << a->name.substr(0, a->name.size() - 3) << '/' << c << '\n' << a->seq.substr(a->hq_start, n) << "\n+\n" << a->qual.substr(a->hq_start, n) << '\n';
		pfputs(c == '1' ? outputs.fd1 : outputs.fd2, output_buffer.str());
		counts.seq_prepped += n;
	}
	counts.reads_prepped += reads.size();
}

static void write_fasta(const Outputs &outputs, const std::vector<Read> &reads, Counts &counts) {
	std::vector<Read>::const_iterator a(reads.begin());
	const std::vector<Read>::const_iterator end_a(reads.end());
	for (; a != end_a; ++a) {
		const size_t n(a->lfpe_end - a->lfpe_start);
		output_buffer.str("");
		output_buffer << '>' << a->name << '\n';
		pfputs(outputs.fd2, output_buffer.str());
		output_buffer << a->seq.substr(a->lfpe_start, n) << '\n';
		pfputs(outputs.fd1, output_buffer.str());
		output_buffer.str("");
		output_buffer << static_cast<int>(a->qual[a->lfpe_start] - 33);
		for (size_t i(a->lfpe_start + 1); i < a->lfpe_end; ++i) {
			output_buffer << ' ' << static_cast<int>(a->qual[i] - 33);
		}
		output_buffer << '\n';
		pfputs(outputs.fd2, output_buffer.str());
		counts.seq_prepped += n;
	}
	counts.reads_prepped += reads.size();
}

static void write_fasta_clipped(const Outputs &outputs, const std::vector<Read> &reads, Counts &counts) {
	std::vector<Read>::const_iterator a(reads.begin());
	const std::vector<Read>::const_iterator end_a(reads.end());
	for (; a != end_a; ++a) {
		const size_t n(a->hq_end - a->hq_start);
		output_buffer.str("");
		output_buffer << '>' << a->name << '\n';
		pfputs(outputs.fd2, output_buffer.str());
		output_buffer << a->seq.substr(a->hq_start, n) << '\n';
		pfputs(outputs.fd1, output_buffer.str());
		output_buffer.str("");
		output_buffer << static_cast<int>(a->qual[a->hq_start] - 33);
		for (size_t i(a->hq_start); i != a->hq_end - 1; ++i) {
			output_buffer << ' ' << static_cast<int>(a->qual[i] - 33);
		}
		output_buffer << '\n';
		pfputs(outputs.fd2, output_buffer.str());
		counts.seq_prepped += n;
	}
	counts.reads_prepped += reads.size();
}

static void prepare_for_writing(const Options &opts, const Library &library, Outputs &outputs) {
	if (opts.print_to_stdout) {
		outputs.fd1 = STDOUT_FILENO;
		if (library.is_rnaseq() || library.is_paired) {
			outputs.write_output = write_fastq_clipped;
		} else {
			outputs.write_output = write_fastq;
		}
	} else if (library.is_rnaseq()) {
		outputs.output_files.push_back(opts.library + ".prepped.R1.fastq.gz");
		outputs.output_files.push_back(opts.library + ".prepped.R2.fastq.gz");
		std::list<std::string> args(1, "gzip");
		outputs.fd1 = write_fork(args, outputs.output_files[0] + ".tmp");
		outputs.fd2 = write_fork(args, outputs.output_files[1] + ".tmp");
		outputs.write_output = write_fastq_split_clipped;
	} else if (library.is_paired) {
		if (opts.output_fasta) {
			outputs.output_files.push_back(opts.library + ".prepped.uncomp.fasta.bz2");
			outputs.output_files.push_back(opts.library + ".prepped.uncomp.qual.bz2");
			std::list<std::string> args(1, "bzip2");
			outputs.fd1 = write_fork(args, outputs.output_files[0] + ".tmp");
			outputs.fd2 = write_fork(args, outputs.output_files[1] + ".tmp");
			outputs.write_output = write_fasta_clipped;
		} else {
			outputs.output_files.push_back(opts.library + ".prepped.uncomp.fastq.bz2");
			std::list<std::string> args(1, "bzip2");
			outputs.fd1 = write_fork(args, outputs.output_files[0] + ".tmp");
			outputs.write_output = write_fastq_clipped;
		}
	} else {
		if (opts.output_fasta) {
			outputs.output_files.push_back(opts.library + ".prepped.fasta.bz2");
			outputs.output_files.push_back(opts.library + ".prepped.qual.bz2");
			std::list<std::string> args(1, "bzip2");
			outputs.fd1 = write_fork(args, outputs.output_files[0] + ".tmp");
			outputs.fd2 = write_fork(args, outputs.output_files[1] + ".tmp");
			outputs.write_output = write_fasta;
		} else {
			outputs.output_files.push_back(opts.library + ".prepped.fastq.bz2");
			std::list<std::string> args(1, "bzip2");
			outputs.fd1 = write_fork(args, outputs.output_files[0] + ".tmp");
			outputs.write_output = write_fastq;
		}
	}
	// prepare auxiliary output files
	if (!opts.contaminant_fasta.empty()) {
		outputs.contaminant_file = opts.library + ".contam.fastq.bz2";
		std::list<std::string> args(1, "bzip2");
		outputs.fd_contaminant = write_fork(args, outputs.contaminant_file + ".tmp");
	}
	if (!opts.no_simple_filter) {
		outputs.simple_file = opts.library + ".simpleReads.fastq.bz2";
		std::list<std::string> args(1, "bzip2");
		outputs.fd_simple = write_fork(args, outputs.simple_file + ".tmp");
	}
	outputs.singleton_file = opts.library + ".singletons.bz2";
	std::list<std::string> args(1, "bzip2");
	outputs.fd_singleton = write_fork(args, outputs.singleton_file + ".tmp");
}

// convert read name into local format
static void convert_header_to_name(std::string &name) {
	if (name.empty()) {
		throw LocalException("blank header line in fastq file");
	}
	size_t i(1);
	name.erase(0, 1);	// remove [@>] prefix
	const size_t n(name.size());
	// translate [:#-] to _
	for (; i != n && !isspace(name[i]); ++i) {
		if (name[i] == ':' || name[i] == '#' || name[i] == '-') {
			name[i] = '_';
		}
	}
	// slash separator
	if (name[n - 2] == '/' && (name[n - 1] == '1' || name[n - 1] == '2')) {
		name.replace(n - 2, 1, "-R");
		return;
	} else if (i != n) {
		const size_t j(i);
		// skip any more whitespace
		for (++i; i != n && isspace(name[i]); ++i) { }
		if ((name[i] == '1' || name[i] == '2') && name[i + 1] == ':') {
			name.erase(i + 1);
			name.replace(j, i - j, "-R");
			return;
		} else if (name[n - 2] == ':' && (name[n - 1] == '1' || name[n - 1] == '2')) {
			name.replace(j, n - 1 - j, "-R");
			return;
		}
	}
	throw LocalException("could not parse read name");
}

static int get_next_read(const int fd, const int qual_offset, std::vector<Read> &reads) {
	std::string name, seq, qual;
	if (pfgets(fd, name) == -1) {		// EOF
		return 0;
	}
	convert_header_to_name(name);
	if (pfgets(fd, seq) == -1) {
		throw LocalException("truncated fastq file");
	}
	pfgets(fd, qual);			// skip qual header
	if (pfgets(fd, qual) == -1) {
		throw LocalException("truncated fastq file");
	}
	if (qual_offset != 0) {			// normalize quals to 33 offset
		for (size_t i(0); i != qual.size(); ++i) {
			qual[i] -= qual_offset;
		}
	}
	reads.push_back(Read(std::move(name), std::move(seq), std::move(qual)));
	return 1;
}

static std::list<Read> next_reads;		// next pairs or partial pairs
static std::map<std::string, Read> waiting_reads; // unpaired reads, [base_read_name] = read

// when one interleaved file finishes, finish up the rest;
// dump any pairs found into next_reads

static void finish_off_pairs(const int fd, const int qual_offset, std::vector<Read> &reads) {
	for (;;) {
		const std::map<std::string, Read>::iterator a(waiting_reads.find(reads[0].name.substr(0, reads[0].name.size() - 3)));
		if (a != waiting_reads.end()) {
			next_reads.push_back(std::move(reads[0]));
			next_reads.push_back(std::move(a->second));
			waiting_reads.erase(a);
			if (waiting_reads.empty()) {	// no more pairs to be had
				return;
			}
		}
		reads.clear();				// either unpaired or moved
		if (!get_next_read(fd, qual_offset, reads)) {
			break;
		}
	}
	waiting_reads.clear();
}

static int get_next_reads(const std::vector<int> &input_fds, const int qual_offset, const int paired_reads, std::vector<Read> &reads) {
	reads.clear();
	if (input_fds.size() == 1) {		// don't have to interleave manually
		if (!get_next_read(input_fds[0], qual_offset, reads)) {
			return 0;
		} else if (paired_reads && !get_next_read(input_fds[0], qual_offset, reads)) {
			return 0;
		} else {
			return 1;
		}
	}
	while (next_reads.empty()) {		// find next pair
		const int got0(get_next_read(input_fds[0], qual_offset, reads));
		const int got1(get_next_read(input_fds[1], qual_offset, reads));
		if (got0 && got1) {
			// check if they're a pair
			if (reads[0].name.compare(0, reads[0].name.size() - 3, reads[1].name, 0, reads[1].name.size() - 3) == 0) {
				return 1;
			}
			const std::map<std::string, Read>::iterator a0(waiting_reads.find(reads[0].name.substr(0, reads[0].name.size() - 3)));
			const std::map<std::string, Read>::iterator a1(waiting_reads.find(reads[1].name.substr(0, reads[1].name.size() - 3)));
			if (a0 != waiting_reads.end()) {
				if (a1 != waiting_reads.end()) {	// save for later
					next_reads.push_back(std::move(reads[1]));
					next_reads.push_back(std::move(a1->second));
					waiting_reads.erase(a1);
				}
				reads[1] = std::move(a0->second);
				waiting_reads.erase(a0);
				return 1;
			} else if (a1 != waiting_reads.end()) {
				reads[0] = std::move(a1->second);
				waiting_reads.erase(a1);
				return 1;
			}
		} else if ((got0 || got1) && !waiting_reads.empty()) {
			// one file ended earlier, so look for unpaired matches
			finish_off_pairs(input_fds[got0 ? 0 : 1], qual_offset, reads);
			break;
		} else {			// end of files
			return 0;
		}
	}
	if (next_reads.empty()) {
		return 0;
	}
	reads.push_back(std::move(next_reads.front()));
	next_reads.pop_front();
	if (paired_reads) {
		reads.push_back(std::move(next_reads.front()));
		next_reads.pop_front();
	}
	return 1;
}

// return if read was not lost to N's; trims N's in place;
// modifies quality to match changes to sequence
static int trim_ns(Read &read) {
	size_t i;
	// trim trailing N's
	i = read.seq.find_last_not_of('N');
	if (i == std::string::npos) {		// all N's
		return 0;
	}
	if (++i != read.seq.size()) {
		read.seq.erase(i);
		read.qual.erase(i);
	}
	// trim leading N's
	i = read.seq.find_first_not_of('N');
	if (i != 0) {
		read.seq.erase(0, i);
		read.qual.erase(0, i);
	}
	// check for a leading basepair following by two or more N's
	if (read.seq.size() > 2 && read.seq[1] == 'N' && read.seq[2] == 'N') {
		return 0;
	}
	// remove singular N (and non-N) at front
	if (read.seq.size() > 1 && read.seq[1] == 'N') {
		read.seq.erase(0, 2);
		read.qual.erase(0, 2);
		return 1;
	}
	i = read.seq.find('N');
	if (i == std::string::npos) {		// no more N's
		return 1;
	} else if (read.seq[i + 1] == 'N') {
		// if there's a multi-N run, just take the starting sequence
		read.seq.erase(i);
		read.qual.erase(i);
		return 1;
	}
	// if the pattern is ...N...N, remove everything from the second N on
	const size_t j(read.seq.find('N', i + 1));
	if (j != std::string::npos) {
		read.seq.erase(j);
		read.qual.erase(j);
	}
	return 1;
}

// add hq_start/hq_end to read, returns if it found a high quality region
static int find_high_quality(Read &read) {
	const size_t window_size(20);
	if (read.qual.size() < window_size) {	// not enough quality
		read.hq_start = read.hq_end;
		return 0;
	}
	const int quality_cutoff(58);		// 25 + 33 offset
	const int window_cutoff(window_size * quality_cutoff);
	// find start of high quality region
	size_t i(0);
	int window_total(0);
	for (; i != window_size; ++i) {
		window_total += read.qual[i];
	}
	for (; i != read.qual.size() && window_total < window_cutoff; ++i) {
		window_total += read.qual[i] - read.qual[i - window_size];
	}
	if (window_total < window_cutoff) {	// no high quality sequence
		read.hq_start = read.hq_end;
		return 0;
	}
	read.hq_start = i - window_size;
	// find end of last high quality region
	window_total = 0;
	for (i = read.qual.size() - 1; i != read.qual.size() - window_size - 1; --i) {
		window_total += read.qual[i];
	}
	for (; window_total < window_cutoff; --i) {
		window_total += read.qual[i] - read.qual[i + window_size];
	}
	read.hq_end = i + window_size + 1;	// exclusive end
	return 1;
}

// add hq_start and hq_end to each read, returns true if enough sequence left
static int hq_clip(std::vector<Read> &reads, const size_t minimum_read_length, const int is_paired, const int no_clipping) {
	if (!no_clipping && !find_high_quality(reads[0])) {
		return 0;			// no high quality sequence
	}
	const size_t n(reads[0].hq_length());
	if (reads.size() == 2) {
		if (!no_clipping && !find_high_quality(reads[1])) {
			return 0;		// no high quality sequence on second read
		}
		const size_t n2(reads[1].hq_length());
		if (is_paired) {
			// pairs require both reads to be long enough
			if (n < minimum_read_length || n2 < minimum_read_length) {
				return 0;
			}
		} else if (n < minimum_read_length && n2 < minimum_read_length) {
			return 0;
		}
	} else if (n < minimum_read_length) {
		return 0;
	}
	return 1;
}

// checks to see if the range is long enough, or the ends are close to hq start/end
static inline bool check_for_early_end_condition(const size_t hq_end, const size_t min_region_length, const size_t hq_region_spacing, const size_t linker_range_start, const size_t linker_range_end) {
	return linker_range_end >= linker_range_start + min_region_length || hq_end <= linker_range_start + hq_region_spacing;
}

static inline bool check_for_late_end_condition(const size_t hq_start, const size_t hq_region_spacing, const size_t linker_range_end) {
	return linker_range_end <= hq_start + hq_region_spacing;
}

// vector clipping - returns if enough read is left after clipping
static int lfpe_clip(Read &read, const std::set<std::string> &linker_kmers, const int unclip_odd_case, const size_t mer_size, const size_t minimum_read_length) {
	const int failed_clipping(read.hq_end < read.hq_start + minimum_read_length);
	const size_t hq_region_spacing(36);
	const size_t min_region_length(hq_region_spacing / 2);
	const size_t collapse_spacing(5);
	// find ranges that are composed of linker kmers
	size_t linker_range_start(-1), linker_range_end(-1);
	size_t i(0);
	size_t end_i(read.seq.size() - mer_size + 1);
	// look for first linker kmer
	for (; i != end_i && linker_kmers.find(read.seq.substr(i, mer_size)) == linker_kmers.end(); ++i) { }
	if (i == end_i) {			// no vector found
		if (failed_clipping) {
			read.hq_end = 0;	// mark as bad read for rnaseq
			return 0;
		}
		return 1;
	}
	for (;;) {
		const size_t start(i);
		// find end of linker range
		for (++i; i != end_i && linker_kmers.find(read.seq.substr(i, mer_size)) != linker_kmers.end(); ++i) { }
		if (linker_range_start == static_cast<size_t>(-1)) {	// first range
			linker_range_start = start;
			linker_range_end = i + mer_size - 1;
		} else if (start <= linker_range_end + collapse_spacing) {
			linker_range_end = i + mer_size - 1;		// extend range
		} else if (start < linker_range_end + min_region_length) {
			// check for N's that might extend the previous region
			for (;;) {
				size_t j(linker_range_end);
				for (; j != start && read.seq[j] != 'N'; ++j) { }
				if (j == start) {	// no more N's
					// see if N's got us close enough to merge
					if (start <= linker_range_end + collapse_spacing) {
						linker_range_end = i + mer_size - 1;
					} else if (check_for_late_end_condition(read.hq_start, hq_region_spacing, linker_range_end)) {
						goto BREAK_OUTER_LOOP;
					} else {	// if not, move on
						linker_range_start = start;
						linker_range_end = i + mer_size - 1;
					}
					break;
				} else if (j <= linker_range_end + collapse_spacing) {
					// extending
					// check for run of N's
					// (don't have to worry about running off end because
					// linkers kmers don't have N's)
					for (++j; read.seq[j] == 'N'; ++j) { }
					linker_range_end = j;
				} else if (check_for_late_end_condition(read.hq_start, hq_region_spacing, linker_range_end)) {
					goto BREAK_OUTER_LOOP;
				} else {	// can't extend, so start with the N's
					linker_range_start = j;
					for (++j; read.seq[j] == 'N'; ++j) { }
					linker_range_end = j;
				}
				if (check_for_early_end_condition(read.hq_end, min_region_length, hq_region_spacing, linker_range_start, linker_range_end)) {
					goto BREAK_OUTER_LOOP;
				}
			}
		} else if (check_for_late_end_condition(read.hq_start, hq_region_spacing, linker_range_end)) {
			break;
		} else {						// start new range
			linker_range_start = start;
			linker_range_end = i + mer_size - 1;
		}
		if (check_for_early_end_condition(read.hq_end, min_region_length, hq_region_spacing, linker_range_start, linker_range_end)) {
			break;
		}
		// advance to next linker kmer
		if (i != end_i) {
			for (++i; i != end_i && linker_kmers.find(read.seq.substr(i, mer_size)) == linker_kmers.end(); ++i) { }
		}
		if (i == end_i) {			// no (effective) vector found
			if (check_for_late_end_condition(read.hq_start, hq_region_spacing, linker_range_end)) {
				break;
			} else if (failed_clipping) {
				read.hq_end = 0;	// mark as bad read for rnaseq
				return 0;
			} else if (unclip_odd_case) {
				// not sure about this, but remove clipping in this case
				read.hq_start = 0;
				read.hq_end = read.seq.size();
			}
			return 1;
		}
	}
	BREAK_OUTER_LOOP:
	// see if vector found too close to beginning
	if (linker_range_start < minimum_read_length) {
		read.hq_end = 0;	// mark as bad read for rnaseq
		return 0;
	}
	// change high quality region to full read prior to vector
	read.hq_start = 0;
	read.hq_end = linker_range_start;
	if (unclip_odd_case) {
		// this clipping extends to printing the ranges, whereas the hq clipping
		// (for this case) is only used for analysis
		read.lfpe_start = 0;
		read.lfpe_end = linker_range_start;
	}
	return 1;
}

// returns if read wasn't lost to poly-a
// (we also trim poly-t, in case it's a complimented poly-a)

static int trim_polya(Read &read, const std::set<std::string> &linker_7mers, Counts &counts) {
	const size_t collapse(5);
	const size_t padding(25);	// area to check for poly-a/t
	const size_t minimum_nonpoly_basepairs(50);
	std::vector<std::pair<size_t, size_t> > poly_ranges;
	// surprisingly, this doesn't get checked by the clip routines
	if (read.hq_end < read.hq_start + minimum_nonpoly_basepairs) {
		++counts.reads_lost_to_polya;
		read.hq_end = 0;	// mark as bad
		return 0;
	}
	// TODO: after verification against the perl, rewrite this to just
	// remove poly-a/t with a foot in the padding area
	// -9 to leave space for a 10 run
	for (size_t start(read.hq_start); start != read.hq_end - 9; ++start) {
		if (read.seq[start] != 'A' && read.seq[start] != 'T') {
			continue;
		}
		// is it a 10+ run?
		size_t end(start + 1);
		const size_t long_enough(start + 10);
		for (; end != long_enough && read.seq[end] == read.seq[start]; ++end) { }
		if (end != long_enough) {
			continue;
		}
		// find end of run
		for (; end != read.hq_end && read.seq[end] == read.seq[start]; ++end) { }
		if (poly_ranges.empty()) {
			poly_ranges.push_back(std::make_pair(start, end));
			continue;
		}
		const size_t spacing(start - poly_ranges.back().second);
		if (start <= poly_ranges.back().second + collapse) {	// collapse nearby ranges
			poly_ranges.back().second = end;
		} else if (spacing < 10) {	// add A's or T's in the middle
			const std::string s(read.seq.substr(poly_ranges.back().second, spacing));
			size_t last_end(0);	// the amount we'll extend the previous run
			size_t at_start, at_end(0);
			while ((at_start = s.find_first_of("AT", at_end)) != std::string::npos) {
				for (at_end = at_start + 1; at_end != s.size() && s[at_end] == s[at_start]; ++at_end) { }
				if (spacing - at_end <= collapse) {
					// AT's run into current range
					start -= spacing - at_start;
					break;
				} else if (at_start <= last_end + collapse) {	// extend last
					last_end = at_end;
				}
			}
			poly_ranges.back().second += last_end;
			// did the AT's fill in the gap?
			if (start <= poly_ranges.back().second + collapse) {	// yep
				poly_ranges.back().second = end;
			} else {						// nope
				poly_ranges.push_back(std::make_pair(start, end));
			}
		} else {
			poly_ranges.push_back(std::make_pair(start, end));
		}
	}
	std::vector<std::pair<size_t, size_t> >::const_iterator a(poly_ranges.begin());
	const std::vector<std::pair<size_t, size_t> >::const_iterator end_a(poly_ranges.end());
	const size_t start_cutoff(read.hq_start + padding);
	for (; a != end_a; ++a) {
		// any range that would satisfy both conditions will
		// get marked as bad down below
		if (a->first <= start_cutoff) {
			read.hq_start = a->second;
		} else if (read.hq_end <= a->second + padding) {
			read.hq_end = a->first;
			break;
		}
	}
	// check for linker 7-mers on ends, and trim
	if (linker_7mers.find(read.seq.substr(read.hq_start, 7)) != linker_7mers.end()) {
		read.hq_start += 7;
	}
	if (linker_7mers.find(read.seq.substr(read.hq_end - 7, 7)) != linker_7mers.end()) {
		read.hq_end -= 7;
	}
	if (read.hq_end < read.hq_start + minimum_nonpoly_basepairs) {
		++counts.reads_lost_to_polya;
		read.hq_end = 0;	// mark as bad
		return 0;
	}
	return 1;
}

static std::vector<int> base_lookup(256, -1);

static void init_base_lookup() {
	base_lookup['A'] = 0;
	base_lookup['C'] = 1;
	base_lookup['G'] = 2;
	base_lookup['T'] = 3;
}

// loop until we find a good triplet or hit the end of high quality sequence
static int init_triplet(const std::string &seq, size_t &i, const size_t end_i, int &j) {
	for (;;) {
		// look for next good basepair
		for (; i < end_i - 2 && base_lookup[seq[i]] == -1; ++i) { }
		if (i >= end_i - 2) {		// ran out of high quality sequence
			return 0;
		}
		// see if next two basepairs are also good; if not, skip past them
		if (base_lookup[seq[i + 2]] == -1) {
			i += 3;
		} else if (base_lookup[seq[i + 1]] == -1) {
			i += 2;
		} else {
			j = (((base_lookup[seq[i]] << 2) | base_lookup[seq[i + 1]]) << 2) | base_lookup[seq[i + 2]];
			i += 3;
			return 1;
		}
	}
}

static int find_next_triplet(const std::string &seq, size_t &i, const size_t end_i, int &j) {
	if (base_lookup[seq[i]] != -1) {
		j = ((j << 2) & 63) | base_lookup[seq[i]];
		++i;
		return 1;
	}
	return init_triplet(seq, ++i, end_i, j);
}

// returns is simple sequence score is too high
static bool find_simple_sequence(const Read &read) {
	// count all ACGT triplets
	std::vector<size_t> triplet_counts(64, 0);
	size_t i(read.hq_start);
	int j(0);
	if (!init_triplet(read.seq, i, read.hq_end, j)) {	// no good sequence
		return 1;					// means no simple sequence :|
	}
	++triplet_counts[j];
	while (i != read.hq_end && find_next_triplet(read.seq, i, read.hq_end, j)) {
		++triplet_counts[j];
	}
	// compute score
	size_t count(0), sigma(0);
	std::vector<size_t>::const_iterator a(triplet_counts.begin());
	const std::vector<size_t>::const_iterator end_a(triplet_counts.end());
	for (; a != end_a; ++a) {
		if (*a > 1) {
			++count;
			sigma += *a * (*a - 1);
		}
	}
	if (count > 2) {
		// cluster the 3mer counts to find errant 3mers arising from sequencing error
		std::sort(triplet_counts.begin(), triplet_counts.end());
		size_t start(0);
		// first, skip singletons
		// (we already know there are at least two values that satisfy this)
		for (; triplet_counts[start] < 2; ++start) { }
		// second, do a simple separation into two clusters
		// (extend the lower cluster until values are closer to the highest value)
		size_t lower_total(triplet_counts[start]);
		// while abs(distance to lower mean) < abs(distance to higher mean)
		for (i = start + 1; triplet_counts[i] - lower_total / (i - start) < triplet_counts.back() - triplet_counts[i]; ++i) {
			lower_total += triplet_counts[i];
		}
		size_t higher_total(triplet_counts[i]);
		for (size_t k(i + 1); k != triplet_counts.size(); ++k) {
			higher_total += triplet_counts[k];
		}
		// if the difference between means is greater than 10, recompute
		if (double(higher_total) / (triplet_counts.size() - i) - double(lower_total) / (i - start) > 10) {
			count = triplet_counts.size() - i;
			sigma = 0;
			for (; i != triplet_counts.size(); ++i) {
				sigma += triplet_counts[i] * (triplet_counts[i] - 1);
			}
		} else if (count > 4) {		// regular sequence, but reset count
			count = 4;
		}
	}
	const double n(read.hq_length() - 2);
	return double(count) * sigma / n / (n - count) >= .9604;
}

// we only check adjacent reads to see if they're a pair, but this means
// we need to keep track of the last read to the last call of this

static void count_singletons(const int fd, const std::vector<Read> &reads, Counts &counts) {
	static std::string last_read;
	if (reads.empty()) {		// clean up after all reads processed
		if (!last_read.empty()) {
			++counts.reads_singleton;
			pfputs(fd, last_read);
			pfputc(fd, '\n');
		}
		return;
	}
	size_t i(0);
	if (!last_read.empty()) {			// deal with leftover read
		if (reads[0].name.compare(0, reads[0].name.size() - 3, last_read) == 0) {
			++i;
		} else {
			++counts.reads_singleton;
			pfputs(fd, last_read);
			pfputc(fd, '\n');
		}
		last_read.clear();
	}
	for (; i < reads.size() - 1; ++i) {
		if (reads[i].name.compare(0, reads[i].name.size() - 3, reads[i + 1].name, 0, reads[i + 1].name.size() - 3) == 0) {
			++i;
		} else {
			++counts.reads_singleton;
			pfwrite(fd, reads[i].name.c_str(), reads[i].name.size() - 3);
			pfputc(fd, '\n');
		}
	}
	if (i != reads.size()) {	// we had one left
		last_read = reads.back().name.substr(0, reads.back().name.size() - 3);
	}
}

// modifies reads in place, removing contaminted reads (and their pairs, if need be);
// in addition to screening, also prints out non-contaminant reads

static void screen_contaminants(std::vector<Read> &reads, const std::string &contaminant_fasta, const int paired_reads, const int fd) {
	(void)reads;
	(void)contaminant_fasta;
	(void)paired_reads;
	(void)fd;
	// XXX
}

static void process_reads(const Options &opts, const Library &library, Outputs &outputs, const int qual_offset, const std::vector<int> &input_fds, Counts &counts, const std::set<std::string> &linker_kmers, const std::set<std::string> &linker_7mers) {
	(void)outputs;
	std::vector<Read> batch_reads;	// for batched contaminant processing
	batch_reads.reserve(50000);
	std::vector<Read> reads;
	while (counts.reads_extracted < opts.max_reads && get_next_reads(input_fds, qual_offset, opts.paired_reads, reads)) {
		counts.reads_extracted += reads.size();
		counts.seq_extracted += reads[0].seq.size();
		if (reads.size() == 2) {
			counts.seq_extracted += reads[1].seq.size();
		}
		if (!trim_ns(reads[0]) || (reads.size() == 2 && !trim_ns(reads[1]))) {
			counts.reads_lost_to_ns += reads.size();
			continue;
		}
		reads[0].set_limits();
		if (reads.size() == 2) {
			reads[1].set_limits();
		}
		// soft clip to high quality regions; diversity simply disables hq clipping
		if (!hq_clip(reads, opts.minimum_read_length, library.is_paired, opts.diversity)) {
			counts.reads_lost_to_lq += reads.size();
			continue;
		}
		const int unclip_odd_case(!library.is_paired && !library.is_rnaseq());
		// check for linker kmers and soft clip
		size_t worked(lfpe_clip(reads[0], linker_kmers, unclip_odd_case, opts.mer_size, opts.minimum_read_length));
		if (reads.size() == 2 && (worked || library.is_rnaseq())) {
			worked += lfpe_clip(reads[1], linker_kmers, unclip_odd_case, opts.mer_size, opts.minimum_read_length);
		}
		if (library.is_rnaseq() && worked != 0) {	// poly-a trimming
			if (reads[0].hq_end != 0 && !trim_polya(reads[0], linker_7mers, counts)) {
				--worked;
			}
			if (reads[1].hq_end != 0 && !trim_polya(reads[1], linker_7mers, counts)) {
				--worked;
			}
		}
		if (worked == reads.size()) {	// all's good
		} else if (worked != 0 && library.is_rnaseq()) {
			// revcomp working one to make other
			const int working(reads[0].hq_end != 0 ? 0 : 1);
			reads[1 - working].create_from_pair(reads[working]);
			++counts.reads_flipped;
		} else {			// lost to vector
			counts.reads_lost_to_vector += reads.size();
			continue;
		}
		// check for too much simple sequence
		if (!opts.no_simple_filter && find_simple_sequence(reads[0]) && (reads.size() == 1 || find_simple_sequence(reads[1]))) {
			write_fastq_no_count(outputs.fd_simple, reads);
			counts.reads_lost_to_simple += reads.size();
			for (size_t i(0); i != reads.size(); ++i) {
				counts.seq_lost_to_simple += reads[i].seq.size();
			}
			continue;
		}
		// check for contaminate
		// as it's inefficient to do a blat run per read, we process
		// them in batches; as a result, they don't get written out one
		// at a time, but in a batch after contaminant screening
		if (opts.contaminant_fasta.empty()) {
			count_singletons(outputs.fd_singleton, reads, counts);
			outputs.write_output(outputs, reads, counts);
		} else {
			batch_reads.insert(batch_reads.end(), std::move(reads.begin()), std::move(reads.end()));
			if (batch_reads.size() >= 50000) {
				screen_contaminants(batch_reads, opts.contaminant_fasta, opts.paired_reads, outputs.fd_contaminant);
				outputs.write_output(outputs, batch_reads, counts);
				batch_reads.clear();
			}
		}
	}
	reads.clear();
	if (!batch_reads.empty()) {
		screen_contaminants(batch_reads, opts.contaminant_fasta, opts.paired_reads, outputs.fd_contaminant);
		outputs.write_output(outputs, batch_reads, counts);
		batch_reads.clear();
	}
	count_singletons(outputs.fd_singleton, reads, counts);
	counts.print_summary(opts.library, library);
}

static void find_het_rate(const Options &opts, const Outputs &outputs) {
	const std::string output_file(opts.library + "_Queryreads.blat.gz");
	std::ostringstream command("/mnt/local/EXBIN/clip -B 1000 -f 25 -L 0 -c ");
	command <<
		(opts.minimum_read_length == 75 ? 251 : 151) << " " <<
		outputs.output_files[0] <<
		" | /home/raid2/LINUXOPT/AnacondaEnv/PREP_ENV/bin/blat " << opts.het_rate <<
		" stdin -out=blast8 -noHead stdout | awk '{if ($1 != $2 && $3 > 95 && $12 > " <<
		(opts.minimum_read_length == 75 ? 450 : 250) << ") {print}} | gzip > " <<
		output_file << ".tmp && mv " << output_file << ".tmp " << output_file;
	// don't actually care what the return value is, but want to prevent warning
	const int whatever(system(command.str().c_str()));
	(void)whatever;
}

int main(int argc, char **argv) {
	int had_error(0);
	init_comp_lookup();
	init_base_lookup();
	try {
		class Options opts;
		if (get_opts(argc, argv, opts)) {
			return 0;
		}
		class Library library;
		read_config_file(opts, library);
		if (library.is_smrna()) {
			throw LocalException("smRNA is not implemented");
		}
		apply_library_defaults(opts, library);
		std::set<std::string> linker_kmers, linker_7mers;
		get_linker_kmers(opts.linker_file, opts.mer_size, library.is_rnaseq(), linker_kmers, linker_7mers);
		const int qual_offset(33 - get_qual_offset(library.input_files[0]));
		std::vector<int> input_fds;
		for (size_t i(0); i != library.input_files.size(); ++i) {
			const int fd(open_compressed(library.input_files[i]));
			if (fd == -1) {
				throw LocalException("could not open input file");
			}
			input_fds.push_back(fd);
		}
		Outputs outputs;
		prepare_for_writing(opts, library, outputs);
		Counts counts;
		process_reads(opts, library, outputs, qual_offset, input_fds, counts, linker_kmers, linker_7mers);
		for (size_t i(0); i != input_fds.size(); ++i) {
			close_compressed(input_fds[i]);
		}
		outputs.close_all();
		outputs.rename_all(counts);
		if (!opts.print_to_stdout) {
			update_config_file(opts);
		}
		if (!opts.het_rate.empty()) {
			find_het_rate(opts, outputs);
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
