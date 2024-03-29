#include "pbbam/BamReader.h"	// BamReader
#include "pbbam/BamRecord.h"	// BamRecord
#include "pbbam/BamWriter.h"	// BamWriter
#include <fstream>	// ifstream, ofstream
#include <getopt.h>	// getopt(), optarg, optind
#include <iostream>	// cerr, cin, istream
#include <limits>	// numeric_limits<>
#include <sstream>	// istringstream, ostringstream
#include <string>	// getline(), string
#include <sys/types.h>	// int32_t, size_t
#include <vector>	// vector<>

// this takes a subreads.bam and an uncompressed ccs.fastq file and splits
// them into N parts, also converting the fastq to fasta; this reads the
// entire fastq (well, header and sequence) into memory to prevent the
// need for rereading (you can pipe a compressed file to it, for example),
// so a touch memory intensive (~32gb on a 25gb fastq.gz file)

// this also screens the bam file as it splits it, filtering out any
// subreads from holes that aren't present in the fastq file

static unsigned int opt_chunks(1024);

static void print_usage() {
	std::cerr <<
		"usage: split_bam [-n ##] <ccs.fastq> [subreads.bam]\n"
		"    -n ##  number of chunks to split into [1024]\n";
}

static int get_opts(const int argc, char ** const argv) {
	unsigned int x;
	int c;
	while ((c = getopt(argc, argv, "n:")) != EOF) {
		switch (c) {
		    case 'n':
			std::istringstream(optarg) >> x;
			if (x) {			// if zero, use default
				opt_chunks = x;
			}
			break;
		    default:
			std::cerr << "bad option: " << char(c) << "\n";
			print_usage();
			return -1;
		}
	}
	if ((optind + 1 != argc && optind + 2 != argc) || !*argv[optind] || (optind + 2 == argc && !*argv[optind + 1])) {
		print_usage();
		return -1;
	}
	return 0;
}

static int read_fastq(const std::string &ccs_fastq, std::vector<std::string> &reads) {
	const std::streamsize max(std::numeric_limits<std::streamsize>::max());
	// bit of a hack to choose between stdin or an actual file - we open
	// whatever, but ignore it if "-" was given
	std::ifstream f_in(ccs_fastq);	// ignore this if we're reading from stdin
	std::istream &f(ccs_fastq == "-" ? std::cin : f_in);
	if (!f) {
		std::cerr << "could not open file: " << ccs_fastq << "\n";
		return -1;
	}
	// to prevent excess memory allocation, use buffer to read in data and then
	// copy to vector (cost is an extra copy, but we do use a lot of memory)
	std::string header, seq;
	for (;;) {
		if (!std::getline(f, header)) {		// get header
			if (f.eof()) {			// done reading
				return 0;
			}
			std::cerr << "failed to read header: " << reads.size() << "\n";
			return -1;
		}
		header[0] = '>';			// replace leading @ with >
		if (!std::getline(f, seq)) {		// get sequence
			std::cerr << "failed to read sequence: " << header << "\n";
			return -1;
		}
		reads.push_back(std::string());
		std::string &s(reads.back());
		s.reserve(header.size() + seq.size() + 2);
		s = header;
		s += '\n';
		s += seq;
		s += '\n';
		if (!f.ignore(max, '\n')) {		// skip quality header
			std::cerr << "failed to ignore quality header: " << header << "\n";
			return -1;
		}
		if (!f.ignore(max, '\n')) {		// skip quality
			std::cerr << "failed to ignore quality: " << header << "\n";
			return -1;
		}
	}
}

// like read_fastq, but reads in as fastq, not as fasta
static int read_fastq_entire(const std::string &ccs_fastq, std::vector<std::string> &reads) {
	// bit of a hack to choose between stdin or an actual file - we open
	// whatever, but ignore it if "-" was given
	std::ifstream f_in(ccs_fastq);	// ignore this if we're reading from stdin
	std::istream &f(ccs_fastq == "-" ? std::cin : f_in);
	if (!f) {
		std::cerr << "could not open file: " << ccs_fastq << "\n";
		return -1;
	}
	// to prevent excess memory allocation, use buffer to read in data and then
	// copy to vector (cost is an extra copy, but we do use a lot of memory)
	std::string header, seq, qual_header, qual;
	for (;;) {
		if (!std::getline(f, header)) {		// get header
			if (f.eof()) {			// done reading
				return 0;
			}
			std::cerr << "failed to read header: " << reads.size() << "\n";
			return -1;
		}
		if (!std::getline(f, seq)) {		// get sequence
			std::cerr << "failed to read sequence: " << header << "\n";
			return -1;
		}
		if (!std::getline(f, qual_header)) {	// get quality header
			std::cerr << "failed to quality header: " << header << "\n";
			return -1;
		}
		if (!std::getline(f, qual)) {		// get quality
			std::cerr << "failed to quality: " << header << "\n";
			return -1;
		}
		reads.push_back(std::string());
		std::string &s(reads.back());
		s.reserve(header.size() + seq.size() + qual_header.size() + qual.size() + 4);
		s = header;
		s += '\n';
		s += seq;
		s += '\n';
		s += qual_header;
		s += '\n';
		s += qual;
		s += '\n';
	}
}

static int32_t hole_number(const std::string &header) {
	size_t i(header.find('/') + 1);
	const size_t j(header.find('/', i));
	// probably faster than using istringstream since it's less general
	int32_t k(0);
	for (; i < j; ++i) {
		k = k * 10 + header[i] - '0';
	}
	return k;
}

static int write_fastas(const std::vector<std::string> &reads, std::vector<int32_t> &holes, std::vector<int32_t> &cutoffs) {
	const size_t base_bin_size(reads.size() / opt_chunks);
	unsigned int remainder(reads.size() % opt_chunks);
	std::ostringstream x;
	std::vector<std::string>::const_iterator a(reads.begin());
	for (unsigned int i(0); i < opt_chunks; ++i) {
		x.str("");
		x << "ccs." << i << ".fasta";
		std::ofstream f(x.str());
		if (!f) {
			std::cerr << "failed to open: " << x.str() << "\n";
			return -1;
		}
		size_t bin_size(base_bin_size);
		if (remainder) {
			--remainder;
			++bin_size;
		}
		do {
			if (!f.write(a->data(), a->size())) {
				std::cerr << "failed to write: " << x.str() << "\n";
				return -1;
			}
			holes.push_back(hole_number(*a));
			++a;
		} while (--bin_size);
		cutoffs.push_back(holes.back());
	}
	return 0;
}

static int write_fastqs(const std::vector<std::string> &reads) {
	const size_t base_bin_size(reads.size() / opt_chunks);
	unsigned int remainder(reads.size() % opt_chunks);
	std::ostringstream x;
	std::vector<std::string>::const_iterator a(reads.begin());
	for (unsigned int i(0); i < opt_chunks; ++i) {
		x.str("");
		x << "ccs." << i << ".fastq";
		std::ofstream f(x.str());
		if (!f) {
			std::cerr << "failed to open: " << x.str() << "\n";
			return -1;
		}
		size_t bin_size(base_bin_size);
		if (remainder) {
			--remainder;
			++bin_size;
		}
		do {
			if (!f.write(a->data(), a->size())) {
				std::cerr << "failed to write: " << x.str() << "\n";
				return -1;
			}
			++a;
		} while (--bin_size);
	}
	return 0;
}

static int split_ccs(const char * const ccs_file, std::vector<int32_t> &holes, std::vector<int32_t> &cutoffs) {
	std::vector<std::string> reads;
	if (read_fastq(ccs_file, reads)) {
		std::cerr << "error reading fastq file\n";
		return -1;
	}
	if (reads.size() < opt_chunks) {
		opt_chunks = reads.size();
	}
	holes.reserve(reads.size());
	cutoffs.reserve(opt_chunks);
	if (write_fastas(reads, holes, cutoffs)) {
		std::cerr << "error splitting fastq into fastas\n";
		return -1;
	}
	return 0;
}

static int split_ccs_fastq(const char * const ccs_file) {
	std::vector<std::string> reads;
	if (read_fastq_entire(ccs_file, reads)) {
		std::cerr << "error reading fastq file\n";
		return -1;
	}
	if (reads.size() < opt_chunks) {
		opt_chunks = reads.size();
	}
	if (write_fastqs(reads)) {
		std::cerr << "error splitting fastq\n";
		return -1;
	}
	return 0;
}

int main(const int argc, char ** const argv) {
	if (get_opts(argc, argv)) {
		return 1;
	}
	if (optind + 1 == argc) {	// just split fastq
		return split_ccs_fastq(argv[optind]);
	}
	std::vector<int32_t> holes, cutoffs;
	if (split_ccs(argv[optind], holes, cutoffs)) {
		return 1;
	}
	PacBio::BAM::BamReader f_in(argv[optind + 1]);
	PacBio::BAM::BamRecord r;
	std::ostringstream x;
	std::vector<int32_t>::const_iterator a(holes.begin());
	const std::vector<int32_t>::const_iterator end_a(holes.end());
	if (!f_in.GetNext(r)) {
		std::cerr << "error: empty bam file\n";
		return 1;
	}
	for (unsigned int i(0); i < opt_chunks; ++i) {
		const int32_t cutoff(cutoffs[i]);
		x.str("");
		x << "subreads." << i << ".bam";
		PacBio::BAM::BamWriter f_out(x.str(), f_in.Header());
		// print out subreads hole by hole
		do {
			// skip holes not in ccs
			while (r.HoleNumber() < *a && f_in.GetNext(r)) { }
			// write out all subreads for hole
			do {
				f_out.Write(r);
			} while (f_in.GetNext(r) && r.HoleNumber() == *a);
			// check a first, as r may be invalid
		} while (++a != end_a && r.HoleNumber() <= cutoff);
	}
	return 0;
}
