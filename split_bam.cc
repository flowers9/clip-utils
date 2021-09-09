#include "pbbam/BamReader.h"	// BamReader
#include "pbbam/BamRecord.h"	// BamRecord
#include "pbbam/BamWriter.h"	// BamWriter
#include <fstream>	// ifstream, ofstream
#include <iostream>	// cerr, cin, istream
#include <limits>	// numeric_limits<>
#include <sstream>	// ostringstream
#include <string>	// getline(), string
#include <sys/types.h>	// int32_t, size_t
#include <vector>	// vector<>

// this takes a subreads.bam and an uncompressed ccs.fastq file and splits
// them into N parts, also converting the fastq to a fasta; this reads
// the entire fastq (well, header and sequence) into memory to prevent
// need for rereading (you can pipe a compressed file to it, for example),
// so a touch memory intensive

#define CHUNKS 64

// appends input stream up to \n to end of string
static int get_line(std::istream &f, std::string &output) {
	char c;
	while (f.get(c)) {
		output += c;
		if (c == '\n') {
			return 0;
		}
	}
	return -1;
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
	for (;;) {
		// we'll read directly into the string to avoid extra copying
		reads.push_back(std::string());
		if (get_line(f, reads.back())) {	// get header
			if (f.eof()) {
				reads.pop_back();
				break;
			}
			std::cerr << "failed to read header: " << reads.size() << "\n";
			return -1;
		}
		reads.back()[0] = '>';			// replace leading @ with >
		if (get_line(f, reads.back())) {	// append sequence
			std::cerr << "failed to read sequence: " << reads.back();
			return -1;
		}
		if (!f.ignore(max, '\n')) {		// skip quality header
			std::cerr << "failed to ignore quality header: " << reads.back();
			return -1;
		}
		if (!f.ignore(max, '\n')) {		// skip quality
			std::cerr << "failed to ignore quality: " << reads.back();
			return -1;
		}
	}
	return 0;
}

static int32_t hole_number(const std::string &header) {
	const size_t i(header.find('/'));
	size_t j(header.find('/', i + 1));
	int32_t k(0);
	for (--j; i < j; --j) {
		k = k * 10 + header[j] - '0';
	}
	return k;
}

static int write_fastas(const std::vector<std::string> &reads, std::vector<int32_t> &holes, std::vector<int32_t> &cutoffs) {
	holes.reserve(reads.size());
	cutoffs.reserve(CHUNKS);
	const size_t base_bin_size(reads.size() / CHUNKS);
	int remainder(reads.size() % CHUNKS);
	std::ostringstream x;
	std::vector<std::string>::const_iterator a(reads.begin());
	for (int i(0); i < CHUNKS; ++i) {
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

int main(const int argc, const char ** const argv) {
	if (argc != 3 || !*argv[1] || !*argv[2]) {
		std::cerr << "usage: split_ccs_bam <subreads.bam> <ccs.fastq>\n";
		return 1;
	}
	std::vector<int32_t> holes, cutoffs;
	{	// allow reads to fall out of scope to recover memory
		std::vector<std::string> reads;
		if (read_fastq(argv[2], reads)) {
			std::cerr << "error reading fastq file\n";
			return 1;
		}
		if (reads.size() < CHUNKS) {
			std::cerr << "too few reads to split: " << reads.size() << "\n";
			return 1;
		}
		if (write_fastas(reads, holes, cutoffs)) {
			std::cerr << "error splitting fastq into fastas\n";
			return 1;
		}
	}
	PacBio::BAM::BamReader f_in(argv[1]);
	PacBio::BAM::BamRecord r;
	std::ostringstream x;
	std::vector<int32_t>::const_iterator a(holes.begin());
	for (int i(0); i < CHUNKS; ++i) {
		const int32_t &cutoff(cutoffs[i]);
		x.str("");
		x << "subreads." << i << ".bam";
		PacBio::BAM::BamWriter f_out(x.str(), f_in.Header());
		if (i != 0 && r.HoleNumber() == *a) {
			f_out.Write(r);
		}
		for (;;) {
			// skip holes not in ccs
			while (f_in.GetNext(r) && r.HoleNumber() < *a) { }
			// write out all subreads for hole
			do {
				f_out.Write(r);
			while (f_in.GetNext(r) && r.HoleNumber() == *a);
			// check against holes.end() first, as r may be invalid
			if (++a == holes.end() || cutoff < r.HoleNumber()) {
				break;
			}
		}
	}
	return 0;
}
