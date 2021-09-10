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
// them into N parts, also converting the fastq to fasta; this reads the
// entire fastq (well, header and sequence) into memory to prevent the
// need for rereading (you can pipe a compressed file to it, for example),
// so a touch memory intensive (~32gb on a 25gb fastq.gz file)

// this also screens the bam file as it splits it, filtering out any
// subreads from holes that aren't present in the fastq file

static unsigned int opt_chunks(64);

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
			std::cerr << "failed to read sequence: " << reads.back();
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
			std::cerr << "failed to ignore quality header: " << reads.back();
			return -1;
		}
		if (!f.ignore(max, '\n')) {		// skip quality
			std::cerr << "failed to ignore quality: " << reads.back();
			return -1;
		}
	}
}

static int32_t hole_number(const std::string &header) {
	size_t i(header.find('/') + 1);
	const size_t j(header.find('/', i));
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

int main(const int argc, const char ** const argv) {
	if (argc != 3 || !*argv[1] || !*argv[2]) {
		std::cerr << "usage: split_ccs_bam <subreads.bam> <ccs.fastq>\n";
		return 1;
	}
	std::vector<int32_t> holes, cutoffs;
	if (split_ccs(argv[2], holes, cutoffs)) {
		return 1;
	}
	PacBio::BAM::BamReader f_in(argv[1]);
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
