#include "open_compressed.h"	// close_compressed(), open_compressed(), pfgets(), skip_next_line()
#include "pbbam/BamReader.h"	// BamReader
#include "pbbam/BamRecord.h"	// BamRecord
#include "pbbam/BamWriter.h"	// BamWriter
#include <iostream>	// cerr
#include <set>		// set<>
#include <sstream>	// istringstream
#include <stdlib.h>	// exit()
#include <string>	// string
#include <sys/types.h>	// size_t, int32_t

// filter out all subreads from zmws present in fasta/q

// assume read name is of form .*/\d+/.*
static int32_t get_zmw(const std::string &s) {
	size_t i(s.find('/'));
	if (i == std::string::npos) {
		std::cerr << "Error: could not parse read name: " << s << "\n";
		exit(2);
	}
	const size_t j(s.find('/', ++i));
	if (j == std::string::npos) {
		std::cerr << "Error: could not parse read name: " << s << "\n";
		exit(3);
	}
	std::istringstream x(s.substr(i, j - i));
	int32_t z;
	x >> z;
	return z;
}

static void get_zmw_list(const char * const filename, std::set<int32_t> &zmws) {
	const int fd(open_compressed(filename));
	if (fd == -1) {
		std::cerr << "Error: could not open " << filename << ": " << strerror(errno) << "\n";
		exit(1);
	}
	std::string line;
	while (pfgets(fd, line) != -1) {
		if (line.empty()) {
		} else if (line[0] == '>') {			// fasta
			zmws.insert(get_zmw(line));
			skip_next_line(fd);	// sequence
		} else if (line[0] == '@') {			// fastq
			zmws.insert(get_zmw(line));
			skip_next_line(fd);	// sequence
			skip_next_line(fd);	// quality header
			skip_next_line(fd);	// quality
		} else {					// plain read names?
			zmws.insert(get_zmw(line));
		}
	}
	close_compressed(fd);
}

int main(const int argc, const char * const * const argv) {
	if (argc != 4 || !*argv[1] || !*argv[2] || !*argv[3]) {
		std::cerr << "usage: filter_bam <ccs.fasta/q|list of read names> <subreads.bam> <output.bam>\n";
		return 1;
	}
	std::set<int32_t> zmws;
	get_zmw_list(argv[1], zmws);
	PacBio::BAM::BamReader f_in(argv[2]);
	// note that BamWriter uses 4 threads for compression by default
	PacBio::BAM::BamWriter f_out(argv[3], f_in.Header());
	PacBio::BAM::BamRecord r;
	while (f_in.GetNext(r)) {
		if (zmws.find(r.HoleNumber()) == zmws.end()) {
			f_out.Write(r);
		}
	}
	return 0;
}
