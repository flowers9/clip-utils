#include "pbbam/BamReader.h"	// BamReader
#include "pbbam/BamRecord.h"	// BamRecord
#include "pbbam/BamWriter.h"	// BamWriter
#include <iostream>	// cerr
#include <string>	// string
#include <sys/types.h>	// size_t

// only print alignments that match subread to same ccs read, by hole number
int main(const int argc, const char ** const argv) {
	if (argc != 3 || !*argv[1] || !*argv[2]) {
		std::cerr << "usage: filter_bam <aligned.subreads.bam> <subreads_to_ccs.bam>\n";
		return 1;
	}
	PacBio::BAM::BamReader f_in(argv[1]);
	// note that BamWriter uses 4 threads for compression by default
	PacBio::BAM::BamWriter f_out(argv[2], f_in.Header());
	PacBio::BAM::BamRecord r;
	while (f_in.GetNext(r)) {
		// the read names contain the hole number: ^[^/]*/hole_number/
		// (note we can get the numeric hole number for the read
		// (r.HoleNumber()), but not for the reference, so probably fastest
		// to compare as strings rather than convert and compare ints)
		const std::string &s1(r.FullName());
		const std::string &s2(r.ReferenceName());
		const size_t i1(s1.find('/') + 1);
		const size_t j1(s1.find('/', i1));
		const size_t i2(s2.find('/') + 1);
		const size_t j2(s2.find('/', i2));
		if (!s1.compare(i1, j1 - i1, s2, i2, j2 - i2)) {
			f_out.Write(r);
		}
	}
	return 0;
}
