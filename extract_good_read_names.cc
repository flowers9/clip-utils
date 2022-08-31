#include "pbbam/BamReader.h"	// BamReader
#include "pbbam/BamRecord.h"	// BamRecord
#include <iostream>	// cerr, cout

// get list of all reads with passes >= 3

int main(const int argc, const char * const * const argv) {
	if (argc < 2) {
		std::cerr << "usage: filter_bam <ccs1.bam> [ccs2.bam ...]\n";
		return 1;
	}
	PacBio::BAM::BamRecord r;
	for (int i(1); i < argc; ++i) {
		PacBio::BAM::BamReader f(argv[i]);
		while (f.GetNext(r)) {
			if (r.NumPasses() >= 3) {
				std::cout << r.FullName() << "\n";
			}
		}
	}
	return 0;
}
