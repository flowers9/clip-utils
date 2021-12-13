// returns a flat array in row/column order of the maximum read size of
// each well normalized to the average + 2 * std_dev; with -a, the total
// read size for each well is use; with -p, the number of passes (or reads
// per well, depending on what information the file contains)
//
// first line is "max_row max_col"

#include "pbbam/BamReader.h"	// BamReader
#include "pbbam/BamRecord.h"	// BamRecord
#include <exception>	// exception
#include <getopt.h>	// getopt(), optind
#include <iostream>	// cerr, cout
#include <map>		// map<>
#include <math.h>	// sqrt()
#include <stdio.h>	// EOF
#include <string>	// string
#include <sys/types.h>	// size_t, uint32_t

// well_id is 16 bits row (high), 16 bits column (low)
#define WELL_TO_ROW(x) ((x) >> 16)
#define WELL_TO_COL(x) ((x) & 0xFFFF)
#define ROW_COL_TO_WELL(x,y) (((x) << 16) | y)

int opt_aggregate, opt_passes;

class LocalException : public std::exception {
    private:
	const std::string error_;
	const int show_usage_;
    public:
	explicit LocalException(const std::string &s) : error_(s), show_usage_(0) { }
	LocalException(const std::string &s, int i) : error_(s), show_usage_(i) { }
	virtual ~LocalException(void) throw() { }
	virtual const char * what() const throw() {
		return error_.c_str();
	}
	const int &show_usage(void) const throw() {
		return show_usage_;
	}
};

static void print_usage(void) {
	std::cerr <<
		"usage: extract_bam_well_sizes [-ap] <pacbio_bam>\n"
		"    -a  add all read lengths together for each well\n"
		"    -p  use passes per well instead of read length\n";
}

static void get_opts(int argc, char **argv) {
	opt_aggregate = 0;
	opt_passes = 0;
	int c;
	while ((c = getopt(argc, argv, "ap")) != EOF) {
		switch (c) {
		    case 'a':
			opt_aggregate = 1;
			break;
		    case 'p':
			opt_passes = 1;
			break;
		    default:
			throw LocalException("bad options: " + static_cast<char>(c), 1);
		}
	}
	if (optind + 1 != argc) {
		throw LocalException("need to specify bam file", 1);
	}
}

static void read_bam(const char * const file, std::map<uint32_t, uint32_t> &well_read_size) {
	// BamReader throws on error
	PacBio::BAM::BamReader bam_file(file);
	PacBio::BAM::BamRecord bam_record;
	while (bam_file.GetNext(bam_record)) {
		const uint32_t well_id(bam_record.HoleNumber());
		const std::map<uint32_t, uint32_t>::iterator a(well_read_size.find(well_id));
		if (opt_passes) {
			if (bam_record.HasNumPasses()) {
				// use += as subreads.bam sets it to one for each read
				well_read_size[well_id] += bam_record.NumPasses();
			} else {	// if not present, just count number of reads
				++well_read_size[well_id];
			}
		} else {	// using read size (max or total)
			uint32_t read_size;
			// in addition to being more correct, faster than Sequence().length()
			// but, not always present for ccs reads
			if (bam_record.HasQueryStart() && bam_record.HasQueryEnd()) {
				const uint32_t start(bam_record.QueryStart());
				const uint32_t end(bam_record.QueryEnd());
				read_size = start < end ? end - start : start - end;
			} else {
				read_size = bam_record.Sequence().length();
			}
			if (a == well_read_size.end()) {
				if (opt_aggregate && bam_record.HasNumPasses()) {
					well_read_size[well_id] = read_size * bam_record.NumPasses();
				} else {
					well_read_size[well_id] = read_size;
				}
			} else if (opt_aggregate) {		// total read sizes
				if (bam_record.HasNumPasses()) {
					// handle aggregating for ccs bams, and safe for
					// subread bams as passes will equal one
					a->second += read_size * bam_record.NumPasses();
				} else {
					a->second += read_size;
				}
			} else if (a->second < read_size) {	// update maximum read size
				a->second = read_size;
			}
		}
	}
}

static void calc_values(const std::map<uint32_t, uint32_t> &well_read_size, uint32_t &max_row, uint32_t &max_col, double &c) {
	max_row = max_col = 0;
	double x(0), x2(0);		// read size avg and std_dev
	std::map<uint32_t, uint32_t>::const_iterator a(well_read_size.begin());
	const std::map<uint32_t, uint32_t>::const_iterator end_a(well_read_size.end());
	for (; a != end_a; ++a) {
		const uint32_t row(WELL_TO_ROW(a->first));
		if (max_row < row) {
			max_row = row;
		}
		const uint32_t col(WELL_TO_COL(a->first));
		if (max_col < col) {
			max_col = col;
		}
		x += a->second;
		x2 += a->second * a->second;
	}
	x /= well_read_size.size();
	x2 /= well_read_size.size();
	const double standard_deviation(sqrt(x2 - x * x));
	c = x + 2 * standard_deviation;			// normalization constant
}

// print normalized max read sizes for all wells (0 for wells without reads)
static void print_well_sizes(const std::map<uint32_t, uint32_t> &well_read_size, const uint32_t max_row, const uint32_t max_col, const double c) {
	std::cout << max_row << " " << max_col << "\n";
	uint32_t next_row(0), next_col(0);
	std::map<uint32_t, uint32_t>::const_iterator a(well_read_size.begin());
	const std::map<uint32_t, uint32_t>::const_iterator end_a(well_read_size.end());
	for (; a != end_a; ++a) {
		const uint32_t row(WELL_TO_ROW(a->first));
		const uint32_t col(WELL_TO_COL(a->first));
		// print filler zeros
		for (; next_row != row; ++next_row) {
			for (; next_col != max_col; ++next_col) {
				std::cout << "0\n";
			}
			next_col = 0;
		}
		for (; next_col != col; ++next_col) {
			std::cout << "0\n";
		}
		// increment next_row & next_col
		if (++next_col == max_col) {
			next_col = 0;
			++next_row;
		}
		if (a->second < c) {
			std::cout << (a->second / c) << "\n";
		} else {
			std::cout << "1\n";
		}
	}
}

int main(int argc, char **argv) {
	int had_error(0);
	try {
		get_opts(argc, argv);
		// use map (rather than vector), as well numbers are likely to be sparse
		std::map<uint32_t, uint32_t> well_read_size;	// {well} = max read size
		read_bam(argv[optind], well_read_size);
		// now that we have max read sizes per well, calculate some other maxes
		// (do this after reading so you don't have multiple passes per well)
		uint32_t max_row, max_col;
		double c;					// normalization constant
		calc_values(well_read_size, max_row, max_col, c);
		// put in ending sentinel so array is fully filled
		const uint32_t max_well_id(ROW_COL_TO_WELL(max_row, max_col));
		if (well_read_size.find(max_well_id) == well_read_size.end()) {
			well_read_size[max_well_id] = 0;
		}
		// +1 to make max values usable as end conditions
		print_well_sizes(well_read_size, max_row + 1, max_col + 1, c);
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
