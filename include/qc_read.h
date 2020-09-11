#ifndef _QC_READ_H
#define _QC_READ_H

#include "read.h"	// Read
#include <map>		// map<>
#include <string>	// string
#include <sys/types.h>	// size_t

class QC_Read : public Read {
    private:
	void print_quality_range(size_t, size_t) const;
    public:
	unsigned int contigs;	// n_runs + 1, except that empty contigs
				// are ignored (such as those made by
				// leading or trailing gap N's)
	unsigned int n1_runs;			// gaps
	unsigned int n2_runs;			// non-gaps
	size_t n1_count;	// gap N's
	size_t n2_count;	// non-gap N's
	size_t lq_count;
	size_t gc_count;
	QC_Read(void) : contigs(0), n1_runs(0), n2_runs(0), n1_count(0), n2_count(0), lq_count(0), gc_count(0) { }
	explicit QC_Read(const std::string &__s) : Read(__s), contigs(0), n1_runs(0), n2_runs(0), n1_count(0), n2_count(0), lq_count(0), gc_count(0) { }
	~QC_Read(void) { }
	void calc_stats(std::map<size_t, unsigned int> &, std::map<size_t, unsigned int> &);
};

extern bool opt_print_n_quality;

#endif // !_QC_READ_H
