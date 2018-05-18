#include "qc_read.h"
#include "read.h"	/* opt_quality_cutoff */
#include <map>		/* map<> */
#include <string>	/* string */
#include <sys/types.h>	/* size_t */
#include <vector>	/* vector<> */

bool opt_print_n_quality = 0;

/* print a range of quality values */

void QC_Read::print_quality_range(size_t i, size_t end) const {
	size_t j = end - i;
	for (; i != end; ++i) {
		printf("%lu	%u\n", j, get_quality(i));
	}
}

/*
 * check the sequence to find the number of N's, and number of runs of them;
 * also, count the number of low quality bases (<40, non-N)
 */

void QC_Read::calc_stats(std::map<size_t, unsigned int> &n_hist, std::map<size_t, unsigned int> &lq_hist) {
	size_t end = size();
	contigs = 1;
	n1_runs = 0;
	n1_count = 0;
	n2_runs = 0;
	n2_count = 0;
	size_t i;
	for (i = 0; i != end && get_sequence(i) != 'N'; ++i) { }
	if (i == 0 && get_quality(0) != 0) {
		contigs = 0;	/* starting gap-N - don't count empty contig */
	}
	while (i != size()) {
		unsigned char x = get_quality(i);
		size_t j = i;
		for (++i; i != end && get_sequence(i) == 'N' && get_quality(i) == x; ++i) { }
		++n_hist[i - j];
		if (x != 0) {	/* an N-run of non-zero quality is a gap */
			++contigs;
			++n1_runs;
			n1_count += i - j;
		} else {
			++n2_runs;
			n2_count += i - j;
		}
		if (opt_print_n_quality) {
			print_quality_range(j, i);
		}
		for (; i != end && get_sequence(i) != 'N'; ++i) { }
	}
	if (get_sequence(end - 1) == 'N' && get_quality(end - 1) != 0) {
		--contigs;	/* ending gap-N - don't count empty contig */
	}
	/* count all low quality base pairs */
	lq_count = 0;
	/* only count large contigs (>= 8k non-gap bp) for lq histogram */
	if (size() - n1_count >= 8000) {
		for (i = 0; i != end; ++i) {
			unsigned char x = get_quality(i);
			if (x < opt_quality_cutoff) {
				++lq_count;
				++lq_hist[x / 5];
			}
		}
		lq_hist[0] -= n1_count;
	} else {
		for (i = 0; i != end; ++i) {
			if (get_quality(i) < opt_quality_cutoff) {
				++lq_count;
			}
		}
	}
	/* subtract out gap N's (which are all quality 1) */
	lq_count -= n1_count;
	gc_count = 0;
	for (i = 0; i != end; ++i) {
		switch (get_sequence(i)) {
		    case 'C':
		    case 'G':
			++gc_count;
			break;
		}
	}
}
