#include "hashz.h"	/* hashz */
#include "hist_lib_hashz.h"
#include "pattern.h"	/* Pattern */
#include "read.h"	/* Read */
#include "time_used.h"	/* elapsed_time(), start_time() */
#include <ctype.h>	/* lowercase() */
#include <gmp.h>	/* mpz_add(), mpz_add_ui(), mpz_clear(), mpz_clrbit(), mpz_cmp(), mpz_fdiv_q_2exp(), mpz_init2(), mpz_init_set(), mpz_mul_2exp(), mpz_realloc2(), mpz_setbit(), mpz_set_ui(), mpz_sizeinbase(), mpz_sub_ui(), mpz_tstbit() */
#include <list>		/* list<> */
#include <map>		/* map<> */
#include <stdio.h>	/* fprintf(), stderr */
#include <string>	/* string */
#include <sys/types.h>	/* size_t */

Pattern opt_include;
bool opt_feedback = 1;
bool opt_mask_lowercase = 0;
hashz::value_type opt_repeat_threshold = 20;
hashz::value_type opt_repeat_threshold_upper = static_cast<hashz::value_type>(-1);
int opt_phred20_anchor = -1;
int opt_repeat_coverage = 1;
size_t opt_skip_size = 0;
std::map<std::string, bool> opt_exclude;

/* these four are constants calculated from the mer length */
static hashz::key_type bp_comp[4];
static unsigned long mer_bits;
static unsigned long mer_length;		/* mer_length - 1 */

#define INCREMENT(key,comp_key,i) \
	mpz_clrbit(key, mer_bits - 1); \
	mpz_clrbit(key, mer_bits - 2); \
	mpz_mul_2exp(key, key, 2); \
	mpz_add_ui(key, key, i); \
	mpz_fdiv_q_2exp(comp_key, comp_key, 2); \
	mpz_add(comp_key, comp_key, bp_comp[i]);

/*
 * given the sequence, create the key and comped key for the first mer
 * length - 1 proper (i.e., ACGT) base pairs, returning the current
 * position in the sequence (or end, if there aren't at least mer length
 * proper base pairs); assumes key and comp_key are already initialized
 */

static size_t preload_keys(const Read &a, size_t s, size_t end, hashz::key_type &key, hashz::key_type &comp_key) {
	a.next_good_sequence(s);
	if (s == a.size()) {	/* no good characters left */
		return end;
	}
	size_t end2 = s + mer_length;
	if (end2 > end) {		/* less than mer length sequence */
		return end;
	}
	mpz_set_ui(key, 0);
	mpz_set_ui(comp_key, 0);
	for (; s != end2; ++s) {
		int i = a.get_seq(s);
		if (i != -1) {
			INCREMENT(key, comp_key, i);
		} else {	// non-base character - advance to
				// the next proper base and start over
			++s;
			a.next_good_sequence(s);
			if (s == a.size()) {	/* no good characters left */
				return end;
			}
			end2 = s + mer_length;
			if (end2 > end) {
				/* less than mer length sequence left */
				return end;
			}
			mpz_set_ui(key, 0);
			mpz_set_ui(comp_key, 0);
			--s;
		}
	}
	return s;
}

/*
 * find n-mer's and count up how many of each there are, looking at both
 * forward and comped versions of the sequence
 */

bool add_sequence_mers(std::list<Read>::const_iterator a, std::list<Read>::const_iterator end_a, hashz &mer_list) {
	if (opt_feedback) {
		start_time();
		fprintf(stderr, "%lu : %10lu entries used (%5.2f%%), %lu overflow\n", time(NULL), mer_list.size(), (double)100 * mer_list.size() / mer_list.capacity(), mer_list.overflow_size());
	}
	hashz::key_type key, comp_key;
	mpz_init2(key, mer_bits);
	mpz_init2(comp_key, mer_bits);
	for (; a != end_a; ++a) {
		/* print feedback every 10 minutes */
		if (opt_feedback && elapsed_time() >= 600) {
			start_time();
			fprintf(stderr, "%lu : %10lu entries used (%5.2f%%), %lu overflow\n", time(NULL), mer_list.size(), (double)100 * mer_list.size() / mer_list.capacity(), mer_list.overflow_size());
		}
		if ((!opt_include.empty() && !opt_include.is_match(a->name())) || opt_exclude.find(a->name()) != opt_exclude.end()) {
			continue;
		}
		if (a->size() < opt_skip_size) {
			continue;
		}
		size_t end = a->quality_stop;
		/* set key with first n-mer - 1 bases */
		size_t s = preload_keys(*a, a->quality_start, end, key, comp_key);
		for (; s != end; ++s) {
			int i = a->get_seq(s);
			if (i == -1) {	/* non-base character - start over */
				s = preload_keys(*a, s, end, key, comp_key);
				--s;
				continue;
			}
			INCREMENT(key, comp_key, i);
			if (!mer_list.increment(mpz_cmp(key, comp_key) < 0 ? key : comp_key)) {
				return 0;
			}
		}
	}
	mpz_clear(key);
	mpz_clear(comp_key);
	if (opt_feedback) {
		fprintf(stderr, "%lu : %10lu entries used (%5.2f%%), %lu overflow\n", time(NULL), mer_list.size(), (double)100 * mer_list.size() / mer_list.capacity(), mer_list.overflow_size());
	}
	return 1;
}

bool add_sequence_mers(std::list<Read>::const_iterator a, std::list<Read>::const_iterator end_a, hashz &mer_list, const std::map<std::string, hashz::offset_type> &opt_readnames_exclude) {
	if (opt_feedback) {
		start_time();
		fprintf(stderr, "%lu : %10lu entries used (%5.2f%%), %lu overflow\n", time(NULL), mer_list.size(), (double)100 * mer_list.size() / mer_list.capacity(), mer_list.overflow_size());
	}
	hashz::key_type key, comp_key;
	mpz_init2(key, mer_bits);
	mpz_init2(comp_key, mer_bits);
	for (; a != end_a; ++a) {
		/* print feedback every 10 minutes */
		if (opt_feedback && elapsed_time() >= 600) {
			start_time();
			fprintf(stderr, "%lu : %10lu entries used (%5.2f%%), %lu overflow\n", time(NULL), mer_list.size(), (double)100 * mer_list.size() / mer_list.capacity(), mer_list.overflow_size());
		}
		if (!opt_include.empty() && !opt_include.is_match(a->name())) {
			continue;
		}
		hashz::offset_type x = 0;
		if (!opt_readnames_exclude.empty()) {
			std::map<std::string, hashz::offset_type>::const_iterator b = opt_readnames_exclude.find(a->name());
			if (b != opt_readnames_exclude.end()) {
				x = b->second;
			}
		}
		size_t end = a->quality_stop;
		/* set key with first n-mer - 1 bases */
		size_t s = preload_keys(*a, a->quality_start, end, key, comp_key);
		for (; s != end; ++s) {
			int i = a->get_seq(s);
			if (i == -1) {	/* non-base character - start over */
				s = preload_keys(*a, s, end, key, comp_key);
				--s;
				continue;
			}
			INCREMENT(key, comp_key, i);
			if (x) {
				if (!mer_list.increment_alt(mpz_cmp(key, comp_key) < 0 ? key : comp_key, x)) {
					return 0;
				}
			} else if (!mer_list.increment(mpz_cmp(key, comp_key) < 0 ? key : comp_key)) {
				return 0;
			}
		}
	}
	mpz_clear(key);
	mpz_clear(comp_key);
	if (opt_feedback) {
		fprintf(stderr, "%lu : %10lu entries used (%5.2f%%), %lu overflow\n", time(NULL), mer_list.size(), (double)100 * mer_list.size() / mer_list.capacity(), mer_list.overflow_size());
	}
	return 1;
}

/* convert key to sequence */

std::string convert_key(const hashz::key_type &key) {
	const char values[4] = { 'A', 'C', 'G', 'T' };
	std::string sequence;
	long i;
	for (i = mer_bits - 1; i != -1; i -= 2) {
		sequence += values[mpz_tstbit(key, i) * 2 + mpz_tstbit(key, i - 1)];
	}
	return sequence;
}

/* initialize mer-related constants */

void init_mer_constants(unsigned long opt_mer_length) {
	mer_length = opt_mer_length - 1;
	mer_bits = opt_mer_length * 2;
	mpz_init2(bp_comp[0], mer_bits);
	mpz_init2(bp_comp[1], mer_bits);
	mpz_init2(bp_comp[2], mer_bits);
	mpz_init2(bp_comp[3], mer_bits);
	mpz_set_ui(bp_comp[0], 3);
	mpz_set_ui(bp_comp[1], 2);
	mpz_set_ui(bp_comp[2], 1);
	mpz_set_ui(bp_comp[3], 0);
	mpz_mul_2exp(bp_comp[0], bp_comp[0], mer_bits - 2);
	mpz_mul_2exp(bp_comp[1], bp_comp[1], mer_bits - 2);
	mpz_mul_2exp(bp_comp[2], bp_comp[2], mer_bits - 2);
}

/* count number of kmers, repetitive kmers, and unique repetitive kmers */

void count_kmers(const Read &a, const hashz &mer_list, size_t &kmers, size_t &r_kmers, size_t &ur_kmers) {
	kmers = r_kmers = 0;
	if (!opt_include.empty() && !opt_include.is_match(a.name())) {
		return;
	}
	hashz::key_type key, comp_key;
        mpz_init2(key, mer_bits + 1);
        mpz_set_ui(key, 1);
        mpz_mul_2exp(key, key, mer_bits);
        mpz_sub_ui(key, key, 1);
        char *key_cstr = new char[mpz_sizeinbase(key, 62) + 2];
	mpz_realloc2(key, mer_bits);
	mpz_init2(comp_key, mer_bits);
	std::map<std::string, bool> r_kmers_list; /* list of repetitive kmers */
	size_t end = a.quality_stop;
	/* set key with first n-mer - 1 bases */
	size_t s = preload_keys(a, a.quality_start, end, key, comp_key);
	for (; s < end; ++s) {
		int i = a.get_seq(s);
		if (i == -1) {
			s = preload_keys(a, s, end, key, comp_key);
			--s;
			continue;
		}
		INCREMENT(key, comp_key, i);
		++kmers;
		/* is this repetitive enough to count as a repeat? */
		hashz::value_type x = mer_list.value(mpz_cmp(key, comp_key) < 0 ? key : comp_key);
		if (opt_repeat_threshold <= x && x < opt_repeat_threshold_upper) {
			++r_kmers;
			mpz_get_str(key_cstr, 62, mpz_cmp(key, comp_key) < 0 ? key : comp_key);
			r_kmers_list[key_cstr] = 1;
		}
	}
	mpz_clear(key);
	mpz_clear(comp_key);
	ur_kmers = r_kmers_list.size();
}

/* check to see if position s should be masked */

static void check_mask(size_t s, const std::list<int> &window, int total, std::string &mask) {
	if (total >= opt_repeat_coverage) {	/* add mask */
		mask[s] = 'X';
	} else if ((size_t)total < window.size()) {
	} else if (s > 0 && mask[s - 1] == 'X') { /* extend existing mask */
		mask[s] = 'X';
		return;
	} else {		/* conditional - will match whatever follows */
		mask[s] = 'R';
		return;
	}
	/* check for conditional end */
	if (s > 0 && mask[s - 1] == 'R') {
		char c = mask[s];
		do {
			mask[--s] = c;
		} while (s != 0 && mask[s - 1] == 'R');
	}
}

/* create mask for highly repetitive regions - 'X's are to be masked out */

static void create_mask(const Read &a, const hashz &mer_list, std::string &mask) {
	mask.resize(a.size(), ' ');
	hashz::key_type key, comp_key;
	mpz_init2(key, mer_bits);
	mpz_init2(comp_key, mer_bits);
	int total = 0;		/* number of repeats in current window */
	std::list<int> window;
	size_t end = a.quality_stop;
	/* set key with first n-mer - 1 bases */
	size_t s = preload_keys(a, a.quality_start, end, key, comp_key);
	for (; s < end; ++s) {
		int i = a.get_seq(s);
		if (i == -1) {	// non-base character - advance to
				// the next proper base and start over
			size_t t = s - mer_length;
			/* fill out window for short sections */
			window.insert(window.begin(), mer_length + 1 - window.size(), 0);
			for (; window.size() > 1; ++t) {
				total -= window.front();
				window.pop_front();
				check_mask(t, window, total, mask);
			}
			total = 0;
			window.clear();
			s = preload_keys(a, s, end, key, comp_key);
			--s;
			continue;
		}
		INCREMENT(key, comp_key, i);
		/* if window is full sized (mer length), pop first value */
		if (window.size() == (size_t)(mer_length + 1)) {
			total -= window.front();
			window.pop_front();
		}
		/* is this repetitive enough to count as a repeat? */
		hashz::value_type x = mer_list.value(mpz_cmp(key, comp_key) < 0 ? key : comp_key);
		int j = (opt_repeat_threshold <= x && x < opt_repeat_threshold_upper) ? 1 : 0;
		total += j;
		window.push_back(j);
		check_mask(s - mer_length, window, total, mask);
	}
	mpz_clear(key);
	mpz_clear(comp_key);
	/* fill out window for short sections */
	window.insert(window.begin(), mer_length + 1 - window.size(), 0);
	for (s -= mer_length; window.size() > 1; ++s) {
		total -= window.front();
		window.pop_front();
		check_mask(s, window, total, mask);
	}
}

/*
 * find the start of the first unmasked region that has at least
 * opt_phred20_anchor phred20's, and the end of the last such region
 * (non-ACGT base pairs count as masking, in that regions can't contain them)
 */

static void find_phred20_anchors(const Read &a, const std::string &mask, size_t *phred20_start, size_t *phred20_stop) {
	int total = 0;		/* number of phred20's since last repeat */
	size_t s = a.quality_start;
	size_t end = a.quality_stop;
	size_t start = end;
	size_t stop = end;
	size_t last = s;
	/* get start */
	for (; s != end; ++s) {
		if (mask[s] == 'X') {
			total = 0;
			last = s + 1;
		} else if (!a.is_good_basepair(s)) {
			total = 0;
			last = s + 1;
		} else if (a.is_high_quality(s)) {
			++total;
			if (total == opt_phred20_anchor) {
				start = last;
				break;
			}
		}
	}
	/* if there's a start, get the end */
	if (s != end) {
		total = 0;
		s = end - 1;
		last = s;
		for (;; --s) {
			if (mask[s] == 'X') {
				total = 0;
				last = s - 1;
			} else if (!a.is_good_basepair(s)) {
				total = 0;
				last = s - 1;
			} else if (a.is_high_quality(s)) {
				++total;
				if (total == opt_phred20_anchor) {
					stop = last;
					break;
				}
			}
		}
	}
	*phred20_start = start;
	*phred20_stop = stop;
}

/*
 * mask out regions of the read's sequence wherever mask has an X, except
 * where anchored by phred20's, if start and stop are set
 */

static void mask_repeats(Read &a, const std::string &mask, size_t phred20_start, size_t phred20_stop) {
	/* don't mask between phred20_start and phred20_stop */
	size_t s;
	for (s = a.quality_start; s < phred20_start; ++s) {
		if (mask[s] == 'X') {
			a.set_sequence(s, 'X');
		}
	}
	for (s = phred20_stop + 1; s < a.quality_stop; ++s) {
		if (mask[s] == 'X') {
			a.set_sequence(s, 'X');
		}
	}
}

/* same as mask_repeats, but changes sequence to lowercase instead of X */

#ifndef COMPRESS_READS
static void mask_repeats_lowercase(Read &a, const std::string &mask, size_t phred20_start, size_t phred20_stop) {
	/* don't mask between phred20_start and phred20_stop */
	size_t s;
	for (s = a.quality_start; s < phred20_start; ++s) {
		if (mask[s] == 'X') {
			a.set_sequence(s, tolower(a.get_sequence(s)));
		}
	}
	for (s = phred20_stop + 1; s < a.quality_stop; ++s) {
		if (mask[s] == 'X') {
			a.set_sequence(s, tolower(a.get_sequence(s)));
		}
	}
}
#endif

/*
 * mask out highly repetitive regions in the read's sequence, unless
 * they're anchored by opt_phred20_anchor phred20's on both sides; the
 * low-quality tails are ignored; returns number of masked basepairs
 */

void screen_repeats(Read &a, const hashz &mer_list) {
	if (!opt_include.empty() && !opt_include.is_match(a.name())) {
		return;
	}
	std::string mask;
	create_mask(a, mer_list, mask);
	size_t phred20_start;
	size_t phred20_stop;
	if (opt_phred20_anchor == -1) {
		phred20_start = phred20_stop = a.quality_stop;
	} else {
		find_phred20_anchors(a, mask, &phred20_start, &phred20_stop);
	}
#ifndef COMPRESS_READS
	if (opt_mask_lowercase) {
		mask_repeats_lowercase(a, mask, phred20_start, phred20_stop);
	} else
#endif
	{
		mask_repeats(a, mask, phred20_start, phred20_stop);
	}
}

/* check to see if base pair is highly repetitive */

static int check_unique(bool is_phred20, const std::list<int> &window, int total, int *state) {
	if (total >= opt_repeat_coverage) {	/* highly repetitive */
		*state = -2;
		return 0;
	} else if ((size_t)total < window.size()) {
		int i = *state > 0 ? *state : 0;	/* not a repeat */
		*state = -1;
		return i + (is_phred20 ? 1 : 0);
	} else if (*state == -2) {		/* extend existing repeat */
		return 0;
	} else if (*state == -1) {		/* start conditional */
		*state = is_phred20 ? 1 : 0;
		return 0;
	} else if (is_phred20) { /* extend conditional (with increase) */
		++(*state);
		return 0;
	} else {		/* extend conditional (without increase) */
		return 0;
	}
}

/*
 * count phred20's of non-highly repetitive base pairs; can also return
 * total number of phred20's for comparison
 */

static unsigned long count_phreds(const Read &a, const hashz &mer_list, unsigned long *total_phreds_out) {
	unsigned long total_phreds = 0;
	unsigned long total_unique_phreds = 0;
	hashz::key_type key, comp_key;
	mpz_init2(key, mer_bits);
	mpz_init2(comp_key, mer_bits);
	int total = 0;		/* number of repeats in current window */
	std::list<int> window;
	size_t end = a.quality_stop;
	/* set key with first n-mer - 1 bases */
	size_t s = preload_keys(a, a.quality_start, end, key, comp_key);
	/*
	 * -2 if last basepair was repetitive, -1 if it was not, non-negative
	 * if it might count (value equal to size of run of maybes)
	 */
	int state = -1;
	for (; s < end; ++s) {
		int i = a.get_seq(s);
		if (i == -1) {	// non-base character - advance to
				// the next proper base and start over
			size_t t = s - mer_length;
			for (; window.size() > 1; ++t) {
				total -= window.front();
				window.pop_front();
				if (a.is_high_quality(s)) {
					++total_phreds;
				}
				total_unique_phreds += check_unique(a.is_high_quality(s), window, total, &state);
			}
			total = 0;
			window.clear();
			s = preload_keys(a, s, end, key, comp_key);
			--s;
			continue;
		}
		INCREMENT(key, comp_key, i);
		/* if window is full sized (mer length), pop first value */
		if (window.size() == (size_t)(mer_length + 1)) {
			total -= window.front();
			window.pop_front();
		}
		/* is this repetitive enough to count as a repeat? */
		hashz::value_type x = mer_list.value(mpz_cmp(key, comp_key) < 0 ? key : comp_key);
		int j = (opt_repeat_threshold <= x && x < opt_repeat_threshold_upper) ? 1 : 0;
		total += j;
		window.push_back(j);
		if (a.is_high_quality(s)) {
			++total_phreds;
		}
		total_unique_phreds += check_unique(a.is_high_quality(s), window, total, &state);
	}
	mpz_clear(key);
	mpz_clear(comp_key);
	for (s -= mer_length; window.size() > 1; ++s) {
		total -= window.front();
		window.pop_front();
		if (a.is_high_quality(s)) {
			++total_phreds;
		}
		total_unique_phreds += check_unique(a.is_high_quality(s), window, total, &state);
	}
	if (state > 0) { /* conditional collapses to non-highly repetitive */
		total_unique_phreds += state;
	}
	if (total_phreds_out != NULL) {
		*total_phreds_out = total_phreds;
	}
	return total_unique_phreds;
}

/*
 * count the number of non-highly repetitive phred20's in each read; also
 * returns the total number of phred20's, for comparison
 */

unsigned long count_unique_phreds(const std::list<Read> &read_list, const hashz &mer_list, unsigned long *total_phreds_out) {
	unsigned long total_phreds = 0;
	unsigned long total_unique_phreds = 0;
	std::list<Read>::const_iterator a = read_list.begin();
	std::list<Read>::const_iterator end = read_list.end();
	for (; a != end; ++a) {
		unsigned long read_phreds;
		total_unique_phreds += count_phreds(*a, mer_list, &read_phreds);
		total_phreds += read_phreds;
	}
	if (total_phreds_out != NULL) {
		*total_phreds_out = total_phreds;
	}
	return total_unique_phreds;
}

/* assumes key_out is already initialized */
void reverse_key(const hashz::key_type &key_in, hashz::key_type &key_out) {
	mpz_set_ui(key_out, 0);
	unsigned long i, j;
	for (i = 0, j = mer_bits - 1; i != mer_bits; i += 2, j -= 2) {
		if (!mpz_tstbit(key_in, i + 1)) {
			mpz_setbit(key_out, j);
		}
		if (!mpz_tstbit(key_in, i)) {
			mpz_setbit(key_out, j - 1);
		}
	}
}
