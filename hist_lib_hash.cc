#include "hash.h"	// hash
#include "hist_lib_hash.h"
#include "pattern.h"	// Pattern
#include "read.h"	// Read
#include "time_used.h"	// elapsed_time(), start_time()
#include <algorithm>	// swap()
#include <ctype.h>	// lowercase()
#include <list>		// list<>
#include <map>		// map<>
#include <stdio.h>	// fprintf(), stderr
#include <stdlib.h>	// exit()
#include <string>	// string
#include <sys/types.h>	// size_t

size_t opt_mer_length;	// this is actually mer length - 1,
			// for convenience of calculations below
Pattern opt_include;
bool opt_feedback = 1;
bool opt_mask_lowercase = 0;
bool opt_reverse_mask = 0;
hash::value_type opt_repeat_threshold = 20;
hash::value_type opt_repeat_threshold_upper = static_cast<hash::value_type>(-1);
int opt_phred20_anchor = -1;
size_t opt_repeat_coverage = 1;
size_t opt_skip_size = 0;
std::map<std::string, bool> opt_exclude;

// these four are constants calculated from the mer length
static hash::key_type bp_comp[4];
static hash::key_type mer_mask;

// given the sequence, create the key and comped key for the first mer
// length - 1 proper (i.e., ACGT) base pairs, returning the current
// position in the sequence (or end, if there aren't at least mer length
// proper base pairs); key and comp_key are presumed to have no set bits
// greater than the mer_mask

static size_t preload_keys(const Read &a, size_t s, size_t end, hash::key_type &key, hash::key_type &comp_key) {
	a.next_good_sequence(s);
	if (s == a.size()) {	// no good characters left
		return end;
	}
	size_t end2 = s + opt_mer_length;
	if (end2 > end) {		// less than mer length sequence
		return end;
	}
	for (; s != end2; ++s) {
		int i = a.get_seq(s);
		if (i != -1) {
			key = ((key << 2) & mer_mask) | i;
			comp_key = (comp_key >> 2) | bp_comp[i];
		} else {	// non-base character - advance to
				// the next proper base and start over
			++s;
			a.next_good_sequence(s);
			if (s == a.size()) {	// no good characters left
				return end;
			}
			end2 = s + opt_mer_length;
			if (end2 > end) {
				// less than mer length sequence left
				return end;
			}
			--s;
		}
	}
	return s;
}

// find n-mer's and count up how many of each there are, looking at both
// forward and comped versions of the sequence

bool add_sequence_mers(std::list<Read>::const_iterator a, const std::list<Read>::const_iterator end_a, hash &mer_list) {
	for (; a != end_a; ++a) {
		// print feedback every 10 minutes
		if (opt_feedback && elapsed_time() >= 600) {
			start_time();
			fprintf(stderr, "%lu : %10lu entries used (%5.2f%%), %lu overflow\n", time(NULL), mer_list.size(), static_cast<double>(100) * mer_list.size() / mer_list.capacity(), mer_list.overflow_size());
		}
		if (a->size() < opt_skip_size) {
			continue;
		}
		if ((!opt_include.empty() && !opt_include.is_match(a->name())) || opt_exclude.find(a->name()) != opt_exclude.end()) {
			continue;
		}
		hash::key_type key(0);
		hash::key_type comp_key(0);
		const size_t end(a->quality_stop);
		// set key with first n-mer - 1 bases
		size_t s(preload_keys(*a, a->quality_start, end, key, comp_key));
		for (; s != end; ++s) {
			const int i(a->get_seq(s));
			if (i == -1) {	// non-base character - start over
				s = preload_keys(*a, s, end, key, comp_key);
				--s;
				continue;
			}
			key = ((key << 2) & mer_mask) | i;
			comp_key = (comp_key >> 2) | bp_comp[i];
			if (!mer_list.increment(key < comp_key ? key : comp_key)) {
				return 0;
			}
		}
	}
	return 1;
}

bool add_sequence_mers(std::list<Read>::const_iterator a, std::list<Read>::const_iterator end_a, hash &mer_list, const std::map<std::string, hash::offset_type> &opt_readnames_exclude) {
	for (; a != end_a; ++a) {
		// print feedback every 10 minutes
		if (opt_feedback && elapsed_time() >= 600) {
			start_time();
			fprintf(stderr, "%lu : %10lu entries used (%5.2f%%), %lu overflow\n", time(NULL), mer_list.size(), static_cast<double>(100) * mer_list.size() / mer_list.capacity(), mer_list.overflow_size());
		}
		if (!opt_include.empty() && !opt_include.is_match(a->name())) {
			continue;
		}
		hash::offset_type x = 0;
		if (!opt_readnames_exclude.empty()) {
			std::map<std::string, hash::offset_type>::const_iterator b = opt_readnames_exclude.find(a->name());
			if (b != opt_readnames_exclude.end()) {
				x = b->second;
			}
		}
		hash::key_type key = 0;
		hash::key_type comp_key = 0;
		size_t end = a->quality_stop;
		// set key with first n-mer - 1 bases
		size_t s = preload_keys(*a, a->quality_start, end, key, comp_key);
		for (; s != end; ++s) {
			int i = a->get_seq(s);
			if (i == -1) {	// non-base character - start over
				s = preload_keys(*a, s, end, key, comp_key);
				--s;
				continue;
			}
			key = ((key << 2) & mer_mask) | i;
			comp_key = (comp_key >> 2) | bp_comp[i];
			if (x) {
				if (!mer_list.increment_alt(key < comp_key ? key : comp_key, x)) {
					return 0;
				}
			} else if (!mer_list.increment(key < comp_key ? key : comp_key)) {
				return 0;
			}
		}
	}
	return 1;
}

// converts key to sequence

std::string convert_key(hash::key_type key) {
	const char values[4] = { 'A', 'C', 'G', 'T' };
	std::string sequence(opt_mer_length + 1, '\0');
	for (size_t i(0); i <= opt_mer_length; ++i, key >>= 2) {
		sequence[opt_mer_length - i] = values[key & 3];
	}
	return sequence;
}

// initialize mer-related constants

void init_mer_constants() {
	// start it once, here, so multiple files and batches don't get
	// extra feedbacks
	if (opt_feedback) {
		start_time();
	}
	if (opt_mer_length > sizeof(hash::key_type) * 4) {
		fprintf(stderr, "Error: mer length too long: %d > %d\n", static_cast<int>(opt_mer_length), static_cast<int>(sizeof(hash::key_type) * 4));
		exit(1);
	} else if (opt_mer_length == sizeof(hash::key_type) * 4) {
		// have to special case this one
		mer_mask = static_cast<hash::key_type>(-1);
	} else {
		mer_mask = (static_cast<hash::key_type>(1) << (2 * opt_mer_length)) - 1;
	}
	--opt_mer_length;
	bp_comp[0] = static_cast<hash::key_type>(3) << (2 * opt_mer_length);
	bp_comp[1] = static_cast<hash::key_type>(2) << (2 * opt_mer_length);
	bp_comp[2] = static_cast<hash::key_type>(1) << (2 * opt_mer_length);
	bp_comp[3] = 0;
}

void print_final_input_feedback(const hash &mer_list) {
	if (opt_feedback && mer_list.size() != 0) {
		fprintf(stderr, "%lu : %10lu entries used (%5.2f%%), %lu overflow\n", time(NULL), mer_list.size(), static_cast<double>(100) * mer_list.size() / mer_list.capacity(), mer_list.overflow_size());
	}
}

void clear_mer_list(hash &mer_list) {
	print_final_input_feedback(mer_list);
	mer_list.clear();
}

// count number of kmers, repetitive kmers, and unique repetitive kmers

void count_kmers(const Read &a, const hash &mer_list, size_t &kmers, size_t &r_kmers, size_t &ur_kmers) {
	kmers = r_kmers = 0;
	if (!opt_include.empty() && !opt_include.is_match(a.name())) {
		return;
	}
	hash::key_type key(0);
	hash::key_type comp_key(0);
	std::map<hash::key_type, char> r_kmers_list; // list of repetitive kmers
	size_t end(a.quality_stop);
	// set key with first n-mer - 1 bases
	size_t s(preload_keys(a, a.quality_start, end, key, comp_key));
	for (; s < end; ++s) {
		const int i(a.get_seq(s));
		if (i == -1) {
			s = preload_keys(a, s, end, key, comp_key);
			--s;
			continue;
		}
		key = ((key << 2) & mer_mask) | i;
		comp_key = (comp_key >> 2) | bp_comp[i];
		++kmers;
		// is this repetitive enough to count as a repeat?
		const hash::value_type x(mer_list.value(key < comp_key ? key : comp_key));
		if (opt_repeat_threshold <= x && x < opt_repeat_threshold_upper) {
			++r_kmers;
			r_kmers_list[key < comp_key ? key : comp_key] = 0;
		}
	}
	ur_kmers = r_kmers_list.size();
}

// check to see if position s should be masked

static void check_mask(size_t s, const std::list<int> &window, size_t total, std::string &mask) {
	if (total >= opt_repeat_coverage) {	// add mask
		mask[s] = 'X';
	} else if (total < window.size()) {
	} else if (s > 0 && mask[s - 1] == 'X') { // extend existing mask
		mask[s] = 'X';
		return;
	} else {		// conditional - will match whatever follows
		mask[s] = 'R';
		return;
	}
	// check for conditional end
	if (s > 0 && mask[s - 1] == 'R') {
		char c = mask[s];
		do {
			mask[--s] = c;
		} while (s != 0 && mask[s - 1] == 'R');
	}
}

// create mask for highly repetitive regions - 'X's are to be masked out

static void create_mask(const Read &a, const hash &mer_list, std::string &mask) {
	mask.resize(a.size(), ' ');
	hash::key_type key = 0;
	hash::key_type comp_key = 0;
	int total = 0;		// number of repeats in current window
	std::list<int> window;
	size_t end = a.quality_stop;
	// set key with first n-mer - 1 bases
	size_t s = preload_keys(a, a.quality_start, end, key, comp_key);
	for (; s < end; ++s) {
		int i = a.get_seq(s);
		if (i == -1) {	// non-base character - advance to
				// the next proper base and start over
			size_t t = s - opt_mer_length;
			// fill out window for short sections
			window.insert(window.begin(), opt_mer_length + 1 - window.size(), 0);
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
		key = ((key << 2) & mer_mask) | i;
		comp_key = (comp_key >> 2) | bp_comp[i];
		// if window is full sized (mer length), pop first value
		if (window.size() == opt_mer_length + 1) {
			total -= window.front();
			window.pop_front();
		}
		// is this repetitive enough to count as a repeat?
		hash::value_type x = mer_list.value(key < comp_key ? key : comp_key);
		int j = opt_repeat_threshold <= x && x < opt_repeat_threshold_upper ? 1 : 0;
		total += j;
		window.push_back(j);
		check_mask(s - opt_mer_length, window, total, mask);
	}
	// fill out window for short sections
	window.insert(window.begin(), opt_mer_length + 1 - window.size(), 0);
	for (s -= opt_mer_length; window.size() > 1; ++s) {
		total -= window.front();
		window.pop_front();
		check_mask(s, window, total, mask);
	}
}

// find the start of the first unmasked region that has at least
// opt_phred20_anchor phred20's, and the end of the last such region
// (non-ACGT base pairs count as masking, in that regions can't contain them)

static void find_phred20_anchors(const Read &a, const std::string &mask, size_t *phred20_start, size_t *phred20_stop) {
	int total = 0;		// number of phred20's since last repeat
	size_t s = a.quality_start;
	size_t end = a.quality_stop;
	size_t start = end;
	size_t stop = end;
	size_t last = s;
	// get start
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
	// if there's a start, get the end
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

// mask out regions of the read's sequence wherever mask has an X, except
// where anchored by phred20's, if start and stop are set

static void mask_repeats(Read &a, const std::string &mask, size_t phred20_start, size_t phred20_stop) {
	// don't mask between phred20_start and phred20_stop
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

// same as mask_repeats, but changes sequence to lowercase instead of X

#ifndef COMPRESS_READS
static void mask_repeats_lowercase(Read &a, const std::string &mask, size_t phred20_start, size_t phred20_stop) {
	// don't mask between phred20_start and phred20_stop
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

// mask out highly repetitive regions in the read's sequence, unless
// they're anchored by opt_phred20_anchor phred20's on both sides; the
// low-quality tails are ignored; returns number of masked basepairs

void screen_repeats(Read &a, const hash &mer_list) {
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
	if (opt_reverse_mask) {
		for (size_t i(0); i != mask.size(); ++i) {
			mask[i] = mask[i] == 'X' ? ' ' : 'X';
		}
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

// check to see if base pair is highly repetitive

static int check_unique(bool is_phred20, const std::list<int> &window, size_t total, int *state) {
	if (total >= opt_repeat_coverage) {	// highly repetitive
		*state = -2;
		return 0;
	} else if ((size_t)total < window.size()) {
		int i = *state > 0 ? *state : 0;	// not a repeat
		*state = -1;
		return i + (is_phred20 ? 1 : 0);
	} else if (*state == -2) {		// extend existing repeat
		return 0;
	} else if (*state == -1) {		// start conditional
		*state = is_phred20 ? 1 : 0;
		return 0;
	} else if (is_phred20) { // extend conditional (with increase)
		++(*state);
		return 0;
	} else {		// extend conditional (without increase)
		return 0;
	}
}

// count phred20's of non-highly repetitive base pairs; can also return
// total number of phred20's for comparison

static unsigned long count_phreds(const Read &a, const hash &mer_list, unsigned long *total_phreds_out) {
	unsigned long total_phreds = 0;
	unsigned long total_unique_phreds = 0;
	hash::key_type key = 0;
	hash::key_type comp_key = 0;
	int total = 0;		// number of repeats in current window
	std::list<int> window;
	size_t end = a.quality_stop;
	// set key with first n-mer - 1 bases
	size_t s = preload_keys(a, a.quality_start, end, key, comp_key);
	// -2 if last basepair was repetitive, -1 if it was not, non-negative
	// if it might count (value equal to size of run of maybes)
	int state = -1;
	for (; s < end; ++s) {
		int i = a.get_seq(s);
		if (i == -1) {	// non-base character - advance to
				// the next proper base and start over
			size_t t = s - opt_mer_length;
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
		key = ((key << 2) & mer_mask) | i;
		comp_key = (comp_key >> 2) | bp_comp[i];
		// if window is full sized (mer length), pop first value
		if (window.size() == opt_mer_length + 1) {
			total -= window.front();
			window.pop_front();
		}
		// is this repetitive enough to count as a repeat?
		hash::value_type x = mer_list.value(key < comp_key ? key : comp_key);
		int j = opt_repeat_threshold <= x && x < opt_repeat_threshold_upper ? 1 : 0;
		total += j;
		window.push_back(j);
		if (a.is_high_quality(s)) {
			++total_phreds;
		}
		total_unique_phreds += check_unique(a.is_high_quality(s), window, total, &state);
	}
	for (s -= opt_mer_length; window.size() > 1; ++s) {
		total -= window.front();
		window.pop_front();
		if (a.is_high_quality(s)) {
			++total_phreds;
		}
		total_unique_phreds += check_unique(a.is_high_quality(s), window, total, &state);
	}
	if (state > 0) { // conditional collapses to non-highly repetitive
		total_unique_phreds += state;
	}
	if (total_phreds_out != NULL) {
		*total_phreds_out = total_phreds;
	}
	return total_unique_phreds;
}

// count the number of non-highly repetitive phred20's in each read; also
// returns the total number of phred20's, for comparison

unsigned long count_unique_phreds(const std::list<Read> &read_list, const hash &mer_list, unsigned long *total_phreds_out) {
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

hash::key_type reverse_key(hash::key_type key) {
	size_t i;
	hash::key_type x;
	for (i = 0, x = 0; i <= opt_mer_length; ++i, key >>= 2) {
		x <<= 2;
		x += 3 - (key & 3);
	}
	return x;
}
