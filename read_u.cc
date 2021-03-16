#include "itoa.h"	// itoa()
#include "pattern.h"	// Pattern
#include "read.h"
#include <algorithm>	// min()
#include <limits.h>	// UCHAR_MAX
#include <map>		// map<>
#include <math.h>	// ceil(), floor()
#include <sstream>	// istringstream, ostringstream
#include <stdio.h>	// FILE, fprintf(), stderr
#include <stdlib.h>	// strtol()
#include <string>	// string
#include <sys/types.h>	// size_t
#include <utility>	// make_pair(), pair<>
#include <vector>	// vector<>

Pattern opt_linker;
bool opt_N_is_vector = 0;	// treat N as an X, for finding vector
bool opt_add_range = 0;
bool opt_all_p20 = 1;		// check all basepairs for p20s
bool opt_clip_quality = 0;
bool opt_clip_vector = 0;
bool opt_pacbio = 0;		// perform pacbio name changes when clipping
bool opt_strict_quality = 0;	// harder ends, no interior low quality windows
bool opt_strip_trailing_zero_qual = 0;
double opt_base_cutoff(0);
double opt_repeat_clip(0);	// repeat average to clip from end
int opt_quality_cutoff = 20;	// for setting vector and quality clipping
size_t opt_line_length = 50;	// wrapping point for lines
size_t opt_minimum_clip = 0;
std::map<std::string, std::string> read_name_translation;

static bool *good_base = NULL;

static char comp_lookup[256];

void init_read_comp() {
	for (int i(0); i != 256; ++i) {
		comp_lookup[i] = i;
	}
	comp_lookup['A'] = 'T';
	comp_lookup['C'] = 'G';
	comp_lookup['G'] = 'C';
	comp_lookup['T'] = 'A';
	comp_lookup['a'] = 't';
	comp_lookup['c'] = 'g';
	comp_lookup['g'] = 'c';
	comp_lookup['t'] = 'a';
}

// set read to be the reverse complement of the given read

void Read::set_comp(const Read &a) {
	quality.assign(a.quality.rbegin(), a.quality.rend());
	const size_t length(a.sequence_.size());
	// XXX - vectors
	header = a.header;
	quality_start = length - a.quality_stop;
	quality_stop = length - a.quality_start;
	vector_start = length - a.vector_stop;
	vector_stop = length - a.vector_start;
	phred_count = a.phred_count;
	// I couldn't think up a better way to assign the space,
	// and it does avoid having to do the reverse myself
	sequence_.assign(a.sequence_.rbegin(), a.sequence_.rend());
	for (size_t i(0); i != length; ++i) {
		sequence_[i] = comp_lookup[static_cast<int>(sequence_[i])];
	}
}

// find the largest continuous sequence of non-vector ('X') and set the
// vector start and stop points for the read

void Read::set_vector_endpoints() {
	// check for any non-vector at all
	if (vectors.empty()) {
		return;
	}
	std::vector<std::pair<size_t, size_t> >::const_iterator a(vectors.begin());
	const std::vector<std::pair<size_t, size_t> >::const_iterator end_a(vectors.end());
	size_t best_count(0);
	for (; a != end_a; ++a) {
	const size_t n(count_quality(*a));
		if (best_count < n) {
			best_count = n;
			vector_start = a->first;
			vector_stop = a->second;
		}
	}
	vectors.clear();
}

// gets starting and stopping location of sequence to print,
// and returns if there is, in fact, any sequence to print

bool Read::get_output_endpoints(size_t &i, size_t &end) const {
	if (opt_clip_quality) {
		i = quality_start;
	} else if (opt_clip_vector) {
		i = vector_start;
	} else {
		i = 0;
	}
	if (opt_clip_quality) {
		end = quality_stop;
	} else if (opt_clip_vector) {
		end = vector_stop;
	} else {
		end = size();
	}
	return i != end;
}

// modify pacbio header to reflect trimming changes; return 0 if the
// original header could not be parsed in pacbio style (^[^/]*/\d*/\d+_\d+$)

static bool make_pacbio_header(const std::string &name, const size_t i, const size_t j, std::string &s) {
	size_t k(name.find('/'));
	if (k == static_cast<size_t>(-1)) {
		return 0;
	}
	k = name.find('/', k + 1);
	if (k == static_cast<size_t>(-1)) {
		return 0;
	}
	size_t l(name.find_first_not_of("0123456789", k + 1));
	if (l == static_cast<size_t>(-1) || l == k || name[l] != '_') {
		return 0;
	}
	if (name.find_first_not_of("0123456789", l + 1) != static_cast<size_t>(-1)) {
		return 0;
	}
	s = name.substr(0, k + 1);
	std::istringstream t(name.substr(k + 1, l - (k + 1)));
	size_t x;
	t >> x;
	x += i;
	std::ostringstream u;
	u << x;
	s += u.str();
	s += '_';
	t.str(name.substr(l + 1));
	t.clear();
	t >> x;
	x -= j;		// might check for wrapping here
	u.str("");
	u << x;
	s += u.str();
	return 1;
}

// print a sequence header, adding clip ranges if specified;
// returns if clipping region is not empty; [i,j] are the endpoints

bool Read::print_header(FILE * const fp, const size_t i, const size_t j) const {
	std::string x(name());
	const size_t n(x.size() + 1);
	const std::map<std::string, std::string>::const_iterator a(read_name_translation.find(x));
	if (!opt_add_range) {
		std::string s;
		if (opt_pacbio && (i != 0 || j != size()) && make_pacbio_header(x, i, size() - j, s)) {
			fprintf(fp, ">%s%s\n", s.c_str(), header.substr(n).c_str());
		} else if (a != read_name_translation.end()) {
			fprintf(fp, ">%s%s\n", a->second.c_str(), header.substr(n).c_str());
		} else {
			fprintf(fp, "%s\n", header.c_str());
		}
	} else if (quality_start == quality_stop) {
		return 0;
	} else {
		if (a != read_name_translation.end()) {
			x = a->second;
		}
		if (opt_clip_quality) {
			fprintf(fp, ">%s %d %ld%s\n", x.c_str(), 1, quality_stop - quality_start, header.substr(n).c_str());
		} else if (opt_clip_vector) {
			fprintf(fp, ">%s %ld %ld%s\n", x.c_str(), quality_start + 1 - vector_start, quality_stop - vector_start, header.substr(n).c_str());
		} else {
			fprintf(fp, ">%s %ld %ld%s\n", x.c_str(), quality_start + 1, quality_stop, header.substr(n).c_str());
		}
	}
	return 1;
}

// initialize good_base, for quick lookup of good bases

static void init_phred_counting() {
	if (good_base == NULL) {
		good_base = new bool[UCHAR_MAX + 1];
		for (size_t i(0); i != UCHAR_MAX + 1; ++i) {
			good_base[i] = 0;
		}
		good_base['A'] = 1;
		good_base['C'] = 1;
		good_base['G'] = 1;
		good_base['T'] = 1;
		good_base['a'] = 1;
		good_base['c'] = 1;
		good_base['g'] = 1;
		good_base['t'] = 1;
	}
}

// set the quality start and stop points for the read

void Read::set_quality_endpoints() {
	if (vector_stop < opt_minimum_clip) {
		quality_start = quality_stop = vector_start;
		return;
	}
	size_t i(vector_start);
	const size_t end(vector_stop);
	std::vector<unsigned char> window(20, 0);
	int total(0);
	// find first spot window total is >= opt_quality_cutoff on average
	for (; total < 20 * opt_quality_cutoff && i < end; ++i) {
		total += quality[i] - window[i % 20];
		window[i % 20] = quality[i];
	}
	// no window satisfied the limit, so quality clip is empty
	if (i == end && total < 20 * opt_quality_cutoff) {
		quality_start = quality_stop = vector_start;
		return;
	}
	// be careful when subtracting from unsigned values not to wrap
	quality_start = (i < 20) ? 0 : (i - 20);
	// need to check for this, as, if the total was reached in less than
	// 20 spaces, i - 20 will be before vector_start
	if (quality_start < vector_start) {
		quality_start = vector_start;
	}
	if (quality_start < opt_minimum_clip) {
		quality_start = opt_minimum_clip;
	}
	// re-initialize window
	window.assign(20, 0);
	total = 0;
	i = end - 1;
	// find last spot window total is >= opt_quality_cutoff on overage;
	// don't have to worry about running off the beginning, as we already
	// know that there exists a window that will satisfy the minimum total
	for (; total < 20 * opt_quality_cutoff; --i) {
		total += quality[i] - window[i % 20];
		window[i % 20] = quality[i];
	}
	quality_stop = i + 21;
	// again, if the total is reached in less than 20 spaces,
	// i + 21 could be beyond vector_stop
	if (quality_stop > vector_stop) {
		quality_stop = vector_stop;
	}
	// in case the minimum clip ended up beyond quality stop
	if (quality_stop < quality_start) {
		quality_stop = quality_start;
	}
	if (opt_repeat_clip < 1) {		// do we clip on repeat strings?
		return;
	}
	// find first position such that the remainder of the string has, on
	// average, runs at least opt_repeat_clip in length
	double total_sequence(0);
	double runs(0);
	char last_bp(0);
	size_t last(static_cast<size_t>(-1));
	for (i = quality_stop - 1; i > quality_start; --i) {
		++total_sequence;
		if (last_bp != sequence_[i]) {
			last_bp = sequence_[i];
			runs += opt_repeat_clip;
		}
		if (total_sequence >= runs) {
			last = i;
		}
	}
	if (last == static_cast<size_t>(-1)) {	// didn't find anything
		return;
	}
	// now find first position that has a repeat of opt_repeat_clip length
	size_t run(0);
	const size_t max_run(size_t(floor(opt_repeat_clip)));
	for (last_bp = 0; ; ++last) {
		if (last_bp != sequence_[last]) {
			last_bp = sequence_[last];
			run = 1;
		} else if (++run == max_run) {
			break;
		}
	}
	quality_stop = last - max_run + 1;
}

// find longest stretch of sequence using stricter rules;
// return if we changed the quality start & stop

bool Read::find_strict_window(const std::pair<size_t, size_t> &x, int &best_score) {
	if (x.second < opt_minimum_clip) {
		return 0;
	}
	bool changed(0);
	std::vector<unsigned char> window(20, 0);
	int total(0);
	size_t i(x.first);
	const size_t end(x.second);
	while (i < end) {
		// find next spot total is >= opt_quality_cutoff on average
		for (; total < 20 * opt_quality_cutoff && i < end; ++i) {
			total += quality[i] - window[i % 20];
			window[i % 20] = quality[i];
		}
		// no window satisfied the limit, so quality clip is empty
		if (i == end && total < 20 * opt_quality_cutoff) {
			break;
		}
		// be careful when subtracting from unsigned values not to wrap
		size_t start = (i < 20) ? 0 : (i - 20);
		// need to check for this, as, if the total was reached in less than
		// 20 spaces, i - 20 will be before the start of the range
		if (start < x.first) {
			start = x.first;
		}
		int run_total(total);
		// get rid of any low quality bases at the leading edge
		for (; quality[start] < opt_quality_cutoff; ++start) {
			run_total -= quality[start];
		}
		// now find end of window
		for (; total >= 20 * opt_quality_cutoff && i < end; ++i) {
			total += quality[i] - window[i % 20];
			window[i % 20] = quality[i];
		}
		size_t stop(i);
		// get rid of any low quality bases at the trailing edge
		for (--stop; quality[stop] < opt_quality_cutoff; --stop) {
			run_total -= quality[stop];
		}
		++stop;
		if (stop <= opt_minimum_clip) {
			continue;
		}
		for (; start < opt_minimum_clip; ++start) {
			run_total -= quality[start];
		}
		if (best_score < run_total && opt_base_cutoff != 0) {
			std::vector<size_t> count(256, 0);
			for (size_t j(start); j != stop; ++j) {
				++count[sequence_[j]];
			}
			const size_t cutoff(size_t(ceil((stop - start) * opt_base_cutoff)));
			std::vector<size_t>::const_iterator b(count.begin());
			const std::vector<size_t>::const_iterator end_b(count.end());
			for (; b != end_b && *b < cutoff; ++b) { }
			if (b != end_b) {
				run_total = 0;
			}
		}
		if (best_score < run_total) {
			best_score = run_total;
			changed = 1;
			quality_start = start;
			quality_stop = stop;
		}
	}
	return changed;
}

// set vector and quality endpoints using strict standards for quality
// windows, and choose vector endpoints based on best quality window,
// rather than just number of quality basepairs

void Read::set_strict_endpoints() {
	// check for any non-vector at all
	if (vectors.empty()) {
		return;
	}
	int best_score(0);
	std::vector<std::pair<size_t, size_t> >::const_iterator a(vectors.begin());
	const std::vector<std::pair<size_t, size_t> >::const_iterator end_a(vectors.end());
	for (; a != end_a; ++a) {
		if (find_strict_window(*a, best_score)) {
			vector_start = a->first;
			vector_stop = a->second;
		}
	}
	vectors.clear();
}

// prior to setting vector endpoints, the possible endpoints must be found

void Read::record_vectors() {
	const char * const vector = opt_N_is_vector ? "NX" : "X";
	std::string::size_type j(sequence_.find_first_not_of(vector, 0));
	while (j != std::string::npos) {
		const std::string::size_type i(sequence_.find_first_of(vector, j + 1));
		if (i == std::string::npos) {
			vectors.push_back(std::make_pair(j, size()));
			break;
		}
		vectors.push_back(std::make_pair(j, i));
		j = sequence_.find_first_not_of(vector, i + 1);
	}
}

// count phred20s in the quality region; non-ACGT basepairs are ignored

void Read::count_phreds() {
	init_phred_counting();
	size_t i;
	for (phred_count = 0, i = quality_start; i != quality_stop; ++i) {
		if (quality[i] >= 20 && (opt_all_p20 || good_base[(int)sequence_[i]])) {
			++phred_count;
		}
	}
}

// check to see if read's sequence and quality match up

void Read::consistency_check(bool opt_warnings) {
	if (quality.empty()) {
		if (opt_warnings) {
			fprintf(stderr, "Warning: sequence with no quality: %s\n", name().c_str());
		}
		set_quality(opt_quality_cutoff);
	} else if (sequence_.size() != quality.size()) {
		if (opt_warnings) {
			fprintf(stderr, "Warning: sequence and quality of different lengths (%lu vs %lu): %s\n", sequence_.size(), quality.size(), name().c_str());
		}
		if (sequence_.size() < quality.size()) {
			quality.resize(sequence_.size());
		} else {
			set_quality(opt_quality_cutoff);
		}
	}
}

/* print read header, followed by sequence */

void Read::print_sequence(FILE *fp) const {
	size_t i, end;
	if (!get_output_endpoints(i, end) || !print_header(fp, i, end)) {
		return;
	}
	if (opt_line_length) {
		size_t next = i + opt_line_length;
		for (; next < end; i = next, next += opt_line_length) {
			fprintf(fp, "%s\n", sequence_.substr(i, opt_line_length).c_str());
		}
	}
	/* last bit might not be a full opt_line_length */
	fprintf(fp, "%s\n", sequence_.substr(i, end - i).c_str());
}

/*
 * print read header and read info, followed by quality; max_qual is used
 * to force qualities into a given range, since some programs (like phrap)
 * aren't happy with even 0-255
 */

void Read::print_quality(FILE * const fp, const unsigned char max_qual) const {
	size_t i, end;
	if (!get_output_endpoints(i, end) || !print_header(fp, i, end)) {
		return;
	}
	std::vector<unsigned char>::const_iterator a(quality.begin() + i);
	const std::vector<unsigned char>::const_iterator end_a(quality.begin() + end);
	while (a != end_a) {
		if (max_qual < *a) {
			fprintf(fp, "%d", max_qual);
		} else {
			fprintf(fp, "%d", *a);
		}
		size_t j(1);
		const size_t end_j(opt_line_length ? opt_line_length : end);
		for (++a; j < end_j && a != end_a; ++j, ++a) {
			if (max_qual < *a) {
				fprintf(fp, " %d", max_qual);
			} else {
				fprintf(fp, " %d", *a);
			}
		}
		fprintf(fp, "\n");
	}
}

/*
 * replace sequence with an 'X' for any position with a phred value lower than
 * the cutoff; make sure to do this after clipping to vector
 */

void Read::mask_by_phred(size_t cutoff) {
	size_t i;
	for (i = 0; i != size(); ++i) {
		if (quality[i] < cutoff) {
			sequence_[i] = 'X';
		}
	}
}

Read Read::operator=(const Read &a) {
	header = a.header;
	quality_start = a.quality_start;
	quality_stop = a.quality_stop;
	vector_start = a.vector_start;
	vector_stop = a.vector_stop;
	phred_count = a.phred_count;
	sequence_ = a.sequence_;
	quality = a.quality;
	vectors = a.vectors;
	return *this;
}

Read Read::operator=(const Read *a) {
	sequence_ = a->sequence_;
	quality = a->quality;
	vectors = a->vectors;
	header = a->header;
	quality_start = a->quality_start;
	quality_stop = a->quality_stop;
	vector_start = a->vector_start;
	vector_stop = a->vector_stop;
	phred_count = a->phred_count;
	return *this;
}

void Read::clip_linker() {
	if (opt_linker.is_match(sequence_)) {
		// trim linker off sequence, quality, and vector list
		const size_t n(opt_linker[0].rm_so);
		sequence_.resize(n);
		quality.resize(n);
		std::vector<std::pair<size_t, size_t> >::iterator a(vectors.begin());
		const std::vector<std::pair<size_t, size_t> >::iterator end_a(vectors.end());
		for (; a != end_a && a->second <= n; ++a) { }
		if (a != end_a) {
			if (a->first < n) {
				a->second = n;
				++a;
			}
			vectors.erase(a, end_a);
		}
	}
}

void Read::add_quality(const std::string &line, const bool opt_warnings) {
	char *t;
	const char *s(line.c_str());
	long i(strtol(s, &t, 10));
	while (s != t) {
		if (i > UCHAR_MAX) {
			quality.push_back(UCHAR_MAX);
		} else if (i < 0) {
			quality.push_back(0);
		} else {
			quality.push_back(i);
		}
		i = strtol(s = t, &t, 10);
	}
	if (opt_strip_trailing_zero_qual && quality.size() == sequence_.size() + 1 && quality.back() == 0) {
		quality.pop_back();
	}
	consistency_check(opt_warnings);
	clip_linker();
	if (opt_strict_quality) {
		set_strict_endpoints();
	} else {
		set_vector_endpoints();
		if (opt_clip_quality) {
			set_quality_endpoints();
		} else {
			quality_start = vector_start;
			quality_stop = vector_stop;
		}
	}
	count_phreds();
}

void Read::add_quality_fastq(const std::string &line, const bool opt_warnings) {
	size_t i(0);
	const size_t end_i(line.size());
	for (; i != end_i; ++i) {
		const int x(line[i] - 33);		// convert to int
		if (x > UCHAR_MAX) {
			quality.push_back(UCHAR_MAX);
		} else if (x < 0) {
			quality.push_back(0);
		} else {
			quality.push_back(x);
		}
	}
	if (opt_strip_trailing_zero_qual && quality.size() == sequence_.size() + 1 && quality.back() == 0) {
		quality.pop_back();
	}
	consistency_check(opt_warnings);
	clip_linker();
	if (opt_strict_quality) {
		set_strict_endpoints();
	} else {
		set_vector_endpoints();
		if (opt_clip_quality) {
			set_quality_endpoints();
		} else {
			quality_start = vector_start;
			quality_stop = vector_stop;
		}
	}
	count_phreds();
}

// set quality to fixed value across entire length

void Read::set_quality(unsigned char x) {
	quality.clear();
	quality.resize(size(), x);
	clip_linker();
	set_vector_endpoints();
	quality_start = vector_start;
	if (opt_clip_quality && (vector_stop < opt_minimum_clip || x < opt_quality_cutoff)) {
		quality_stop = quality_start;
	} else {
		quality_stop = vector_stop;
	}
}

// return a subseq of this read, [start, stop)

Read Read::subseq(size_t start, size_t stop) const {
	Read a;
	a.header = '>' + name() + '_' + itoa(start + 1) + ' ' + itoa(stop - start);
	a.sequence_ = sequence_.substr(start, stop - start);
	a.quality.assign(quality.begin() + start, quality.begin() + stop);
	if (vector_start < stop && start < vector_stop) {
		a.vector_start = vector_start > start ? vector_start - start : 0;
		a.vector_stop = std::min(vector_stop - start, a.size());
	} else {
		a.vector_start = a.vector_stop = 0;
	}
	if (quality_start < stop && start < quality_stop) {
		a.quality_start = quality_start > start ? quality_start - start : 0;
		a.quality_stop = std::min(quality_stop - start, a.size());
	} else {
		a.quality_start = a.quality_stop = 0;
	}
	return a;
}

// count number of qualities above cutoff in given region
// (if no quality is available, just return sequence length);
// start is inclusive, end is not

size_t Read::count_quality(const std::pair<size_t, size_t> &a) const {
	if (opt_base_cutoff == 0) {
		size_t n(0);
		for (size_t i(a.first); i != a.second; ++i) {
			if (quality[i] >= opt_quality_cutoff) {
				++n;
			}
		}
		return n;
	} else {
		std::vector<size_t> count(256, 0);
		size_t n(0);
		for (size_t i(a.first); i != a.second; ++i) {
			if (quality[i] >= opt_quality_cutoff) {
				++n;
			}
			++count[sequence_[i]];
		}
		const size_t cutoff(size_t(ceil((a.second - a.first) * opt_base_cutoff)));
		std::vector<size_t>::const_iterator b(count.begin());
		const std::vector<size_t>::const_iterator end_b(count.end());
		for (; b != end_b && *b < cutoff; ++b) { }
		return b == end_b ? n : 0;
	}
}

size_t Read::count_masked() const {
	size_t i, n;
	for (i = n = 0; i != sequence_.size(); ++i) {
		if (sequence_[i] == 'X') {
			++n;
		}
	}
	return n;
}
