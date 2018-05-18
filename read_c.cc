#include "itoa.h"	// itoa()
#include "pattern.h"	// Pattern
#include "read.h"
#include <algorithm>	// min()
#include <limits.h>	// UCHAR_MAX
#include <map>		// map<>
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
bool opt_strip_trailing_zero_qual = 0;
int opt_quality_cutoff = 20;	// for setting vector and quality clipping
size_t opt_line_length = 50;	// wrapping point for lines
size_t opt_minimum_clip = 0;
std::map<std::string, std::string> read_name_translation;

// find the largest continuous sequence of non-vector ('X') and set the
// vector start and stop points for the read

void Read::set_vector_endpoints() {
	// check for any vector at all
	if (vectors.empty()) {
		vector_start = 0;
		vector_stop = size();
		return;
	}
	std::vector<std::pair<size_t, size_t> >::const_iterator a(vectors.begin());
	const std::vector<std::pair<size_t, size_t> >::const_iterator end_a(vectors.end());
	std::vector<std::pair<size_t, size_t> >::const_iterator best(a);
	size_t best_count(count_quality(*a));
	for (++a; a != end_a; ++a) {
		const size_t n(count_quality(*a));
		if (best_count < n || (best_count == n && best->second - best->first < a->second - a->first)) {
			best_count = n;
			best = a;
		}
	}
	vector_start = best->first;
	vector_stop = best->second;
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

// print a sequence header, adding clip ranges if specified;
// returns if clipping region is not empty

bool Read::print_header(FILE * const fp) const {
	std::string x(name());
	const size_t i(x.size() + 1);
	const std::map<std::string, std::string>::const_iterator a(read_name_translation.find(x));
	if (!opt_add_range) {
		if (a != read_name_translation.end()) {
			fprintf(fp, ">%s%s\n", a->second.c_str(), header.substr(i).c_str());
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
			fprintf(fp, ">%s %d %ld%s\n", x.c_str(), 1, quality_stop - quality_start, header.substr(i).c_str());
		} else if (opt_clip_vector) {
			fprintf(fp, ">%s %ld %ld%s\n", x.c_str(), quality_start + 1 - vector_start, quality_stop - vector_start, header.substr(i).c_str());
		} else {
			fprintf(fp, ">%s %ld %ld%s\n", x.c_str(), quality_start + 1, quality_stop, header.substr(i).c_str());
		}
	}
	return 1;
}

/* set the quality start and stop points for the read */

void Read::set_quality_endpoints(const std::vector<unsigned char> &quality_in) {
	if (vector_stop < opt_minimum_clip) {
		quality_start = quality_stop = vector_start;
		return;
	}
	size_t i = vector_start;
	size_t end = vector_stop;
	std::vector<unsigned char> window;
	/* initialize 20-space window */
	window.resize(20, 0);
	int total = 0;
	/* find first spot window total is >= opt_quality_cutoff on average */
	for (; total < 20 * opt_quality_cutoff && i < end; ++i) {
		total += quality_in[i] - window[i % 20];
		window[i % 20] = quality_in[i];
	}
	/* no window satisfied the limit, so quality clip is empty */
	if (i == end && total < 20 * opt_quality_cutoff) {
		quality_start = quality_stop = vector_start;
		return;
	}
	/* be careful when subtracting from unsigned values not to wrap */
	if (i < 20) {
		quality_start = 0;
	} else {
		quality_start = i - 20;
	}
	/*
	 * need to check for this, as, if the total was reached in less than
	 * 20 spaces, i - 20 will be before vector_start
	 */
	if (quality_start < vector_start) {
		quality_start = vector_start;
	}
	if (quality_start < opt_minimum_clip) {
		quality_start = opt_minimum_clip;
	}
	/* re-initialize window */
	window.assign(20, 0);
	total = 0;
	i = end - 1;
	/* find last spot window total is >= opt_quality_cutoff on overage */
	/*
	 * don't have to worry about running off the beginning, as we already
	 * know that there exists a window that will satisfy the minimum total
	 */
	for (; total < 20 * opt_quality_cutoff; --i) {
		total += quality_in[i] - window[i % 20];
		window[i % 20] = quality_in[i];
	}
	quality_stop = i + 21;
	/*
	 * again, if the total is reached in less than 20 spaces,
	 * i + 21 could be beyond vector_stop
	 */
	if (quality_stop > vector_stop) {
		quality_stop = vector_stop;
	}
	/* in case the minimum clip ended up beyond quality stop */
	if (quality_stop < quality_start) {
		quality_stop = quality_start;
	}
}

// prior to setting vector endpoints, the possible endpoints must be found

void Read::record_vectors() {
	const unsigned char cutoff = opt_N_is_vector ? 1 : 0;
	size_t i;
	for (i = 0; i != size() && get_quality_raw(i) > cutoff; ++i) { }
	if (i == size()) {
		return;
	}
	vectors.push_back(std::make_pair(0, i));
	for (;;) {
		size_t j = i + 1;
		for (; j != size() && get_quality_raw(j) <= cutoff; ++j) { }
		if (j == size()) {
			break;
		}
		for (i = j + 1; i != size() && get_quality_raw(i) > cutoff; ++i) { }
		vectors.push_back(std::make_pair(j, i));
		if (i == size()) {
			break;
		}
	}
}

/* count phred20s in the quality region; non-ACGT basepairs may be ignored */

void Read::count_phreds(const std::vector<unsigned char> &quality_in) {
	size_t i;
	for (phred_count = 0, i = quality_start; i != quality_stop; ++i) {
		if (quality_in[i] >= 20 && (opt_all_p20 || get_quality_raw(i) > 1)) {
			++phred_count;
		}
	}
}

/* print read header, followed by sequence */

void Read::print_sequence(FILE *fp) const {
	size_t i, end;
	if (!get_output_endpoints(i, end) || !print_header(fp)) {
		return;
	}
	if (opt_line_length) {
		size_t next = i + opt_line_length;
		for (; next < end; next += opt_line_length) {
			for (; i != next; ++i) {
				fputc(get_sequence(i), fp);
			}
			fputc('\n', fp);
		}
	}
	if (i != end) {
		/* last bit might not be a full opt_line_length */
		for (; i != end; ++i) {
			fputc(get_sequence(i), fp);
		}
		fputc('\n', fp);
	}
}

/*
 * print read header and read info, followed by quality; max_qual is used
 * to force qualities into a given range, since some programs (like phrap)
 * aren't happy with even 0-255
 */

void Read::print_quality(FILE *fp, unsigned char max_qual) const {
	size_t i, end;
	if (!get_output_endpoints(i, end) || !print_header(fp)) {
		return;
	}
	if (opt_line_length) {
		size_t next = i + opt_line_length;
		for (; next < end; next += opt_line_length) {
			for (; i != next; ++i) {
				/* the trailing space on the line is expected */
				fprintf(fp, "%d ", get_quality(i));
			}
			fputc('\n', fp);
		}
	}
	if (i != end) {
		/* last bit might not be a full opt_line_length */
		for (; i != end; ++i) {
			fprintf(fp, "%d ", get_quality(i));
		}
		fputc('\n', fp);
	}
}

/*
 * replace sequence with an 'X' for any position with a phred value lower than
 * the cutoff; make sure to do this after clipping to vector
 */

void Read::mask_by_phred(size_t cutoff) {
	size_t i;
	/* just in case the sizes are different, use the smaller one */
	size_t end = sequence.size() < quality.size() ? sequence.size() : quality.size();
	for (i = 0; i != end; ++i) {
		if (quality[i] < cutoff) {
			sequence[i] = 'X';
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
	sequence = a.sequence;
	quality = a.quality;
	vectors = a.vectors;
	_size = a._size;
	return *this;
}

Read Read::operator=(const Read *a) {
	header = a->header;
	quality_start = a->quality_start;
	quality_stop = a->quality_stop;
	vector_start = a->vector_start;
	vector_stop = a->vector_stop;
	phred_count = a->phred_count;
	sequence = a->sequence;
	quality = a->quality;
	vectors = a->vectors;
	_size = a->_size;
	return *this;
}

void Read::add_sequence(const std::string &line) {
	_size = line.size();
	sequence.resize((size() + 3) / 4);
	quality.resize(sequence.size());
	size_t i;
	for (i = 0; i != line.size(); ++i) {
		set_sequence(i, line[i]);
	}
	if (opt_clip_vector) {
		record_vectors();
	}
}

void Read::add_quality(const std::string &line, bool opt_warnings) {
	std::vector<unsigned char> quality_tmp;
	quality_tmp.reserve(size());
	size_t n = 0;
	char *t;
	const char *s = line.c_str();
	long i = strtol(s, &t, 10);
	while (s != t && n != size()) {
		quality_tmp[n] = i;
		if (get_quality_raw(n) < 2) {
			if (i < 1) {		// default value of 0
			} else if (i == 1) {
				set_sequence_raw(n, 1);
			} else if (i < opt_quality_cutoff) {
				set_sequence_raw(n, 2);
			} else {
				set_sequence_raw(n, 3);
			}
		} else if (i < opt_quality_cutoff && get_quality_raw(n) == 3) {
			set_quality_raw(n, 2);
		}
		++n;
		i = strtol(s = t, &t, 10);
	}
	if (opt_warnings && (s != t || n != size())) {
		while (s != t) {	// finish counting n
			++n;
			strtol(s = t, &t, 10);
		}
		if (!opt_strip_trailing_zero_qual || n != size() + 1) {
			fprintf(stderr, "Warning: sequence and quality of different lengths (%lu vs %lu): %s\n", size(), n, name().c_str());
			for (; n < size(); ++n) {	// < required
				quality_tmp[n] = 0;
			}
		}
	}
	set_vector_endpoints();
	if (opt_clip_quality) {
		set_quality_endpoints(quality_tmp);
	} else {
		quality_start = vector_start;
		quality_stop = vector_stop;
	}
	count_phreds(quality_tmp);
	qual_set = 1;
}

// return a subseq of this read, [start, stop)

Read Read::subseq(size_t start, size_t stop) const {
	Read a;
	a.header = '>' + name() + '_' + itoa(start + 1) + ' ' + itoa(stop - start);
	size_t i;
	for (i = start; i != stop; ++i) {
		a.set_sequence_raw(i - start, get_sequence_raw(i));
		a.set_quality_raw(i - start, get_quality_raw(i));
	}
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

/*
 * count number of qualities above cutoff in given region;
 * start is inclusive, end is not
 */

size_t Read::count_quality(const std::pair<size_t, size_t> &a) const {
	size_t n = 0;
	size_t i = a.first;
	for (; i != a.second; ++i) {
		if (get_quality_raw(i) == 3) {
			++n;
		}
	}
	return n;
}

size_t Read::count_masked() const {
	size_t i, n;
	for (i = n = 0; i != _size; ++i) {
		if (get_quality_raw(i) == 0) {
			++n;
		}
	}
	return n;
}
