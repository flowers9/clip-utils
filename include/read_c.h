#ifndef _READ_C_H
#define _READ_C_H

#include "get_name.h"	/* get_name() */
#include "pattern.h"	/* Pattern */
#include <limits.h>	/* UCHAR_MAX */
#include <map>		/* map<> */
#include <stdio.h>	/* FILE, stdout */
#include <string>	/* string */
#include <sys/types.h>	/* size_t */
#include <utility>	/* pair<> */
#include <vector>	/* vector<> */

extern Pattern opt_linker;
extern bool opt_N_is_vector;
extern bool opt_add_range;
extern bool opt_all_p20;
extern bool opt_clip_quality;
extern bool opt_clip_vector;
extern bool opt_strip_trailing_zero_qual;
extern int opt_quality_cutoff;
extern size_t opt_line_length;
extern size_t opt_minimum_clip;
extern std::map<std::string, std::string> read_name_translation;

class Read {
    private:
	size_t count_quality(const std::pair<size_t, size_t> &) const;
	void record_vectors(void);
	void set_vector_endpoints(void);
	bool get_output_endpoints(size_t &, size_t &) const;
	bool print_header(FILE *fp) const;
	void clip_linker(void);
    private:
	std::vector<unsigned char> quality;
	std::vector<std::pair<size_t, size_t> > vectors;
    public:
	std::string header;
	/* sequence sans low quality tails; start is inclusive, stop is not */
	size_t quality_start, quality_stop;
	/* longest continuous sequence without vector */
	size_t vector_start, vector_stop;
	unsigned int phred_count;	/* phred20s in quality region */
	~Read(void) { }
	Read operator=(const Read &);
	Read operator=(const Read *);
	std::string name(void) const {
		return get_name(header);
	}
	void add_quality(const std::string &, bool);
	Read subseq(size_t, size_t) const;
	void next_good_sequence(size_t &__s) const {
		for (; __s != size() && !is_good_basepair(__s); ++__s) { }
	}
	void print_sequence(FILE *fp = stdout) const;
	void print_quality(FILE *fp = stdout, unsigned char = UCHAR_MAX) const;
	void mask_by_phred(size_t);
	size_t count_masked(void) const;
    private:
	unsigned char get_sequence_raw(size_t __x) const {
		return (sequence[__x >> 2] >> (2 * (__x & 3))) & 3;
	}
	unsigned char get_quality_raw(size_t __x) const {
		return (quality[__x >> 2] >> (2 * (__x & 3))) & 3;
	}
	void set_sequence_raw(size_t __x, char __c) {
		unsigned char &__i = sequence[__x >> 2];
		__x = 2 * (__x & 3);
		__i = (__i & (~(3 << __x))) | (__c << __x);
	}
	void set_quality_raw(size_t __x, unsigned char __c) {
		unsigned char &__i = quality[__x >> 2];
		__x = 2 * (__x & 3);
		__i = (__i & (~(3 << __x))) | (__c << __x);
	}
	void set_quality_endpoints(const std::vector<unsigned char> &);
	void count_phreds(const std::vector<unsigned char> &);
    private:
	bool qual_set;
	size_t _size;
	std::vector<unsigned char> sequence;
    public:
	Read(void) : quality_start(0), quality_stop(0), vector_start(0), vector_stop(0), phred_count(0), qual_set(0), _size(0) { }
	explicit Read(const std::string &__s) : header(__s), quality_start(0), quality_stop(0), vector_start(0), vector_stop(0), phred_count(0), qual_set(0), _size(0) { }
	explicit Read(const std::string &__s, const std:string &__t) : header(__s), quality_start(0), quality_stop(0), vector_start(0), vector_stop(0), phred_count(0), qual_set(0), _size(0) {
		add_sequence(__t);
	}
	size_t size(void) const {
		return _size;
	}
	char get_sequence(size_t __x) const {
		switch (get_quality_raw(__x)) {
		    case 0:
			return 'X';
		    case 1:
			return 'N';
		    default:
			return "ACGT"[get_sequence_raw(__x)];
		}
	}
	unsigned char get_quality(size_t __x) const {
		unsigned char __i = get_quality_raw(__x);
		if (__i < 2) {
			__i = get_sequence_raw(__x);
		}
		switch (__i) {
		    case 0:
			return 0;
		    case 1:
			return 1;
		    case 2:
			return opt_quality_cutoff / 2;
		    default:
			return opt_quality_cutoff;
		}
	}
	int get_seq(size_t __x) const {	// returns 0-3, -1 for bad sequence
		return get_quality_raw(__x) < 2 ? -1 : get_sequence_raw(__x);
	}
	void set_sequence(size_t __x, char __c) {
		switch (__c) {
		    case 'A':
		    case 'a':
			set_sequence_raw(__x, 0);
			break;
		    case 'C':
		    case 'c':
			set_sequence_raw(__x, 1);
			break;
		    case 'G':
		    case 'g':
			set_sequence_raw(__x, 2);
			break;
		    case 'T':
		    case 't':
			set_sequence_raw(__x, 3);
			break;
		    case 'N':
		    case 'n':
			set_sequence_raw(__x, 0);
			set_quality_raw(__x, 1);
			return;
		    default:	// everything else becomes an X
			set_sequence_raw(__x, 0);
			set_quality_raw(__x, 0);
			return;
		}
		set_quality_raw(__x, 3);
	}
	void set_quality(size_t __x, unsigned char __c) {
		if (get_quality_raw(__x) > 1) {
			set_quality_raw(__x, __c < opt_quality_cutoff ? 2 : 3);
		} else if (__c > 1) {
			set_sequence_raw(__x, __c < opt_quality_cutoff ? 2 : 3);
		} else {
			set_sequence_raw(__x, __c);
		}
	}
	void add_sequence(const std::string &);
	void set_quality(unsigned char __c) {
		set_vector_endpoints();
		quality_start = vector_start;
		if (opt_clip_quality && (vector_stop < opt_minimum_clip || __c < opt_quality_cutoff)) {
			quality_stop = quality_start;
		} else {
			quality_stop = vector_stop;
		}
		qual_set = 1;
	}
	bool is_good_basepair(size_t __s) const {
		return get_quality_raw(__s) > 1;
	}
	bool is_high_quality(size_t __s) const {
		return get_quality_raw(__s) == 3;
	}
	bool has_quality(void) const {
		return qual_set;
	}

};

#endif /* !_READ_C_H */
