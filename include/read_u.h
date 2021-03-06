#ifndef _READ_U_H
#define _READ_U_H

#include "get_name.h"	// get_name()
#include "pattern.h"	// Pattern
#include <limits.h>	// UCHAR_MAX
#include <map>		// map<>
#include <stdio.h>	// FILE, stdout
#include <string>	// string
#include <sys/types.h>	// size_t
#include <utility>	// pair<>
#include <vector>	// vector<>

extern Pattern opt_linker;
extern bool opt_N_is_vector;
extern bool opt_add_range;
extern bool opt_all_p20;
extern bool opt_clip_quality;
extern bool opt_clip_vector;
extern bool opt_pacbio;
extern bool opt_strict_quality;
extern bool opt_strip_trailing_zero_qual;
extern double opt_base_cutoff;
extern double opt_repeat_clip;
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
	bool print_header(FILE *fp, size_t, size_t) const;
	void clip_linker(void);
    private:
	std::vector<unsigned char> quality;
	std::vector<std::pair<size_t, size_t> > vectors;
    public:
	std::string header;
	// sequence sans low quality tails; start is inclusive, stop is not
	size_t quality_start, quality_stop;
	// longest continuous sequence without vector
	size_t vector_start, vector_stop;
	unsigned int phred_count;	// phred20s in quality region
	~Read(void) { }
	Read operator=(const Read &);
	Read operator=(const Read *);
	std::string name(void) const {
		return get_name(header);
	}
	void add_quality(const std::string &, bool);
	void add_quality_fastq(const std::string &, bool);
	Read subseq(size_t, size_t) const;
	void next_good_sequence(size_t &__s) const {
		for (; __s != size() && !is_good_basepair(__s); ++__s) { }
	}
	void print_sequence(FILE *fp = stdout) const;
	void print_quality(FILE *fp = stdout, unsigned char = UCHAR_MAX) const;
	void mask_by_phred(size_t);
	size_t count_masked(void) const;
    private:
	bool find_strict_window(const std::pair<size_t, size_t> &, int &);
	void set_strict_endpoints(void);
	void consistency_check(bool);
	void set_quality_endpoints(void);
	void count_phreds(void);
    private:
	std::string sequence;
    public:
	Read(void) : quality_start(0), quality_stop(0), vector_start(0), vector_stop(0), phred_count(0) { }
	explicit Read(const std::string &__s) : header(__s), quality_start(0), quality_stop(0), vector_start(0), vector_stop(0), phred_count(0) { }
	explicit Read(const std::string &__s, const std::string &__t) : header(__s), quality_start(0), quality_stop(0), vector_start(0), vector_stop(0), phred_count(0) {
		add_sequence(__t);
	}
	explicit Read(const std::string &__s, const std::string &__t, const std::string &__u, const bool __b) : header(__s), quality_start(0), quality_stop(0), vector_start(0), vector_stop(0), phred_count(0) {
		add_sequence(__t);
		add_quality_fastq(__u, __b);
	}
	size_t size(void) const {
		return sequence.size();
	}
	const char &get_sequence(size_t __x) const {
		return sequence[__x];
	}
	const unsigned char &get_quality(size_t __x) const {
		return quality[__x];
	}
	int get_seq(size_t __x) const {	// return 0-3, -1 for bad sequence
		switch (sequence[__x]) {
		    case 'A':
		    case 'a':
			return 0;
		    case 'C':
		    case 'c':
			return 1;
		    case 'G':
		    case 'g':
			return 2;
		    case 'T':
		    case 't':
			return 3;
		    default:
			return -1;
		}
	}
	void set_sequence(size_t __x, char __c) {
		sequence[__x] = __c;
	}
	void set_quality(size_t __x, unsigned char __c) {
		quality[__x] = __c;
	}
	void add_sequence(const std::string &__s) {
		sequence = __s;
		if (opt_clip_vector || opt_strict_quality) {
			record_vectors();
		} else {
			vector_start = 0;
			vector_stop = size();
		}
	}
	void set_quality(unsigned char);
	bool is_good_basepair(size_t __s) const {
		const std::string __good_bps("ACGTacgt");
		return __good_bps.find(get_sequence(__s)) != __good_bps.npos;
	}
	bool is_high_quality(size_t __s) const {
		return quality[__s] >= opt_quality_cutoff;
	}
	bool has_quality(void) const {
		return !quality.empty();
	}

};

#endif // !_READ_U_H
