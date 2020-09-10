#ifndef _HIST_LIB_HASHN_H
#define _HIST_LIB_HASHN_H

#include "hashn.h"	// hashn
#include "pattern.h"	// Pattern
#include "read.h"	// Read
#include <map>		// map<>
#include <string>	// string
#include <sys/types.h>	// size_t

extern Pattern opt_include;
extern bool opt_feedback;
extern bool opt_mask_lowercase;
extern bool opt_reverse_mask;
extern hashn::value_type opt_repeat_threshold;
extern hashn::value_type opt_repeat_threshold_upper;
extern int opt_phred20_anchor;
extern int opt_repeat_coverage;
extern size_t opt_skip_size;
extern std::map<std::string, bool> opt_exclude;

extern std::string convert_key(const hashn::key_type_base &);
// must be called before any of the following can be used
extern void init_mer_constants(unsigned long);
extern void clear_mer_list(hashn &);
extern void print_final_input_feedback(const hashn &);
extern bool add_sequence_mers(std::list<Read>::const_iterator, std::list<Read>::const_iterator, hashn &, size_t);
extern bool add_sequence_mers(std::list<Read>::const_iterator, std::list<Read>::const_iterator, hashn &, const std::map<std::string, hashn::offset_type> &, size_t);
extern void count_kmers(const Read &, const hashn &, size_t &, size_t &, size_t &);
extern void screen_repeats(Read &, const hashn &);
extern unsigned long count_unique_phreds(const std::list<Read> &, const hashn &, unsigned long * = NULL);
extern void reverse_key(const hashn::key_type_base &, hashn::key_type &);

#endif // !_HIST_LIB_HASHN_H
