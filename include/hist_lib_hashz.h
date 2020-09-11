#ifndef _HIST_LIB_HASHZ_H
#define _HIST_LIB_HASHZ_H

#include "hashz.h"	// hashz, hashz::key_type, hashz::value_type
#include "pattern.h"	// Pattern
#include <map>		// map<>
#include <string>	// string
#include <sys/types.h>	// size_t

class Read;

extern Pattern opt_include;
extern bool opt_feedback;
extern bool opt_mask_lowercase;
extern hashz::value_type opt_repeat_threshold;
extern hashz::value_type opt_repeat_threshold_upper;
extern int opt_phred20_anchor;
extern int opt_repeat_coverage;
extern size_t opt_skip_size;
extern std::map<std::string, bool> opt_exclude;

extern std::string convert_key(const hashz::key_type &);
// must be called before any of the following can be used
extern void init_mer_constants(unsigned long);
extern bool add_sequence_mers(std::list<Read>::const_iterator, std::list<Read>::const_iterator, hashz &);
extern bool add_sequence_mers(std::list<Read>::const_iterator, std::list<Read>::const_iterator, hashz &, const std::map<std::string, hashz::offset_type> &);
extern void count_kmers(const Read &, const hashz &, size_t &, size_t &, size_t &);
extern void screen_repeats(Read &, const hashz &);
extern unsigned long count_unique_phreds(const std::list<Read> &, const hashz &, unsigned long * = NULL);
extern void reverse_key(const hashz::key_type &, hashz::key_type &);

#endif // !_HIST_LIB_HASHZ_H
