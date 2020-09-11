#ifndef _HIST_LIB_HASH_H
#define _HIST_LIB_HASH_H

#include "hash.h"	// hash
#include "hash_read_hits.h"	// hash_read_hits::read_type
#include "pattern.h"	// Pattern
#include <map>		// map<>
#include <string>	// string
#include <sys/types.h>	// size_t

class KmerLookupInfo;
class Read;

extern Pattern opt_include;
extern bool opt_feedback;
extern bool opt_mask_lowercase;
extern bool opt_reverse_mask;
extern hash::value_type opt_repeat_threshold;
extern hash::value_type opt_repeat_threshold_upper;
extern int opt_phred20_anchor;
extern size_t opt_mer_length;
extern size_t opt_repeat_coverage;
extern size_t opt_skip_size;
extern std::map<std::string, bool> opt_exclude;

extern std::string convert_key(hash::key_type);
extern std::string convert_key_hp(hash::key_type, int = 0);
extern void init_mer_constants(void);	// must be called before any of the
					// following can be used
extern void clear_mer_list(hash &);
extern void print_final_input_feedback(const hash &);
extern bool add_sequence_mers(std::list<Read>::const_iterator, const std::list<Read>::const_iterator, hash &, size_t);
extern void add_sequence_mers_index(std::list<Read>::const_iterator, const std::list<Read>::const_iterator, KmerLookupInfo &, size_t, size_t);
extern bool add_sequence_mers_hp(std::list<Read>::const_iterator, const std::list<Read>::const_iterator, hash &, size_t);
extern bool add_sequence_mers(std::list<Read>::const_iterator, std::list<Read>::const_iterator, hash &, const std::map<std::string, hash::offset_type> &, size_t);
extern void count_read_hits(const std::string &, const KmerLookupInfo &, std::map<hash_read_hits::read_type, int> &, hash_read_hits::value_type);
extern void count_kmers(const Read &, const hash &, size_t &, size_t &, size_t &);
extern void screen_repeats(Read &, const hash &);
extern unsigned long count_unique_phreds(const std::list<Read> &, const hash &, unsigned long * = NULL);
extern hash::key_type reverse_key(hash::key_type);

#endif // !_HIST_LIB_HASH_H
