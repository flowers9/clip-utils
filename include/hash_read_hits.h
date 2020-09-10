#ifndef _HASH_READ_HITS_H
#define _HASH_READ_HITS_H

// This hash is designed to hold a list of reads for each kmer; it's missing
// some of hash features (such as handling OOM conditions) as all memory
// usage is pre-allocated - it's designed to be run after a counting pass
// to determine how many read hits there will be for each kmer.

// Currently it generates a 90% full hash; possibly that may become a run-time
// variable in the future.

#include "hash.h"	// hash
#include <limits.h>	// UCHAR_MAX, ULONG_MAX
#include <map>		// map<>
#include <stdint.h>	// uint32_t, uint64_t
#include <string>	// string
#include <utility>	// pair<>

// a non-palindrome random pattern with one of the two highest bits set,
// and such that the complement is of lower value (so it only collides
// when using 32-mers, and not even then if you only store the "smaller"
// key of key/key-complement); if the given n-mer is found, the repeat
// count for it will probably be low and it will be skipped by the iterator;
// the current pattern corresponds to TTAGCTGGGAAGGCTTATTGTGTCGTCGGATG

// choose size of constant by what will hold it
#if ULONG_MAX >= 0xffffffffffffffffULL
#define INVALID_KEY 0xf27a829f3eedb68eUL
#else
#define INVALID_KEY 0xf27a829f3eedb68eULL
#endif

class hash_read_hits {
    public:	// type declarations
	typedef uint64_t key_type;
	typedef unsigned char small_value_type;
	typedef size_t value_type;
	typedef size_t offset_type;
	typedef uint32_t read_type;		// no more than 2^32 reads
	typedef size_t read_offset_type;
	enum { max_small_value = UCHAR_MAX };
    private:
	offset_type used_elements;
	offset_type modulus;
	offset_type collision_modulus;
	offset_type read_list_size;
	key_type *key_list;
	small_value_type *value_list;
	read_offset_type *read_offset_list;	// offset into read_list for given kmer
	read_type *read_list;			// list of read hits for each kmer
	std::map<offset_type, value_type> value_map;	// for overflow
    private:
	std::string boilerplate(void) const;
	offset_type insert_key(offset_type, key_type);
	offset_type insert_offset(key_type);
	offset_type find_offset(key_type) const;
    public:
	explicit hash_read_hits() : used_elements(0), modulus(0), collision_modulus(0), read_list_size(0), key_list(0), value_list(0), read_offset_list(0), read_list(0) { }
	explicit hash_read_hits(hash &mer_list, double hash_usage = 0.9);
	~hash_read_hits() {
		delete[] key_list;
		delete[] value_list;
		delete[] read_offset_list;
		delete[] read_list;
	}
	void add_read(key_type, read_type);
	void get_reads(key_type, std::map<read_type, int> &reads) const;
	// the -1s are to allow for having to keep at least one INVALID_KEY
	// in the array for lookup termination purposes
	offset_type size(void) const {
		return used_elements - 1;
	}
	offset_type capacity(void) const {
		return modulus - 1;
	}
	bool empty(void) const {
		return used_elements == 1;
	}
	offset_type overflow_size(void) const {
		return value_map.size();
	}
	void save(int) const;
	void restore(int);
};

#endif // !_HASH_READ_HITS_H
