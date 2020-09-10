#ifndef _KMER_LOOKUP_INFO_H
#define _KMER_LOOKUP_INFO_H

#include "hash.h"	// hash
#include "hash_read_hits.h"	// hash_read_hits, hash_read_hits::read_type
#include <new>		// new[]
#include <string.h>	// memcpy()
#include <string>	// string
#include <sys/types.h>	// size_t

// we're using hash_read_hits::read_type instead of uint32_t so in case we ever
// feel the need to bump that up to uint64_t we only have to change one typedef

// class with all the information to do kmer -> read lookups
class KmerLookupInfo {
    private:
	// read name stuff
	hash_read_hits::read_type count;	// number of reads
	hash_read_hits::read_type data_size;
	hash_read_hits::read_type *list;	// offsets into data
	char *data;
    public:
	hash_read_hits kmer_hash;
	// total_length doesn't include ending nulls
	explicit KmerLookupInfo() : count(0), data_size(0), list(0), data(0) { };
	explicit KmerLookupInfo(const size_t total_reads, const size_t total_name_size, hash &mer_list, const double hash_usage = 0.9) : count(0), data_size(0), list(new hash_read_hits::read_type[total_reads]), data(new char[total_name_size + total_reads]), kmer_hash(mer_list, hash_usage) { };
	~KmerLookupInfo() {
		if (list) {
			delete[] list;
		}
		if (data) {
			delete[] data;
		}
	};
	void add_read_name(const std::string &name) {
		list[count] = data_size;
		// +1 to size() to include ending null
		memcpy(&data[data_size], name.c_str(), name.size() + 1);
		data_size += name.size() + 1;
		++count;
	}
	void save(int) const;
	void restore(int);
};

#endif // !_KMER_LOOKUP_INFO_H
