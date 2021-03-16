#ifndef _KMER_LOOKUP_INFO_H
#define _KMER_LOOKUP_INFO_H

#include "hash_read_hits.h"	// hash_read_hits, hash_read_hits::read_type
#include <new>		// new[]
#include <stdint.h>	// uint32_t
#include <string.h>	// memcpy()
#include <string>	// string
#include <sys/types.h>	// size_t

class hash;

// we're using hash_read_hits::read_type instead of uint32_t so in case we ever
// feel the need to bump that up to uint64_t we only have to change one typedef

// class with all the information to do kmer -> read lookups
class KmerLookupInfo {
    private:
	size_t mer_length_;
	// read name stuff
	hash_read_hits::read_type count;	// number of reads
	hash_read_hits::read_type data_size;
	hash_read_hits::read_type *list;	// offsets to given read name in data
	uint32_t *read_kmers_;			// number of kmers in read
	char *data;				// read names
    public:
	hash_read_hits kmer_hash;
	// total_length doesn't include ending nulls
	explicit KmerLookupInfo() : mer_length_(0), count(0), data_size(0), list(0), read_kmers_(0), data(0) { }
	explicit KmerLookupInfo(const size_t mer_length_in, const size_t total_reads, const size_t total_name_size, hash &mer_list, const double hash_usage = 0.9) : mer_length_(mer_length_in), count(0), data_size(0), list(new hash_read_hits::read_type[total_reads]), read_kmers_(new uint32_t[total_reads]), data(new char[total_name_size + total_reads]), kmer_hash(mer_list, hash_usage) { }
	~KmerLookupInfo() {
		if (list) {
			delete[] list;
		}
		if (read_kmers_) {
			delete[] read_kmers_;
		}
		if (data) {
			delete[] data;
		}
	}
	size_t mer_length() const {
		return mer_length_;
	}
	hash_read_hits::read_type read_count() const {
		return count;
	}
	void add_read_name(const std::string &name) {
		list[count] = data_size;
		read_kmers_[count] = 0;
		// +1 to size() to include ending null
		memcpy(&data[data_size], name.c_str(), name.size() + 1);
		data_size += name.size() + 1;
		++count;
	}
	void set_kmer_count(const uint32_t kmer_count) {
		read_kmers_[count - 1] = kmer_count;
	}
	const char *read_name(const hash_read_hits::read_type i) const {
		return &data[list[i]];
	}
	uint32_t read_kmers(const hash_read_hits::read_type i) const {
		return read_kmers_[i];
	}
	void save(int) const;
	void restore(int);
};

#endif // !_KMER_LOOKUP_INFO_H
