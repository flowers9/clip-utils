#ifndef _HASHL_INDEX_H
#define _HASHL_INDEX_H

// this class is for finding kmers among a subset of all kmers in a
// given reference, without taking up gobs of memory but none the less
// being fairly speedy.
//
// In large part, it has much of the same internal structure as hashl,
// except that it only stores the kmer offsets, and those are stored
// sorted by kmer.  Thus, all lookups require an O(log(n)) search, but
// it's a sparse search that can be performed by leaving most of the index
// untouched on disk.
//
// key values are offsets into an internal array; a metadata blob is also stored

#include "hashl.h"
#include "hashl_key_type.h"	// hashl_key_type<>
#include <stdint.h>	// uint64_t
#include <string>	// string
#include <sys/types.h>	// size_t
#include <vector>	// vector<>

class hashl_index {
    public:	// type declarations
	typedef unsigned long data_offset_type;
	typedef uint64_t base_type;
	typedef std::vector<base_type> vector_key_type;
	typedef hashl_key_type<hashl_index> key_type;

	class const_iterator {
	    private:
		const hashl_index &list;
		size_t offset_;
	    public:
		explicit const_iterator(const hashl_index &a, const size_t i) : list(a) {
			offset_ = i < list.size() ? i : list.size();
		}
		const_iterator(const const_iterator &a) : list(a.list), offset_(a.offset_) { }
		~const_iterator() { }
		const size_t &offset() const {
			return offset_;
		}
		void key(key_type &key) const {
			key.copy_in(list.data, list.key_list[offset_]);
		}
		bool operator==(const const_iterator &__a) const {
			return &list == &__a.list && offset_ == __a.offset_;
		}
		bool operator!=(const const_iterator &__a) const {
			return !(*this == __a);
		}
		const_iterator operator++() {
			++offset_;
			return *this;
		}
		const_iterator operator++(int) {
			const const_iterator __tmp(*this);
			++offset_;
			return __tmp;
		}
	};

    protected:
	//std::vector<data_offset_type> key_list;
	std::vector<base_type> data;
	std::vector<char> metadata;
	size_t bit_width;
	size_t word_width;
	int fd;
    protected:
	std::string boilerplate() const;
    public:
	// can only be initialized from an uncompressed file
	explicit hashl_index(int);
	~hashl_index() { }
	size_t size() const {
		return key_list.size();
	}
	bool empty() const {
		return key_list.empty();
	}
	size_t bits() const {
		return bit_width;
	}
	size_t words() const {
		return word_width;
	}
	const_iterator cbegin() const {
		return const_iterator(*this, 0);
	}
	const_iterator cend() const {
		return const_iterator(*this, key_list.size());
	}
	const_iterator find(const key_type &key) const;
	const_iterator find(const key_type &key, const key_type &comp_key) const;
	const std::vector<char> &get_metadata() const {
		return metadata;
	}
	const std::vector<base_type> &get_data() const {
		return data;
	}
	// start and length are in bits, not basepairs
	void get_sequence(data_offset_type start, data_offset_type length, std::string &) const;
	void print() const;
};

#endif // !_HASHL_INDEX_H
