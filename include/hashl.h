#ifndef _HASHL_H
#define _HASHL_H

// This is a class designed to be a memory efficient storage for counting
// n-mers; keys and values are stored in separate arrays to avoid
// alignment issues
//
// key values are offsets into an internal array; a metadata blob is also stored

#include "hashl_key_type.h"	// hashl_key_type<>
#include <limits.h>	// UCHAR_MAX, ULONG_MAX
#include <stdint.h>	// uint64_t
#include <string>	// string
#include <utility>	// pair<>
#include <vector>	// vector<>

class hashl {
    public:	// type declarations
	typedef unsigned char small_value_type;
	typedef unsigned long hash_offset_type;
	typedef unsigned long data_offset_type;
	typedef uint64_t base_type;
	typedef hashl_key_type<hashl> key_type;
	// invalid_value must be greater than max_small_value
	enum { max_small_value = UCHAR_MAX - 1, invalid_value = UCHAR_MAX, invalid_key = ULONG_MAX };

	class iterator {
	    private:
		hashl &list;
		hash_offset_type offset_;
	    public:
		explicit iterator(hashl &a, const hash_offset_type i) : list(a), offset_(i) {
			if (i >= list.modulus || list.empty()) {
				offset_ = list.modulus;
			} else {
				offset_ = i;
				if (list.key_list[offset_] == invalid_key) {
					increment();
				}
			}
		}
		~iterator() { }
		// value()/key() undefined if called when pointing to end()
		small_value_type &operator*() {
			return list.value_list[offset_];
		}
		const hash_offset_type &offset() const {
			return offset_;
		}
		void key(key_type &key) const {
			key.copy_in(list.data, list.key_list[offset_]);
		}
		bool operator==(const iterator &__a) const {
			return &list == &__a.list && offset_ == __a.offset_;
		}
		bool operator!=(const iterator &__a) const {
			return !(*this == __a);
		}
		void increment() {
			if (offset_ < list.modulus) {
				for (++offset_; offset_ < list.modulus && list.key_list[offset_] == invalid_key; ++offset_) { }
			}
		}
		iterator operator++() {
			increment();
			return *this;
		}
		iterator operator++(int) {
			const iterator __tmp(*this);
			increment();
			return __tmp;
		}
	};

	class const_iterator {
	    private:
		const hashl &list;
		hash_offset_type offset_;
	    public:
		explicit const_iterator(const hashl &a, const hash_offset_type i) : list(a) {
			if (i >= list.modulus || list.empty()) {
				offset_ = list.modulus;
			} else {
				offset_ = i;
				if (list.key_list[offset_] == invalid_key) {
					increment();
				}
			}
		}
		const_iterator(const const_iterator &a) : list(a.list), offset_(a.offset_) { }
		~const_iterator() { }
		// value()/key() undefined if called when pointing to end()
		const small_value_type &operator*() const {
			return list.value_list[offset_];
		}
		const hash_offset_type &offset() const {
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
		void increment() {
			if (offset_ < list.modulus) {
				for (++offset_; offset_ < list.modulus && list.key_list[offset_] == invalid_key; ++offset_) { }
			}
		}
		const_iterator operator++() {
			increment();
			return *this;
		}
		const_iterator operator++(int) {
			const const_iterator __tmp(*this);
			increment();
			return __tmp;
		}
	};

    protected:
	std::vector<data_offset_type> key_list;
	std::vector<small_value_type> value_list;
	std::vector<small_value_type> value_list_backup;	// only used for filtering
	std::vector<base_type> data;
	std::vector<char> metadata;
	hash_offset_type used_elements;
	hash_offset_type modulus;
	hash_offset_type collision_modulus;
	size_t bit_width;
	size_t word_width;
    protected:
	std::string boilerplate() const;
	hash_offset_type find_offset(const key_type &key) const;
	hash_offset_type find_offset(const key_type &key, const key_type &comp_key) const;
	hash_offset_type insert_offset(const key_type &key, const key_type &comp_key, data_offset_type);
    private:
	hash_offset_type insert_key(hash_offset_type, data_offset_type);
    public:
	explicit hashl() : used_elements(0), modulus(0), collision_modulus(0), bit_width(0), word_width(0) { }
	// size of hash, bit size of key_type, sequence data
	explicit hashl(const hash_offset_type size, const size_t bits, std::vector<base_type> &data_in) {
		init(size, bits, data_in);
	}
	~hashl() { }
	// this trashes data_in by swapping it into the hash
	void init(hash_offset_type size, size_t bits, std::vector<base_type> &data_in);
	void init_from_file(int);
	// will not insert new key
	void increment(const key_type &key, const key_type &comp_key);
	// will insert new key if missing, returns false if insertion failed
	bool increment(const key_type &key, const key_type &comp_key, data_offset_type);
	// will insert new key if missing, or change an existing one to invalid
	bool insert_unique(const key_type &key, const key_type &comp_key, data_offset_type);
	// will insert a key with an invalid value, or convert existing value to invalid
	bool insert_invalid(const key_type &key, const key_type &comp_key, data_offset_type);
	small_value_type value(const key_type &) const;
	std::pair<data_offset_type, small_value_type> entry(const key_type &) const;
	hash_offset_type size() const {
		return used_elements;
	}
	hash_offset_type capacity() const {
		return modulus;
	}
	bool empty() const {
		return used_elements == 0;
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
		return const_iterator(*this, modulus);
	}
	iterator begin() {
		return iterator(*this, 0);
	}
	iterator end() {
		return iterator(*this, modulus);
	}
	const_iterator find(const key_type &key) const {
		return const_iterator(*this, find_offset(key));
	}
	const_iterator find(const key_type &key, const key_type &comp_key) const {
		return const_iterator(*this, find_offset(key, comp_key));
	}
	void save(int) const;
	void set_metadata(std::vector<char> &metadata_in) {
		metadata.swap(metadata_in);
	}
	const std::vector<char> &get_metadata() const {
		return metadata;
	}
	const std::vector<base_type> &get_data() const {
		return data;
	}
	// start and length are in bits, not basepairs
	void get_sequence(data_offset_type start, data_offset_type length, std::string &) const;
	void resize(hash_offset_type new_size);
	void purge_invalid_values();
	// add in new hashl - add new data, add or modify values
	// (<min_cutoff => ignored, <=cutoff => ++, >cutoff => invalid)
	bool add(const hashl &, small_value_type min_cutoff = 0, small_value_type max_cutoff = 1);
	void print() const;
	// save and zero value_list
	void filtering_prep(bool backup_values = 1);
	// restore value_list for values in min-max, set to invalid_value for rest
	void filtering_finish(small_value_type min, small_value_type max);
	// this destroy the hash, as it involves sorting key_list in place
	//void save_index(int);
};

#endif // !_HASHL_H
