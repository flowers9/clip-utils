#ifndef _HASHL_H
#define _HASHL_H

// This is a class designed to be a memory efficient storage for counting
// n-mers; keys and values are stored in separate arrays to avoid
// alignment issues
//
// key values are offsets into an internal array; a metadata blob is also stored

#include <limits.h>	// UCHAR_MAX, ULONG_MAX
#include <stdint.h>	// uint64_t
#include <string>	// string
#include <utility>	// pair<>
#include <vector>	// vector<>
#include <iostream>

class hashl {
    public:	// type declarations
	typedef unsigned char small_value_type;
	typedef unsigned long value_type;
	typedef unsigned long hash_offset_type;
	typedef unsigned long data_offset_type;
	typedef uint64_t base_type;
	// invalid_value must be greater than max_small_value
	enum { max_small_value = UCHAR_MAX - 1, invalid_value = UCHAR_MAX, invalid_key = ULONG_MAX };

	class key_type {
	    private:
		std::vector<base_type> k;	// stored in reverse - high word in [0]
		const size_t word_width;
	    public:
		void copy_in(const std::vector<base_type> &, const data_offset_type);
		bool operator==(const key_type &__a) const {
			for (size_t __i(0); __i < word_width; ++__i) {
				if (k[__i] != __a.k[__i]) {
					return 0;
				}
			}
			return 1;
		}
		bool operator!=(const key_type &__a) const {
			return !(*this == __a);
		}
		base_type hash(void) const {
			base_type __x(k[0]);
			for (size_t __i(1); __i < word_width; ++__i) {
				__x ^= k[__i];
			}
			return __x;
		}
		int basepair(const size_t __i) const {
			const size_t __n(__i / (sizeof(base_type) * 8));
			return (k[word_width - 1 - __n] >> (__i - __n * sizeof(base_type) * 8)) & 3;
		}
	    private:
		const size_t bit_shift;		// precalc for push_front
		const base_type high_mask;	// precalc for push_back
	    public:
		explicit key_type(const hashl &__a) : k(__a.words(), 0), word_width(__a.words()), bit_shift((__a.bits() - 2) % (sizeof(base_type) * 8)), high_mask(static_cast<base_type>(-1) >> (sizeof(base_type) * 8 - __a.bits() % (sizeof(base_type) * 8))) { }
		~key_type(void) { }
		bool operator<(const key_type &__a) const {
			for (size_t __i(0); __i != word_width; ++__i) {
				if (k[__i] != __a.k[__i]) {
					return k[__i] < __a.k[__i];
				}
			}
			return 0;
		}
		void push_back(const base_type __x) {
			const size_t __n(sizeof(base_type) * 8 - 2);
			size_t __i(0);
			for (; __i < word_width - 1; ++__i) {
				k[__i] = (k[__i] << 2) | (k[__i + 1] >> __n);
			}
			k[__i] = (k[__i] << 2) | __x;
			k[0] &= high_mask;
		}
		void push_front(const base_type __x) {
			const size_t __n(sizeof(base_type) * 8 - 2);
			for (size_t __i(word_width - 1); __i > 0; --__i) {
				k[__i] = (k[__i - 1] << __n) | (k[__i] >> 2);
			}
			k[0] = (__x << bit_shift) | (k[0] >> 2);
		}
		void make_complement(const key_type &);
		void convert_to_string(std::string &) const;
		bool equal(const std::vector<base_type> &, const data_offset_type) const;
	};

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
		~iterator(void) { }
		// value()/key() undefined if called when pointing to end()
		small_value_type &operator*(void) {
			return list.value_list[offset_];
		}
		const hash_offset_type &offset(void) const {
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
		void increment(void) {
			if (offset_ < list.modulus) {
				for (++offset_; offset_ < list.modulus && list.key_list[offset_] == invalid_key; ++offset_) { }
			}
		}
		iterator operator++(void) {
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
		~const_iterator(void) { }
		// value()/key() undefined if called when pointing to end()
		const small_value_type &operator*(void) const {
			return list.value_list[offset_];
		}
		const hash_offset_type &offset(void) const {
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
		void increment(void) {
			if (offset_ < list.modulus) {
				for (++offset_; offset_ < list.modulus && list.key_list[offset_] == invalid_key; ++offset_) { }
			}
		}
		const_iterator operator++(void) {
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
	size_t bit_width;			// only used by key_type
	size_t word_width;
    protected:
	std::string boilerplate(void) const;
	hash_offset_type find_offset(const key_type &key) const;
	hash_offset_type find_offset(const key_type &key, const key_type &comp_key) const;
	hash_offset_type insert_offset(const key_type &key, const key_type &comp_key, data_offset_type);
    private:
	hash_offset_type insert_key(hash_offset_type, data_offset_type);
    public:
	explicit hashl(void) : used_elements(0), modulus(0), collision_modulus(0), bit_width(0), word_width(0) { }
	// size of hash, bit size of key_type, sequence data
	explicit hashl(const hash_offset_type size, const size_t bits, std::vector<base_type> &data_in) {
		init(size, bits, data_in);
	}
	~hashl(void) { }
	// this trashes data_in by swapping it into the hash
	void init(hash_offset_type size, size_t bits, std::vector<base_type> &data_in);
	void init_from_file(int);
	// will not insert new key
	void increment(const key_type &key, const key_type &comp_key);
	// will insert new key if missing, returns false if insertion failed
	bool increment(const key_type &key, const key_type &comp_key, data_offset_type);
	value_type value(const key_type &) const;
	std::pair<data_offset_type, value_type> entry(const key_type &) const;
	hash_offset_type size(void) const {
		return used_elements;
	}
	hash_offset_type capacity(void) const {
		return modulus;
	}
	bool empty(void) const {
		return used_elements == 0;
	}
	size_t bits(void) const {
		return bit_width;
	}
	size_t words(void) const {
		return word_width;
	}
	const_iterator cbegin(void) const {
		return const_iterator(*this, 0);
	}
	const_iterator cend(void) const {
		return const_iterator(*this, modulus);
	}
	iterator begin(void) {
		return iterator(*this, 0);
	}
	iterator end(void) {
		return iterator(*this, modulus);
	}
	const_iterator find(const key_type &key) const {
		return const_iterator(*this, find_offset(key));
	}
	void save(int) const;
	void set_metadata(std::vector<char> &metadata_in) {
		metadata.swap(metadata_in);
	}
	const std::vector<char> &get_metadata(void) const {
		return metadata;
	}
	const std::vector<base_type> &get_data(void) const {
		return data;
	}
	// start and length are in bits, not basepairs
	void get_sequence(data_offset_type start, data_offset_type length, std::string &) const;
	void resize(hash_offset_type new_size);
	// set values <= cutoff to 1, values > cutoff to invalid_value
	void normalize(small_value_type min_cutoff = 0, small_value_type max_cutoff = 1);
	// add in new hashl - add new data, add or modify values
	// (<min_cutoff => ignored, <=cutoff => ++, >cutoff => invalid)
	bool add(const hashl &, small_value_type min_cutoff = 0, small_value_type max_cutoff = 1);
	void print(void) const;
	// save and zero value_list
	void filtering_prep(void);
	// restore value_list for values in min-max, set to invalid_value for rest
	void filtering_finish(small_value_type min, small_value_type max);
};

#endif // !_HASHL_H
