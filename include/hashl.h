#ifndef _HASHL_H
#define _HASHL_H

// This is a class designed to be a memory efficient storage for counting
// n-mers; keys and values are stored in separate arrays to avoid
// alignment issues; a map is used to handle the infrequent case of values
// beyond the size of the small value type
//
// key values are lookups into an internal array; a metadata blob is also stored

#include "refcount_array.h"	// recount_array
#include <limits.h>	// UCHAR_MAX, ULONG_MAX
#include <map>		// map<>
#include <stdint.h>	// uint64_t
#include <string>	// string

class hashl {
    public:	// type declarations
	typedef unsigned char small_value_type;
	typedef unsigned long value_type;
	typedef unsigned long offset_type;
	typedef uint64_t base_type;
	enum { max_small_value = UCHAR_MAX, invalid_key = ULONG_MAX };

	class key_type {
	    private:
		size_t word_width;
		base_type *k;		// stored in reverse - high word in [0]
		friend class hashl;	// read access to k needed for key_equal()
	    public:
		void copy_in(const base_type *, const offset_type);
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
		explicit key_type(const hashl &__a) : word_width(__a.words()), k(new base_type[__a.words()]), bit_shift((__a.bits() - 2) % (sizeof(base_type) * 8)), high_mask(static_cast<base_type>(-1) >> (sizeof(base_type) * 8 - __a.bits() % (sizeof(base_type) * 8))) {
			for (size_t __i(0); __i != word_width; ++__i) {
				k[__i] = 0;
			}
		}
		~key_type(void) {
			delete[] k;
		}
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
			for (; __i != word_width - 1; ++__i) {
				k[__i] = (k[__i] << 2) | (k[__i + 1] >> __n);
			}
			k[__i] = (k[__i] << 2) | __x;
			k[0] &= high_mask;
		}
		void push_front(const base_type __x) {
			const size_t __n(sizeof(base_type) * 8 - 2);
			for (size_t __i(word_width - 1); __i != 0; --__i) {
				k[__i] = (k[__i - 1] << __n) | (k[__i] >> 2);
			}
			k[0] = (__x << bit_shift) | (k[0] >> 2);
		}
		void make_complement(const key_type &);
	};

	class const_iterator {	// only useful for pulling out data
	    private:
		const hashl * const list;
		offset_type offset;
	    private:
		void get_value(void);
	    public:
		value_type value;
		explicit const_iterator(const hashl * const a, const offset_type i) : list(a), offset(i) {
			get_value();
		}
		~const_iterator(void) { }
		void get_key(key_type &key) const {
			key.copy_in(list->data, offset);
		}
		bool operator==(const const_iterator &__a) const {
			return list == __a.list && offset == __a.offset;
		}
		bool operator!=(const const_iterator &__a) const {
			return !(*this == __a);
		}
		void increment(void) {
			if (offset < list->modulus) {
				for (++offset; offset < list->modulus && list->key_list[offset] == invalid_key; ++offset) { }
				get_value();
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
	offset_type used_elements;
	offset_type modulus;
	offset_type collision_modulus;
	offset_type data_size;
	size_t metadata_size;
	size_t bit_width;			// only used by key_type
	size_t word_width;
	offset_type *key_list;
	small_value_type *value_list;
	const base_type *data;
	const void *metadata;
	std::map<offset_type, value_type> value_map;	// for overflow
    protected:
	std::string boilerplate(void) const;
	offset_type find_offset(const key_type &) const;
	offset_type insert_offset(const key_type &key, const key_type &comp_key, offset_type);
    private:
	offset_type insert_key(offset_type, offset_type);
	void get_key(offset_type, key_type &) const;
	bool key_equal(offset_type, const key_type &) const;
    public:
	explicit hashl(void) : used_elements(0), modulus(0), collision_modulus(0), data_size(0), metadata_size(0), bit_width(0), word_width(0), key_list(0), value_list(0), data(0), metadata(0) { }
	// size of hash, bit size of key_type, sequence data
	explicit hashl(const offset_type i, const size_t j, const base_type * const d, const offset_type d_size) : metadata_size(0), metadata(0) {
		init(i, j, d, d_size);
	}
	~hashl(void);
	void init(offset_type, size_t, const base_type *, offset_type);
	void init_from_file(int);
	// will not insert new key
	bool increment(const key_type &);
	// will insert new key if missing
	bool increment(const key_type &key, const key_type &comp_key, offset_type);
	value_type value(const key_type &) const;
	offset_type size(void) const {
		return used_elements;
	}
	offset_type capacity(void) const {
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
	offset_type overflow_size(void) const {
		return value_map.size();
	}
	const_iterator begin(void) const;
	const_iterator end(void) const;
	const_iterator find(const key_type &key) const {
		return const_iterator(this, find_offset(key));
	}
	void save(int) const;
	void set_metadata(const void *metadata_in, size_t metadata_size_in);
	void get_metadata(const void * &metadata_out, size_t &metadata_size_out) const;
	// throw out current hash, but not data or metadata
	void resize(offset_type);
};

#endif // !_HASHL_H
