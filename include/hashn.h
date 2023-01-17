#ifndef _HASHN_H
#define _HASHN_H

// This is a class designed to be a memory efficient storage for counting
// n-mers; keys and values are stored in separate arrays to avoid
// alignment issues; a map is used to handle the infrequent case of values
// beyond the size of the small value type
//
// The alt_list/alt_map arrays are available for storing extra information
// associated with each element in an efficient manner.

#include "refcount_array.h"	// refcount_array
#include <algorithm>	// swap()
#include <limits.h>	// UCHAR_MAX, ULONG_MAX
#include <list>		// list<>
#include <map>		// map<>
#include <stdint.h>	// uint64_t
#include <string>	// string
#include <utility>	// pair<>

// a non-palindrome random pattern with one of the two highest bits set,
// and such that the complement is of lower value (so it only collides
// when using multiple of 32 mers, and then only for 64+ and only for
// palindromic sequences leading with the key); the current pattern
// corresponds to TTAGCTGGGAAGGCTTATTGTGTCGTCGGATG

// choose size of constant by what will hold it
#if ULONG_MAX >= 0xffffffffffffffffULL
#define INVALID_KEY 0xf27a829f3eedb68eUL
#else
#define INVALID_KEY 0xf27a829f3eedb68eULL
#endif

class hashn {
    public:	// type declarations
	typedef unsigned char small_value_type;
	typedef unsigned long value_type;
	typedef unsigned long offset_type;
	typedef uint64_t base_type;
	enum { max_small_value = UCHAR_MAX };
	enum no_space_response_t { CLEAN_HASH = 1, TMP_FILE = 2 };

	class key_type_base {
	    protected:
		size_t word_width;
		base_type *k;		// stored in reverse - high word in [0]
	    public:
		explicit key_type_base(void) { }
		explicit key_type_base(size_t __i, base_type * const __j) : word_width(__i), k(__j) { }
		~key_type_base(void) { }
		void copy_out(base_type * const __x) const {
			for (size_t __i(0); __i != word_width; ++__i) {
				__x[__i] = k[__i];
			}
		}
		void swap(base_type * const __x) const {
			for (size_t __i(0); __i != word_width; ++__i) {
				std::swap(__x[__i], k[__i]);
			}
		}
		bool equal(const base_type * const __x) const {
			for (size_t __i(0); __i != word_width; ++__i) {
				if (k[__i] != __x[__i]) {
					return 0;
				}
			}
			return 1;
		}
		bool operator==(const key_type_base &__a) const {
			return equal(__a.k);
		}
		bool operator!=(const key_type_base &__a) const {
			return !equal(__a.k);
		}
		key_type_base &operator=(const key_type_base &__a) {
			word_width = __a.word_width;
			k = __a.k;
			return *this;
		}
		base_type hash(void) const {
			base_type __x(k[0]);
			for (size_t __i(1); __i != word_width; ++__i) {
				__x ^= k[__i];
			}
			return __x;
		}
		int basepair(const size_t __i) const {
			const size_t __n(__i / (sizeof(base_type) * 8));
			return (k[word_width - 1 - __n] >> (__i - __n * sizeof(base_type) * 8)) & 3;
		}
		std::string string(void) const;
	};

	// key_type with data kept in an allocated array
	class key_type : public key_type_base {
	    private:
		const size_t bit_shift;		// precalc for push_front
		const base_type high_mask;	// precalc for push_back
	    public:
		explicit key_type(void) : bit_shift(0), high_mask(0) { }
		explicit key_type(const hashn &__a) : key_type_base(__a.words(), new base_type[__a.words()]), bit_shift((__a.bits() - 2) % (sizeof(base_type) * 8)), high_mask(static_cast<base_type>(-1) >> (sizeof(base_type) * 8 - __a.bits() % (sizeof(base_type) * 8))) {
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
		void copy_in(const base_type * const __x) {
			for (size_t __i(0); __i != word_width; ++__i) {
				k[__i] = __x[__i];
			}
		}
	};

	// key_type with data internal to the array
	class key_type_internal : public key_type_base {
	    public:
		explicit key_type_internal(void) { }
		explicit key_type_internal(const hashn &__a, base_type * const __j) : key_type_base(__a.words(), __j) { }
		~key_type_internal(void) { }
		void assign(const hashn &__a, base_type * const __i) {
			word_width = __a.words();
			k = __i;
		}
		void assign(const hashn &__a, const offset_type __i) {
			assign(__a, __a.key_list + __i * __a.words());
		}
	};

	// pointers into buffer for reading from files, or into hash for memory
	class sort_key {
	    private:
		const size_t word_width;
	    public:
		base_type *k;
		explicit sort_key(const int i, base_type * const j) : word_width(i), k(j) { }
		explicit sort_key(const sort_key &i) : word_width(i.word_width), k(i.k) { }
		bool operator<(const sort_key &__a) const {
			for (size_t __i(0); __i != word_width; ++__i) {
				if (k[__i] != __a.k[__i]) {
					return k[__i] < __a.k[__i];
				}
			}
			return 0;
		}
		bool operator==(const sort_key &__x) const {
			for (size_t __i(0); __i != word_width; ++__i) {
				if (k[__i] != __x.k[__i]) {
					return 0;
				}
			}
			return 1;
		}
	};

	class const_iterator {	// only useful for pulling out data
	    private:
		const hashn * const list;
		offset_type offset;
		refcount_array<base_type> key_buffer;
		std::map<sort_key, std::pair<value_type, int> > next_keys;
	    public:
		key_type_internal key;
		value_type value;
		explicit const_iterator(const hashn *, offset_type);
		explicit const_iterator(const hashn *, offset_type, const std::map<sort_key, std::pair<value_type, int> > &, refcount_array<base_type> &);
		~const_iterator(void) { }
		bool operator==(const const_iterator &__a) const {
			return list == __a.list && offset == __a.offset;
		}
		bool operator!=(const const_iterator &__a) const {
			return !(*this == __a);
		}
		void increment(void);
		const_iterator operator++(void) {
			increment();
			return *this;
		}
		const_iterator operator++(int) {
			const const_iterator __tmp(*this);
			increment();
			return __tmp;
		}
		void get_alt_values(value_type []) const;
	};

    private:
	int no_space_response;
    protected:
	// INVALID_KEY repeated to fill an entire key size
	key_type_internal invalid_key;
	offset_type used_elements;
	offset_type modulus;
	offset_type collision_modulus;
	size_t bit_width;			// only used by key_type
	size_t word_width;
	offset_type alt_size;
	base_type *key_list;
	small_value_type *value_list;
	// alt_list is a two dimensional array (modulus * alt_size) declared
	// as a single to maintain locality of values to reduce cache misses
	small_value_type *alt_list;
	std::map<std::string, value_type> value_map;	// for overflow
	std::map<std::string, value_type> *alt_map;	// for alt overflows
    private:
	std::string tmp_file_prefix;			// for TMP_FILE response
	std::list<std::string> state_files;		// for TMP_FILE response
    protected:
	std::string boilerplate(void) const;
	offset_type find_offset(const key_type_base &) const;
	offset_type insert_offset(const key_type_base &);
    private:
	offset_type find_empty_offset(const base_type *) const;
	void rehash(void);
	void rehash_alt(void);
	bool clean_hash(void);
	offset_type insert_key(const offset_type, base_type *, const key_type_base &);
	void shell_sort(offset_type, offset_type);
	void calculate_offsets(offset_type, offset_type, offset_type *, int, int) const;
	void radix_sort_internal(offset_type, offset_type, offset_type *, int);
	void radix_sort(offset_type);
	void save_state(void);
	void squash_hash(void);
	bool get_next_entry(int, sort_key &, value_type &, offset_type &) const;
	void prep_for_readback(offset_type &, std::map<sort_key, std::pair<value_type, int> > &, refcount_array<base_type> &);
	base_type hash(const base_type * const __k) const {
		base_type __x(__k[0]);
		for (size_t __i(1); __i != word_width; ++__i) {
			__x ^= __k[__i];
		}
		return __x;
	}
	void copy(base_type * const __x, const base_type * const __y) {
		for (size_t __i(0); __i != word_width; ++__i) {
			__x[__i] = __y[__i];
		}
	}
	void swap(base_type * const __x, base_type * const __y) {
		for (size_t __i(0); __i != word_width; ++__i) {
			std::swap(__x[__i], __y[__i]);
		}
	}
	bool equal(const base_type * const __x, const base_type * const __y) const {
		for (size_t __i(0); __i != word_width; ++__i) {
			if (__x[__i] != __y[__i]) {
				return 0;
			}
		}
		return 1;
	}
	bool less_than(const base_type * const __x, const base_type * const __y) const {
		for (size_t __i(0); __i != word_width; ++__i) {
			if (__x[__i] != __y[__i]) {
				return __x[__i] < __y[__i];
			}
		}
		return 0;
	}
    public:
	explicit hashn(void) : no_space_response(0), used_elements(0), modulus(0), collision_modulus(0), bit_width(0), word_width(0), alt_size(0), key_list(0), value_list(0), alt_list(0), alt_map(0) { }
	// size of hash, bit size of key_type, size of alt array
	explicit hashn(offset_type i, const unsigned long j, const offset_type k = 0) : no_space_response(0) {
		init(i, j, k);
	}
	~hashn(void);
	void init(offset_type, unsigned long, offset_type = 0);
	void init_from_file(int);
	bool increment(const key_type &);
	bool increment_alt(const key_type &, offset_type);
	bool assign(const key_type &, value_type);
	value_type value(const key_type &) const;
	value_type value(const key_type &, value_type []) const;
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
	size_t bits(void) const {
		return bit_width;
	}
	size_t words(void) const {
		return word_width;
	}
	offset_type overflow_size(void) const {
		return value_map.size();
	}
	void clear(bool = 0);
	const_iterator begin(void);
	const_iterator end(void) const {
		return const_iterator(this, modulus);
	}
	const_iterator find(key_type key) const {
		return const_iterator(this, find_offset(key));
	}
	void save(int) const;
	void set_no_space_response(int, const std::string & = "NONE");
};

#endif // !_HASHN_H
