#ifndef _HASH_H
#define _HASH_H

// This is a class designed to be a memory efficient storage for counting
// n-mers; keys and values are stored in separate arrays to avoid
// alignment issues; a map is used to handle the infrequent case of values
// beyond the size of the small value type

// The alt_list/alt_map arrays are available for storing extra information
// associated with each element in an efficient manner.

#include <limits.h>	// UCHAR_MAX, ULONG_MAX
#include <list>		// list<>
#include <map>		// map<>
#include <stdint.h>	// uint64_t
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

class hash {
    public:	// type declarations
	typedef uint64_t key_type;
	typedef unsigned char small_value_type;
	typedef size_t value_type;
	typedef size_t offset_type;
	enum { max_small_value = UCHAR_MAX };
	enum no_space_response_t { CLEAN_HASH = 1, TMP_FILE = 2 };

	class const_iterator {	// only useful for pulling out data
	    private:
		const hash *list;
		offset_type offset;
		std::map<key_type, std::pair<value_type, int> > next_keys;
	    public:
		key_type key;
		value_type value;
		explicit const_iterator(void) : list(0), offset(0), key(INVALID_KEY), value(0) { }
		// this can't be explicit
		const_iterator(const const_iterator &__a) : list(__a.list), offset(__a.offset), next_keys(__a.next_keys), key(__a.key), value(__a.value) { }
		explicit const_iterator(const hash *, offset_type);
		explicit const_iterator(const hash *, offset_type, const std::map<key_type, std::pair<value_type, int> > &);
		~const_iterator(void) { }
		const_iterator &operator=(const const_iterator &);
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
	bool can_overflow;		// allow values greater than max_small_value
	int no_space_response;
	int max_key_size;		// for radix sorting
	offset_type used_elements;
	offset_type modulus;
	offset_type collision_modulus;
	offset_type alt_size;
	key_type *key_list;
	small_value_type *value_list;
	// alt_list is a two dimensional array (modulus * alt_size) declared
	// as a single to maintain locality of values to reduce cache misses
	small_value_type *alt_list;
	std::string tmp_file_prefix;			// for TMP_FILE response
	std::map<key_type, value_type> value_map;	// for overflow
	std::map<key_type, value_type> *alt_map;	// for alt overflows
	std::list<std::string> state_files;		// for TMP_FILE response
    private:
	std::string boilerplate(void) const;
	offset_type find_empty_offset(key_type) const;
	void rehash(void);
	void rehash_alt(void);
	bool clean_hash(void);
	offset_type insert_key(offset_type, key_type);
	offset_type insert_offset(key_type);
	offset_type find_offset(key_type) const;
	bool add(key_type, value_type);
	bool add_alt(key_type, value_type, const value_type []);
	void shell_sort(offset_type, offset_type);
	void calculate_offsets(offset_type, offset_type, offset_type *, int) const;
	void radix_sort_internal(offset_type, offset_type, offset_type *, int);
	void radix_sort(offset_type);
	void save_state(void);
	void squash_hash(void);
	bool get_next_entry(int, key_type &, value_type &, offset_type &) const;
	void prep_for_readback(offset_type &, std::map<key_type, std::pair<value_type, int> > &);
    public:
	explicit hash(void) : can_overflow(1), no_space_response(0), max_key_size(sizeof(key_type) * 8), used_elements(0), modulus(0), collision_modulus(0), alt_size(0), key_list(0), value_list(0), alt_list(0), tmp_file_prefix(""), alt_map(0), state_files() { }
	explicit hash(const offset_type size_in, const offset_type alt_size_in = 0) : can_overflow(1), no_space_response(0), max_key_size(sizeof(key_type) * 8), tmp_file_prefix("") {
		init(size_in, alt_size_in);
	}
	~hash(void);
	void init(offset_type, offset_type = 0);
	void init_from_file(int);
	bool increment(key_type);
	bool increment_alt(key_type, offset_type);
	value_type value(key_type) const;
	value_type value(key_type, value_type []) const;
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
	void clear(bool = 0);
	// when using the TMP_FILE no_space_response, you should call
	// clear() after iterating over the hash, and you shouldn't make
	// changes to the hash during iteration
	const_iterator begin(void);
	const_iterator end(void) const {
		return const_iterator(this, modulus);
	}
	const_iterator find(const key_type key) const {
		return const_iterator(this, find_offset(key));
	}
	void save(int) const;
	void clean_hash(value_type, value_type);
	bool add_hash(hash &);
	void set_overflow(const bool i) {
		can_overflow = i;
	}
	// max_key_size must be a multiple of 8, and at least 8
	void set_max_key_size(const unsigned int i) {
		if (i < 8) {
			max_key_size = 8;
		} else if (i > sizeof(key_type) * 8) {
			max_key_size = sizeof(key_type) * 8;
		} else {
			max_key_size = ((i + 7) / 8) * 8;
		}
	}
	void set_no_space_response(int, const std::string & = "NONE");
};

#endif // !_HASH_H
