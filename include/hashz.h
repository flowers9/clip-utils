#ifndef _HASHZ_H
#define _HASHZ_H

/*
 * This is a class designed to be a memory efficient storage for counting
 * n-mers; keys and values are stored in separate arrays to avoid
 * alignment issues; a map is used to handle the infrequent case of values
 * beyond the size of the small value type
 *
 * The alt_list/alt_map arrays are available for storing extra information
 * associated with each element in an efficient manner.
 */

#include <gmp.h>	/* mpz_clear(), mpz_t */
#include <limits.h>	/* UCHAR_MAX */
#include <map>		/* map<> */
#include <string>	/* string */

class hashz {
    public:	/* type declarations */
	typedef mpz_t key_type;
	typedef unsigned char small_value_type;
	typedef unsigned long value_type;
	typedef unsigned long offset_type;
	enum { max_small_value = UCHAR_MAX };

	class const_iterator {	/* only useful for pulling out data */
	    private:
		const hashz *list;	/* hope *list doesn't move... */
		offset_type offset;
	    public:
		key_type key;
		value_type value;
		const_iterator(void) : list(NULL), offset(0), value(0) { }
		const_iterator(const const_iterator &);
		const_iterator(const hashz *, offset_type);
		~const_iterator(void) {
			if (list != NULL) {
				mpz_clear(key);
			}
		}
		const_iterator &operator=(const const_iterator &);
		bool operator==(const const_iterator &) const;
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
	/*
	 * a non-palindrome random pattern with one of the two highest
	 * bits set, and such that the complement is of lower value (so
	 * it  doesn't collide if you only store the "smaller" key of
	 * key/key-complement); if the given n-mer is found, the repeat
	 * count for it will probably be low and it will be skipped by
	 * the iterator;
	 */
	key_type invalid_key;
	char *key_cstr;			/* dedicated space for key conversion */
	offset_type used_elements;
	offset_type modulus;
	offset_type collision_modulus;
	key_type *key_list;
	small_value_type alt_size;
	small_value_type *value_list;
	std::map<std::string, value_type> value_map;	/* for overflow */
	small_value_type **alt_list;
	std::map<std::string, value_type> *alt_map;	/* for alt overflows */
    private:
	offset_type insert_key(offset_type, const key_type &);
	offset_type insert_offset(const key_type &);
	offset_type find_offset(const key_type &) const;
    public:
	/* size of hash, bit size of key_type, size of alt array */
	explicit hashz(offset_type, unsigned long, small_value_type = 0);
	~hashz(void);
	bool increment(const key_type &);
	bool increment_alt(const key_type &, offset_type);
	value_type value(const key_type &) const;
	value_type value(const key_type &, value_type []) const;
	const offset_type &size(void) const {
		return used_elements;
	}
	const offset_type &capacity(void) const {
		return modulus;
	}
	const key_type &get_invalid_key(void) const {
		return invalid_key;
	}
	std::map<std::string, value_type>::size_type overflow_size(void) const {
		return value_map.size();
	}
	void clear(void);
	const_iterator begin(void) const;
	const_iterator end(void) const;
};

#endif /* !_HASHZ_H */
