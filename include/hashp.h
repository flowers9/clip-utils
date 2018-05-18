#ifndef _HASHP_H
#define _HASHP_H

// This is a class designed to be a memory efficient storage for
// n-mer lookup (n <= 32); two array lookup values (start and stop)
// are associated with each entry, for secondary checking against the
// full sequence

#include <limits.h>	// ULONG_MAX
#include <new>		// delete[]
#include <stdint.h>	// uint64_t
#include <sys/types.h>	// size_t

// default invalid key - it will get changed if a collision with a valid
// key happens

// choose size of constant by what will hold it
#if ULONG_MAX >= 0xffffffffffffffffULL
#define INVALID_KEY 0xf27a829f3eedb68eUL
#else
#define INVALID_KEY 0xf27a829f3eedb68eULL
#endif

class hashp {
    public:	// type declarations
	typedef uint64_t key_type;
	typedef unsigned long offset_type;
	typedef size_t value_type;

	class const_iterator {	// only useful for pulling out data
	    private:
		const hashp *list_;
		offset_type offset_;
	    public:
		key_type key;
		value_type v1_out, v2_out;
		explicit const_iterator(void) : list_(0), offset_(0), key(INVALID_KEY), v1_out(-1), v2_out(-1) { }
		// this can't be explicit
		const_iterator(const const_iterator &__a) : list_(__a.list_), offset_(__a.offset_), key(__a.key), v1_out(__a.v1_out), v2_out(__a.v2_out) { }
		explicit const_iterator(const hashp *, offset_type);
		~const_iterator(void) { }
		const_iterator &operator=(const const_iterator &);
		bool operator==(const const_iterator &__a) const {
			return list_ == __a.list_ && offset_ == __a.offset_;
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
	};

    private:
	offset_type used_elements;
	offset_type modulus;
	offset_type collision_modulus;
	key_type invalid_key;
	key_type *key_list;
	value_type *v1, *v2;
    private:
	offset_type insert_key(offset_type, key_type);
	offset_type insert_offset(key_type);
	offset_type find_offset(key_type) const;
    public:
	explicit hashp(void) : used_elements(1), modulus(1), collision_modulus(0), invalid_key(INVALID_KEY), key_list(0), v1(0), v2(0) { }
	explicit hashp(const offset_type size_in) : used_elements(1), invalid_key(INVALID_KEY) {
		init(size_in);
	}
	~hashp(void);
	void init(offset_type);
	bool add(key_type, value_type, value_type);
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
	bool has_key(const key_type key, value_type &v1_out, value_type &v2_out) const {
		const offset_type i(find_offset(key));
		if (i == modulus) {
			return 0;
		}
		v1_out = v1[i];
		v2_out = v2[i];
		return 1;
	}
	const_iterator begin(void);
	const_iterator end(void) const {
		return const_iterator(this, modulus);
	}
};

#endif // !_HASHP_H
