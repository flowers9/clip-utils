#include "hashp.h"
#include "next_prime.h"	// next_prime()
#include <new>		// new
#include <sys/types.h>	// size_y

hashp::~hashp(void) {
	delete[] key_list;
	delete[] v1;
	delete[] v2;
}

void hashp::init(offset_type size_asked) {
	++size_asked;		// to account for minimum of one invalid_key 
	if (size_asked < 3) {	// to avoid collision_modulus == modulus
		size_asked = 3;
	}
	modulus = next_prime(size_asked);
	// collision_modulus just needs to be relatively prime with modulus;
	// since modulus is prime, any value will do - I made it prime for fun
	collision_modulus = next_prime(size_asked / 2);
	key_list = new key_type[modulus];
	v1 = new value_type[modulus];
	v2 = new value_type[modulus];
	// initialize keys
	for (offset_type i = 0; i != modulus; ++i) {
		key_list[i] = invalid_key;
	}
}

// insert a key at a particular location

hashp::offset_type hashp::insert_key(offset_type i, key_type key) {
	if (used_elements == modulus) {
		return modulus;		// hash table is full
	}
	++used_elements;
	key_list[i] = key;
	return i;
}

// find a key, or insert it if it doesn't exist; return modulus if hash is full

hashp::offset_type hashp::insert_offset(key_type key) {
	offset_type i(key % modulus);
	if (key_list[i] == invalid_key) {
		return insert_key(i, key);
	} else if (key_list[i] == key) {
		return i;
	}
	const offset_type j(collision_modulus - (key % collision_modulus));
	for (;;) {	// search over all elements
		i = (i + j) % modulus;
		if (key_list[i] == invalid_key) {
			return insert_key(i, key);
		} else if (key_list[i] == key) {
			return i;
		}
	}
}

// find a key; return modulus as offset if not found

hashp::offset_type hashp::find_offset(key_type key) const {
	offset_type i(key % modulus);
	if (key_list[i] == key) {
		return i;
	} else if (key_list[i] == invalid_key) {
		return modulus;
	}
	const offset_type j(collision_modulus - (key % collision_modulus));
	for (;;) {	// search over all elements
		i = (i + j) % modulus;
		if (key_list[i] == key) {
			return i;
		} else if (key_list[i] == invalid_key) {
			return modulus;
		}
	}
}

bool hashp::add(const key_type key, const value_type j, const value_type k) {
	if (key == invalid_key) {	// check if new key matches invalid_key
		// find new invalid key
		key_type x(invalid_key);
		for (++x; find_offset(x) != modulus; ++x) { }
		// change invalid_key to new invalid key over hash
		for (offset_type i = 0; i != modulus; ++i) {
			if (key_list[i] == invalid_key) {
				key_list[i] = x;
			}
		}
		invalid_key = x;
	}
	const offset_type i(insert_offset(key));
	if (i == modulus) {
		return 0;
	}
	v1[i] = j;
	v2[i] = k;
	return 1;
}

hashp::const_iterator hashp::begin() {
	if (used_elements == 1) {	// only key is INVALID_KEY
		return end();
	}
	const_iterator a(this, 0);
	if (a.key == INVALID_KEY) {	// advance to first valid value
		++a;
	}
	return a;
}

hashp::const_iterator::const_iterator(const hashp * const b, const offset_type i) : list_(b), offset_(i) {
	if (list_ == 0 || offset_ == list_->modulus) {
		key = INVALID_KEY;
		v1_out = v2_out = static_cast<value_type>(-1);
	} else {
		key = list_->key_list[offset_];
		v1_out = list_->v1[offset_];
		v2_out = list_->v2[offset_];
	}
}

hashp::const_iterator &hashp::const_iterator::operator=(const const_iterator &a) {
	list_ = a.list_;
	offset_ = a.offset_;
	key = a.key;
	v1_out = a.v1_out;
	v2_out = a.v2_out;
	return *this;
}

// advance to next used element

void hashp::const_iterator::increment() {
	if (offset_ == list_->modulus) {
		return;
	}
	for (++offset_; offset_ != list_->modulus && list_->key_list[offset_] == INVALID_KEY; ++offset_) { }
	if (offset_ != list_->modulus) {
		key = list_->key_list[offset_];
		v1_out = list_->v1[offset_];
		v2_out = list_->v2[offset_];
	} else {
		key = INVALID_KEY;
		v1_out = v2_out = static_cast<value_type>(-1);
	}
}
