#include "hashz.h"
#include "next_prime.h"	/* next_prime() */
#include <gmp.h>	/* gmp_randclear(), gmp_randinit_default(), gmp_randstate_t, mpz_clear(), mpz_comp(), mpz_fdiv_ui(), mpz_get_str(), mpz_init2(), mpz_init_set(), mpz_mul_2exp(), mpz_realloc(), mpz_set(), mpz_setbit(), mpz_set_ui(), mpz_sizeinbase(), mpz_sub_ui(), mpz_urandomb() */
#include <map>		/* map<> */
#include <stdio.h>	/* fprintf(), stderr */
#include <stdlib.h>	/* exit() */

hashz::hashz(offset_type size_asked, unsigned long bits, small_value_type alt_size_in) {
	if (alt_size_in > 8 * sizeof(offset_type)) {
		fprintf(stderr, "Error: hash alt size too large: %d > %ld\n", alt_size_in, 8 * sizeof(offset_type));
		exit(1);
	}
	alt_size = alt_size_in;
	used_elements = 0;
	if (size_asked < 3) {	/* to avoid collision_modulus == modulus */
		size_asked = 3;
	}
	modulus = next_prime(size_asked);
	/*
	 * collision_modulus just needs to be relatively prime with modulus;
	 * since modulus is prime, any value will do - I made it prime for fun
	 */
	collision_modulus = next_prime(size_asked / 2);
	key_list = new key_type[modulus];
	if (key_list == NULL) {
		fprintf(stderr, "Error: unable to allocate hash key array\n");
		exit(1);
	}
	value_list = new small_value_type[modulus];
	if (value_list == NULL) {
		fprintf(stderr, "Error: unable to allocate hash value array\n");
		exit(1);
	}
	if (alt_size == 0) {
		alt_list = NULL;
		alt_map = NULL;
	} else {
		alt_list = new small_value_type *[alt_size];
		if (alt_list == NULL) {
			fprintf(stderr, "Error: unable to allocate hash alt list array\n");
			exit(1);
		}
		small_value_type i;
		for (i = 0; i != alt_size; ++i) {
			alt_list[i] = new small_value_type[modulus];
			if (alt_list[i] == NULL) {
				fprintf(stderr, "Error: unable to allocate hash alt list %d array\n", i);
				exit(1);
			}
		}
		alt_map = new std::map<std::string, value_type>[alt_size];
		if (alt_map == NULL) {
			fprintf(stderr, "Error: unable to allocate hash alt map array\n");
			exit(1);
		}
	}
	/* allocate space for key conversion */
	mpz_init2(invalid_key, bits + 1);
	mpz_set_ui(invalid_key, 1);
	mpz_mul_2exp(invalid_key, invalid_key, bits);
	mpz_sub_ui(invalid_key, invalid_key, 1);
	key_cstr = new char[mpz_sizeinbase(invalid_key, 62) + 2];
	/* initialize invalid_key - random with high bit set */
	mpz_realloc2(invalid_key, bits);
	gmp_randstate_t state;
	gmp_randinit_default(state);
	mpz_urandomb(invalid_key, state, bits);
	gmp_randclear(state);
	mpz_setbit(invalid_key, bits - 1);
	mpz_setbit(invalid_key, 1);			/* non-palindrome */
	/* initialize keys */
	offset_type i;
	for (i = 0; i != modulus; ++i) {
		mpz_init_set(key_list[i], invalid_key);
	}
}

hashz::~hashz() {
	delete[] key_cstr;
	mpz_clear(invalid_key);
	offset_type i;
	for (i = 0; i != modulus; ++i) {
		mpz_clear(key_list[i]);
	}
	delete[] key_list;
	delete[] value_list;
	if (alt_list != NULL) {
		for (i = 0; i != alt_size; ++i) {
			delete[] alt_list[i];
		}
		delete[] alt_list;
	}
	if (alt_map != NULL) {
		delete[] alt_map;
	}
}

/* insert a key at a particular location */

hashz::offset_type hashz::insert_key(offset_type i, const key_type &key) {
	if (used_elements + 1 == modulus) {	/* hash is full */
		return modulus;
	}
	++used_elements;
	mpz_set(key_list[i], key);
	value_list[i] = 0;
	offset_type j;
	for (j = 0; j != alt_size; ++j) {
		alt_list[j][i] = 0;
	}
	return i;
}

/*
 * find a key, or insert it if it doesn't exist; return modulus if hash is full
 */

hashz::offset_type hashz::insert_offset(const key_type &key) {
	offset_type i = mpz_fdiv_ui(key, modulus);
	if (mpz_cmp(key_list[i], invalid_key) == 0) {
		return insert_key(i, key);
	} else if (mpz_cmp(key_list[i], key) == 0) {
		return i;
	}
	offset_type j = collision_modulus - mpz_fdiv_ui(key, collision_modulus);
	for (;;) {	/* search over all elements */
		i = (i + j) % modulus;
		if (mpz_cmp(key_list[i], invalid_key) == 0) {
			return insert_key(i, key);
		} else if (mpz_cmp(key_list[i], key) == 0) {
			return i;
		}
	}
}

/* find a key; return modulus as offset if not found */

hashz::offset_type hashz::find_offset(const key_type &key) const {
	offset_type i = mpz_fdiv_ui(key, modulus);
	if (mpz_cmp(key_list[i], key) == 0) {
		return i;
	} else if (mpz_cmp(key_list[i], invalid_key) == 0) {
		return modulus;
	}
	offset_type j = collision_modulus - mpz_fdiv_ui(key, collision_modulus);
	for (;;) {	/* search over all elements */
		i = (i + j) % modulus;
		if (mpz_cmp(key_list[i], key) == 0) {
			return i;
		} else if (mpz_cmp(key_list[i], invalid_key) == 0) {
			return modulus;
		}
	}
}

/* increment the count for a key */

bool hashz::increment(const key_type &key) {
	offset_type i = insert_offset(key);
	if (i == modulus) {	/* insert failed */
		return 0;
	}
	if (value_list[i] != max_small_value) {
		++value_list[i];
	} else {
		mpz_get_str(key_cstr, 62, key);
		++value_map[key_cstr];
	}
	return 1;
}

/* increment only the alt values, using x as a bit flag to mark which ones */

bool hashz::increment_alt(const key_type &key, offset_type x) {
	offset_type i = insert_offset(key);
	if (i == modulus) {	/* insert failed */
		return 0;
	}
	offset_type z;
	small_value_type j;
	mpz_get_str(key_cstr, 62, key);
	std::string key_str(key_cstr);
	for (j = 0, z = 1; j != alt_size; ++j, z <<= 1) {
		if (x & z) {
			if (alt_list[j][i] != max_small_value) {
				++alt_list[j][i];
			} else {
				++alt_map[j][key_str];
			}
		}
	}
	return 1;
}

/* return the value associated with a key */

hashz::value_type hashz::value(const key_type &key) const {
	offset_type i = find_offset(key);
	if (i == modulus) {	/* key not found */
		return 0;
	} else if (value_list[i] != max_small_value) {
		return value_list[i];
	} else {
		mpz_get_str(key_cstr, 62, key);
		/* use find() to avoid inserting a value into value_map */
		std::map<std::string, value_type>::const_iterator a = value_map.find(key_cstr);
		if (a == value_map.end()) {
			return max_small_value;
		} else {
			return a->second + max_small_value;
		}
	}
}

hashz::value_type hashz::value(const key_type &key, value_type x[]) const {
	offset_type i = find_offset(key);
	if (i == modulus) {	/* key not found */
		return 0;
	}
	mpz_get_str(key_cstr, 62, key);
	std::string key_str(key_cstr);
	small_value_type j;
	for (j = 0; j != alt_size; ++j) {
		if (alt_list[j][i] != max_small_value) {
			x[j] = alt_list[j][i];
		} else {
			/* use find() to avoid inserting a value into value_map */
			std::map<std::string, value_type>::const_iterator a = alt_map[j].find(key_str);
			if (a == alt_map[j].end()) {
				x[j] = max_small_value;
			} else {
				x[j] = a->second + max_small_value;
			}
		}
	}
	if (value_list[i] != max_small_value) {
		return value_list[i];
	} else {
		/* use find() to avoid inserting a value into value_map */
		std::map<std::string, value_type>::const_iterator a = value_map.find(key_str);
		if (a == value_map.end()) {
			return max_small_value;
		} else {
			return a->second + max_small_value;
		}
	}
}

/* reset hash to empty state */

void hashz::clear() {
	used_elements = 0;
	/* initialize keys */
	offset_type i;
	for (i = 0; i != modulus; ++i) {
		mpz_set(key_list[i], invalid_key);
	}
	value_map.clear();
	for (i = 0; i != alt_size; ++i) {
		alt_map[i].clear();
	}
}

hashz::const_iterator hashz::begin() const {
	if (used_elements == 0) {
		return end();
	}
	const_iterator a(this, 0);
	/* advance to first valid value */
	if (mpz_cmp(a.key, invalid_key) == 0) {
		++a;
	}
	return a;
}

hashz::const_iterator hashz::end() const {
	return const_iterator(this, modulus);
}

hashz::const_iterator::const_iterator(const const_iterator &a) {
	list = a.list;
	offset = a.offset;
	value = a.value;
	if (list != NULL) {
		mpz_init_set(key, a.key);
	}
}

hashz::const_iterator::const_iterator(const hashz *b, offset_type i) {
	list = b;
	offset = i;
	if (list == NULL) {
		value = 0;
	} else if (offset == list->modulus) {
		value = 0;
		mpz_init_set(key, list->invalid_key);
	} else {
		mpz_init_set(key, list->key_list[i]);
		if (list->value_list[offset] != max_small_value) {
			value = list->value_list[offset];
		} else {
			mpz_get_str(list->key_cstr, 62, key);
			std::map<std::string, value_type>::const_iterator a = list->value_map.find(list->key_cstr);
			if (a == list->value_map.end()) {
				value = max_small_value;
			} else {
				value = a->second + max_small_value;
			}
		}
	}
}

hashz::const_iterator &hashz::const_iterator::operator=(const const_iterator &a) {
	if (list != NULL) {
		if (a.list != NULL) {
			mpz_set(key, a.key);
		} else {
			mpz_clear(key);
		}
	} else if (a.list != NULL) {
		mpz_init_set(key, a.key);
	}
	list = a.list;
	offset = a.offset;
	value = a.value;
	return *this;
}

/* advance to next used element */

void hashz::const_iterator::increment() {
	if (offset == list->modulus) {
		return;
	}
	for (;;) {
		++offset;
		if (offset == list->modulus) {
			mpz_set(key, list->invalid_key);
			value = 0;
			return;
		} else if (mpz_cmp(list->key_list[offset], list->invalid_key) != 0) {
			mpz_set(key, list->key_list[offset]);
			if (list->value_list[offset] != max_small_value) {
				value = list->value_list[offset];
			} else {
				mpz_get_str(list->key_cstr, 62, key);
				std::map<std::string, value_type>::const_iterator a = list->value_map.find(list->key_cstr);
				if (a == list->value_map.end()) {
					value = max_small_value;
				} else {
					value = a->second + max_small_value;
				}
			}
			return;
		}
	}
}

/* iterator has to point to the same list and offset to match */

bool hashz::const_iterator::operator==(const const_iterator &a) const {
	return list == a.list && offset == a.offset;
}

/* extract list of alt_values associated with current offset */

void hashz::const_iterator::get_alt_values(value_type x[]) const {
	mpz_get_str(list->key_cstr, 62, key);
	std::string key_str(list->key_cstr);
	offset_type i;
	for (i = 0; i != list->alt_size; ++i) {
		if (list->alt_list[i][offset] != max_small_value) {
			x[i] = list->alt_list[i][offset];
		} else {
			std::map<std::string, value_type>::const_iterator a = list->alt_map[i].find(key_str);
			if (a == list->alt_map[i].end()) {
				x[i] = max_small_value;
			} else {
				x[i] = a->second + max_small_value;
			}
		}
	}
}
