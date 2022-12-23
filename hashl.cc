#include "hashl.h"
#include "itoa.h"	// itoa()
#include "local_endian.h"	// big_endian
#include "next_prime.h"	// next_prime()
#include "open_compressed.h"	// pfread()
#include "write_fork.h"	// pfwrite()
#include <map>		// map<>
#include <new>		// new
#include <stdio.h>	// fprintf(), stderr
#include <stdlib.h>	// exit()
#include <string.h>	// memcmp(), memcpy()
#include <string>	// string

hashl::~hashl() {
	delete[] key_list;
	delete[] value_list;
}

// description beginning of saved file

std::string hashl::boilerplate() const {
	std::string s("hashl\n");
	s += itoa(sizeof(base_type));
	s += " bytes\n";
#ifdef big_endian
	s += "big endian\n";
#else
	s += "little endian\n";
#endif
	return s;
}

void hashl::init(const hash_offset_type size_asked, const size_t bits_in, const base_type *data_in, const data_offset_type data_size_in) {
	bit_width = bits_in;
	data = data_in;
	data_size = data_size_in;
	word_width = (bit_width + 8 * sizeof(base_type) - 1) / (8 * sizeof(base_type));
	resize(size_asked);
}

void hashl::init_from_file(const int fd) {
	const std::string s(boilerplate());
	char t[s.size()];
	pfread(fd, t, s.size());
	if (memcmp(t, s.c_str(), s.size()) != 0) {
		fprintf(stderr, "Error: could not read hash from file: header mismatch\n");
		exit(1);
	}
	pfread(fd, &modulus, sizeof(modulus));
	pfread(fd, &collision_modulus, sizeof(collision_modulus));
	pfread(fd, &used_elements, sizeof(used_elements));
	pfread(fd, &bit_width, sizeof(bit_width));
	word_width = (bit_width + 8 * sizeof(base_type) - 1) / (8 * sizeof(base_type));
	pfread(fd, &metadata_size, sizeof(metadata_size));
	// read in in two steps, as metadata is const (as is data, below)
	void * const metadata_tmp(new char[metadata_size]);
	pfread(fd, metadata_tmp, metadata_size);
	metadata = metadata_tmp;
	pfread(fd, &data_size, sizeof(data_size));
	base_type * const data_tmp(new base_type[data_size]);
	pfread(fd, data_tmp, sizeof(base_type) * data_size);
	data = data_tmp;
	value_list = new small_value_type[modulus];
	pfread(fd, value_list, sizeof(small_value_type) * modulus);
	key_list = new data_offset_type[modulus];
	for (hash_offset_type i(0); i < modulus; ++i) {
		if (value_list[i] == 0) {
			key_list[i] = invalid_key;
		} else {
			pfread(fd, &key_list[i], sizeof(data_offset_type));
		}
	}
	// read in overflow map
	size_t x;
	pfread(fd, &x, sizeof(x));
	for (; x; --x) {
		hash_offset_type i;
		value_type j;
		pfread(fd, &i, sizeof(i));
		pfread(fd, &j, sizeof(j));
		value_map[i] = j;
	}
}

// insert a key at a particular location

hashl::hash_offset_type hashl::insert_key(const hash_offset_type i, const data_offset_type offset) {
	if (used_elements == modulus) {
		return modulus;		// hash table is full
	}
	++used_elements;
	key_list[i] = offset;
	value_list[i] = 0;	// init counts
	return i;
}

// create key from bit offset into data

void hashl::key_type::copy_in(const base_type *data, const data_offset_type i) {
	// move to start of sequence in data
	data += i / (sizeof(base_type) * 8);
	// how many bits we have in the first word
	const base_type starting_bits(sizeof(base_type) * 8 - i % (sizeof(base_type) * 8));
	// how many bits the first word is supposed to have for a key
	const base_type high_bits(bit_shift + 2);
	if (starting_bits == high_bits) {
		k[0] = data[0] & high_mask;
		for (size_t j(1); j < word_width; ++j) {
			k[j] = data[j];
		}
	} else if (starting_bits < high_bits) {		// shift left to fill up first word
		const int shift_left(high_bits - starting_bits);
		const int shift_right(sizeof(base_type) * 8 - shift_left);
		k[0] = ((data[0] << shift_left) | (data[1] >> shift_right)) & high_mask;
		for (size_t j(1); j < word_width; ++j) {
			k[j] = (data[j] << shift_left) | (data[j + 1] >> shift_right);
		}
	} else {					// shift right to empty out first word
		const int shift_right(starting_bits - high_bits);
		const int shift_left(sizeof(base_type) * 8 - shift_right);
		k[0] = (data[0] >> shift_right) & high_mask;
		for (size_t j(1); j < word_width; ++j) {
			k[j] = (data[j - 1] << shift_left) | (data[j] >> shift_right);
		}
	}
}

// generate internal key from bit offset into data, compare to key;
// equivalent to above subroutine, just with more breakpoints

bool hashl::key_type::equal(const base_type *data, const data_offset_type i) const {
	// move to start of sequence in data
	data += i / (sizeof(base_type) * 8);
	// how many bits we have in the first word
	const base_type starting_bits(sizeof(base_type) * 8 - i % (sizeof(base_type) * 8));
	// how many bits the first word is supposed to have for a key
	const base_type high_offset(bit_shift + 2);
	if (starting_bits == high_offset) {
		if (k[0] != (data[0] & high_mask)) {
			return 0;
		}
		for (size_t j(1); j < word_width; ++j) {
			if (k[j] != data[j]) {
				return 0;
			}
		}
	} else if (starting_bits < high_offset) {	// shift left to fill up first word
		const int shift_left(high_offset - starting_bits);
		const int shift_right(sizeof(base_type) * 8 - shift_left);
		if (k[0] != (((data[0] << shift_left) | (data[1] >> shift_right)) & high_mask)) {
			return 0;
		}
		for (size_t j(1); j < word_width; ++j) {
			if (k[j] != ((data[j] << shift_left) | (data[j + 1] >> shift_right))) {
				return 0;
			}
		}
	} else {					// shift right to empty out first word
		const int shift_right(starting_bits - high_offset);
		const int shift_left(sizeof(base_type) * 8 - shift_right);
		if (k[0] != ((data[0] >> shift_right) & high_mask)) {
			return 0;
		}
		for (size_t j(1); j < word_width; ++j) {
			if (k[j] != ((data[j - 1] << shift_left) | (data[j] >> shift_right))) {
				return 0;
			}
		}
	}
	return 1;
}

// find a key, or insert it if it doesn't exist; return modulus if hash is full

hashl::hash_offset_type hashl::insert_offset(const key_type &key, const key_type &comp_key, const data_offset_type offset) {
	const base_type key_hash(key < comp_key ? key.hash() : comp_key.hash());
	hash_offset_type i(key_hash % modulus);
	if (key_list[i] == invalid_key) {		// insert
		return insert_key(i, offset);
	} else if (key.equal(data, key_list[i]) || comp_key.equal(data, key_list[i])) {
		return i;				// already present
	}
	const hash_offset_type j(collision_modulus - key_hash % collision_modulus);
	for (;;) {	// search over all elements
		i = (i + j) % modulus;
		if (key_list[i] == invalid_key) {
			return insert_key(i, offset);
		} else if (key.equal(data, key_list[i]) || comp_key.equal(data, key_list[i])) {
			return i;
		}
	}
}

// make reverse complement of given key
// TODO: could be more efficient

void hashl::key_type::make_complement(const hashl::key_type &key) {
	const size_t bit_width(bit_shift + 2 + (word_width - 1) * sizeof(base_type) * 8);
	for (size_t i(0); i < bit_width; i += 2) {
		push_back(3 - key.basepair(i));
	}
}

// find a key; return modulus as offset if not found

hashl::hash_offset_type hashl::find_offset(const key_type &key) const {
	key_type comp_key(*this);
	comp_key.make_complement(key);
	const base_type key_hash(key < comp_key ? key.hash() : comp_key.hash());
	hash_offset_type i(key_hash % modulus);
	if (key_list[i] == invalid_key) {
		return modulus;
	} else if (key.equal(data, key_list[i]) || comp_key.equal(data, key_list[i])) {
		return i;
	}
	const hash_offset_type j(collision_modulus - key_hash % collision_modulus);
	for (;;) {	// search over all elements
		i = (i + j) % modulus;
		if (key_list[i] == invalid_key) {
			return modulus;
		} else if (key.equal(data, key_list[i]) || comp_key.equal(data, key_list[i])) {
			return i;
		}
	}
}

// increment the count for a key (but don't create a new entry if it doesn't exist)

bool hashl::increment(const key_type &key) {
	const hash_offset_type i(find_offset(key));
	if (i == modulus) {	// insert failed
		return 0;
	}
	if (value_list[i] < max_small_value) {
		++value_list[i];
	} else {
		++value_map[i];
	}
	return 1;
}

bool hashl::increment(const key_type &key, const key_type &comp_key, const data_offset_type offset) {
	const hash_offset_type i(insert_offset(key, comp_key, offset));
	if (i == modulus) {	// insert failed
		return 0;
	}
	if (value_list[i] < max_small_value) {
		++value_list[i];
	} else {
		++value_map[i];
	}
	return 1;
}

// return the value associated with a key

hashl::value_type hashl::value(const key_type &key) const {
	const hash_offset_type i(find_offset(key));
	if (i == modulus) {	// key not found
		return 0;
	} else if (value_list[i] < max_small_value) {
		return value_list[i];
	} else {
		// use find() to avoid inserting a value into value_map
		const std::map<hash_offset_type, value_type>::const_iterator a(value_map.find(i));
		if (a == value_map.end()) {
			return max_small_value;
		} else {
			return a->second + max_small_value;
		}
	}
}

hashl::const_iterator hashl::begin() const {
	if (used_elements == 0) {
		return end();
	}
	const_iterator a(this, 0);
	// advance to first valid value
	if (key_list[0] == invalid_key) {
		++a;
	}
	return a;
}

hashl::const_iterator hashl::end() const {
	return const_iterator(this, modulus);
}

void hashl::const_iterator::get_value() {
	if (offset < list->modulus) {
		value = list->value_list[offset];
		if (value == max_small_value) {
			const std::map<hash_offset_type, value_type>::const_iterator b(list->value_map.find(offset));
			if (b != list->value_map.end()) {
				value += b->second;
			}
		}
	} else {
		value = 0;
	}
}

void hashl::save(const int fd) const {
	const std::string s(boilerplate());
	pfwrite(fd, s.c_str(), s.size());
	pfwrite(fd, &modulus, sizeof(modulus));
	pfwrite(fd, &collision_modulus, sizeof(collision_modulus));
	pfwrite(fd, &used_elements, sizeof(used_elements));
	pfwrite(fd, &bit_width, sizeof(bit_width));
	pfwrite(fd, &metadata_size, sizeof(metadata_size));
	pfwrite(fd, metadata, metadata_size);
	pfwrite(fd, &data_size, sizeof(data_size));
	pfwrite(fd, data, sizeof(base_type) * data_size);
	pfwrite(fd, value_list, sizeof(small_value_type) * modulus);
	for (hash_offset_type i(0); i < modulus; ++i) {
		if (key_list[i] != invalid_key) {
			pfwrite(fd, &key_list[i], sizeof(data_offset_type));
		}
	}
	const size_t x(value_map.size());
	pfwrite(fd, &x, sizeof(x));
	std::map<hash_offset_type, value_type>::const_iterator a(value_map.begin());
	const std::map<hash_offset_type, value_type>::const_iterator end_a(value_map.end());
	for (; a != end_a; ++a) {
		pfwrite(fd, &a->first, sizeof(hash_offset_type));
		pfwrite(fd, &a->second, sizeof(value_type));
	}
}

void hashl::set_metadata(const void * const metadata_in, const size_t metadata_size_in) {
	metadata = metadata_in;
	metadata_size = metadata_size_in;
}

void hashl::get_metadata(const void * &metadata_out, size_t &metadata_size_out) const {
	metadata_out = metadata;
	metadata_size_out = metadata_size;
}

void hashl::resize(hash_offset_type size_asked) {
	const size_t old_modulus(modulus);
	const data_offset_type * const old_key_list(key_list);
	const small_value_type * const old_value_list(value_list);
	// key for map is offset into hash, which will change, so this needs updating, too
	std::map<hash_offset_type, value_type> old_value_map;
	old_value_map.swap(value_map);
	if (size_asked < 3) {	// to avoid collision_modulus == modulus
		size_asked = 3;
	}
	modulus = next_prime(size_asked);
	// collision_modulus just needs to be relatively prime with modulus;
	// since modulus is prime, any value will do - I made it prime for fun
	collision_modulus = next_prime(size_asked / 2);
	key_list = new data_offset_type[modulus];
	value_list = new small_value_type[modulus];
	// initialize keys and values
	for (hash_offset_type i(0); i < modulus; ++i) {
		key_list[i] = invalid_key;
	}
	for (hash_offset_type i(0); i < modulus; ++i) {
		value_list[i] = 0;
	}
	// copy over old hash keys and values
	key_type key(*this), comp_key(*this);
	for (hash_offset_type i(0); i < old_modulus; ++i) {
		if (old_key_list[i] != invalid_key) {
			key.copy_in(data, old_key_list[i]);
			comp_key.make_complement(key);
			const base_type key_hash(key < comp_key ? key.hash() : comp_key.hash());
			hash_offset_type new_i(key_hash % modulus);
			if (key_list[new_i] != invalid_key) {
				const hash_offset_type j(collision_modulus - key_hash % collision_modulus);
				do {
					new_i = (new_i + j) % modulus;
				} while (key_list[new_i] != invalid_key);
			}
			key_list[new_i] = old_key_list[i];
			value_list[new_i] = old_value_list[i];
			if (old_value_list[i] == max_small_value) {
				// use find() to avoid inserting a value into old_value_map
				const std::map<hash_offset_type, value_type>::const_iterator a(old_value_map.find(i));
				if (a != old_value_map.end()) {
					value_map[new_i] = a->second;
				}
			}
		}
	}
	delete[] old_key_list;
	delete[] old_value_list;
}
