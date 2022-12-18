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

void hashl::init(offset_type size_asked, const size_t bits_in, const base_type *data_in, const offset_type data_size_in) {
	bit_width = bits_in;
	data = data_in;
	data_size = data_size_in;
	word_width = (bit_width + 8 * sizeof(base_type) - 1) / (8 * sizeof(base_type));
	used_elements = 0;
	if (size_asked < 3) {	// to avoid collision_modulus == modulus
		size_asked = 3;
	}
	modulus = next_prime(size_asked);
	// collision_modulus just needs to be relatively prime with modulus;
	// since modulus is prime, any value will do - I made it prime for fun
	collision_modulus = next_prime(size_asked / 2);
	key_list = new offset_type[modulus];
	value_list = new small_value_type[modulus];
	// initialize keys and values
	for (offset_type i = 0; i != modulus; ++i) {
		key_list[i] = invalid_key;
	}
	for (offset_type i = 0; i != modulus; ++i) {
		value_list[i] = 0;
	}
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
	void *metadata_tmp = new char[metadata_size];
	pfread(fd, metadata_tmp, metadata_size);
	metadata = metadata_tmp;
	pfread(fd, &data_size, sizeof(data_size));
	base_type *data_tmp = new base_type[data_size];
	pfread(fd, data_tmp, sizeof(base_type) * data_size);
	data = data_tmp;
	value_list = new small_value_type[modulus];
	pfread(fd, value_list, sizeof(small_value_type) * modulus);
	key_list = new offset_type[modulus];
	for (offset_type i(0); i != modulus; ++i) {
		if (value_list[i] == 0) {
			key_list[i] = invalid_key;
		} else {
			pfread(fd, &key_list[i], sizeof(offset_type));
		}
	}
	// read in overflow map
	size_t x;
	pfread(fd, &x, sizeof(x));
	for (; x != 0; --x) {
		offset_type i;
		value_type j;
		pfread(fd, &i, sizeof(i));
		pfread(fd, &j, sizeof(j));
		value_map[i] = j;
	}
}

// insert a key at a particular location

hashl::offset_type hashl::insert_key(const offset_type i, const offset_type offset) {
	if (used_elements == modulus) {
		return modulus;		// hash table is full
	}
	++used_elements;
	key_list[i] = offset;
	value_list[i] = 0;	// init counts
	return i;
}

// create key from offset

void hashl::key_type::copy_in(const base_type * const data, const offset_type i) {
	offset_type j(i / (sizeof(base_type) * 8));
	// how many bits we have in the first word
	const base_type starting_bits(sizeof(base_type) * 8 - i % (sizeof(base_type) * 8));
	// how many bits the first word is supposed to have for a key
	const base_type high_offset(bit_width % (sizeof(base_type) * 8));
	const base_type high_mask(static_cast<base_type>(-1) >> (sizeof(base_type) * 8 - high_offset));
	if (starting_bits == high_offset) {
		size_t m(0);
		k[m] = data[j] & high_mask;
		for (++m, ++j; m < word_width; ++m, ++j) {
			k[m] = data[j];
		}
	} else if (starting_bits < high_offset) {	// shift left to fill up first word
		const int shift_left(high_offset - starting_bits);
		const int shift_right(sizeof(base_type) * 8 - shift_left);
		size_t m(0);
		k[m] = ((data[j] << shift_left) | (data[j + 1] >> shift_right)) & high_mask;
		for (++m, ++j; m < word_width; ++m, ++j) {
			k[m] = (data[j] << shift_left) | (data[j + 1] >> shift_right);
		}
	} else {			// shift right to empty out first word
		const int shift_right(starting_bits - high_offset);
		const int shift_left(sizeof(base_type) * 8 - shift_right);
		size_t m(0);
		k[m] = (data[j] >> shift_right) & high_mask;
		for (++m; m < word_width; ++m, ++j) {
			k[m] = (data[j] << shift_left) | (data[j + 1] >> shift_right);
		}
	}
}

// generate internal key from offset, compare to given key;
// equivalent to above subroutine, just with less storage and more breakpoints

bool hashl::key_equal(const offset_type i, const key_type &key) const {
	offset_type j(i / (sizeof(base_type) * 8));
	// how many bits we have in the first word
	const base_type starting_bits(sizeof(base_type) * 8 - i % (sizeof(base_type) * 8));
	// how many bits the first word is supposed to have for a key
	const base_type high_offset(bit_width % (sizeof(base_type) * 8));
	const base_type high_mask(static_cast<base_type>(-1) >> (sizeof(base_type) * 8 - high_offset));
	if (starting_bits == high_offset) {
		size_t m(0);
		if (key.k[m] != (data[j] & high_mask)) {
			return 0;
		}
		for (++m, ++j; m < word_width; ++m, ++j) {
			if (key.k[m] != data[j]) {
				return 0;
			}
		}
	} else if (starting_bits < high_offset) {	// shift left to fill up first word
		const int shift_left(high_offset - starting_bits);
		const int shift_right(sizeof(base_type) * 8 - shift_left);
		size_t m(0);
		if (key.k[m] != (((data[j] << shift_left) | (data[j + 1] >> shift_right)) & high_mask)) {
			return 0;
		}
		for (++m, ++j; m < word_width; ++m, ++j) {
			if (key.k[m] != ((data[j] << shift_left) | (data[j + 1] >> shift_right))) {
				return 0;
			}
		}
	} else {			// shift right to empty out first word
		const int shift_right(starting_bits - high_offset);
		const int shift_left(sizeof(base_type) * 8 - shift_right);
		size_t m(0);
		if (key.k[m] != ((data[j] >> shift_right) & high_mask)) {
			return 0;
		}
		for (++m; m < word_width; ++m, ++j) {
			if (key.k[m] != ((data[j] << shift_left) | (data[j + 1] >> shift_right))) {
				return 0;
			}
		}
	}
	return 1;
}

// find a key, or insert it if it doesn't exist; return modulus if hash is full

hashl::offset_type hashl::insert_offset(const key_type &key, const offset_type offset) {
	const base_type key_hash(key.hash());
	offset_type i(key_hash % modulus);
	if (key_list[i] == invalid_key) {	// insert
		return insert_key(i, offset);
	} else if (key_equal(i, key)) {		// already present
		return i;
	}
	const offset_type j(collision_modulus - key_hash % collision_modulus);
	for (;;) {	// search over all elements
		i = (i + j) % modulus;
		if (key_list[i] == invalid_key) {
			return insert_key(i, offset);
		} else if (key_equal(i, key)) {
			return i;
		}
	}
}

// find a key; return modulus as offset if not found

hashl::offset_type hashl::find_offset(const key_type &key) const {
	const base_type key_hash(key.hash());
	offset_type i(key_hash % modulus);
	if (key_list[i] == invalid_key) {
		return modulus;
	} else if (key_equal(i, key)) {
		return i;
	}
	const offset_type j(collision_modulus - key_hash % collision_modulus);
	for (;;) {	// search over all elements
		i = (i + j) % modulus;
		if (key_list[i] == invalid_key) {
			return modulus;
		} else if (key_equal(i, key)) {
			return i;
		}
	}
}

// increment the count for a key (but don't create a new entry if it doesn't exist)

bool hashl::increment(const key_type &key) {
	const offset_type i(find_offset(key));
	if (i == modulus) {	// insert failed
		return 0;
	}
	if (value_list[i] != max_small_value) {
		++value_list[i];
	} else {
		++value_map[i];
	}
	return 1;
}

bool hashl::increment(const key_type &key, const offset_type offset) {
	const offset_type i(insert_offset(key, offset));
	if (i == modulus) {	// insert failed
		return 0;
	}
	if (value_list[i] != max_small_value) {
		++value_list[i];
	} else {
		++value_map[i];
	}
	return 1;
}

// return the value associated with a key

hashl::value_type hashl::value(const key_type &key) const {
	const offset_type i(find_offset(key));
	if (i == modulus) {	// key not found
		return 0;
	} else if (value_list[i] != max_small_value) {
		return value_list[i];
	} else {
		// use find() to avoid inserting a value into value_map
		const std::map<offset_type, value_type>::const_iterator a(value_map.find(i));
		if (a == value_map.end()) {
			return max_small_value;
		} else {
			return a->second + max_small_value;
		}
	}
}

hashl::const_iterator hashl::begin() {
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

void hashl::const_iterator::get_value() {
	if (offset < list->modulus) {
		if (list->value_list[offset] != max_small_value) {
			value = list->value_list[offset];
		} else {
			const std::map<offset_type, value_type>::const_iterator b(list->value_map.find(offset));
			if (b == list->value_map.end()) {
				value = max_small_value;
			} else {
				value = b->second + max_small_value;
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
	pfwrite(fd, data, data_size);
	pfwrite(fd, value_list, sizeof(small_value_type) * modulus);
	for (offset_type i(0); i != modulus; ++i) {
		if (key_list[i] != invalid_key) {
			pfwrite(fd, &key_list[i], sizeof(offset_type));
		}
	}
	const size_t x(value_map.size());
	pfwrite(fd, &x, sizeof(x));
	std::map<offset_type, value_type>::const_iterator a(value_map.begin());
	const std::map<offset_type, value_type>::const_iterator end_a(value_map.end());
	for (; a != end_a; ++a) {
		pfwrite(fd, &a->first, sizeof(offset_type));
		pfwrite(fd, &a->second, sizeof(value_type));
	}
}

void hashl::set_metadata(const void * const data, const size_t data_size) {
	metadata = data;
	metadata_size = data_size;
}

void hashl::get_metadata(const void * &data, size_t &data_size) const {
	data = metadata;
	data_size = metadata_size;
}
