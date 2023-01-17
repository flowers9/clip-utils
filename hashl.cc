#include "hashl.h"
#include "itoa.h"	// itoa()
#include "local_endian.h"	// big_endian
#include "next_prime.h"	// next_prime()
#include "open_compressed.h"	// pfread()
#include "write_fork.h"	// pfwrite()
#include <map>		// map<>
#include <stdio.h>	// fprintf(), stderr
#include <stdlib.h>	// exit()
#include <string.h>	// memcmp(), memcpy()
#include <string>	// string
#include <vector>	// vector<>

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

void hashl::init(const hash_offset_type size_asked, const size_t bits_in, std::vector<base_type> &data_in) {
	bit_width = bits_in;
	word_width = (bit_width + 8 * sizeof(base_type) - 1) / (8 * sizeof(base_type));
	data.swap(data_in);
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
	uint64_t metadata_size;
	pfread(fd, &metadata_size, sizeof(metadata_size));
	metadata.assign(metadata_size, 0);
	pfread(fd, &metadata[0], metadata_size);
	uint64_t data_size;
	pfread(fd, &data_size, sizeof(data_size));
	data.assign(data_size, 0);
	pfread(fd, &data[0], sizeof(base_type) * data_size);
	value_list.assign(modulus, 0);
	pfread(fd, &value_list[0], sizeof(small_value_type) * modulus);
	key_list.assign(modulus, invalid_key);
	for (hash_offset_type i(0); i < modulus; ++i) {
		if (value_list[i]) {
			pfread(fd, &key_list[i], sizeof(data_offset_type));
		}
	}
}

// insert a key at a particular location

hashl::hash_offset_type hashl::insert_key(const hash_offset_type i, const data_offset_type offset) {
	if (used_elements == modulus) {
		return modulus;		// hash table is full
	}
	++used_elements;
	key_list[i] = offset;
	value_list[i] = 0;
	return i;
}

// create key from bit offset into data

void hashl::key_type::copy_in(const std::vector<base_type> &data, const data_offset_type offset) {
	// start of sequence in data
	const size_t i(offset / (sizeof(base_type) * 8));
	// how many bits we have in the first word
	const base_type starting_bits(sizeof(base_type) * 8 - offset % (sizeof(base_type) * 8));
	// how many bits the first word is supposed to have for a key
	const base_type high_bits(bit_shift + 2);
	if (starting_bits == high_bits) {
		k[0] = data[i] & high_mask;
		for (size_t j(1); j < word_width; ++j) {
			k[j] = data[i + j];
		}
	} else if (starting_bits < high_bits) {		// shift left to fill up first word
		const int shift_left(high_bits - starting_bits);
		const int shift_right(sizeof(base_type) * 8 - shift_left);
		k[0] = ((data[i] << shift_left) | (data[i + 1] >> shift_right)) & high_mask;
		for (size_t j(1); j < word_width; ++j) {
			k[j] = (data[i + j] << shift_left) | (data[i + j + 1] >> shift_right);
		}
	} else {					// shift right to empty out first word
		const int shift_right(starting_bits - high_bits);
		const int shift_left(sizeof(base_type) * 8 - shift_right);
		k[0] = (data[i] >> shift_right) & high_mask;
		for (size_t j(1); j < word_width; ++j) {
			k[j] = (data[i + j - 1] << shift_left) | (data[i + j] >> shift_right);
		}
	}
}

// generate internal key from bit offset into data, compare to key;
// equivalent to above subroutine, just with more breakpoints

bool hashl::key_type::equal(const std::vector<base_type> &data, const data_offset_type offset) const {
	// start of sequence in data
	const size_t i(offset / (sizeof(base_type) * 8));
	// how many bits we have in the first word
	const base_type starting_bits(sizeof(base_type) * 8 - offset % (sizeof(base_type) * 8));
	// how many bits the first word is supposed to have for a key
	const base_type high_offset(bit_shift + 2);
	if (starting_bits == high_offset) {
		if (k[0] != (data[i] & high_mask)) {
			return 0;
		}
		for (size_t j(1); j < word_width; ++j) {
			if (k[j] != data[i + j]) {
				return 0;
			}
		}
	} else if (starting_bits < high_offset) {	// shift left to fill up first word
		const int shift_left(high_offset - starting_bits);
		const int shift_right(sizeof(base_type) * 8 - shift_left);
		if (k[0] != (((data[i] << shift_left) | (data[i + 1] >> shift_right)) & high_mask)) {
			return 0;
		}
		for (size_t j(1); j < word_width; ++j) {
			if (k[j] != ((data[i + j] << shift_left) | (data[i + j + 1] >> shift_right))) {
				return 0;
			}
		}
	} else {					// shift right to empty out first word
		const int shift_right(starting_bits - high_offset);
		const int shift_left(sizeof(base_type) * 8 - shift_right);
		if (k[0] != ((data[i] >> shift_right) & high_mask)) {
			return 0;
		}
		for (size_t j(1); j < word_width; ++j) {
			if (k[j] != ((data[i + j - 1] << shift_left) | (data[i + j] >> shift_right))) {
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
	}
	return 1;
}

// return the value associated with a key (or zero if key not found)

hashl::value_type hashl::value(const key_type &key) const {
	const hash_offset_type i(find_offset(key));
	return i < modulus ? value_list[i] : 0;
}

hashl::const_iterator hashl::begin() const {
	if (used_elements == 0) {
		return end();
	}
	const_iterator a(*this, 0);
	// advance to first valid value
	if (key_list[0] == invalid_key) {
		++a;
	}
	return a;
}

hashl::const_iterator hashl::end() const {
	return const_iterator(*this, modulus);
}

void hashl::const_iterator::get_value() {
	value = offset < list.modulus ? list.value_list[offset] : 0;
}

void hashl::save(const int fd) const {
	const std::string s(boilerplate());
	pfwrite(fd, s.c_str(), s.size());
	pfwrite(fd, &modulus, sizeof(modulus));
	pfwrite(fd, &collision_modulus, sizeof(collision_modulus));
	pfwrite(fd, &used_elements, sizeof(used_elements));
	pfwrite(fd, &bit_width, sizeof(bit_width));
	uint64_t tmp;
	pfwrite(fd, &(tmp = metadata.size()), sizeof(tmp));
	pfwrite(fd, &metadata[0], metadata.size());
	pfwrite(fd, &(tmp = data.size()), sizeof(tmp));
	pfwrite(fd, &data[0], sizeof(base_type) * data.size());
	pfwrite(fd, &value_list[0], sizeof(small_value_type) * modulus);
	for (hash_offset_type i(0); i < modulus; ++i) {
		if (key_list[i] != invalid_key) {
			pfwrite(fd, &key_list[i], sizeof(data_offset_type));
		}
	}
}

void hashl::resize(hash_offset_type size_asked) {
	const size_t old_modulus(modulus);
	if (size_asked < 3) {	// to avoid collision_modulus == modulus
		size_asked = 3;
	}
	modulus = next_prime(size_asked);
	// collision_modulus just needs to be relatively prime with modulus;
	// since modulus is prime, any value will do - I made it prime for fun
	collision_modulus = next_prime(size_asked / 2);
	// initialize keys and values
	std::vector<data_offset_type> old_key_list(modulus, invalid_key);
	key_list.swap(old_key_list);
	std::vector<small_value_type> old_value_list(modulus, 0);
	value_list.swap(old_value_list);
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
		}
	}
}

// identical in effect to add()ing this hash to an empty hash
// (except that values <min_cutoff are set to zero rather than removed)

void hashl::normalize(const small_value_type min_cutoff, const small_value_type max_cutoff) {
	for (hash_offset_type i(0); i < modulus; ++i) {
		if (value_list[i] < min_cutoff) {
			value_list[i] = 0;
		} else if (value_list[i] > max_cutoff) {
			value_list[i] = invalid_value;
		} else {
			value_list[i] = 1;
		}
	}
}

// add in new hash: any values <min_cutoff ignored, <=max_cutoff increment existing value,
// values of >max_cutoff set to invalid_value

bool hashl::add(const hashl &a, const small_value_type min_cutoff, const small_value_type max_cutoff) {
	if (used_elements + a.used_elements > a.modulus) {
		resize(used_elements + a.used_elements);
	}
	// copy over data (and make sure to update offsets when adding new entries)
	const size_t offset(data.size() * sizeof(base_type) * 8);
	// this pads out the existing data so we don't have to shift all the new data
	// TODO: see if shifting would actually slow things down much
	data.reserve(data.size() + a.data.size());
	data.insert(data.end(), a.data.begin(), a.data.end());
	// loop over incoming hash and increment or invalidate entries as needed
	key_type key(a), comp_key(a);
	for (size_t i(0); i < a.modulus; ++i) {
		if (a.key_list[i] != invalid_key) {
			key.copy_in(a.data, a.key_list[i]);
			comp_key.make_complement(key);
			const hash_offset_type new_i(insert_offset(key, comp_key, a.key_list[i] + offset));
			if (new_i == modulus) {
				return 0;
			} else if (a.value_list[i] < min_cutoff) {
			} else if (a.value_list[i] > max_cutoff) {
				value_list[new_i] = invalid_value;
			} else if (value_list[new_i] < max_small_value) {
				++value_list[new_i];
			}
		}
	}
	return 1;
}
