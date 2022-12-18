#include "hash.h"
#include "itoa.h"	// itoa()
#include "local_endian.h"	// big_endian
#include "next_prime.h"	// next_prime()
#include "open_compressed.h"	// close_compressed(), open_compressed(), pfread()
#include "write_fork.h"	// close_fork(), close_fork_wait(), pfwrite(), write_fork()
#include <algorithm>	// swap()
#include <cassert>	// assert()
#include <errno.h>	// errno
#include <limits.h>	// UINT64_MAX
#include <list>		// list<>
#include <map>		// map<>
#include <new>		// new
#include <sstream>	// ostringstream
#include <stdio.h>	// fprintf(), stderr
#include <stdlib.h>	// exit()
#include <string.h>	// memcmp(), strerror()
#include <string>	// string
#include <sys/stat.h>	// S_ISDIR(), stat(), struct stat
#include <unistd.h>	// unlink()
#include <utility>	// make_pair(), pair<>

hash::~hash() {
	delete[] key_list;
	delete[] value_list;
	delete[] alt_list;
	delete[] alt_map;
	std::list<std::string>::const_iterator a(state_files.begin());
	const std::list<std::string>::const_iterator end_a(state_files.end());
	for (; a != end_a; ++a) {
		unlink(a->c_str());
	}
}

// description beginning of saved file

std::string hash::boilerplate() const {
	std::string s("hash\n");
	s += itoa(sizeof(key_type));
	s += " bytes\n";
#ifdef big_endian
	s += "big endian\n";
#else
	s += "little endian\n";
#endif
	return s;
}

void hash::init(offset_type size_asked, offset_type alt_size_in) {
	if (alt_size_in > 8 * sizeof(offset_type)) {
		fprintf(stderr, "Error: hash alt size too large: %lu > %lu\n", alt_size_in, 8 * sizeof(offset_type));
		exit(1);
	}
	alt_size = alt_size_in;
	used_elements = 1;	// to account for minimum of one INVALID_KEYs
	++size_asked;
	if (size_asked < 3) {	// to avoid collision_modulus == modulus
		size_asked = 3;
	}
	modulus = next_prime(size_asked);
	// collision_modulus just needs to be relatively prime with modulus;
	// since modulus is prime, any value will do - I made it prime for fun
	collision_modulus = next_prime(size_asked / 2);
	key_list = new key_type[modulus];
	value_list = new small_value_type[modulus];
	if (alt_size == 0) {
		alt_list = NULL;
		alt_map = NULL;
	} else {
		alt_list = new small_value_type[modulus * alt_size];
		alt_map = new std::map<key_type, value_type>[alt_size];
	}
	// initialize keys; values are initialized as keys are entered
	for (offset_type i = 0; i != modulus; ++i) {
		key_list[i] = INVALID_KEY;
	}
}

void hash::init_from_file(const int fd) {
	const std::string s(boilerplate());
	char t[s.size()];
	pfread(fd, t, s.size());
	if (memcmp(s.c_str(), t, s.size()) != 0) {
		fprintf(stderr, "Error: could not read hash from file: header mismatch\n");
		exit(1);
	}
	pfread(fd, &modulus, sizeof(modulus));
	pfread(fd, &collision_modulus, sizeof(collision_modulus));
	pfread(fd, &used_elements, sizeof(used_elements));
	pfread(fd, &alt_size, sizeof(alt_size));
	key_list = new key_type[modulus];
	value_list = new small_value_type[modulus];
	if (alt_size == 0) {
		alt_list = NULL;
		alt_map = NULL;
	} else {
		alt_list = new small_value_type[modulus * alt_size];
		alt_map = new std::map<key_type, value_type>[alt_size];
	}
	// read in values (they're the smallest size)
	pfread(fd, &value_list[0], sizeof(small_value_type) * modulus);
	// read in keys for non-zero values
	for (offset_type i(0); i != modulus; ++i) {
		if (value_list[i] == 0) {
			key_list[i] = INVALID_KEY;
		} else {
			pfread(fd, &key_list[i], sizeof(key_type));
		}
	}
	// read in overflow map
	offset_type x;
	pfread(fd, &x, sizeof(x));
	for (; x != 0; --x) {
		key_type i;
		value_type j;
		pfread(fd, &i, sizeof(i));
		pfread(fd, &j, sizeof(j));
		value_map[i] = j;
	}
	if (alt_size != 0) {
		for (offset_type i(0), j(0); i != modulus; ++i) {
			if (value_list[i] == 0) {
				const offset_type end_j(j + alt_size);
				for (; j != end_j; ++j) {
					alt_list[j] = 0;
				}
			} else {
				pfread(fd, &alt_list[j], sizeof(small_value_type) * alt_size);
				j += alt_size;
			}
		}
		// alt map overflows
		for (offset_type k(0); k != alt_size; ++k) {
			std::map<key_type, value_type> &z = alt_map[k];
			pfread(fd, &x, sizeof(x));
			for (; x != 0; --x) {
				key_type i;
				value_type j;
				pfread(fd, &i, sizeof(i));
				pfread(fd, &j, sizeof(j));
				z[i] = j;
			}
		}
	}
}

// returns next empty spot found, or modulus if it spots the key first

hash::offset_type hash::find_empty_offset(key_type key) const {
	offset_type i(key % modulus);
	if (key_list[i] == INVALID_KEY) {
		return i;
	} else if (key_list[i] == key) {
		return modulus;
	}
	const offset_type j(collision_modulus - (key % collision_modulus));
	for (;;) {	// search over all elements
		i = (i + j) % modulus;
		if (key_list[i] == INVALID_KEY) {
			return i;
		} else if (key_list[i] == key) {
			return modulus;
		}
	}
}

void hash::rehash() {
	// first pass, fill all non-collision slots
	for (offset_type i(0); i != modulus; ++i) {
		if (key_list[i] != INVALID_KEY) {
			const offset_type j(key_list[i] % modulus);
			// swap if non-collision slot has a collision fill
			if (i != j && (key_list[j] == INVALID_KEY || key_list[j] % modulus != j)) {
				std::swap(key_list[j], key_list[i]);
				std::swap(value_list[j], value_list[i]);
				--i;
			}
		}
	}
	// subsequent passes, shift collisions to empty slots; the shifts
	// are entropy lowering operations, meaning there is no risk of
	// an infinite loop
	bool changed;
	do {
		changed = 0;
		for (offset_type i(0); i != modulus; ++i) {
			if (key_list[i] != INVALID_KEY) {
				const offset_type j(find_empty_offset(key_list[i]));
				if (j != modulus) {
					changed = 1;
					key_list[j] = key_list[i];
					key_list[i] = INVALID_KEY;
					value_list[j] = value_list[i];
				}
			}
		}
	} while (changed);
}

// same as rehash(), but moves alt_list as well

void hash::rehash_alt() {
	// first pass, fill all non-collision slots
	offset_type i(0);
	for (; i != modulus; ++i) {
		if (key_list[i] != INVALID_KEY) {
			const offset_type j(key_list[i] % modulus);
			// swap if non-collision slot has a collision fill
			if (i != j && (key_list[j] == INVALID_KEY || key_list[j] % modulus != j)) {
				std::swap(key_list[j], key_list[i]);
				std::swap(value_list[j], value_list[i]);
				offset_type x(j * alt_size);
				offset_type y(i * alt_size);
				const offset_type end_x(x + alt_size);
				for (; x != end_x; ++x, ++y) {
					std::swap(alt_list[x], alt_list[y]);
				}
				--i;
			}
		}
	}
	// subsequent passes, shift collisions to empty slots; the shifts
	// are entropy lowering operations, meaning there is no risk of
	// an infinite loop
	bool changed;
	do {
		changed = 0;
		for (i = 0; i != modulus; ++i) {
			if (key_list[i] != INVALID_KEY) {
				const offset_type j(find_empty_offset(key_list[i]));
				if (j != modulus) {
					changed = 1;
					key_list[j] = key_list[i];
					key_list[i] = INVALID_KEY;
					value_list[j] = value_list[i];
					offset_type x(j * alt_size);
					offset_type y(i * alt_size);
					const offset_type end_x(x + alt_size);
					for (; x != end_x; ++x, ++y) {
						alt_list[x] = alt_list[y];
					}
				}
			}
		}
	} while (changed);
}

// remove keys with a value of 1

bool hash::clean_hash(void) {
	for (offset_type i(0); i != modulus; ++i) {
		if (key_list[i] != INVALID_KEY && value_list[i] == 1) {
			key_list[i] = INVALID_KEY;
			--used_elements;
		}
	}
	if (used_elements == modulus) {
		return 0;
	} else if (alt_list == NULL) {
		rehash();
	} else {
		rehash_alt();
	}
	return 1;
}

// remove all keys with value < min || max < value; if min or max is zero,
// they're ignored for the comparison

void hash::clean_hash(value_type min, value_type max) {
	if (min == 0 && max == 0) {
		return;
	} else if (min != 0 && max != 0 && max < min) {
		clear();
		return;
	}
	const offset_type starting_used_elements(used_elements);
	const bool big_max(max > max_small_value - 1);
	small_value_type s_min, s_max;
	if (max == 0) {
		if (min < max_small_value) {
			s_max = max_small_value;
		} else {	// to trigger large value comparisons
			s_max = max_small_value - 1;
		}
	} else if (big_max) {
		s_max = max_small_value - 1;
		max -= max_small_value;
	} else {
		s_max = max;
		max = 0;
	}
	if (min > max_small_value) {
		s_min = max_small_value;
		min -= max_small_value;
	} else {
		s_min = min;
		min = 0;
	}
	for (offset_type i(0); i != modulus; ++i) {
		if (key_list[i] == INVALID_KEY) {
		} else if (s_min <= value_list[i] && value_list[i] <= s_max) {
		} else if (value_list[i] != max_small_value) {
			key_list[i] = INVALID_KEY;
			--used_elements;
		} else {
			const std::map<key_type, value_type>::iterator a(value_map.find(key_list[i]));
			const value_type big_value(a == value_map.end() ? 0 : a->second);
			if (big_value < min || (big_max && max < big_value)) {
				if (a != value_map.end()) {
					value_map.erase(a);
				}
				key_list[i] = INVALID_KEY;
				--used_elements;
			}
		}
	}
	if (used_elements == starting_used_elements) {
		// no changes, so no need to rehash
	} else if (alt_list == NULL) {
		rehash();
	} else {
		rehash_alt();
	}
}

// insert a key at a particular location

hash::offset_type hash::insert_key(const offset_type i, const key_type key) {
	if (used_elements == modulus) {
		if ((no_space_response & CLEAN_HASH) && clean_hash()) {
			// have to redo positioning after clean_hash()
			return insert_offset(key);
		} else if (no_space_response & TMP_FILE) {
			radix_sort(modulus);
			save_state();
			clear(1);
			// have to redo positioning after clear()
			return insert_offset(key);
		} else {
			return modulus;		// hash table is full
		}
	}
	++used_elements;
	key_list[i] = key;
	value_list[i] = 0;
	offset_type j(i * alt_size);
	const offset_type end_j(j + alt_size);
	for (; j != end_j; ++j) {
		alt_list[j] = 0;
	}
	return i;
}

// find a key, or insert it if it doesn't exist; return modulus if hash is full

hash::offset_type hash::insert_offset(key_type key) {
	offset_type i(key % modulus);
	if (key_list[i] == INVALID_KEY) {
		return insert_key(i, key);
	} else if (key_list[i] == key) {
		return i;
	}
	const offset_type j(collision_modulus - (key % collision_modulus));
	for (;;) {	// search over all elements
		i = (i + j) % modulus;
		if (key_list[i] == INVALID_KEY) {
			return insert_key(i, key);
		} else if (key_list[i] == key) {
			return i;
		}
	}
}

// find a key; return modulus if not found

hash::offset_type hash::find_offset(const key_type key) const {
	offset_type i(key % modulus);
	if (key_list[i] == key) {
		return i;
	} else if (key_list[i] == INVALID_KEY) {
		return modulus;
	}
	const offset_type j(collision_modulus - (key % collision_modulus));
	for (;;) {	// search over all elements
		i = (i + j) % modulus;
		if (key_list[i] == key) {
			return i;
		} else if (key_list[i] == INVALID_KEY) {
			return modulus;
		}
	}
}

// increment the count for a key

bool hash::increment(key_type key) {
	const offset_type i(insert_offset(key));
	if (i == modulus) {	// insert failed
		return 0;
	}
	if (value_list[i] != max_small_value) {
		++value_list[i];
	} else if (can_overflow) {
		++value_map[key];
	}
	return 1;
}

// increment only the alt values, using x as a bit flag to mark which ones

bool hash::increment_alt(key_type key, offset_type x) {
	const offset_type i(insert_offset(key));
	if (i == modulus) {	// insert failed
		return 0;
	}
	const offset_type start_j(i * alt_size);
	const offset_type end_j(start_j + alt_size);
	for (offset_type j(start_j), z(1); j != end_j; ++j, z <<= 1) {
		if (x & z) {
			if (alt_list[j] != max_small_value) {
				++alt_list[j];
			} else if (can_overflow) {
				++alt_map[j - start_j][key];
			}
		}
	}
	return 1;
}

// return the value associated with a key

hash::value_type hash::value(key_type key) const {
	const offset_type i(find_offset(key));
	if (i == modulus) {	// key not found
		return 0;
	} else if (value_list[i] != max_small_value) {
		return value_list[i];
	} else {
		// use find() to avoid inserting a value into value_map
		const std::map<key_type, value_type>::const_iterator a(value_map.find(key));
		if (a == value_map.end()) {
			return max_small_value;
		} else {
			return a->second + max_small_value;
		}
	}
}

hash::value_type hash::value(key_type key, value_type x[]) const {
	const offset_type i(find_offset(key));
	if (i == modulus) {	// key not found
		return 0;
	}
	offset_type j(0);
	const offset_type alt_offset(i * alt_size);
	for (; j != alt_size; ++j) {
		if (alt_list[j + alt_offset] != max_small_value) {
			x[j] = alt_list[j + alt_offset];
		} else {
			const std::map<key_type, value_type> &b(alt_map[j]);
			// use find() to avoid inserting a value into alt_map
			const std::map<key_type, value_type>::const_iterator a(b.find(key));
			if (a == b.end()) {
				x[j] = max_small_value;
			} else {
				x[j] = a->second + max_small_value;
			}
		}
	}
	if (value_list[i] != max_small_value) {
		return value_list[i];
	} else {
		// use find() to avoid inserting a value into value_map
		const std::map<key_type, value_type>::const_iterator a(value_map.find(key));
		if (a == value_map.end()) {
			return max_small_value;
		} else {
			return a->second + max_small_value;
		}
	}
}

// reset hash to empty state

void hash::clear(const bool mostly_clear) {
	used_elements = 1;	// to account for minimum of one INVALID_KEYs
	// initialize keys
	for (offset_type i(0); i != modulus; ++i) {
		key_list[i] = INVALID_KEY;
	}
	value_map.clear();
	for (offset_type i(0); i != alt_size; ++i) {
		alt_map[i].clear();
	}
	if (!mostly_clear) {
		std::list<std::string>::const_iterator a(state_files.begin());
		const std::list<std::string>::const_iterator end_a(state_files.end());
		for (; a != end_a; ++a) {
			unlink(a->c_str());
		}
		state_files.clear();
	}
}

hash::const_iterator hash::begin() {
	if (state_files.empty()) {
		if (used_elements == 1) {	// only key is INVALID_KEY
			return end();
		}
		const_iterator a(this, 0);
		if (a.key == INVALID_KEY) {	// advance to first valid value
			++a;
		}
		return a;
	} else {			// set up readback from tmp files
		offset_type index(static_cast<offset_type>(-1));
		std::map<key_type, std::pair<value_type, int> > next_keys;
		prep_for_readback(index, next_keys);
		return const_iterator(this, index, next_keys);
	}
}

hash::const_iterator::const_iterator(const hash * const b, const offset_type i) : list(b), offset(i) {
	if (list == NULL || offset == list->modulus) {
		key = INVALID_KEY;
		value = 0;
	} else {
		key = list->key_list[offset];
		if (list->value_list[offset] != max_small_value) {
			value = list->value_list[offset];
		} else {
			const std::map<key_type, value_type>::const_iterator a(list->value_map.find(key));
			if (a == list->value_map.end()) {
				value = max_small_value;
			} else {
				value = a->second + max_small_value;
			}
		}
	}
}

hash::const_iterator::const_iterator(const hash * const b, const offset_type i, const std::map<key_type, std::pair<value_type, int> > &c) : list(b), offset(i), next_keys(c) {
	if (list == NULL || next_keys.empty()) {
		key = INVALID_KEY;
		value = 0;
	} else {
		const std::map<key_type, std::pair<value_type, int> >::iterator a(next_keys.begin());
		key = a->first;
		value = a->second.first;
	}
}

hash::const_iterator &hash::const_iterator::operator=(const const_iterator &a) {
	list = a.list;
	offset = a.offset;
	key = a.key;
	value = a.value;
	next_keys = a.next_keys;
	return *this;
}

// advance to next used element

void hash::const_iterator::increment() {
	if (offset == list->modulus) {
		return;
	}
	if (next_keys.empty()) {
		for (++offset; offset != list->modulus && list->key_list[offset] == INVALID_KEY; ++offset) { }
		if (offset < list->modulus) {
			key = list->key_list[offset];
			value = list->value_list[offset];
			if (value == max_small_value) {
				const std::map<key_type, value_type>::const_iterator a(list->value_map.find(key));
				if (a != list->value_map.end()) {
					value += a->second;
				}
			}
		} else {
			key = INVALID_KEY;
			value = 0;
		}
	} else {
		const int fd(next_keys.begin()->second.second);
		next_keys.erase(next_keys.begin());
		// get next unique entry from that file
		for (;;) {
			// check if file is out of entries
			if (!list->get_next_entry(fd, key, value, offset)) {
				if (fd != -1) {
					close_compressed(fd);
				}
				break;
			}
			std::map<key_type, std::pair<value_type, int> >::iterator c(next_keys.find(key));
			if (c == next_keys.end()) {
				next_keys[key] = std::make_pair(value, fd);
				break;
			} else { // keep going until we get a unique entry
				c->second.first += value;
			}
		}
		if (!next_keys.empty()) {
			const std::map<key_type, std::pair<value_type, int> >::iterator a(next_keys.begin());
			key = a->first;
			value = a->second.first;
		} else {
			key = INVALID_KEY;
			value = 0;
			// mark end of iteration; can't happen through
			// get_next_entry() because that maxes offset out
			// at list->used_elements, which is guaranteed to
			// be at least one less than list->modulus, thanks
			// to squash_hash()
			offset = list->modulus;
		}
	}
}

// extract list of alt_values associated with current offset

void hash::const_iterator::get_alt_values(value_type x[]) const {
	const offset_type alt_offset(offset * list->alt_size);
	for (offset_type i(0); i != list->alt_size; ++i) {
		if (list->alt_list[alt_offset + i] != max_small_value) {
			x[i] = list->alt_list[alt_offset + i];
		} else {
			// use find() to avoid inserting a value into alt_map
			const std::map<key_type, value_type>::const_iterator a(list->alt_map[i].find(key));
			if (a == list->alt_map[i].end()) {
				x[i] = max_small_value;
			} else {
				x[i] = a->second + max_small_value;
			}
		}
	}
}

void hash::save(const int fd) const {
	const std::string s(boilerplate());
	pfwrite(fd, s.c_str(), s.size());
	pfwrite(fd, &modulus, sizeof(modulus));
	pfwrite(fd, &collision_modulus, sizeof(collision_modulus));
	pfwrite(fd, &used_elements, sizeof(used_elements));
	pfwrite(fd, &alt_size, sizeof(alt_size));
	// normally, unused entries have uninitialized data, but here we
	// store them with zero values so we can use them to mark entries
	// in the other arrays that we don't bother to write to file
	const small_value_type zero(0);
	for (offset_type i(0); i != modulus; ++i) {
		pfwrite(fd, key_list[i] == INVALID_KEY ? &zero : &value_list[i], sizeof(small_value_type));
	}
	for (offset_type i(0); i != modulus; ++i) {
		if (key_list[i] != INVALID_KEY) {
			pfwrite(fd, &key_list[i], sizeof(key_type));
		}
	}
	offset_type x;
	x = value_map.size();
	pfwrite(fd, &x, sizeof(x));
	std::map<key_type, value_type>::const_iterator a(value_map.begin());
	const std::map<key_type, value_type>::const_iterator end_a(value_map.end());
	for (; a != end_a; ++a) {
		pfwrite(fd, &a->first, sizeof(key_type));
		pfwrite(fd, &a->second, sizeof(value_type));
	}
	if (alt_size != 0) {
		for (offset_type i(0); i != modulus; ++i) {
			if (key_list[i] != INVALID_KEY) {
				pfwrite(fd, &alt_list[i * alt_size], sizeof(small_value_type) * alt_size);
			}
		}
		// alt map overflows
		for (offset_type j(0); j != alt_size; ++j) {
			std::map<key_type, value_type> &z = alt_map[j];
			x = z.size();
			pfwrite(fd, &x, sizeof(x));
			std::map<key_type, value_type>::const_iterator b(z.begin());
			const std::map<key_type, value_type>::const_iterator end_b(z.end());
			for (; b != end_b; ++b) {
				pfwrite(fd, &b->first, sizeof(key_type));
				pfwrite(fd, &b->second, sizeof(value_type));
			}
		}
	}
}

// add value to possibly existing entry

bool hash::add(const key_type key, const value_type new_value) {
	const offset_type i(insert_offset(key));
	if (i == modulus) {	// insert failed
		return 0;
	}
	if (value_list[i] + new_value <= max_small_value) {
		value_list[i] += new_value;
	} else if (!can_overflow) {
		value_list[i] = max_small_value;
	} else if (value_list[i] != max_small_value) {
		value_map[key] = value_list[i] + new_value - max_small_value;
		value_list[i] = max_small_value;
	} else {
		value_map[key] += new_value;
	}
	return 1;
}

// add value and alt values to possibly existing entry

bool hash::add_alt(const key_type key, const value_type new_value, const value_type alt_values[]) {
	const offset_type i(insert_offset(key));
	if (i == modulus) {	// insert failed
		return 0;
	}
	if (value_list[i] + new_value <= max_small_value) {
		value_list[i] += new_value;
	} else if (!can_overflow) {
		value_list[i] = max_small_value;
	} else if (value_list[i] != max_small_value) {
		value_map[key] = value_list[i] + new_value - max_small_value;
		value_list[i] = max_small_value;
	} else {
		value_map[key] += new_value;
	}
	const offset_type start_j(i * alt_size);
	const offset_type end_j(start_j + alt_size);
	for (offset_type j(start_j); j != end_j; ++j) {
		if (alt_list[j] + alt_values[j - start_j] <= max_small_value) {
			alt_list[j] += alt_values[j - start_j];
		} else if (!can_overflow) {
			alt_list[j] = max_small_value;
		} else if (alt_list[j] != max_small_value) {
			alt_map[j - start_j][key] = alt_list[j] + alt_values[j - start_j] - max_small_value;
			alt_list[j] = max_small_value;
		} else {
			alt_map[j - start_j][key] += alt_values[j - start_j];
		}
	}
	return 1;
}

// add values from hash table to this hash; returns if hashes are
// successfully added (failure is generally because of lack of space)

bool hash::add_hash(hash &h) {
	// see if hashes are compatible
	if (alt_size != 0 && h.alt_size != 0 && alt_size != h.alt_size) {
		fprintf(stderr, "Error: cannot add hashes: different size alt arrays\n");
		exit(1);
	}
	if (alt_size == 0) {
		const_iterator a(h.begin());
		const const_iterator end_a(h.end());
		for (; a != end_a; ++a) {
			if (!add(a.key, a.value)) {
				return 0;
			}
		}
	} else {
		value_type alt_values[alt_size];
		const_iterator a(h.begin());
		const const_iterator end_a(h.end());
		for (; a != end_a; ++a) {
			a.get_alt_values(alt_values);
			if (!add_alt(a.key, a.value, alt_values)) {
				return 0;
			}
		}
	}
	return 1;
}

// a shell sort does an insertion sort on an array over a range of gap sizes;
// insert sorts are inefficient as they can involve shifting many values,
// but they can be good for small arrays and mostly ordered arrays; the
// shell sort does insert sorts on spaced sets of increasing size so that
// as the arrays get larger, they are also more ordered

// optimized (-O2) shell_sort() beats qsort() for n < ~32k (on random data)

void hash::shell_sort(const offset_type start_index, const offset_type stop_index) {
	// if you plan to only use arrays shorter than 701 (or some
	// smaller value), reduce the array to only the ones that cover
	// your size range (701, 301, 132)
	const static offset_type gaps[] = { 57, 23, 10, 4, 1 };
	const static int gap_count(sizeof(gaps) / sizeof(offset_type));
	for (int i(0); i != gap_count; ++i) {
		const offset_type gap(gaps[i]);
		const offset_type start_index_gap(start_index + gap);
		// do an insertion sort with the given spacing (this will
		// skip gap sizes larger than stop_index - start_index)
		for (offset_type j(start_index_gap); j < stop_index; ++j) {
			if (key_list[j] < key_list[j - gap]) {
				const key_type my_key(key_list[j]);
				const small_value_type my_value(value_list[j]);
				offset_type k(j);
				for (; k >= start_index_gap && my_key < key_list[k - gap]; k -= gap) {
					key_list[k] = key_list[k - gap];
					value_list[k] = value_list[k - gap];
				}
				key_list[k] = my_key;
				value_list[k] = my_value;
			}
		}
	}
}

void hash::calculate_offsets(const offset_type start_index, const offset_type stop_index, offset_type * offsets, const int shift) const {
	// we actually count in offsets[1-256], but we don't care about
	// the value in offsets[256], so don't bother initializing it
	offsets[0] = start_index;
	for (int i(1); i != 256; ++i) {		// initialize array
		offsets[i] = 0;
	}
	++offsets;	// offset[0] will always be start_index, and we always
			// have extra space at the end of the array,
			// so we won't overrun
	for (offset_type i(start_index); i != stop_index; ++i) {	// count
		++offsets[(key_list[i] >> shift) & 255];
	}
	--offsets;
	// as what we actually want are starting offsets into the
	// destination array, accumulate the totals, so offset[n] points
	// to the starting location for bin n
	for (int i(1); i != 256; ++i) {
		offsets[i] += offsets[i - 1];
	}
}

void hash::radix_sort_internal(const offset_type start_index, const offset_type stop_index, offset_type * offsets, const int shift) {
	// use shell sort for terminal sorting;
	// this seems to give very similar results when set between 128 and 4k
	if (stop_index - start_index < 512) {
		shell_sort(start_index, stop_index);
		return;
	}
	calculate_offsets(start_index, stop_index, offsets, shift);
	// swap unbinned elements into bins, growing sorted bins as we go
	offset_type * const unbinned_start(offsets + 256);
	for (int i(0); i != 256; ++i) {
		// -1 here so we can use prefix ++ instead of suffix ++ below
		unbinned_start[i] = offsets[i] - 1;
	}
	++offsets;	// shift offsets to be end of bins, instead of start;
	// we never sort the last bin, as it's guaranteed to be full by the
	// time we get to it, so no worries about it never having been set
	for (int i(0); i != 255; ++i) {
		// iterate over the unsorted portion of each bin
		offset_type j(unbinned_start[i] + 1);
		while (j != offsets[i]) {
			const int my_bin((key_list[j] >> shift) & 255);
			if (my_bin != i) {
				std::swap(key_list[++unbinned_start[my_bin]], key_list[j]);
				std::swap(value_list[unbinned_start[my_bin]], value_list[j]);
			} else {	// already in place!
				++j;
				// don't need to update unbinned_start[i],
				// as we won't use it again
			}
		}
	}
	if (shift == 0) {	// we're done!
		return;
	}
	--offsets;
	// recurse down to the next set of bins
	for (int i(0); i != 255; ++i) {
		if (offsets[i] != offsets[i + 1]) {
			radix_sort_internal(offsets[i], offsets[i + 1], offsets + 256, shift - 8);
		}
	}
	if (offsets[255] != stop_index) {
		radix_sort_internal(offsets[255], stop_index, offsets + 256, shift - 8);
	}
}

// msb in-place radix sort with one-byte radix (256);
// appears to be about twice as fast as qsort on 1g random data (although that
// was for just the key_value array, not with the small_values_array as well);
// of course, the hash is unusable as a hash after this operation

void hash::radix_sort(const offset_type elements) {
	// this does not permute alt values,
	// so don't use it for hashes with alts
	assert(alt_size == 0);
	// 256 entries for each pass, max_key_size / 8 passes, with 256
	// extra for holding the unbinned positions of bins while swapping
	offset_type offsets[max_key_size * 32 + 256];
	radix_sort_internal(0, elements, offsets, max_key_size - 8);
}

void hash::set_no_space_response(int i, const std::string &s) {
	if (alt_size != 0 && (i & TMP_FILE)) {
		fprintf(stderr, "Warning: cannot use TMP_FILE strategy with alt_values; TMP_FILE disabled\n");
		i &= ~TMP_FILE;
	}
	no_space_response = i;
	if (s != "NONE") {
		tmp_file_prefix = s;
		// if a directory, make sure it ends in a / so files go in it
		if (s != "") {
			struct stat buf;
			if (stat(s.c_str(), &buf) == 0 && S_ISDIR(buf.st_mode) && *(s.end()) != '/') {
				tmp_file_prefix += "/";
			}
		}
	}
}

// save sorted hash to disk for later processing
void hash::save_state(void) {
	static int count(-1);
	std::ostringstream s;
	s << tmp_file_prefix << "hash." << ++count << ".gz";
	const std::string file(s.str());
	std::list<std::string> args;
	args.push_back("gzip");
	args.push_back("-c");
	const int fd(write_fork(args, file));
	if (fd == -1) {
		exit(1);
	}
	state_files.push_back(file);
	for (offset_type i(0); i != modulus; ++i) {
		if (key_list[i] != INVALID_KEY) {
			pfwrite(fd, &key_list[i], sizeof(key_type));
			value_type x(value_list[i]);
			if (x == max_small_value) {
				const std::map<key_type, value_type>::const_iterator a(value_map.find(key_list[i]));
				if (a != value_map.end()) {
					x += a->second;
				}
			}
			pfwrite(fd, &x, sizeof(value_type));
		}
	}
	close_fork(fd);
}

// move all non-INVALID_KEY entries to front of hash table;
// the hash table is, of course, broken as a hash after this operation
void hash::squash_hash(void) {
	// used_elements will get reset at next clear()
	if (--used_elements == 0) {	// don't count INVALID_KEY
		return;
	}
	offset_type i(static_cast<offset_type>(-1));
	offset_type j(modulus);
	for (;;) {
		for (++i; i != used_elements && key_list[i] != INVALID_KEY; ++i) { }
		if (i == used_elements) {
			break;
		}
		for (--j; key_list[j] == INVALID_KEY; --j) { }
		key_list[i] = key_list[j];
		value_list[i] = value_list[j];
	}
}

bool hash::get_next_entry(const int fd, key_type &i, value_type &j, offset_type &offset) const {
	if (fd != -1) {	
		if (pfread(fd, &i, sizeof(key_type)) == -1) {
			return 0;
		}
		if (pfread(fd, &j, sizeof(value_type)) == -1) {
			fprintf(stderr, "Error: short read on state file %d\n", fd);
			exit(1);
		}
	} else {		// fd == -1 for in-memory values
		// by construction, this routine will not be called with
		// fd == -1 after it returns zero for it once, so you don't
		// have to worry about incrementing past used_elements
		if (++offset == used_elements) {
			return 0;
		}
		i = key_list[offset];
		j = value_list[offset];
		if (j == max_small_value) {
			const std::map<key_type, value_type>::const_iterator a(value_map.find(i));
			if (a != value_map.end()) {
				j += a->second;
			}
		}
	}
	return 1;
}

// set up for reading back values from saved state files
void hash::prep_for_readback(offset_type &offset, std::map<key_type, std::pair<value_type, int> > &next_keys) {
	// make sure all state files are finished being written to
	close_fork_wait(-1);
	// get the in-memory values into an ordered form
	squash_hash();	// so we don't have to sort INVALID_KEY entries
	radix_sort(used_elements);
	// open the files and read in the first values
	key_type i;
	value_type j;
	if (get_next_entry(-1, i, j, offset)) {	// in-memory first, though
		next_keys[i] = std::make_pair(j, -1);
	}
	std::list<std::string>::const_iterator a(state_files.begin());
	const std::list<std::string>::const_iterator end_a(state_files.end());
	for (; a != end_a; ++a) {
		const int fd(open_compressed(*a));
		if (fd == -1) {
			exit(1);
		}
		for (;;) {	// keep going until we get a unique entry
			// check if file is out of entries
			if (!get_next_entry(fd, i, j, offset)) {
				close_compressed(fd);
				break;
			}
			std::map<key_type, std::pair<value_type, int> >::iterator c(next_keys.find(i));
			if (c == next_keys.end()) {
				next_keys[i] = std::make_pair(j, fd);
				break;
			} else { // keep going until we get a unique entry
				c->second.first += j;
			}
		}
	}
	// by here, next_keys will have a key-sorted list of exactly
	// one key per file
}

bool hash::set_value(const key_type key, const value_type value) {
	const offset_type i(insert_offset(key));
	if (i == modulus) {	// insert failed
		return 0;
	}
	if (value <= max_small_value) {
		value_list[i] = value;
		value_map.erase(key);
	} else {
		value_list[i] = max_small_value;
		if (can_overflow) {
			value_map[key] = value - max_small_value;
		}
	}
	return 1;
}
