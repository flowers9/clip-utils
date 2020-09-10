#include "hashn.h"
#include "itoa.h"	// itoa()
#include "local_endian.h"	// big_endian
#include "next_prime.h"	// next_prime()
#include "open_compressed.h"	// close_compressed(), open_compressed(), pfread()
#include "write_fork.h"	// close_fork(), close_fork_wait(), pfwrite(), write_fork()
#include <algorithm>	// swap()
#include <cassert>	// assert()
#include <list>		// list<>
#include <map>		// map<>
#include <new>		// new
#include <sstream>	// ostringstream
#include <stdio.h>	// fprintf(), stderr
#include <stdlib.h>	// exit()
#include <string.h>	// memcmp(), memcpy()
#include <string>	// string
#include <sys/stat.h>	// S_ISDIR(), stat(), struct stat
#include <unistd.h>	// unlink()
#include <utility>	// make_pair(), pair<>
#include <vector>	// vector<>

hashn::~hashn() {
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

std::string hashn::boilerplate() const {
	std::string s("hashn\n");
	s += itoa(sizeof(base_type));
	s += " bytes\n";
#ifdef big_endian
	s += "big endian\n";
#else
	s += "little endian\n";
#endif
	return s;
}

void hashn::init(offset_type size_asked, const unsigned long bits_in, const offset_type alt_size_in) {
	if (alt_size_in > 8 * sizeof(offset_type)) {
		fprintf(stderr, "Error: hash alt size too large: %lu > %lu\n", alt_size_in, 8 * sizeof(offset_type));
		exit(1);
	}
	bit_width = bits_in;
	word_width = (bit_width + 8 * sizeof(base_type) - 1) / (8 * sizeof(base_type));
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
	// array is modulus + 1 because invalid_key lives at [modulus]
	const offset_type n((modulus + 1) * word_width);
	key_list = new base_type[n];
	value_list = new small_value_type[modulus];
	if (alt_size == 0) {
		alt_list = NULL;
		alt_map = NULL;
	} else {
		alt_list = new small_value_type[modulus * alt_size];
		alt_map = new std::map<std::string, value_type>[alt_size];
	}
	// initialize keys; values are initialized as keys are entered
	invalid_key.assign(*this, modulus);
	for (offset_type i = 0; i != n; ++i) {
		key_list[i] = INVALID_KEY;
	}
}

void hashn::init_from_file(const int fd) {
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
	pfread(fd, &alt_size, sizeof(alt_size));
	pfread(fd, &bit_width, sizeof(bit_width));
	word_width = (bit_width + 8 * sizeof(base_type) - 1) / (8 * sizeof(base_type));
	const offset_type n((modulus + 1) * word_width);
	key_list = new base_type[n];
	value_list = new small_value_type[modulus];
	if (alt_size == 0) {
		alt_list = NULL;
		alt_map = NULL;
	} else {
		alt_list = new small_value_type[modulus * alt_size];
		alt_map = new std::map<std::string, value_type>[alt_size];
	}
	pfread(fd, &value_list[0], sizeof(small_value_type) * modulus);
	base_type *a(key_list);
	for (offset_type i(0); i != modulus; ++i) {
		if (value_list[i] == 0) {
			base_type * const end_a(a + word_width);
			for (; a != end_a; ++a) {
				*a = INVALID_KEY;
			}
		} else {
			pfread(fd, a, sizeof(base_type) * word_width);
			a += word_width;
		}
	}
	invalid_key.assign(*this, modulus);
	pfread(fd, a, sizeof(base_type) * word_width);		// invalid_key
	offset_type x;
	pfread(fd, &x, sizeof(x));
	char buf[sizeof(base_type) * word_width];
	for (; x != 0; --x) {
		value_type j;
		pfread(fd, buf, sizeof(buf));
		pfread(fd, &j, sizeof(j));
		value_map[buf] = j;
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
			std::map<std::string, value_type> &z = alt_map[k];
			pfread(fd, &x, sizeof(x));
			for (; x != 0; --x) {
				value_type j;
				pfread(fd, buf, sizeof(buf));
				pfread(fd, &j, sizeof(j));
				z[buf] = j;
			}
		}
	}
}

// returns next empty spot found, or modulus if it spots the key first

hashn::offset_type hashn::find_empty_offset(const base_type * const z) const {
	const base_type key_hash(hash(z));
	offset_type i(key_hash % modulus);
	const base_type * const k(key_list + i * word_width);
	if (invalid_key.equal(k)) {
		return i;
	} else if (equal(z, k)) {
		return modulus;
	}
	const offset_type j(collision_modulus - key_hash % collision_modulus);
	for (;;) {	// search over all elements
		i = (i + j) % modulus;
		const base_type * const l(key_list + i * word_width);
		if (invalid_key.equal(l)) {
			return i;
		} else if (equal(z, l)) {
			return modulus;
		}
	}
}

void hashn::rehash() {
	// first pass, fill all non-collision slots
	offset_type i(0);
	base_type *k(key_list);
	for (; i != modulus; ++i, k += word_width) {
		if (!invalid_key.equal(k)) {
			const offset_type j(hash(k) % modulus);
			// swap if non-collision slot has a collision fill
			if (i != j) {
				base_type * const l(key_list + j * word_width);
				if (invalid_key.equal(l) || hash(l) % modulus != j) {
					swap(l, k);
					std::swap(value_list[j], value_list[i]);
					--i;
					k -= word_width;
				}
			}
		}
	}
	// subsequent passes, shift collisions to empty slots; the shifts
	// are entropy lowering operations, meaning there is no risk of
	// an infinite loop
	bool changed;
	do {
		changed = 0;
		for (i = 0, k = key_list; i != modulus; ++i, k += word_width) {
			if (!invalid_key.equal(k)) {
				const offset_type j(find_empty_offset(k));
				if (j != modulus) {
					changed = 1;
					copy(key_list + j * word_width, k);
					invalid_key.copy_out(k);
					value_list[j] = value_list[i];
				}
			}
		}
	} while (changed);
}

// same as rehash(), but moves alt_list as well

void hashn::rehash_alt() {
	// first pass, fill all non-collision slots
	offset_type i(0);
	base_type *k(key_list);
	for (; i != modulus; ++i, k += word_width) {
		if (!invalid_key.equal(k)) {
			const offset_type j(hash(k) % modulus);
			// swap if non-collision slot has a collision fill
			if (i != j) {
				base_type * const l(key_list + j * word_width);
				if (invalid_key.equal(l) || hash(l) % modulus != j) {
					swap(l, k);
					std::swap(value_list[j], value_list[i]);
					offset_type x(j * alt_size);
					offset_type y(i * alt_size);
					const offset_type end_x(x + alt_size);
					for (; x != end_x; ++x, ++y) {
						std::swap(alt_list[x], alt_list[y]);
					}
					--i;
					k -= word_width;
				}
			}
		}
	}
	// subsequent passes, shift collisions to empty slots; the shifts
	// are entropy lowering operations, meaning there is no risk of
	// an infinite loop
	bool changed;
	do {
		changed = 0;
		for (i = 0, k = key_list; i != modulus; ++i, k += word_width) {
			if (!invalid_key.equal(k)) {
				const offset_type j(find_empty_offset(k));
				if (j != modulus) {
					changed = 1;
					copy(key_list + j * word_width, k);
					invalid_key.copy_out(k);
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

// remove all keys with value == 1; returns if hash has space left

bool hashn::clean_hash() {
	offset_type i(0);
	base_type * k(key_list);
	for (; i != modulus; ++i, k += word_width) {
		if (!invalid_key.equal(k) && value_list[i] == 1) {
			invalid_key.copy_out(k);
			--used_elements;
		}
	}
	// as you can't decrement elements, once hash is filled, that's it
	if (used_elements == modulus) {
		return 0;
	} else if (alt_list == NULL) {
		rehash();
	} else {
		rehash_alt();
	}
	return 1;
}

// insert a key at a particular location

hashn::offset_type hashn::insert_key(const offset_type i, base_type * const k, const key_type &key) {
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
	key.copy_out(k);	// copy key into table
	value_list[i] = 0;	// init counts
	offset_type j(i * alt_size);
	const offset_type end_j(j + alt_size);
	for (; j != end_j; ++j) {
		alt_list[j] = 0;
	}
	return i;
}

// find a key, or insert it if it doesn't exist; return modulus if hash is full

hashn::offset_type hashn::insert_offset(const key_type &key) {
	const base_type key_hash(key.hash());
	offset_type i(key_hash % modulus);
	base_type * const k(key_list + i * word_width);
	if (invalid_key.equal(k)) {	// insert
		return insert_key(i, k, key);
	} else if (key.equal(k)) {	// already present
		return i;
	}
	const offset_type j(collision_modulus - key_hash % collision_modulus);
	for (;;) {	// search over all elements
		i = (i + j) % modulus;
		base_type * const l(key_list + i * word_width);
		if (invalid_key.equal(l)) {
			return insert_key(i, l, key);
		} else if (key.equal(l)) {
			return i;
		}
	}
}

// find a key; return modulus as offset if not found

hashn::offset_type hashn::find_offset(const key_type &key) const {
	const base_type key_hash(key.hash());
	offset_type i(key_hash % modulus);
	const base_type * const k(key_list + i * word_width);
	if (invalid_key.equal(k)) {
		return modulus;
	} else if (key.equal(k)) {
		return i;
	}
	const offset_type j(collision_modulus - key_hash % collision_modulus);
	for (;;) {	// search over all elements
		i = (i + j) % modulus;
		const base_type * const l(key_list + i * word_width);
		if (invalid_key.equal(l)) {
			return modulus;
		} else if (key.equal(l)) {
			return i;
		}
	}
}

// increment the count for a key

bool hashn::increment(const key_type &key) {
	const offset_type i(insert_offset(key));
	if (i == modulus) {	// insert failed
		return 0;
	}
	if (value_list[i] != max_small_value) {
		++value_list[i];
	} else {
		++value_map[key.string()];
	}
	return 1;
}

// increment only the alt values, using x as a bit flag to mark which ones

bool hashn::increment_alt(const key_type &key, const offset_type x) {
	const offset_type i(insert_offset(key));
	if (i == modulus) {	// insert failed
		return 0;
	}
	const std::string key_str(key.string());
	const offset_type start_j(i * alt_size);
	const offset_type end_j(start_j + alt_size);
	for (offset_type j(start_j), z(1); j != end_j; ++j, z <<= 1) {
		if (x & z) {
			if (alt_list[j] != max_small_value) {
				++alt_list[j];
			} else {
				++alt_map[j - start_j][key_str];
			}
		}
	}
	return 1;
}

// assign the count for a key

bool hashn::assign(const key_type &key, const value_type x) {
	const offset_type i(insert_offset(key));
	if (i == modulus) {	// insert failed
		return 0;
	}
	if (x <= max_small_value) {
		value_list[i] = x;
		const std::map<std::string, value_type>::iterator a(value_map.find(key.string()));
		if (a != value_map.end()) {
			value_map.erase(a);
		}
	} else {
		value_list[i] = max_small_value;
		value_map[key.string()] = x - max_small_value;
	}
	return 1;
}

// return the value associated with a key

hashn::value_type hashn::value(const key_type &key) const {
	const offset_type i(find_offset(key));
	if (i == modulus) {	// key not found
		return 0;
	} else if (value_list[i] != max_small_value) {
		return value_list[i];
	} else {
		// use find() to avoid inserting a value into value_map
		const std::map<std::string, value_type>::const_iterator a(value_map.find(key.string()));
		if (a == value_map.end()) {
			return max_small_value;
		} else {
			return a->second + max_small_value;
		}
	}
}

hashn::value_type hashn::value(const key_type &key, value_type x[]) const {
	const offset_type i(find_offset(key));
	if (i == modulus) {	// key not found
		return 0;
	}
	const std::string key_str(key.string());
	offset_type j(0);
	const offset_type alt_offset(i * alt_size);
	for (; j != alt_size; ++j) {
		if (alt_list[j + alt_offset] != max_small_value) {
			x[j] = alt_list[j + alt_offset];
		} else {
			const std::map<std::string, value_type> &b(alt_map[j]);
			// use find() to avoid inserting a value into value_map
			const std::map<std::string, value_type>::const_iterator a(b.find(key_str));
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
		const std::map<std::string, value_type>::const_iterator a(value_map.find(key_str));
		if (a == value_map.end()) {
			return max_small_value;
		} else {
			return a->second + max_small_value;
		}
	}
}

// reset hash to empty state

void hashn::clear(const bool mostly_clear) {
	used_elements = 1;	// to account for minimum of one INVALID_KEYs
	// initialize keys
	const offset_type n(modulus * word_width);
	offset_type i(0);
	for (; i != n; ++i) {
		key_list[i] = INVALID_KEY;
	}
	value_map.clear();
	for (i = 0; i != alt_size; ++i) {
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

hashn::const_iterator hashn::begin() {
	if (state_files.empty()) {
		if (used_elements == 1) {	// only key is INVALID_KEY
			return end();
		}
		const_iterator a(this, 0);
		// advance to first valid value
		if (invalid_key.equal(key_list)) {
			++a;
		}
		return a;
	} else {
		offset_type index(static_cast<offset_type>(-1));
		std::map<sort_key, std::pair<value_type, int> > next_keys;
		refcount_array<base_type> key_buffer;
		prep_for_readback(index, next_keys, key_buffer);
		return const_iterator(this, index, next_keys, key_buffer);
	}
}

// return a string uniquely associated with this key; to avoid byte layout
// issues, string is padded to full word_width

std::string hashn::key_type_base::string() const {
	const int n(word_width * sizeof(base_type));
	std::string s(n, 0);	// 0 values are valid inside a std::string
	// can't make this a static_cast<>
	const char * const t((char *)k);
	for (int i(0); i != n; ++i) {
		s[i] = t[i];
	}
	return s;
}

hashn::const_iterator::const_iterator(const hashn * const a, const offset_type i) : list(a), offset(i), key(*a, a->key_list + i * a->word_width) {
	if (offset != list->modulus) {
		if (list->value_list[offset] != max_small_value) {
			value = list->value_list[offset];
		} else {
			const std::map<std::string, value_type>::const_iterator b(list->value_map.find(key.string()));
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

hashn::const_iterator::const_iterator(const hashn * const a, const offset_type i, const std::map<sort_key, std::pair<value_type, int> > &c, refcount_array<base_type> &d) : list(a), offset(i), key_buffer(d), next_keys(c) {
	if (offset != list->modulus) {
		const std::map<sort_key, std::pair<value_type, int> >::iterator b(next_keys.begin());
		key.assign(*list, b->first.k);
		value = b->second.first;
	} else {
		key.assign(*list, list->key_list + offset * list->word_width);
		value = 0;
	}
}

// advance to next used element

void hashn::const_iterator::increment() {
	if (offset == list->modulus) {
		return;
	}
	if (next_keys.empty()) {
		for (++offset; offset != list->modulus && list->invalid_key.equal(list->key_list + offset * list->word_width); ++offset) { }
		key.assign(*list, offset);
		if (offset != list->modulus) {
			if (list->value_list[offset] != max_small_value) {
				value = list->value_list[offset];
			} else {
				const std::map<std::string, value_type>::const_iterator a(list->value_map.find(key.string()));
				if (a == list->value_map.end()) {
					value = max_small_value;
				} else {
					value = a->second + max_small_value;
				}
			}
		} else {
			value = 0;
		}
	} else {
		// have to delete and reinsert identical key to get it to
		// sort properly, as it sorts on the value pointed to, and
		// there's no other way to signal when that changes
		const std::map<sort_key, std::pair<value_type, int> >::iterator a(next_keys.begin());
		const int fd(a->second.second);
		sort_key tmp_key(a->first);
		next_keys.erase(a);
		// get next unique entry from that file
		for (;;) {
			// check if file is out of entries
			if (!list->get_next_entry(fd, tmp_key, value, offset)) {
				if (fd != -1) {
					close_compressed(fd);
				}
				break;
			}
			std::map<sort_key, std::pair<value_type, int> >::iterator c(next_keys.find(tmp_key));
			if (c == next_keys.end()) {
				next_keys[tmp_key] = std::make_pair(value, fd);
				break;
			} else { // keep going until we get a unique entry
				c->second.first += value;
			}
		}
		if (!next_keys.empty()) {
			const std::map<sort_key, std::pair<value_type, int> >::iterator b(next_keys.begin());
			key.assign(*list, b->first.k);
			value = b->second.first;
		} else {
			// mark end of iteration; can't happen through
			// get_next_entry() because that maxes offset out
			// at list->used_elements, which is guaranteed to
			// be at least one less than list->modulus, thanks
			// to squash_hash()
			offset = list->modulus;
			key.assign(*list, offset);
			value = 0;
		}
	}
}

// extract list of alt_values associated with current offset

void hashn::const_iterator::get_alt_values(value_type x[]) const {
	const std::string key_str(key.string());
	const offset_type alt_offset(offset * list->alt_size);
	for (offset_type i(0); i != list->alt_size; ++i) {
		if (list->alt_list[alt_offset + i] != max_small_value) {
			x[i] = list->alt_list[alt_offset + i];
		} else {
			// use find() to avoid inserting a value into alt_map
			const std::map<std::string, value_type>::const_iterator a(list->alt_map[i].find(key_str));
			if (a == list->alt_map[i].end()) {
				x[i] = max_small_value;
			} else {
				x[i] = a->second + max_small_value;
			}
		}
	}
}

void hashn::save(const int fd) const {
	const std::string s(boilerplate());
	pfwrite(fd, s.c_str(), s.size());
	pfwrite(fd, &modulus, sizeof(modulus));
	pfwrite(fd, &collision_modulus, sizeof(collision_modulus));
	pfwrite(fd, &used_elements, sizeof(used_elements));
	pfwrite(fd, &alt_size, sizeof(alt_size));
	pfwrite(fd, &bit_width, sizeof(bit_width));
	// normally, unused entries have uninitialized data, but here we
	// store them with zero values so we can use them to mark entries
	// in the other arrays that we don't bother to write to file
	const small_value_type zero(0);
	const base_type *j(key_list);
	for (offset_type i(0); i != modulus; ++i, j += word_width) {
		pfwrite(fd, invalid_key.equal(j) ? &zero : &value_list[i], sizeof(small_value_type));
	}
	j = key_list;
	const base_type * const end_j(j + modulus * word_width);
	for (; j != end_j; j += word_width) {
		if (!invalid_key.equal(j)) {
			pfwrite(fd, j, sizeof(base_type) * word_width);
		}
	}
	pfwrite(fd, j, sizeof(base_type) * word_width);		// invalid_key
	offset_type x;
	x = value_map.size();
	pfwrite(fd, &x, sizeof(x));
	std::map<std::string, value_type>::const_iterator a(value_map.begin());
	const std::map<std::string, value_type>::const_iterator end_a(value_map.end());
	for (; a != end_a; ++a) {
		pfwrite(fd, a->first.c_str(), a->first.size());
		pfwrite(fd, &a->second, sizeof(value_type));
	}
	if (alt_size != 0) {
		j = key_list;
		for (offset_type i(0); j != end_j; j += word_width, ++i) {
			if (!invalid_key.equal(j)) {
				pfwrite(fd, &alt_list[i * alt_size], sizeof(small_value_type) * alt_size);
			}
		}
		// alt map overflows
		for (offset_type k(0); k != alt_size; ++k) {
			std::map<std::string, value_type> &z = alt_map[k];
			x = z.size();
			pfwrite(fd, &x, sizeof(x));
			std::map<std::string, value_type>::const_iterator b(z.begin());
			const std::map<std::string, value_type>::const_iterator end_b(z.end());
			for (; b != end_b; ++b) {
				pfwrite(fd, b->first.c_str(), b->first.size());
				pfwrite(fd, &b->second, sizeof(value_type));
			}
		}
	}
}

// a shell sort does an insertion sort on an array over a range of gap sizes;
// insert sorts are inefficient as they can involve shifting many values,
// but they can be good for small arrays and mostly ordered arrays; the
// shell sort does insert sorts on spaced sets of increasing size so that
// as the arrays get larger, they are also more ordered

// optimized (-O2) shell_sort() beats qsort() for n < ~32k (on random data)

void hashn::shell_sort(const offset_type start_index, const offset_type stop_index) {
	// if you plan to only use arrays shorter than 701 (or some
	// smaller value), reduce the array to only the ones that cover
	// your size range (701, 301, 132)
	const static offset_type gaps[] = { 57, 23, 10, 4, 1 };
	const static int gap_count(sizeof(gaps) / sizeof(offset_type));
	base_type my_key[word_width];	// tmp buffer
	for (int i(0); i != gap_count; ++i) {
		const offset_type gap(gaps[i]);
		const offset_type start_index_gap(start_index + gap);
		const offset_type xgap(gap * word_width);
		// do an insertion sort with the given spacing (this will
		// skip gap sizes larger than stop_index - start_index)
		offset_type j(start_index_gap);
		base_type *x(key_list + start_index * word_width + xgap);
		for (; j < stop_index; ++j, x += word_width) {
			if (less_than(x, x - xgap)) {
				copy(my_key, x);
				const small_value_type my_value(value_list[j]);
				offset_type k(j);
				base_type *y(x);
				for (; k >= start_index_gap && less_than(my_key, y - xgap); k -= gap, y -= xgap) {
					copy(y, y - xgap);
					value_list[k] = value_list[k - gap];
				}
				copy(y, my_key);
				value_list[k] = my_value;
			}
		}
	}
}

void hashn::calculate_offsets(const offset_type start_index, const offset_type stop_index, offset_type * offsets, const int bit_shift, const int word_offset) const {
	// we actually count in offsets[1-256], but we don't care about
	// the value in offsets[256], so don't bother initializing it
	offsets[0] = start_index;
	for (int i(1); i != 256; ++i) {		// initialize array
		offsets[i] = 0;
	}
	++offsets;	// offset[0] will always be start_index, and we always
			// have extra space at the end of the array,
			// so we won't overrun
	const base_type *a(key_list + start_index * word_width + word_offset);
	const base_type * const end_a(key_list + stop_index * word_width + word_offset);
	for (; a != end_a; a += word_width) {	// count
		++offsets[((*a) >> bit_shift) & 255];
	}
	--offsets;
	// as what we actually want are starting offsets into the
	// destination array, accumulate the totals, so offset[n] points
	// to the starting location for bin n
	for (int i(1); i != 256; ++i) {
		offsets[i] += offsets[i - 1];
	}
}

void hashn::radix_sort_internal(const offset_type start_index, const offset_type stop_index, offset_type * offsets, const int shift) {
	const int word_offset(word_width - (shift + 8 + 8 * sizeof(base_type) - 1) / (8 * sizeof(base_type)));
	const int bit_shift(shift % (8 * sizeof(base_type)));
	// use shell sort for terminal sorting;
	// this seems to give very similar results when set between 128 and 4k
	if (stop_index - start_index < 512) {
		shell_sort(start_index, stop_index);
		return;
	}
	calculate_offsets(start_index, stop_index, offsets, bit_shift, word_offset);
	// swap unbinned elements into bins, growing sorted bins as we go
	offset_type * const unbinned_start(offsets + 256);
	memcpy(unbinned_start, offsets, 256 * sizeof(offset_type));
	++offsets;	// shift offsets to be end of bins, instead of start;
	// we never sort the last bin, as it's guaranteed to be full by the
	// time we get to it, so no worries about it never having been set
	for (int i(0); i != 255; ++i) {
		// iterate over the unsorted portion of each bin
		offset_type j(unbinned_start[i]);
		while (j != offsets[i]) {
			const int my_bin((*(key_list + j * word_width + word_offset) >> bit_shift) & 255);
			if (my_bin != i) {
				// I think it's more efficient to do extra
				// swaps here, and not check to see if
				// the target location is already in the
				// correct spot; not sure, though
				swap(key_list + unbinned_start[my_bin] * word_width, key_list + j * word_width);
				std::swap(value_list[unbinned_start[my_bin]], value_list[j]);
				++unbinned_start[my_bin];
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

void hashn::radix_sort(const offset_type elements) {
	// this does not permute alt values,
	// so don't use it for hashes with alts
	assert(alt_size == 0);
	const int shift(bit_width < 8 ? 8 : ((bit_width + 7) / 8) * 8);
	// 256 entries for each pass, shift / 8 passes, with 256 extra for
	// holding the unbinned positions of bins while swapping
	offset_type offsets[shift * 32 + 256];
	radix_sort_internal(0, elements, offsets, shift - 8);
}

void hashn::set_no_space_response(int i, const std::string &s) {
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

// sort hash and save to disk for later processing
void hashn::save_state(void) {
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
	base_type *j(key_list);
	const base_type * const end_j(j + modulus * word_width);
	for (offset_type i(0); j != end_j; j += word_width, ++i) {
		if (!invalid_key.equal(j)) {
			pfwrite(fd, j, sizeof(base_type) * word_width);
			value_type x(value_list[i]);
			if (x == max_small_value) {
				const key_type_internal tmp_key(*this, j);
				const std::map<std::string, value_type>::const_iterator a(value_map.find(tmp_key.string()));
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
void hashn::squash_hash(void) {
	// used_elements will get reset at next clear()
	if (--used_elements == 0) {	// don't count INVALID_KEY
		return;
	}
	offset_type i(static_cast<offset_type>(-1));
	offset_type j(modulus);
	base_type *a(key_list - word_width);
	base_type *b(key_list + modulus * word_width);
	const base_type * const end_a(key_list + used_elements * word_width);
	for (;;) {
		for (a += word_width, ++i; a != end_a && !invalid_key.equal(a); a += word_width, ++i) { }
		if (a == end_a) {
			break;
		}
		for (b -= word_width, --j; invalid_key.equal(b); b -= word_width, --j) { }
		copy(a, b);
		value_list[i] = value_list[j];
	}
}

bool hashn::get_next_entry(const int fd, sort_key &i, value_type &j, offset_type &offset) const {
	if (fd != -1) {	
		if (pfread(fd, i.k, sizeof(base_type) * word_width) == -1) {
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
		i.k += word_width;
		j = value_list[offset];
		if (j == max_small_value) {
			const std::map<std::string, value_type>::const_iterator a(value_map.find(key_type_base(word_width, i.k).string()));
			if (a != value_map.end()) {
				j += a->second;
			}
		}
	}
	return 1;
}

// set up for reading back values from saved state files
void hashn::prep_for_readback(offset_type &offset, std::map<sort_key, std::pair<value_type, int> > &next_keys, refcount_array<base_type> &key_buffer) {
	// make sure all state files are finished being written to
	// (do this first so if there could be any buffering issues, they're
	// hopefully less likely to come up)
	close_fork_wait(-1);
	// get the in-memory values into an ordered form
	squash_hash();	// so we don't have to sort INVALID_KEY entries
	radix_sort(used_elements);
	sort_key i(word_width, key_list - word_width);
	value_type j;
	if (get_next_entry(-1, i, j, offset)) {		// in-memory first
		next_keys[i] = std::make_pair(j, -1);
	}
	// open the files and read in the first values
	// have to count fd's first in order to allocate key_buffer space
	int highest_fd(0);
	std::vector<int> fd_list;
	fd_list.reserve(state_files.size());
	std::list<std::string>::const_iterator a(state_files.begin());
	const std::list<std::string>::const_iterator end_a(state_files.end());
	for (; a != end_a; ++a) {
		const int fd(open_compressed(*a));
		if (fd == -1) {
			exit(1);
		}
		fd_list.push_back(fd);
		if (highest_fd < fd) {
			highest_fd = fd;
		}
	}
	key_buffer.resize((highest_fd + 1) * word_width);
	std::vector<int>::const_iterator b(fd_list.begin());
	const std::vector<int>::const_iterator end_b(fd_list.end());
	for (; b != end_b; ++b) {
		const int fd(*b);
		i.k = &key_buffer[fd * word_width];
		for (;;) {	// keep going until we get a unique entry
			// check if file is out of entries
			if (!get_next_entry(fd, i, j, offset)) {
				close_compressed(fd);
				break;
			}
			std::map<sort_key, std::pair<value_type, int> >::iterator c(next_keys.find(i));
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
