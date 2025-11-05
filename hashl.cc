#include "hashl.h"
#include "hashl_index.h"	// hashl_index
#include "hashl_less.h"	// hashl_less<>
#include "hashl_metadata.h"	// hashl_metadata
#include "itoa.h"	// itoa()
#include "local_endian.h"	// big_endian
#include "next_prime.h"	// next_prime()
#include "open_compressed.h"	// pfread()
#include "write_fork.h"	// pfwrite()
#include <algorithm>	// lower_bound(), sort(), swap()
#include <iomanip>	// setw()
#include <iostream>	// cerr, cout
#include <iterator>	// distance()
#include <map>		// map<>
#include <stdint.h>	// uint64_t
#include <stdlib.h>	// exit()
#include <string.h>	// memcmp(), memcpy()
#include <string>	// string
#include <utility>	// make_pair(), pair<>
#include <vector>	// vector<>

// description beginning of saved file

std::string hashl::boilerplate(void) const {
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

void hashl::init(const hash_offset_type size_asked, const size_type bits_in, std::vector<base_type> &data_in) {
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
		std::cerr << "Error: could not read hash from file: header mismatch\n";
		exit(1);
	}
	pfread(fd, &modulus, sizeof(modulus));
	pfread(fd, &collision_modulus, sizeof(collision_modulus));
	// this value is no longer used, read in here for backward compatibility
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
	// set used_elements based on actual use
	used_elements = 0;
	for (hash_offset_type i(0); i < modulus; ++i) {
		if (value_list[i]) {
			pfread(fd, &key_list[i], sizeof(size_type));
			++used_elements;
		}
	}
}

// insert a key at a particular location

hashl::hash_offset_type hashl::insert_key(const hash_offset_type i, const size_type offset) {
	if (++used_elements == modulus) {	// hash table is full
		return used_elements--;		// (always leave one empty value to mark end)
	}
	key_list[i] = offset;
	value_list[i] = 0;
	return i;
}

// find a key, or insert it if it doesn't exist; return modulus if hash is full

hashl::hash_offset_type hashl::insert_offset(const key_type &key, const key_type &comp_key, const size_type offset) {
	const base_type key_hash(key < comp_key ? key.hash() : comp_key.hash());
	hash_offset_type i(key_hash % modulus);
	if (key_list[i] == invalid_key) {		// insert
		return insert_key(i, offset);
	} else if (key.equal_to(data, key_list[i]) || comp_key.equal_to(data, key_list[i])) {
		return i;				// already present
	}
	const hash_offset_type j(collision_modulus - key_hash % collision_modulus);
	for (;;) {	// search over all elements
		i = (i + j) % modulus;
		if (key_list[i] == invalid_key) {
			return insert_key(i, offset);
		} else if (key.equal_to(data, key_list[i]) || comp_key.equal_to(data, key_list[i])) {
			return i;
		}
	}
}

// find a key; returns modulus if not found

hashl::hash_offset_type hashl::find_offset(const key_type &key, const key_type &comp_key) const {
	const base_type key_hash(key < comp_key ? key.hash() : comp_key.hash());
	hash_offset_type i(key_hash % modulus);
	if (key_list[i] == invalid_key) {
		return modulus;
	} else if (key.equal_to(data, key_list[i]) || comp_key.equal_to(data, key_list[i])) {
		return i;
	}
	const hash_offset_type j(collision_modulus - key_hash % collision_modulus);
	for (;;) {	// search over all elements
		i = (i + j) % modulus;
		if (key_list[i] == invalid_key) {
			return modulus;
		} else if (key.equal_to(data, key_list[i]) || comp_key.equal_to(data, key_list[i])) {
			return i;
		}
	}
}

hashl::hash_offset_type hashl::find_offset(const key_type &key) const {
	key_type comp_key(bit_width, word_width);
	comp_key.make_complement(key);
	return find_offset(key, comp_key);
}

// increment the count for a key (but don't create a new entry if it doesn't exist)

void hashl::increment(const key_type &key, const key_type &comp_key) {
	const hash_offset_type i(find_offset(key, comp_key));
	if (i == modulus) {	// couldn't find it
		return;
	}
	if (value_list[i] < max_small_value) {
		++value_list[i];
	}
}

bool hashl::increment(const key_type &key, const key_type &comp_key, const size_type offset) {
	const hash_offset_type i(insert_offset(key, comp_key, offset));
	if (i == modulus) {	// insert failed
		return 0;
	}
	if (value_list[i] < max_small_value) {
		++value_list[i];
	}
	return 1;
}

// insert key, but if it already exists, mark it as invalid
bool hashl::insert_unique(const key_type &key, const key_type &comp_key, const size_type offset) {
	const hash_offset_type i(insert_offset(key, comp_key, offset));
	if (i == modulus) {	// insert failed
		return 0;
	}
	if (!value_list[i]) {
		value_list[i] = 1;
	} else {
		value_list[i] = invalid_value;
	}
	return 1;
}

bool hashl::insert_invalid(const key_type &key, const key_type &comp_key, const size_type offset) {
	const hash_offset_type i(insert_offset(key, comp_key, offset));
	if (i == modulus) {	// insert failed
		return 0;
	}
	value_list[i] = invalid_value;
	return 1;
}

// return the value associated with a key (or zero if key not found)

hashl::small_value_type hashl::value(const key_type &key) const {
	const hash_offset_type i(find_offset(key));
	return i < modulus ? value_list[i] : 0;
}

// same as above, but also return the data offset

std::pair<hashl::size_type, hashl::small_value_type> hashl::entry(const key_type &key) const {
	const hash_offset_type i(find_offset(key));
	if (i < modulus) {
		return std::make_pair(key_list[i], value_list[i]);
	} else {
		return std::make_pair(0, 0);
	}
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
		// have to remove keys with zero values as they don't get read in
		if (value_list[i] && key_list[i] != invalid_key) {
			pfwrite(fd, &key_list[i], sizeof(size_type));
		}
	}
}

// regenerate key and values tables with new size - holds both new and
// old key and value lists in memory while copying

void hashl::resize(hash_offset_type size_asked) {
	if (size_asked < used_elements) {
		return;
	}
	const size_type old_modulus(modulus);
	if (size_asked < 3) {	// to avoid collision_modulus == modulus
		size_asked = 3;
	}
	modulus = next_prime(size_asked);
	// collision_modulus just needs to be relatively prime with modulus;
	// since modulus is prime, any value will do - I made it prime for fun
	collision_modulus = next_prime(size_asked / 2);
	// initialize keys and values (and save old ones)
	std::vector<size_type> old_key_list(modulus, invalid_key);
	key_list.swap(old_key_list);
	std::vector<small_value_type> old_value_list(modulus, 0);
	value_list.swap(old_value_list);
	// copy over old hash keys and values
	key_type key(bit_width, word_width), comp_key(bit_width, word_width);
	for (hash_offset_type i(0); i < old_modulus; ++i) {
		if (old_key_list[i] != invalid_key && old_value_list[i]) {
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

// recreate the hash without keys having invalid values
void hashl::purge_invalid_values() {
	// remove invalid values
	for (hash_offset_type i(0); i < modulus; ++i) {
		if (value_list[i] == invalid_value) {
			key_list[i] = invalid_key;
			--used_elements;
		}
	}
	// resize the hash for 50% load
	resize(2 * used_elements);
}

// add in new hash: any values <min_cutoff ignored, <=max_cutoff increment
// existing value by 1, values of >max_cutoff set to invalid_value

bool hashl::add(const hashl &a, const small_value_type min_cutoff, const small_value_type max_cutoff) {
	// total possible elements, actually, as there may be duplicates
	const size_type total_elements(used_elements + a.used_elements);
	if (total_elements > modulus * .7) {	// max load of 70%
		resize(total_elements * 2);	// aim for 50% load
	}
	// copy over data (and make sure to update offsets when adding new entries)
	const size_type offset(data.size() * sizeof(base_type) * 8);
	// this pads out the existing data so we don't have to shift all the new data
	// TODO: see if shifting would actually slow things down much
	//       and change metadata to remove padding there, as well
	data.reserve(data.size() + a.data.size());
	data.insert(data.end(), a.data.begin(), a.data.end());
	// loop over incoming hash and increment or invalidate entries as needed
	key_type key(a.bits(), a.words()), comp_key(a.bits(), a.words());
	for (size_type i(0); i < a.modulus; ++i) {
		if (a.key_list[i] != invalid_key && a.value_list[i] >= min_cutoff) {
			key.copy_in(a.data, a.key_list[i]);
			comp_key.make_complement(key);
			const hash_offset_type new_i(insert_offset(key, comp_key, a.key_list[i] + offset));
			if (new_i == modulus) {		// ran out of memory in hash (shouldn't happen)
				return 0;
			} else if (a.value_list[i] > max_cutoff) {	// too high - invalid (also catches invalid_value)
				value_list[new_i] = invalid_value;
			} else if (value_list[new_i] < max_small_value) {
				++value_list[new_i];
			}
		}
	}
	// combine metadata
	hashl_metadata our_md, a_md;
	if (!metadata.empty() && !a.metadata.empty()) {
		// extract metadata from blobs (both ours and a's)
		our_md.unpack(metadata);
		const size_type padding(offset / 2 - our_md.sequence_length());
		a_md.unpack(a.metadata);
		our_md.add(a_md, padding);
		our_md.pack(metadata);
	} else if (!a.metadata.empty()) {			// pad ours, then add a's
		if (offset) {
			// add dummy entry for current hash
			our_md.add_filename("unknown");
			our_md.add_readname("padding");
			our_md.add_read_range(0, offset / 2);
		}
		a_md.unpack(a.metadata);
		our_md.add(a_md);
		our_md.pack(metadata);
	} else if (!metadata.empty() && !a.data.empty()) {	// pad a's, then add to ours
		our_md.unpack(metadata);
		const size_type padding(offset / 2 - our_md.sequence_length());
		// add dummy entry for a
		const size_type a_offset(a.data.size() * sizeof(base_type) * 8);
		a_md.add_filename("unknown");
		a_md.add_readname("padding");
		a_md.add_read_range(0, a_offset / 2);
		our_md.add(a_md, padding);
		our_md.pack(metadata);
	}
	return 1;
}

// for debugging

void hashl::print(const int flags) const {
	int max_offset_width(1), max_key_width(1);
	for (size_type i(10); i < modulus; i *= 10, ++max_offset_width) { }
	for (size_type i(10); i < data.size() * sizeof(base_type) * 8; i *= 10, ++max_key_width) { }
	if (flags & print_hash_header) {
		std::cout << "modulus: " << modulus << "\n"
			<< "collision modulus: " << collision_modulus << "\n"
			<< "used elements: " << used_elements << "\n"
			<< "bit width: " << bit_width << "\n"
			<< "metadata size: " << metadata.size() << "\n"
			<< "data size: " << data.size() * sizeof(base_type) << "\n"
			<< "offset/value/key pairs:\n";
	}
	if (!(flags & (print_hash_index | print_data_offset | print_value | print_keys))) {
		return;		// nothing left to print
	}
	std::string s;
	key_type k(bit_width, word_width);
	for (hash_offset_type i(0); i < modulus; ++i) {
		if (key_list[i] != invalid_key) {
			k.copy_in(data, key_list[i]);
			k.get_sequence(s);
			if (flags & print_hash_index) {
				std::cout << std::setw(max_offset_width) << i << ' ';
			}
			if (flags & print_data_offset) {
				std::cout << std::setw(max_key_width) << key_list[i] << ' ';
			}
			if (flags & print_value) {
				std::cout << std::setw(3) << static_cast<unsigned int>(value_list[i]) << ' ';
			}
			if (flags & print_keys) {
				std::cout << s;
			}
			std::cout << '\n';
		}
	}
}

void hashl::print_sequence(const size_type start, size_type length) const {
	if (start > data.size() * sizeof(base_type) * 8) {
		return;
	} else if (length > data.size() * sizeof(base_type) * 8 - start) {
		length = data.size() * sizeof(base_type) * 8 - start;
	}
	const char values[4] = { 'A', 'C', 'G', 'T' };
	size_type word_offset = start / (sizeof(base_type) * 8);
	size_type bit_offset = sizeof(base_type) * 8 - start % (sizeof(base_type) * 8);
	for (size_type i = 0; i < length; i += 2) {
		if (bit_offset) {
			bit_offset -= 2;
		} else {
			bit_offset = sizeof(base_type) * 8 - 2;
			++word_offset;
		}
		std::cout << static_cast<char>(values[(data[word_offset] >> bit_offset) & 3]);
	}
	std::cout << '\n';
}

void hashl::get_sequence(const size_type start, size_type length, std::string &seq) const {
	seq.clear();
	if (start > data.size() * sizeof(base_type) * 8) {
		return;
	} else if (length > data.size() * sizeof(base_type) * 8 - start) {
		length = data.size() * sizeof(base_type) * 8 - start;
	}
	const char values[4] = { 'A', 'C', 'G', 'T' };
	size_type word_offset = start / (sizeof(base_type) * 8);
	size_type bit_offset = sizeof(base_type) * 8 - start % (sizeof(base_type) * 8);
	for (size_type i = 0; i < length; i += 2) {
		if (bit_offset) {
			bit_offset -= 2;
		} else {
			bit_offset = sizeof(base_type) * 8 - 2;
			++word_offset;
		}
		seq += values[(data[word_offset] >> bit_offset) & 3];
	}
}

void hashl::filtering_prep(const bool backup_values) {
	if (backup_values) {
		// make backup of value_list and zero all values
		value_list_backup.assign(modulus, 0);
		value_list.swap(value_list_backup);
	} else {
		// zero values, but keep invalid_value
		for (hash_offset_type i(0); i < modulus; ++i) {
			if (value_list[i] && value_list[i] != invalid_value) {
				value_list[i] = 0;
			}
		}
	}
}

// with no backup, set values to invalid_value if the new values don't fall between min and max
// with backup, remove values that don't fall between min and max (and resize)
// TODO: add restore_values parameter and double check versus existence of backup vector

void hashl::filtering_finish(const hashl::small_value_type min, const hashl::small_value_type max) {
	if (value_list_backup.empty()) {
		for (hash_offset_type i(0); i < modulus; ++i) {
			if (key_list[i] != invalid_key && (value_list[i] < min || max < value_list[i])) {
				value_list[i] = invalid_value;
			}
		}
	} else {
		value_list.swap(value_list_backup);
		for (hash_offset_type i(0); i < modulus; ++i) {
			if (key_list[i] != invalid_key) {
				// remove entries not present in the filter
				if (!value_list_backup[i]) {
					value_list[i] = 0;
					--used_elements;
				} else if (value_list_backup[i] < min || max < value_list_backup[i]) {
					value_list[i] = invalid_value;
				}
			}
		}
		value_list_backup = std::vector<small_value_type>();
		resize(2 * used_elements);
	}
}

void hashl::save_index(const int fd) {
	// invalidate keys for any invalid_values
	for (hash_offset_type i(0); i < modulus; ++i) {
		if (value_list[i] == invalid_value) {
			key_list[i] = invalid_key;
			--used_elements;
		}
	}
	// we no longer need the values, so free memory
	value_list = std::vector<small_value_type>();
	// shift valid key_list entries to bottom of array
	auto a = key_list.begin();
	auto end_a = a + used_elements;
	// we are guaranteed at least one invalid_key, so this will terminate
	for (; *a != invalid_key; ++a) { }
	auto b = end_a;
	// have to swap() first pair to ensure *end_a == invalid_key
	for (; *b == invalid_key; ++b) { }
	std::swap(*a, *b);
	for (++a; *a != invalid_key; ++a) { }
	// always check a first, so b won't go past end of array
	while (a != end_a) {
		for (++b; *b == invalid_key; ++b) { }
		*a = *b;
		for (++a; *a != invalid_key; ++a) { }
	}
	// remove invalid entries
	key_list.resize(used_elements);
	// sort key_list by *kmer* (not kmer position ;)
	std::sort(key_list.begin(), key_list.end(), [this](const hash_offset_type __a, const hash_offset_type __b) {return hashl_less<hashl>()(*this, __a, __b);});
	// save everything to an index
	hashl_index::save(key_list, data, metadata, bit_width, fd);
	// finish resetting this to pre-initted state
	key_list = std::vector<size_type>();
	data = std::vector<base_type>();
	metadata = std::vector<char>();
	used_elements = 0;
	modulus = 0;
	collision_modulus = 0;
	bit_width = 0;
	word_width = 0;
}

// remove data that is not referenced in the hash, and update hash references to match
void hashl::squash_data() {
	// generate ranges
	std::vector<std::pair<size_type, size_type> > offsets;	// key_list[second] = first
	for (size_type i = 0; i < key_list.size(); ++i) {
		offsets.push_back(std::make_pair(key_list[i], i));
	}
	std::sort(offsets.begin(), offsets.end());
	// make list of used ranges (the ones we'll keep)
	const auto effective_end = std::lower_bound(offsets.begin(), offsets.end(), std::make_pair<size_type, size_type>(invalid_key, 0));
	offsets.resize(std::distance(offsets.begin(), effective_end));
	std::vector<std::pair<size_type, size_type> > ranges;			// (start, stop)
	ranges.push_back(std::make_pair(offsets.front().first, offsets.front().first + bit_width));
	for (const auto &a : offsets) {
		if (ranges.back().second >= a.first) {
			ranges.back().second = a.first + bit_width;
		} else {
			ranges.push_back(std::make_pair(a.first, a.first + bit_width));
		}
	}
	// update hash references
	size_type j = 0, offset = ranges[0].first;
	for (const auto &a : offsets) {
		if (a.first > ranges[j].second) {
			offset += ranges[j + 1].first - ranges[j].second;
			++j;
		}
		key_list[a.second] -= offset;
	}
	offsets = std::vector<std::pair<size_type, size_type> >();
	// move kept data sections into place
	const size_type base_type_bits = sizeof(base_type) * 8;
	size_type current_offset = 0, current_bit = 0;
	// note: one lookout here is preserving data[current_offset]
	// on the tail end of copying in case the next range uses it
	for (const auto &a : ranges) {
		size_type i = a.first / base_type_bits;
		const size_type starting_bit = a.first % base_type_bits;
		size_type end_i = a.second / base_type_bits;
		const size_type ending_bit = a.second % base_type_bits;
		// the bits we're keeping from the head of data[current_offset]
		const base_type start_mask = current_bit ? static_cast<base_type>(-1) << (base_type_bits - current_bit) : 0;
		if (i == end_i) {	// handle small kmers
			const size_type width = ending_bit - starting_bit;
			// doesn't cross a base_type boundary
			if (current_bit + width <= base_type_bits) {
				const base_type mask = (static_cast<base_type>(-1) << (base_type_bits - width)) >> current_bit;
				if (current_bit > starting_bit) {
					data[current_offset] = (data[current_offset] & ~mask) | ((data[i] >> (current_bit - starting_bit)) & mask);
				} else {
					data[current_offset] = (data[current_offset] & ~mask) | ((data[i] << (starting_bit - current_bit)) & mask);
				}
			} else {
				data[current_offset] = (data[current_offset] & start_mask) | ((data[i] >> (current_bit - starting_bit)) & ~start_mask);
				++current_offset;
				const base_type end_mask = static_cast<base_type>(-1) << current_bit;
				data[current_offset] = ((data[i] << (starting_bit + base_type_bits - current_bit)) & end_mask) | (data[current_offset] & ~end_mask);
			}
			current_bit = (current_bit + width) % base_type_bits;
		} else if (starting_bit == current_bit) {
			data[current_offset] = (data[current_offset] & start_mask) | (data[i] & ~start_mask);
			for (++i, ++current_offset; i < end_i; ++i, ++current_offset) {
				data[current_offset] = data[i];
			}
			if (ending_bit) {
				const base_type end_mask = static_cast<base_type>(-1) << (base_type_bits - ending_bit);
				data[current_offset] = (data[i] & end_mask) | (data[current_offset] & ~end_mask);
			}
			current_bit = ending_bit;
		} else if (starting_bit < current_bit) {
			const size_type shift_right = current_bit - starting_bit;
			const size_type shift_left = base_type_bits - shift_right;
			size_type final_bits = ending_bit + shift_right;
			if (final_bits < base_type_bits) {
				--end_i;
			} else {
				final_bits -= base_type_bits;
			}
			data[current_offset] = (data[current_offset] & start_mask) | ((data[i] >> shift_right) & ~start_mask);
			for (++current_offset; i < end_i; ++i, ++current_offset) {
				data[current_offset] = (data[i] << shift_left) | (data[i + 1] >> shift_right);
			}
			if (final_bits) {
				const base_type end_mask = static_cast<base_type>(-1) << (base_type_bits - final_bits);
				data[current_offset] = ((data[i] << shift_left) & end_mask) | (data[current_offset] & ~end_mask);
			}
			current_bit = final_bits;
		} else {	// starting_bit > current_bit
			const size_type shift_left = starting_bit - current_bit;
			const size_type shift_right = base_type_bits - shift_left;
			size_type final_bits;
			if (ending_bit < shift_left) {
				final_bits = base_type_bits + ending_bit - shift_left;
				--end_i;
			} else {
				final_bits = ending_bit - shift_left;
			}
			data[current_offset] = (data[current_offset] & start_mask) | (((data[i] << shift_left) | (data[i + 1] >> shift_right)) & ~start_mask);
			for (++i, ++current_offset; i < end_i; ++i, ++current_offset) {
				data[current_offset] = (data[i] << shift_left) | (data[i + 1] >> shift_right);
			}
			if (ending_bit < shift_left) {
				const base_type end_mask = static_cast<base_type>(-1) << (base_type_bits - final_bits);
				data[current_offset] = (((data[i] << shift_left) | (data[i + 1] >> shift_right)) & end_mask) | (data[current_offset] & ~end_mask);
			} else if (final_bits) {
				const base_type end_mask = static_cast<base_type>(-1) << (base_type_bits - final_bits);
				data[current_offset] = ((data[i] << shift_left) & end_mask) | (data[current_offset] & ~end_mask);
			}
			current_bit = final_bits;
		}
	}
	// zero out leftover data
	if (current_bit) {
		const base_type mask = static_cast<base_type>(-1) << (base_type_bits - current_bit);
		data[current_offset] &= mask;
		++current_offset;
	}
	data.resize(current_offset);
	// convert ranges from bit positions to basepair positions
	for (auto &a : ranges) {
		a.first /= 2;
		a.second /= 2;
	}
	// update metadata
	hashl_metadata md;
	md.unpack(metadata);
	md.update_ranges(ranges);
	md.pack(metadata);
}
