#include "hash.h"	// hash
#include "hash_read_hits.h"
#include "hist_lib_hash.h"	// convert_key()
#include "itoa.h"	// itoa()
#include "kmer_lookup_info.h"	// KmerLookupInfo
#include "local_endian.h"	// big_endian
#include "next_prime.h"	// next_prime()
#include "open_compressed.h"	// pfread()
#include "write_fork.h"	// pfwrite()
#include <cassert>	// assert()
#include <iostream>	// cout
#include <limits.h>	// UINT64_MAX
#include <map>		// map<>
#include <new>		// new[]
#include <stdio.h>	// fprintf(), stderr
#include <stdlib.h>	// exit()
#include <string.h>	// memcmp()
#include <string>	// string

// description beginning of saved file

std::string hash_read_hits::boilerplate() const {
	std::string s("hash_read_hits\n");
	s += itoa(sizeof(key_type));
	s += " bytes\n";
#ifdef big_endian
	s += "big endian\n";
#else
	s += "little endian\n";
#endif
	return s;
}

// allocate arrays and copy in kmers from mer_list
// can't make mer_list const because of the way tmp file readback happens in hash
hash_read_hits::hash_read_hits(hash &mer_list, const double hash_usage) {
	used_elements = 1;	// to account for minimum of one INVALID_KEYs
	assert(0 < hash_usage && hash_usage <= 1);
	size_t size_asked(mer_list.size() / hash_usage + 1);
	if (size_asked < 3) {	// to avoid collision_modulus == modulus
		size_asked = 3;
	}
	modulus = next_prime(size_asked);
	// collision_modulus just needs to be relatively prime with modulus;
	// since modulus is prime, any value will do - I made it prime for fun
	collision_modulus = next_prime(size_asked / 2);
	key_list = new key_type[modulus];
	value_list = new small_value_type[modulus];
	read_offset_list = new read_offset_type[modulus];
	read_list_size = 0;
	hash::const_iterator a(mer_list.begin());
	const hash::const_iterator end_a(mer_list.end());
	for (; a != end_a; ++a) {
		read_list_size += a.value;
	}
	read_list = new read_type[read_list_size];
	// initialize keys; values are initialized as keys are entered
	for (offset_type i = 0; i != modulus; ++i) {
		key_list[i] = INVALID_KEY;
	}
	// now copy over keys from mer_list, and initialize offsets into read_list
	read_offset_type offset(0);
	for (a = mer_list.begin(); a != end_a; ++a) {
		const offset_type i(insert_offset(a.key));
		assert(i != modulus);	// insert should not fail
		read_offset_list[i] = offset;
		offset += a.value;
	}
}

// insert a key at a particular location

hash_read_hits::offset_type hash_read_hits::insert_key(const offset_type i, const key_type key) {
	assert(used_elements != modulus);
	++used_elements;
	key_list[i] = key;
	value_list[i] = 0;
	return i;
}

// find a key, or insert it if it doesn't exist; return modulus if hash is full

hash_read_hits::offset_type hash_read_hits::insert_offset(const key_type key) {
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

// find a key; return modulus as offset if not found

hash_read_hits::offset_type hash_read_hits::find_offset(const key_type key) const {
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

// add a read match to an existing key

void hash_read_hits::add_read(const key_type key, const read_type read) {
	const offset_type i(find_offset(key));
	assert(i != modulus);
	value_type j;
	if (value_list[i] != max_small_value) {
		j = ++value_list[i];
	} else {
		j = ++value_map[i] + max_small_value;
	}
	// -1 because we just incremented it
	read_list[read_offset_list[i] + j - 1] = read;
}

// add reads associated with kmer (if any) to list

void hash_read_hits::get_reads(const key_type key, std::map<read_type, int> &reads, const value_type max_hits) const {
	const offset_type i(find_offset(key));
	if (i == modulus) {	// key not found
		return;
	}
	value_type n(value_list[i]);
	if (n == max_small_value) {
		// use find() to avoid inserting a value into value_map
		const std::map<offset_type, value_type>::const_iterator a(value_map.find(i));
		if (a != value_map.end()) {
			n += a->second;
		}
	}
	if (n <= max_hits) {
		size_t j(read_offset_list[i]);
		const size_t end_j(j + n);
		for (; j != end_j; ++j) {
			++reads[read_list[j]];
		}
	}
}

void hash_read_hits::save(const int fd) const {
	const std::string s(boilerplate());
	pfwrite(fd, s.c_str(), s.size());
	pfwrite(fd, &used_elements, sizeof(used_elements));
	pfwrite(fd, &modulus, sizeof(modulus));
	pfwrite(fd, &collision_modulus, sizeof(collision_modulus));
	pfwrite(fd, &read_list_size, sizeof(read_list_size));
	pfwrite(fd, key_list, modulus * sizeof(key_type));
	pfwrite(fd, value_list, modulus * sizeof(small_value_type));
	pfwrite(fd, read_offset_list, modulus * sizeof(read_offset_type));
	pfwrite(fd, read_list, read_list_size * sizeof(read_type));
	const offset_type x(value_map.size());
	pfwrite(fd, &x, sizeof(x));
	std::map<key_type, value_type>::const_iterator a(value_map.begin());
	const std::map<key_type, value_type>::const_iterator end_a(value_map.end());
	for (; a != end_a; ++a) {
		pfwrite(fd, &a->first, sizeof(key_type));
		pfwrite(fd, &a->second, sizeof(value_type));
	}
}

void hash_read_hits::restore(const int fd) {
	const std::string s(boilerplate());
	char t[s.size()];
	pfread(fd, t, s.size());
	if (memcmp(s.c_str(), t, s.size()) != 0) {
		fprintf(stderr, "Error: could not read hash from file: header mismatch\n");
		exit(1);
	}
	pfread(fd, &used_elements, sizeof(used_elements));
	pfread(fd, &modulus, sizeof(modulus));
	pfread(fd, &collision_modulus, sizeof(collision_modulus));
	pfread(fd, &read_list_size, sizeof(read_list_size));
	key_list = new key_type[modulus];
	value_list = new small_value_type[modulus];
	read_offset_list = new read_offset_type[modulus];
	read_list = new read_type[read_list_size];
	pfread(fd, key_list, modulus * sizeof(key_type));
	pfread(fd, value_list, modulus * sizeof(small_value_type));
	pfread(fd, read_offset_list, modulus * sizeof(read_offset_type));
	pfread(fd, read_list, read_list_size * sizeof(read_type));
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
}

// for debugging - print entire hash
void hash_read_hits::print_hash(const KmerLookupInfo &kmers) const {
	offset_type i(0);
	for (; i != modulus; ++i) {
		if (key_list[i] != INVALID_KEY) {
			value_type n(value_list[i]);
			if (n == max_small_value) {
				// use find() to avoid inserting a value into value_map
				const std::map<offset_type, value_type>::const_iterator a(value_map.find(i));
				if (a != value_map.end()) {
					n += a->second;
				}
			}
			std::cout << convert_key(key_list[i]) << ' ' << n << '\n';
			read_offset_type j(read_offset_list[i]);
			for (n += j; j != n; ++j) {
				std::cout << '\t' << kmers.read_name(read_list[j]) << '\n';
			}
		}
	}
}
