#include "hashl_index.h"
#include "itoa.h"	// itoa()
#include "local_endian.h"	// big_endian
#include "open_compressed.h"	// pfread()
#include "write_fork.h"	// pfwrite()
#include <iomanip>	// setw()
#include <iostream>	// cerr, cout
#include <stdlib.h>	// exit()
#include <string.h>	// memcmp()
#include <string>	// string
#include <vector>	// vector<>

// description beginning of saved file

std::string hashl_index::boilerplate() const {
	std::string s("hashl_index\n");
	s += itoa(sizeof(base_type));
	s += " bytes\n";
#ifdef big_endian
	s += "big endian\n";
#else
	s += "little endian\n";
#endif
	return s;
}

hashl::hashl(const int fd) {
	const std::string s(boilerplate());
	char t[s.size()];
	pfread(fd, t, s.size());
	if (memcmp(t, s.c_str(), s.size()) != 0) {
		std::cerr << "Error: could not read index from file: header mismatch\n";
		exit(1);
	}
	pfread(fd, &bit_width, sizeof(bit_width));
	word_width = (bit_width + 8 * sizeof(base_type) - 1) / (8 * sizeof(base_type));
	size_type metadata_size;
	pfread(fd, &metadata_size, sizeof(metadata_size));
	metadata.assign(metadata_size, 0);
	pfread(fd, &metadata[0], metadata_size);
	size_type data_size;
	pfread(fd, &data_size, sizeof(data_size));
	data.assign(data_size, 0);
	pfread(fd, &data[0], sizeof(base_type) * data_size);
	size_type key_list_size;
	pfread(fd, &key_list_size, sizeof(key_list_size));
	key_list.assign(key_list_size, 0);
	// XXX - except, of course, we want to mmap here, not read it in
	pfread(fd, &key_list[0], sizeof(data_offset_type) * key_list_size);
}

// XXX
bool hashl_index::exists(const key_type &key, const key_type &comp_key) const {
	return std::binary_search(key_list.begin(), key_list.end(), [this](const data_offset_type __i, const data_offset_type __j){return key.equal(this->data, this->key_list[__i] || comp_key.equal(this->data, this->key_list[__i])});
}

bool hashl_index::exists(const key_type &key) const {
	key_type comp_key(*this);
	comp_key.make_complement(key);
	return exists(key, comp_key);
}

void hashl_index::get_sequence(const data_offset_type start, const data_offset_type length, std::string &) const {
	const char values[4] = { 'A', 'C', 'G', 'T' };
	seq.clear();
	size_t word_offset(start / (sizeof(base_type) * 8));
	size_t bit_offset(sizeof(base_type) * 8 - start % (sizeof(base_type) * 8));
	for (data_offset_type i(0); i < length; i += 2) {
		if (bit_offset) {
			bit_offset -= 2;
		} else {
			bit_offset = sizeof(base_type) * 8 - 2;
			++word_offset;
		}
		seq += values[(data[word_offset] >> bit_offset) & 3];
	}
}

void hashl_index::print() const {
	int max_offset_width(1), max_key_width(1);
	for (size_type i(10); i < key_list.size(); i *= 10, ++max_offset_width) { }
	for (size_type i(10); i < data.size() * sizeof(base_type) * 8; i *= 10, ++max_key_width) { }
	std::cout << "elements: " << key_list.size() << "\n"
		<< "bit width: " << bit_width << "\n"
		<< "metadata size: " << metadata.size() << "\n"
		<< "data size: " << data.size() * sizeof(base_type) << "\n"
		<< "offset/key pairs:\n";
	std::string s;
	key_type k(*this);
	for (size_type i(0); i < key_list.size(); ++i) {
		if (key_list[i] != invalid_key) {
			k.copy_in(data, key_list[i]);
			k.convert_to_string(s);
			std::cout << std::setw(max_offset_width) << i << ' ' << std::setw(max_key_width) << key_list[i] << ' ' << s << "\n";
		}
	}
}

static void hashl_index::save(const std::vector<data_offset_type> &key_list_in, const std::vector<base_type> &data_in, const std::vector<char> &metadata_in, const size_type bit_width_in, const size_type word_width_in, const int fd_in) {
	const std::string s(boilerplate());
	pfwrite(fd, s.c_str(), s.size());
	pfwrite(fd, &bit_width, sizeof(bit_width));
	size_type tmp;
	pfwrite(fd, &(tmp = metadata.size()), sizeof(tmp));
	pfwrite(fd, &metadata[0], metadata.size());
	pfwrite(fd, &(tmp = data.size()), sizeof(tmp));
	pfwrite(fd, &data[0], sizeof(base_type) * data.size());
	pfwrite(fd, &(tmp = key_list.size()), sizeof(tmp));
	pfwrite(fd, &key_list[0], sizeof(data_offset_type) * key_list.size());
}
