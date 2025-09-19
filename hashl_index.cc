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

std::string hashl_index::boilerplate() {
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

hashl_index::hashl_index(const int fd) {
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
	pfread(fd, &key_list[0], sizeof(size_type) * key_list_size);
}

// checks for the existence of key or its reverse complement
// (while hashl only stores hashes for key values of key < comp_key,
// our vector can only record the values that actually show up in the
// sequence, whichever one it is)

// returns -1 if kmer is not found

// note: might be faster with a tri-value comparison (-1, 0, 1)?

hashl_index::size_type hashl_index::position(const key_type &key) const {
	// do a binary search, comparing key to kmers at data offsets given by array values
	size_type i = 0, j = key_list.size();
	while (i + 1 < j) {
		const size_type m = (i + j) / 2;
		(key.less_than(data, key_list[m]) ? j : i) = m;
	}
	if (key.equal_to(data, key_list[i])) {
		return key_list[i];
	}
	// and now check the reverse complement
	key_type comp_key(bit_width, word_width);
	comp_key.make_complement(key);
	for (i = 0, j = key_list.size(); i + 1 < j;) {
		const size_type m = (i + j) / 2;
		(comp_key.less_than(data, key_list[m]) ? j : i) = m;
	}
	return comp_key.equal_to(data, key_list[i]) ? key_list[i] : -1;
}

void hashl_index::get_sequence(const size_type start, const size_type length, std::string &seq) const {
	const char values[4] = { 'A', 'C', 'G', 'T' };
	seq.clear();
	size_t word_offset(start / (sizeof(base_type) * 8));
	size_t bit_offset(sizeof(base_type) * 8 - start % (sizeof(base_type) * 8));
	for (size_type i(0); i < length; i += 2) {
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
	for (size_type i = 10; i < key_list.size(); i *= 10, ++max_offset_width) { }
	for (size_type i = 10; i < data.size() * sizeof(base_type) * 8; i *= 10, ++max_key_width) { }
	std::cout << "elements: " << key_list.size() << "\n"
		<< "bit width: " << bit_width << "\n"
		<< "metadata size: " << metadata.size() << "\n"
		<< "data size: " << data.size() * sizeof(base_type) << "\n"
		<< "offset/key pairs:\n";
	std::string s;
	key_type k(bit_width, word_width);
	for (size_type i = 0; i < key_list.size(); ++i) {
		k.copy_in(data, key_list[i]);
		k.get_sequence(s);
		std::cout << std::setw(max_offset_width) << i << ' ' << std::setw(max_key_width) << key_list[i] << ' ' << s << "\n";
	}
}

void hashl_index::save(const std::vector<size_type> &key_list_in, const std::vector<base_type> &data_in, const std::vector<char> &metadata_in, const size_type bit_width_in, const int fd) {
	const std::string s(boilerplate());
	pfwrite(fd, s.c_str(), s.size());
	pfwrite(fd, &bit_width_in, sizeof(bit_width_in));
	size_type tmp;
	pfwrite(fd, &(tmp = metadata_in.size()), sizeof(tmp));
	pfwrite(fd, &metadata_in[0], metadata_in.size());
	pfwrite(fd, &(tmp = data_in.size()), sizeof(tmp));
	pfwrite(fd, &data_in[0], sizeof(base_type) * data_in.size());
	pfwrite(fd, &(tmp = key_list_in.size()), sizeof(tmp));
	pfwrite(fd, &key_list_in[0], sizeof(size_type) * key_list_in.size());
}
