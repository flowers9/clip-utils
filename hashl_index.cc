#include "hashl_index.h"
#include "itoa.h"	// itoa()
#include "local_endian.h"	// big_endian
#include "open_compressed.h"	// pfread()
#include "write_fork.h"	// pfwrite()
#include <errno.h>	// errno
#include <fcntl.h>	// posix_fadvise(), POSIX_FADV_RANDOM
#include <iomanip>	// setw()
#include <iostream>	// cerr, cout
#include <stdlib.h>	// exit()
#include <string.h>	// memcmp(), strerror()
#include <string>	// string
#include <sys/mman.h>	// mmap(), munmap(), MAP_FAILED, MAP_PRIVATE, PROT_READ
#include <unistd.h>	// sysconf(), _SC_PAGE_SIZE
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

// need to initialize key_list up front to prep destructor
hashl_index::hashl_index(const int fd) : key_list(0) {
	const std::string s(boilerplate());
	char t[s.size()];
	size_type key_list_offset = pfread(fd, t, s.size());
	if (memcmp(t, s.c_str(), s.size()) != 0) {
		std::cerr << "Error: could not read index from file: header mismatch\n";
		exit(1);
	}
	key_list_offset += pfread(fd, &bit_width, sizeof(bit_width));
	word_width = (bit_width + 8 * sizeof(base_type) - 1) / (8 * sizeof(base_type));
	size_type metadata_size;
	key_list_offset += pfread(fd, &metadata_size, sizeof(metadata_size));
	metadata.assign(metadata_size, 0);
	key_list_offset += pfread(fd, &metadata[0], metadata_size);
	size_type data_size;
	key_list_offset += pfread(fd, &data_size, sizeof(data_size));
	data.assign(data_size, 0);
	key_list_offset += pfread(fd, &data[0], sizeof(base_type) * data_size);
	key_list_offset += pfread(fd, &key_list_size, sizeof(key_list_size));
	size_type padding_size;
	key_list_offset += pfread(fd, &padding_size, sizeof(padding_size));
	// we don't have to actually read in the padding ;)
	key_list_offset += padding_size;
	// key_list needs to start on a page boundary, and as created it will, but if
	// it's read on a machine with a different page size we might need to offset
	page_offset = key_list_offset % sysconf(_SC_PAGE_SIZE);
	if (page_offset % sizeof(size_type)) {
		// page_offset has to be a multiple of size_type, but it this would
		// only fail on a really weird machine
		std::cerr << "Error: could not align to start of key_list: " << page_offset << " not a multiple of " << sizeof(size_type) << '\n';
		exit(1);
	}
	// inform system we'll be accessing randomly, so don't be clever about reading ahead;
	// this is only advisory anyway, so don't worry about it failing
	posix_fadvise(fd, key_list_offset, 0, POSIX_FADV_RANDOM);
	key_list = static_cast<const size_type *>(mmap(0, key_list_size * sizeof(size_type) + page_offset, PROT_READ, MAP_PRIVATE, fd, key_list_offset - page_offset));
	if (key_list == MAP_FAILED) {
		std::cerr << "Error: mmap(" << errno << "): " << std::string(strerror(errno)) << '\n';
		key_list = 0;	// prevent munmap() on destruction
		exit(1);
	}
	key_list += page_offset / sizeof(size_type);
}

hashl_index::~hashl_index() {
	if (key_list) {
		munmap(key_list - page_offset / sizeof(size_type), key_list_size * sizeof(size_type) + page_offset);
	}
}

// checks for the existence of key or its reverse complement
// (while hashl only stores hashes for key values of key < comp_key,
// our vector can only record the values that actually show up in the
// sequence, whichever one it is)

// returns -1 if kmer is not found

// note: might be faster with a tri-value comparison (-1, 0, 1)?

hashl_index::size_type hashl_index::position(const key_type &key) const {
	// do a binary search, comparing key to kmers at data offsets given by array values
	size_type i = 0, j = key_list_size;
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
	for (i = 0, j = key_list_size; i + 1 < j;) {
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
	for (size_type i = 10; i < key_list_size; i *= 10, ++max_offset_width) { }
	for (size_type i = 10; i < data.size() * sizeof(base_type) * 8; i *= 10, ++max_key_width) { }
	std::cout << "elements: " << key_list_size << "\n"
		<< "bit width: " << bit_width << "\n"
		<< "metadata size: " << metadata.size() << "\n"
		<< "data size: " << data.size() * sizeof(base_type) << "\n"
		<< "offset/key pairs:\n";
	std::string s;
	key_type k(bit_width, word_width);
	for (size_type i = 0; i < key_list_size; ++i) {
		k.copy_in(data, key_list[i]);
		k.get_sequence(s);
		std::cout << std::setw(max_offset_width) << i << ' ' << std::setw(max_key_width) << key_list[i] << ' ' << s << "\n";
	}
}

void hashl_index::save(const std::vector<size_type> &key_list_in, const std::vector<base_type> &data_in, const std::vector<char> &metadata_in, const size_type bit_width_in, const int fd) {
	const std::string s(boilerplate());
	size_type written = pfwrite(fd, s.c_str(), s.size());
	written += pfwrite(fd, &bit_width_in, sizeof(bit_width_in));
	size_type tmp;
	written += pfwrite(fd, &(tmp = metadata_in.size()), sizeof(tmp));
	written += pfwrite(fd, &metadata_in[0], metadata_in.size());
	written += pfwrite(fd, &(tmp = data_in.size()), sizeof(tmp));
	written += pfwrite(fd, &data_in[0], sizeof(base_type) * data_in.size());
	written += pfwrite(fd, &(tmp = key_list_in.size()), sizeof(tmp));
	// now page align the start of key_list
	written += sizeof(tmp);
	// calculate amount of padding we need
	tmp = (sysconf(_SC_PAGE_SIZE) - written % sysconf(_SC_PAGE_SIZE)) % sysconf(_SC_PAGE_SIZE);
	pfwrite(fd, &tmp, sizeof(tmp));
	char buf[tmp] = {0};
	pfwrite(fd, buf, tmp);
	pfwrite(fd, &key_list_in[0], sizeof(size_type) * key_list_in.size());
}
