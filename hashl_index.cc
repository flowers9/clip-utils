#include "hashl_index.h"
#include <string>	// string

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
}

const_iterator hashl_index::find(const key_type &key) const {
}

const_iterator hashl_index::find(const key_type &key, const key_type &comp_key) const {
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
}
