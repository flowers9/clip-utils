#ifndef _HASHL_LESS_H
#define _HASHL_LESS_H

#include <vector>	// vector<>

// given an array and bit-length, compare two bit offsets into it

// this is a little tricky, as the offsets generally don't line up with the data storage
// (c.f., hashl_key_type::equal(), except in that case the internal k vector does align)

// note: It's possible this would be more efficient with only one path that
// shifted both values rather than three paths that minimize the shifting
// (and it would certainly be shorter ;)

template<class T>
class hashl_less {
    private:
	typedef typename T::base_type base_type;
	typedef typename T::hash_offset_type hash_offset_type;
	typedef typename std::vector<base_type>::size_type size_type;
    public:
	bool operator()(const T &blob, const hash_offset_type &a, const hash_offset_type &b) const {
		const size_type bit_width = blob.bits();
		const std::vector<base_type> &data = blob.get_data();
		const size_type a_i = a / (sizeof(base_type) * 8);
		const unsigned int a_starting_bit = sizeof(base_type) * 8 - a % (sizeof(base_type) * 8);
		const size_type b_i = b / (sizeof(base_type) * 8);
		const unsigned int b_starting_bit = sizeof(base_type) * 8 - b % (sizeof(base_type) * 8);
		// breaking out "==" avoids a few bit shifts at the cost of more code
		if (a_starting_bit == b_starting_bit) {
			// compare starting bits
			base_type mask = static_cast<base_type>(-1) >> (sizeof(base_type) * 8 - a_starting_bit);
			if (a_starting_bit == bit_width) {
				return (data[a_i] & mask) < (data[b_i] & mask);
			} else if (a_starting_bit > bit_width) {
				// if the kmer doesn't reach the right side of word, adjust mask
				mask ^= static_cast<base_type>(-1) >> (sizeof(base_type) * 8 - a_starting_bit + bit_width);
				return (data[a_i] & mask) < (data[b_i] & mask);
			} else if ((data[a_i] & mask) != (data[b_i] & mask)) {
				return (data[a_i] & mask) < (data[b_i] & mask);
			}
			// compare all full words
			const size_type words = (bit_width - a_starting_bit) / (sizeof(base_type) * 8) + 1;
			for (size_type j(1); j < words; ++j) {
				if (data[a_i + j] != data[b_i + j]) {
					return data[a_i + j] < data[b_i + j];
				}
			}
			// compare any trailing bits
			const size_type trailing_bit = (bit_width - a_starting_bit) % (sizeof(base_type) * 8);
			if (trailing_bit == 0) {
				return 0;
			}
			mask = static_cast<base_type>(-1) << (sizeof(base_type) * 8 - trailing_bit);
			return (data[a_i + words] & mask) < (data[b_i + words] & mask);
		} else if (a_starting_bit < b_starting_bit) {		// shift b right to align with a
			// compare starting bits
			const unsigned int shift_right = b_starting_bit - a_starting_bit;
			base_type mask = static_cast<base_type>(-1) >> (sizeof(base_type) * 8 - a_starting_bit);
			if (a_starting_bit == bit_width) {
				return (data[a_i] & mask) < ((data[b_i] >> shift_right) & mask);
			} else if (a_starting_bit > bit_width) {
				// if the kmer doesn't reach the right side of word, adjust mask
				mask ^= static_cast<base_type>(-1) >> (sizeof(base_type) * 8 - a_starting_bit + bit_width);
				return (data[a_i] & mask) < ((data[b_i] >> shift_right) & mask);
			} else if ((data[a_i] & mask) != ((data[b_i] >> shift_right) & mask)) {
				return (data[a_i] & mask) < ((data[b_i] >> shift_right) & mask);
			}
			// compare all full words
			const unsigned int shift_left = sizeof(base_type) * 8 - shift_right;
			const size_type words = (bit_width - a_starting_bit) / (sizeof(base_type) * 8) + 1;
			for (size_type j(1); j < words; ++j) {
				if (data[a_i + j] != ((data[b_i + j - 1] << shift_left) | (data[b_i + j] >> shift_right))) {
					return data[a_i + j] < ((data[b_i + j - 1] << shift_left) | (data[b_i + j] >> shift_right));
				}
			}
			// compare any trailing bits
			const size_type trailing_bit = (bit_width - a_starting_bit) % (sizeof(base_type) * 8);
			if (trailing_bit == 0) {
				return 0;
			}
			mask = static_cast<base_type>(-1) << (sizeof(base_type) * 8 - trailing_bit);
			return data[a_i + words] < ((data[b_i + words - 1] << shift_left) | (data[b_i + words] >> shift_right));
		} else {						// shift a right to align with b
			// compare starting bits
			const unsigned int shift_right = a_starting_bit - b_starting_bit;
			base_type mask = static_cast<base_type>(-1) >> (sizeof(base_type) * 8 - b_starting_bit);
			if (b_starting_bit == bit_width) {
				return ((data[a_i] >> shift_right) & mask) < (data[b_i] & mask);
			} else if (b_starting_bit > bit_width) {
				// if the kmer doesn't reach the right side of word, adjust mask
				mask ^= static_cast<base_type>(-1) >> (sizeof(base_type) * 8 - b_starting_bit + bit_width);
				return ((data[a_i] >> shift_right) & mask) < (data[b_i] & mask);
			} else if ((data[b_i] & mask) != ((data[a_i] >> shift_right) & mask)) {
				return ((data[a_i] >> shift_right) & mask) < (data[b_i] & mask);
			}
			// compare all full words
			const unsigned int shift_left = sizeof(base_type) * 8 - shift_right;
			const size_type words = (bit_width - b_starting_bit) / (sizeof(base_type) * 8) + 1;
			for (size_type j(1); j < words; ++j) {
				if (data[b_i + j] != ((data[a_i + j - 1] << shift_left) | (data[a_i + j] >> shift_right))) {
					return ((data[a_i + j - 1] << shift_left) | (data[a_i + j] >> shift_right)) < data[b_i + j];
				}
			}
			// compare any trailing bits
			const size_type trailing_bit = (bit_width - b_starting_bit) % (sizeof(base_type) * 8);
			if (trailing_bit == 0) {
				return 0;
			}
			mask = static_cast<base_type>(-1) << (sizeof(base_type) * 8 - trailing_bit);
			return ((data[a_i + words - 1] << shift_left) | (data[a_i + words] >> shift_right)) < data[b_i + words];
		}
	}
};

#endif // !_HASHL_LESS_H
