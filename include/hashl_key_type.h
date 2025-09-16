#ifndef _HASHL_KEY_TYPE_H
#define _HASHL_KEY_TYPE_H

#include <string>	// string
#include <vector>	// vector<>

template<class T>
class hashl_key_type {
    private:
	typedef typename T::base_type base_type;
	typedef typename T::data_offset_type data_offset_type;
	typedef typename std::vector<base_type>::size_type size_type;
    private:
	std::vector<base_type> k;		// stored in reverse - high word in [0]
    public:
	const std::vector<base_type> &value() const {
		return k;
	}
	bool operator==(const hashl_key_type &__a) const {
		for (size_type __i(0); __i < k.size(); ++__i) {
			if (k[__i] != __a.k[__i]) {
				return 0;
			}
		}
		return 1;
	}
	bool operator!=(const hashl_key_type &__a) const {
		return !(*this == __a);
	}
	base_type hash() const noexcept {
		base_type __x(k[0]);
		for (size_type __i(1); __i < k.size(); ++__i) {
			__x ^= k[__i];
		}
		return __x;
	}
	int basepair(const size_type __i) const {
		const size_type __n(__i / (sizeof(base_type) * 8));
		return (k[k.size() - 1 - __n] >> (__i - __n * sizeof(base_type) * 8)) & 3;
	}
    private:
	const size_type bit_shift;	// precalc for push_front
	const base_type high_mask;	// precalc for push_back
    public:
	explicit hashl_key_type(const T &__a) : k(__a.words(), 0), bit_shift((__a.bits() - 2) % (sizeof(base_type) * 8)), high_mask(static_cast<base_type>(-1) >> (sizeof(base_type) * 8 - __a.bits() % (sizeof(base_type) * 8))) { }
	~hashl_key_type() { }
	bool operator<(const hashl_key_type &__a) const {
		for (size_type __i(0); __i != k.size(); ++__i) {
			if (k[__i] != __a.k[__i]) {
				return k[__i] < __a.k[__i];
			}
		}
		return 0;
	}
	void push_back(const base_type __x) {
		const size_type __n(sizeof(base_type) * 8 - 2);
		size_type __i(0);
		for (; __i < k.size() - 1; ++__i) {
			k[__i] = (k[__i] << 2) | (k[__i + 1] >> __n);
		}
		k[__i] = (k[__i] << 2) | __x;
		k[0] &= high_mask;
	}
	void push_front(const base_type __x) {
		const size_type __n(sizeof(base_type) * 8 - 2);
		for (size_type __i(k.size() - 1); __i > 0; --__i) {
			k[__i] = (k[__i - 1] << __n) | (k[__i] >> 2);
		}
		k[0] = (__x << bit_shift) | (k[0] >> 2);
	}

	// make reverse complement of given key (TODO: could be more efficient)
	void make_complement(const hashl_key_type &key) {
		const size_type bit_width(bit_shift + 2 + (k.size() - 1) * sizeof(base_type) * 8);
		for (size_type i(0); i < bit_width; i += 2) {
			push_back(3 - key.basepair(i));
		}
	}

	void convert_to_string(std::string &sequence) const {
		const char values[4] = { 'A', 'C', 'G', 'T' };
		sequence.clear();
		const size_type bit_width(bit_shift + 2 + (k.size() - 1) * sizeof(base_type) * 8);
		// relies on wrap-around for termination
		for (size_type i(bit_width - 2); i < bit_width; i -= 2) {
			sequence += values[basepair(i)];
		}
	}

	// create key from bit offset into data
	void copy_in(const std::vector<base_type> &data, const data_offset_type offset) {
		// start of sequence in data
		const size_type i(offset / (sizeof(base_type) * 8));
		// how many bits we have in the first word
		const base_type starting_bits(sizeof(base_type) * 8 - offset % (sizeof(base_type) * 8));
		// how many bits the first word is supposed to have for a key
		const base_type high_bits(bit_shift + 2);
		if (starting_bits == high_bits) {
			k[0] = data[i] & high_mask;
			for (size_type j(1); j < k.size(); ++j) {
				k[j] = data[i + j];
			}
		} else if (starting_bits < high_bits) {		// shift left to fill up first word
			const int shift_left(high_bits - starting_bits);
			const int shift_right(sizeof(base_type) * 8 - shift_left);
			k[0] = ((data[i] << shift_left) | (data[i + 1] >> shift_right)) & high_mask;
			for (size_type j(1); j < k.size(); ++j) {
				k[j] = (data[i + j] << shift_left) | (data[i + j + 1] >> shift_right);
			}
		} else {					// shift right to empty out first word
			const int shift_right(starting_bits - high_bits);
			const int shift_left(sizeof(base_type) * 8 - shift_right);
			k[0] = (data[i] >> shift_right) & high_mask;
			for (size_type j(1); j < k.size(); ++j) {
				k[j] = (data[i + j - 1] << shift_left) | (data[i + j] >> shift_right);
			}
		}
	}

	// same as copy_in(), but with breakpoints and not saving the generated key
	bool equal(const std::vector<base_type> &data, const data_offset_type offset) const {
		const size_type i(offset / (sizeof(base_type) * 8));
		const base_type starting_bits(sizeof(base_type) * 8 - offset % (sizeof(base_type) * 8));
		const base_type high_offset(bit_shift + 2);
		if (starting_bits == high_offset) {
			if (k[0] != (data[i] & high_mask)) {
				return 0;
			}
			for (size_type j(1); j < k.size(); ++j) {
				if (k[j] != data[i + j]) {
					return 0;
				}
			}
		} else if (starting_bits < high_offset) {
			const int shift_left(high_offset - starting_bits);
			const int shift_right(sizeof(base_type) * 8 - shift_left);
			if (k[0] != (((data[i] << shift_left) | (data[i + 1] >> shift_right)) & high_mask)) {
				return 0;
			}
			for (size_type j(1); j < k.size(); ++j) {
				if (k[j] != ((data[i + j] << shift_left) | (data[i + j + 1] >> shift_right))) {
					return 0;
				}
			}
		} else {
			const int shift_right(starting_bits - high_offset);
			const int shift_left(sizeof(base_type) * 8 - shift_right);
			if (k[0] != ((data[i] >> shift_right) & high_mask)) {
				return 0;
			}
			for (size_type j(1); j < k.size(); ++j) {
				if (k[j] != ((data[i + j - 1] << shift_left) | (data[i + j] >> shift_right))) {
					return 0;
				}
			}
		}
		return 1;
	}
};

// a hash for broken-out keys (i.e., just the plain vector)

template<class T>
class hashl_key_hash {
    private:
	typedef typename T::base_type base_type;
	typedef typename std::vector<base_type>::size_type size_type;
    public:
	size_type operator()(const std::vector<base_type> &k) const noexcept {
		base_type __x(k[0]);
		for (size_type __i(1); __i < k.size(); ++__i) {
			__x ^= k[__i];
		}
		return __x;
	}
};

#endif // !_HASHL_KEY_TYPE_H
