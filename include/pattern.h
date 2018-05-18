#ifndef _PATTERN_H
#define _PATTERN_H

#include "refcount_array.h"	// refcount_array<>
#include <regex.h>	// regex_t, regmatch_t
#include <string>	// string
#include <sys/types.h>	// size_t

class Pattern {
    public:
	typedef refcount_array<regmatch_t>::const_reference value_type;
    private:
	std::string regex_;
	int is_valid_;	// 0 - invalid, 1 - normal re, 2 - empty re
	regex_t pattern_;
	int cflags_;
	refcount_array<regmatch_t> pmatch_;
    public:
	Pattern(void) : is_valid_(0) { }
	explicit Pattern(const std::string &, size_t = 0, int = 0);
	Pattern(const Pattern &);
	~Pattern(void);
	Pattern &operator=(const Pattern &);
	bool initialize(const std::string &, size_t = 0, int = 0);
	bool is_match(const std::string &, int = 0);
	value_type operator[](const size_t __n) const {
		return pmatch_[__n];
	}
	bool operator<(const Pattern &a) const {
		return regex_ < a.regex_;
	}
	bool empty(void) const {
		return is_valid_ == 0;
	}
	size_t subs(void) const {
		return pmatch_.size();
	}
	const std::string &regex(void) const {
		return regex_;
	};
};

#endif // !_PATTERN_H
