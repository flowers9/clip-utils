#include "pattern.h"
#include "refcount_array.h"	// refcount_array<>
#include <regex.h>	// regexec(), regerror(), regfree(), regmatch_t
#include <stdio.h>	// fprintf(), printf(), stderr
#include <string>	// string
#include <sys/types.h>	// size_t

static int allocate_pattern(regex_t * const pattern, const std::string &regex_in, const int cflags) {
	if (regex_in.empty()) {
		return 2;
	}
	const int i(regcomp(pattern, regex_in.c_str(), cflags));
	if (i == 0) {
		return 1;
	}
	const size_t j(regerror(i, pattern, 0, 0));
	if (j == 0) {
		fprintf(stderr, "Error: regcomp: regex error %d: %s\n", i, regex_in.c_str());
	} else {
		refcount_array<char> err(j);
		regerror(i, pattern, err.array(), j);
		fprintf(stderr, "Error: regcomp: regex error %s: %s\n", err.array(), regex_in.c_str());
	}
	return 0;
}

Pattern::Pattern(const std::string &regex_in, const size_t subs_in, const int cflags_in) {
	if ((is_valid_ = allocate_pattern(&pattern_, regex_in, cflags_in))) {
		regex_ = regex_in;
		cflags_ = cflags_in;
		pmatch_.resize(subs_in);
	}
}

Pattern::Pattern(const Pattern &a) {
	if (a.is_valid_ == 0) {
		is_valid_ = 0;
	} else if ((is_valid_ = allocate_pattern(&pattern_, a.regex_, a.cflags_))) {
		regex_ = a.regex_;
		cflags_ = a.cflags_;
		pmatch_.resize(a.pmatch_.size());
	}
}

Pattern::~Pattern() {
	if (is_valid_ == 1) {
		regfree(&pattern_);
	}
}

Pattern &Pattern::operator=(const Pattern &a) {
	if (is_valid_ == 1) {
		regfree(&pattern_);
	}
	if (a.is_valid_ == 0) {
		is_valid_ = 0;
	} else if ((is_valid_ = allocate_pattern(&pattern_, a.regex_, a.cflags_))) {
		regex_ = a.regex_;
		cflags_ = a.cflags_;
		pmatch_.resize(a.pmatch_.size());
	}
	return *this;
}

bool Pattern::initialize(const std::string &regex_in, const size_t subs_in, const int cflags_in) {
	if (is_valid_ == 1) {
		regfree(&pattern_);
	}
	if ((is_valid_ = allocate_pattern(&pattern_, regex_in, cflags_in))) {
		regex_ = regex_in;
		cflags_ = cflags_in;
		pmatch_.resize(subs_in);
	}
	return is_valid_ != 0;
}

bool Pattern::is_match(const std::string &s, const int eflags) {
	return is_valid_ == 2 || (is_valid_ == 1 && regexec(&pattern_, s.c_str(), pmatch_.size(), pmatch_.array(), eflags) == 0);
}
