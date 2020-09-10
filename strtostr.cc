#include "strtostr.h"
#include <ctype.h>	// isspace()
#include <stdio.h>	// fprintf(), stderr
#include <stdlib.h>	// NULL
#include <string>	// string

// return next (whitespace delimited) word of string, starting at index,
// and update index to point past end of returned word

std::string strtostr(const std::string &s, std::string::size_type * const index) {
	std::string::size_type i(index == NULL ? 0 : *index);
	if (i > s.size()) {		// perform input validation
		return "";
	}
	for (; i != s.size() && isspace(s[i]); ++i) { }
	if (i == s.size()) {		// all whitespace
		if (index != NULL) {
			*index = i;
		}
		return "";
	}
	std::string::size_type j(i + 1);
	for (; j != s.size() && !isspace(s[j]); ++j) { }
	if (index != NULL) {
		*index = j;
	}
	return s.substr(i, j - i);
}

// same as strtostr, but with a delimiter specified instead of whitespace;
// trim leading and trailing whitespace if specified
std::string strtostr(const std::string &s, std::string::size_type * const index, const char delim, const int trim_whitespace) {
	std::string::size_type i(index == NULL ? 0 : *index);
	if (i > s.size()) {		// perform input validation
		return "";
	}
	for (; i != s.size() && (s[i] == delim || (trim_whitespace && isspace(s[i]))); ++i) { }
	if (i == s.size()) {		// all whitespace
		if (index != NULL) {
			*index = i;
		}
		return "";
	}
	std::string::size_type j(i + 1);
	for (; j != s.size() && s[j] != delim; ++j) { }
	if (index != NULL) {
		*index = j;
	}
	if (trim_whitespace) {		// trim trailing whitespace
		for (--j; isspace(s[j]); --j) { }
		++j;
	}
	return s.substr(i, j - i);
}

/*
 * same as strtostr(), except takes a list of whitespace characters and
 * allows for empty elements
 */

std::string strtostr_exact(const std::string &s, const std::string &whitespace, std::string::size_type * const index) {
	std::string::size_type i(index == NULL ? 0 : *index);
	if (i > s.size()) {		// perform input validation
		return "";
	}
	if (i != s.size() && whitespace.find(s[i]) != std::string::npos) {
		++i;
	}
	// check for null entry
	if (i == s.size() || whitespace.find(s[i]) != std::string::npos) {
		if (index != NULL) {
			*index = i;
		}
		return "";
	}
	std::string::size_type j(i + 1);
	for (; j != s.size() && whitespace.find(s[j]) == std::string::npos; ++j) { }
	if (index != NULL) {
		*index = j;
	}
	return s.substr(i, j - i);
}

// like strtostr(), except, if a string begins with a ", it continues until
// the next " (thus allowing whitespace to be included); "'s and \'s may be
// included inside a quoted string by escaping them with \.

std::string strtostr_quoted(const std::string &s, std::string::size_type * const index) {
	std::string::size_type i(index == NULL ? 0 : *index);
	if (i > s.size()) {		// perform input validation
		return "";
	}
	for (; i != s.size() && isspace(s[i]); ++i) { }
	if (i == s.size()) {		// all whitespace
		if (index != NULL) {
			*index = i;
		}
		return "";
	}
	std::string::size_type j(i + 1);
	if (s[i] == '"') {		// go to next (unescaped) quote
		++i;			// skip leading "
		std::string t;
		int escaped(0);
		for (; j != s.size(); ++j) {
			if (escaped) {
				escaped = 0;
			} else if (s[j] == '\\') {
				t += s.substr(i, j - i);
				i = j + 1;		// skip escaping \ ;
				escaped = 1;
			} else if (s[j] == '"') {	// end of string
				break;
			}
		}
		if (i != j) {
			t += s.substr(i, j - i);
		}
		if (j != s.size()) {
			++j;				// skip trailing "
		} else {
			fprintf(stderr, "Warning: unterminated string %s\n", s.substr(index == NULL ? 0 : *index).c_str());
		}
		if (index != NULL) {
			*index = j;
		}
		return t;
	} else {
		for (; j != s.size() && !isspace(s[j]); ++j) { }
		if (index != NULL) {
			*index = j;
		}
		return s.substr(i, j - i);
	}
}
