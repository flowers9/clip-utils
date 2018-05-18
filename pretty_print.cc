#include "pretty_print.h"
#include <string>	// string
#include <sstream>	// ostringstream

// adds commas to a string containing a number

std::string pretty_print(const std::string &s) {
	const std::string::size_type n(s.find('.'));
	std::string::size_type digits(n != std::string::npos ? n : s.size());
	if (digits != 0 && s[0] == '-') {
		--digits;
	}
	if (digits < 4) {
		return s;
	}
	std::string::size_type a(s.size() - 1);
	std::string::size_type b(a + (digits - 1) / 3);
	std::string t;
	if (n != std::string::npos) {
		t.resize(b + 1 - s.size() + n);
		t += s.substr(n);
	} else {
		t.resize(b + 1);
	}
	int i(0);
	for (; a != 0; --a, --b, ++i) {
		t[b] = s[a];
		if (i == 2) {
			i = -1;
			t[--b] = ',';
		}
	}
	t[0] = s[0];	// to avoid possible extraneous , at front
	return t;
}

std::string pretty_print(const std::ostringstream &s) {
	return pretty_print(s.str());
}
