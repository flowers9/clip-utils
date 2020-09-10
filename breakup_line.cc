#include "breakup_line.h"
#include "strtostr.h"	// strtostr(), strtostr_quoted()
#include <string>	// string
#include <vector>	// vector<>

// break up a line into whitespace separated elements

void breakup_line(const std::string &line, std::vector<std::string> &list) {
	std::string::size_type i(0);
	std::string::size_type j(i);
	std::string s = strtostr(line, &i);
	while (j != i) {
		j = i;
		list.push_back(s);
		s = strtostr(line, &i);
	}
}

// like breakup_line, but using a given delimiter instead of whitespace
void breakup_line(const std::string &line, std::vector<std::string> &list, const char delim, const int trim_whitespace) {
	std::string::size_type i(0);
	std::string::size_type j(i);
	std::string s = strtostr(line, &i, delim, trim_whitespace);
	while (j != i) {
		j = i;
		list.push_back(s);
		s = strtostr(line, &i, delim, trim_whitespace);
	}
}

// same as above, but takes a list of whitespace characters and allows
// for empty elements

void breakup_line_exact(const std::string &line, const std::string &whitespace, std::vector<std::string> &list) {
	std::string::size_type i(0);
	std::string::size_type j(i);
	std::string s = strtostr_exact(line, whitespace, &i);
	while (j != i) {
		j = i;
		list.push_back(s);
		s = strtostr_exact(line, whitespace, &i);
	}
}

// same as above, but allow whitespace to be quoted

void breakup_line_quoted(const std::string &line, std::vector<std::string> &list) {
	std::string::size_type i(0);
	std::string::size_type j(i);
	std::string s = strtostr_quoted(line, &i);
	while (j != i) {
		j = i;
		list.push_back(s);
		s = strtostr_quoted(line, &i);
	}
}
