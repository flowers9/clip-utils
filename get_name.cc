#include "get_name.h"
#include <string>	// string

// get the read name from a header line - everything following the > at the
// beginning of the line to the first space (or the end of the line)

std::string get_name(const std::string &header) {
	size_t i(header.find_first_of(" \n", 1));
	if (i == std::string::npos) {
		i = header.length();
	}
	return header.substr(1, i - 1);		// skip leading >
}
