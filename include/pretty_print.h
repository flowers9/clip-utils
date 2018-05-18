#ifndef _PRETTY_PRINT_H
#define _PRETTY_PRINT_H

#include <string>	// string
#include <sstream>	// ostringstream

extern std::string pretty_print(const std::string &);
extern std::string pretty_print(const std::ostringstream &);

// convenience conversion template

template<typename _T>
inline std::string pretty_print(_T x) {
	std::ostringstream s;
	s << x;
	return pretty_print(s.str());
}

#endif // _PRETTY_PRINT_H
