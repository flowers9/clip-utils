#ifndef _ITOA_H
#define _ITOA_H

#include <string>	// string
#include <sstream>	// ostringstream

template<typename _T>
inline std::string itoa(_T x) {
	std::ostringstream s;
	s << x;
	return s.str();
}

#endif // !_ITOA_H
