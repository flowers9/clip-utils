#ifndef _STRTOSTR_H
#define _STRTOSTR_H

#include <stdlib.h>	/* NULL */
#include <string>	/* string */

extern std::string strtostr(const std::string &, std::string::size_type * const = NULL);
extern std::string strtostr(const std::string &, std::string::size_type * const, char, int = 1);
extern std::string strtostr_exact(const std::string &, const std::string &, std::string::size_type * const = NULL);
extern std::string strtostr_quoted(const std::string &, std::string::size_type * const = NULL);

#endif /* !_STRTOSTR_H */
