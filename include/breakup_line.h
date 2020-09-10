#ifndef _BREAKUP_LINE_H
#define _BREAKUP_LINE_H

#include <string>	// string
#include <vector>	// vector<>

extern void breakup_line(const std::string &, std::vector<std::string> &);
extern void breakup_line(const std::string &, std::vector<std::string> &, char, int = 1);
extern void breakup_line_exact(const std::string &, const std::string &, std::vector<std::string> &);
extern void breakup_line_quoted(const std::string &, std::vector<std::string> &);

#endif // !_BREAKUP_LINE_H
