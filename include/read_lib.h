#ifndef _READ_LIB_H
#define _READ_LIB_H

#include "read.h"	/* Read */
#include <list>		/* list<> */
#include <map>		/* map<> */
#include <string>	/* string */
#include <vector>	/* vector<> */

extern bool opt_strip_tracename;
extern std::map<std::string, bool> opt_readname_match;

extern std::string make_read_name(std::string &);
extern std::string make_qual_filename(const char *, bool = 0);
extern void mask_by_phred(std::list<Read> &, unsigned int);
extern int read_sequence(const char *, std::list<Read> &, bool);

#endif /* !_READ_LIB_H */
