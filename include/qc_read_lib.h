#ifndef _QC_READ_LIB_H
#define _QC_READ_LIB_H

#include "qc_read.h"	/* QC_Read */
#include <list>		/* list<> */
#include <map>		/* map<> */
#include <string>	/* string */
#include <sys/types.h>	/* size_t */

extern void qc_calc_stats(std::list<QC_Read> &, std::map<size_t, unsigned int> &, std::map<size_t, unsigned int> &);
extern int qc_read_sequence(const char *, std::list<QC_Read> &, bool);

#endif /* !_QC_READ_LIB_H */
