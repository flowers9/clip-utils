#ifndef _READ_H
#define _READ_H

#ifdef COMPRESS_READS
#include "read_c.h"
#else
#include "read_u.h"
#endif

#endif // !_READ_H
