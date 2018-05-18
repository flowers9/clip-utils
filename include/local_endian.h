#ifndef _LOCAL_ENDIAN_H
#define _LOCAL_ENDIAN_H

#include <sys/types.h>
#include <sys/param.h>
#ifdef _BIG_ENDIAN
     #define big_endian
#endif
#ifdef _LITTLE_ENDIAN
     #define little_endian
#endif
#ifdef __BYTE_ORDER
     #if __BYTE_ORDER == __BIG_ENDIAN
          #define big_endian
     #endif
     #if __BYTE_ORDER == __LITTLE_ENDIAN
          #define little_endian
     #endif
#endif
#ifdef BYTE_ORDER
     #if BYTE_ORDER == BIG_ENDIAN
          #define big_endian
     #endif
     #if BYTE_ORDER == LITTLE_ENDIAN
          #define little_endian
     #endif
#endif

#endif	// !_LOCAL_ENDIAN_H
