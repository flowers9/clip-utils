// Open a filedes to SERVICE at ADDRESS.  If SERVICE is the name of a
// service, then it must exist on the local machine.  SERVICE can also
// be the ASCII representation of a decimal number, in which case it is
// interpreted as the port number to connect to.  Returns a valid file
// descriptor if successful, or -1 if not.

#ifndef _OPEN_TCP_SERVICE_H
#define _OPEN_TCP_SERVICE_H

#include <sys/types.h>	// required: netinet/in.h
#include <netinet/in.h>	// struct sockaddr_in

extern int open_tcp_service(const char *, const char *, struct sockaddr_in *);

#endif // !_OPEN_TCP_SERVICE_H
