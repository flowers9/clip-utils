#include "open_tcp_service.h"
#include <ctype.h>	// isdigit()
#include <netdb.h>	// struct servent, getservbyname()
#include <netinet/in.h>	// struct sockaddr_in, htons()
#include <socket.h>	// socket()
#include <stdio.h>	// fprintf(), perror(), stderr
#include <stdlib.h>	// atoi()
#include <string.h>	// memcpy(), memset()
#include <sys/socket.h>	// AF_INET, connet(), struct sockaddr
#include <unistd.h>	// close()

int open_tcp_service(const char * const service, const char * const address, struct sockaddr_in *host) {
	struct sockaddr_in name;
	if (!host) {
		host = &name;
	}
	// prepare the socket name for binding
	memset(host, 0, sizeof(struct sockaddr_in));
	host->sin_family = AF_INET;
	memcpy(&host->sin_addr.s_addr, address, 4);
	if (isdigit(*service)) {
		host->sin_port = htons(atoi(service));
	} else {
		struct servent * const server(getservbyname(service, "tcp"));
		if (!server) {
			fprintf(stderr, "getservbyname: failed\n");
			return -1;
		}
		host->sin_port = server->s_port;
	}
	const int connection(socket(PF_INET, SOCK_STREAM, 0));
	if (connection == -1) {
		perror("socket");
	} else if (connect(connection, (struct sockaddr *)host, sizeof(struct sockaddr_in)) == -1) {
		close(connection);
		perror("connect");
		return -1;
	}
	return connection;
}
