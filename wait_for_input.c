#include <sys/time.h>

int main() {
	fd_set readfds;
	FD_SET(0, &readfds);
	select(1, &readfds, 0, 0, 0);
	return 0;
}
