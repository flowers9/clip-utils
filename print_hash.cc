#include "hash.h"
#include "open_compressed.h"
#include <stdio.h>

int main(int argc, char **argv) {
	if (argc != 2) {
		fprintf(stderr, "usage: print_array <array_file>\n");
		return 1;
	}
	const int fd(open_compressed(argv[1]));
	hash x;
	x.init_from_file(fd);
	close_compressed(fd);
	printf("%lu %lu\n", x.size(),x.capacity());
	hash::const_iterator a(x.begin());
	const hash::const_iterator end_a(x.end());
	for (; a != end_a; ++a) {
		printf("%lu %lu\n", a.key, a.value);
	}
	return 0;
}
