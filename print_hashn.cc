#include "hashn.h"
#include "hist_lib_hashn.h"
#include "open_compressed.h"
#include <stdio.h>

int main(int argc, char **argv) {
	if (argc != 2) {
		fprintf(stderr, "usage: print_array <array_file>\n");
		return 1;
	}
	const int fd(open_compressed(argv[1]));
	hashn x;
	x.init_from_file(fd);
	close_compressed(fd);
	init_mer_constants(x.bits() / 2);
	printf("%lu %lu\n", x.size(),x.capacity());
	hashn::const_iterator a(x.begin());
	const hashn::const_iterator end_a(x.end());
	for (; a != end_a; ++a) {
		printf("%s %lu\n", convert_key(a.key).c_str(), a.value);
	}
	return 0;
}
