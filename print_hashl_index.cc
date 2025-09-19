#include "hashl_index.h"	// hashl_index
#include "hashl_metadata.h"	// hashl_metadata
#include "open_compressed.h"	// close_compressed(), open_compressed()
#include <getopt.h>	// getopt(), optind
#include <iostream>	// cerr, cout
#include <stdlib.h>	// exit()
#include <string>	// string

static bool opt_just_metadata;

static void print_usage() {
	std::cerr << "usage: print_hashl_index <hashl_index_file>\n"
		"    -h  print this help\n"
		"    -m  only print metadata\n";
	exit(1);
}

static void get_opts(const int argc, char * const * const argv) {
	opt_just_metadata = 0;
	int c;
	while ((c = getopt(argc, argv, "hm")) != EOF) {
		switch (c) {
		    case 'h':
			print_usage();
			break;
		    case 'm':
			opt_just_metadata = 1;
			break;
		    default:
			std::cerr << "Error: unknown option " << char(c) << "\n";
			print_usage();
		}
	}
	if (optind + 1 != argc) {
		print_usage();
	}
}

int main(int argc, char **argv) {
	get_opts(argc, argv);
	const int fd(open_compressed(argv[optind]));
	if (fd == -1) {
			std::cerr << "Error: open: " << argv[optind] << '\n';
			return 1;
	}
	hashl_index x(fd);
	hashl_metadata md;
	md.unpack(x.get_metadata());
	md.print();
	if (opt_just_metadata) {
		close_compressed(fd);
		return 0;
	}
	x.print();
	close_compressed(fd);
	return 0;
}
