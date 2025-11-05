#include "hashl.h"	// hashl
#include "hashl_metadata.h"	// hashl_metadata
#include "open_compressed.h"	// close_compressed(), open_compressed()
#include <getopt.h>	// getopt(), optind
#include <iostream>	// cerr, cout
#include <sstream>	// istringstream
#include <stdlib.h>	// exit()
#include <string>	// string

static bool opt_just_metadata;
static bool opt_no_metadata;
static bool opt_print_sequence;
static int opt_debug_flags;

static void print_usage() {
	std::cerr << "usage: print_hashl <hashl_file>\n"
		"    -h    print this help\n"
		"    -d ## specify hash data fields to print (header = 1,\n"
		"              hash index = 2, data offset = 4, value = 8, key = 16) [31]\n"
		"    -M    don't print metadata\n"
		"    -m    only print metadata\n"
		"    -s    print stored sequence\n";
	exit(1);
}

static void get_opts(const int argc, char * const * const argv) {
	opt_debug_flags = 31;	// print everything
	opt_just_metadata = 0;
	opt_no_metadata = 0;
	opt_print_sequence = 0;
	int c;
	while ((c = getopt(argc, argv, "d:hMms")) != EOF) {
		switch (c) {
		    case 'd':
			std::istringstream(optarg) >> opt_debug_flags;
			break;
		    case 'h':
			print_usage();
			break;
		    case 'M':
			opt_no_metadata = 1;
			break;
		    case 'm':
			opt_just_metadata = 1;
			break;
		    case 's':
			opt_print_sequence = 1;
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
	hashl x;
	x.init_from_file(fd);
	close_compressed(fd);
	if (!opt_no_metadata) {
		hashl_metadata md;
		md.unpack(x.get_metadata());
		md.print();
	}
	if (opt_just_metadata) {
		return 0;
	}
	x.print(opt_debug_flags);
	if (opt_print_sequence) {
		x.print_sequence();
	}
	return 0;
}
