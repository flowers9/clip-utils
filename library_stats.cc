#include "find_library.h"	/* find_library(), init_library_patterns() */
#include "library_read.h"	/* LibraryRead */
#include "library_read_lib.h"	/* library_read_sequence() */
#include "parse_read.h"	/* init_read_patterns(), parse_read_name() */
#include "pretty_print.h"	/* pretty_print() */
#include "read.h"	/* opt_clip_quality, opt_clip_vector */
#include <getopt.h>	// getopt(), optarg, optind
#include <list>		/* list<> */
#include <map>		/* map<> */
#include <stdio.h>	/* EOF, fprintf(), printf(), stderr */
#include <stdlib.h>	// atoi(), exit()
#include <string>	/* string */
#include <sys/types.h>	/* size_t */

static bool opt_warnings;

/* parse read names */

static void process_new_reads(std::list<LibraryRead> &read_list) {
	std::list<LibraryRead>::iterator a = read_list.begin();
	std::list<LibraryRead>::iterator end = read_list.end();
	for (; a != end; ++a) {
		parse_read_name(*a);
	}
}

/* count the total reads and phred 20s in a given fasta file */

static int process_fasta(const std::string &file, std::list<LibraryRead> &read_list) {
	if (library_read_sequence(file.c_str(), read_list, opt_warnings) == -1) {
		return 1;
	}
	process_new_reads(read_list);
	return 0;
}

/* match reads with their pairs */

static void pair_reads(std::list<LibraryRead> &read_list) {
	std::map<std::string, LibraryRead *> read_index;
	std::list<LibraryRead>::iterator a = read_list.begin();
	std::list<LibraryRead>::iterator end = read_list.end();
	/* create index name hash */
	for (; a != end; ++a) {
		std::string index = make_index_name(*a);
		if (!index.empty()) {
			read_index[index] = &(*a);
		}
	}
	/* match pairs, using hash */
	for (a = read_list.begin(); a != end; ++a) {
		if (a->pair == NULL) {
			std::string index = make_index_pair_name(*a);
			if (!index.empty()) {
				std::map<std::string, LibraryRead *>::iterator b = read_index.find(index);
				if (b != read_index.end() && b->second->pair == NULL) {
					a->pair = b->second;
					b->second->pair = &(*a);
				}
			}
		}
	}
}

class LibraryData {
    public:
	size_t reads;
	size_t good_reads;
	size_t good_pairs;
	size_t phred_count;
	LibraryData() : reads(0), good_reads(0), good_pairs(0), phred_count(0) { };
	~LibraryData() { };
};

static void collect_library_stats(const std::list<LibraryRead> &read_list) {
	std::map<std::string, LibraryData> library_list;
	std::list<LibraryRead>::const_iterator a = read_list.begin();
	std::list<LibraryRead>::const_iterator end_a = read_list.end();
	for (; a != end_a; ++a) {
		std::string library = find_library(*a);
		if (!library.empty()) {
			LibraryData *b = &library_list[library];
			++b->reads;
			if (a->phred_count >= 400) {
				++b->good_reads;
				b->phred_count += a->phred_count;
				if (a->pair != NULL) {
					++b->good_pairs;
				}
			}
		}
	}
	size_t total_reads = 0;
	size_t total_good_reads = 0;
	size_t total_good_pairs = 0;
	size_t total_phred_count = 0;
	printf("Library    Reads    Good Reads  Percent Good  Good Pairs   Phred20s\n");
	printf("-------  ---------  ----------  ------------  ----------  ----------\n");
	std::map<std::string, LibraryData>::const_iterator b = library_list.begin();
	std::map<std::string, LibraryData>::const_iterator end_b = library_list.end();
	for (; b != end_b; ++b) {
		total_reads += b->second.reads;
		total_good_reads += b->second.good_reads;
		total_good_pairs += b->second.good_pairs / 2;
		total_phred_count += b->second.phred_count;
		printf("%-7s  %9s  %10s     ", b->first.c_str(), pretty_print(b->second.reads).c_str(), pretty_print(b->second.good_reads).c_str());
		if (b->second.good_reads == 0) {
			printf(" -0-  ");
		} else {
			printf("%5.1f%%", (double)100 * b->second.good_reads / b->second.reads);
		}
		printf("     %10s  %10s\n", pretty_print(b->second.good_pairs / 2).c_str(), pretty_print(b->second.phred_count).c_str());
	}
	printf("-------  ---------  ----------  ------------  ----------  ----------\n");
	printf("Totals   %9s  %10s     ", pretty_print(total_reads).c_str(), pretty_print(total_good_reads).c_str());
	if (total_good_reads == 0) {
		printf(" -0-  ");
	} else {
		printf("%5.1f%%", (double)100 * total_good_reads / total_reads);
	}
	printf("     %10s  %10s\n", pretty_print(total_good_pairs).c_str(), pretty_print(total_phred_count).c_str());
}

/* print usage statement and exit */

static void print_usage() {
	fprintf(stderr, "usage: library_stats [options] file1 [file2] ...\n");
	fprintf(stderr, "    -c  do not clip low quality\n");
	fprintf(stderr, "    -q  turn off all warnings\n");
	fprintf(stderr, "    -v  do not clip vector\n");
	exit(1);
}

int main(int argc, char **argv) {
	opt_clip_quality = 1;
	opt_clip_vector = 1;
	opt_warnings = 1;
	int c;
	while ((c = getopt(argc, argv, "cqv")) != EOF) {
		switch (c) {
		    case 'c':
			opt_clip_quality = 0;
			break;
		    case 'q':
			opt_warnings = 0;
			break;
		    case 'v':
			opt_clip_vector = 0;
			break;
		    default:
			print_usage();
		}
	}
	if (optind == argc) {
		print_usage();
	}
	init_read_patterns();
	init_library_patterns();
	int err = 0;
	std::list<LibraryRead> read_list;
	/* count frequency of n-mers */
	for (; optind < argc; ++optind) {
		err += process_fasta(argv[optind], read_list);
	}
	pair_reads(read_list);
	collect_library_stats(read_list);
	return err;
}
