#include "open_compressed.h"	/* close_compressed(), open_compressed(), pfgets() */
#include "pretty_print.h"	/* pretty_print() */
#include "read.h"	/* Read */
#include "read_lib.h"	/* read_sequence() */
#include <getopt.h>	// getopt(), optind
#include <list>		/* list<> */
#include <map>		/* map<> */
#include <stdio.h>	/* EOF, fprintf(), printf(), stderr */
#include <stdlib.h>	// exit()
#include <string.h>	/* strrchr() */
#include <string>	/* string */

static bool opt_warnings;

/* counts phreds over each read */

static unsigned long count_phreds(std::list<Read> &read_list) {
	unsigned long phred_count = 0;
	std::list<Read>::iterator a = read_list.begin();
	std::list<Read>::iterator end = read_list.end();
	for (; a != end; ++a) {
		phred_count += a->phred_count;
	}
	return phred_count;
}

/* count the total reads and phred 20s in a given fasta file */

static int process_fasta(const std::string &file, unsigned long *total_reads, unsigned long *total_phred_count) {
	std::list<Read> read_list;
	if (read_sequence(file.c_str(), read_list, opt_warnings) == -1) {
		return 1;
	}
	(*total_reads) += read_list.size();
	(*total_phred_count) += count_phreds(read_list);
	return 0;
}

/* go through the list of read sets and process the fasta for each one */

static int process_read_list(char *filename, unsigned long *total_reads, unsigned long *total_phred_count) {
	int fd = open_compressed(filename);
	if (fd == -1) {
		return 1;
	}
	std::list<std::string> read_file_list;
	std::string line;
	while (pfgets(fd, line) != -1) {
		read_file_list.push_back(line);
	}
	close_compressed(fd);
	std::list<std::string>::iterator a = read_file_list.begin();
	std::list<std::string>::iterator end = read_file_list.end();
	char *s = strrchr(filename, '/');
	std::string dir;
	if (s == NULL) {
		dir = "./";
	} else {
		s[1] = '\0';
		dir = filename;
	}
	dir += "traces/";
	int err = 0;
	for (; a != end; ++a) {
		fprintf(stderr, "%s\n", a->c_str());
		std::string file = dir + *a + "/edit_dir/" + *a + ".fasta.screen";
		err += process_fasta(file, total_reads, total_phred_count);
	}
	return err;
}

/* print usage statement and exit */

static void print_usage() {
	fprintf(stderr, "usage: qc_stats2 [-q] file1 [file2] ...\n");
	exit(1);
}

int main(int argc, char **argv) {
	opt_warnings = 1;
	int c;
	while ((c = getopt(argc, argv, "q")) != EOF) {
		switch (c) {
		    case 'q':
			opt_warnings = 0;
			break;
		    default:
			print_usage();
		}
	}
	if (optind == argc) {
		print_usage();
	}
	int err = 0;
	unsigned long total_reads = 0;
	unsigned long total_phred_count = 0;
	for (; optind < argc; ++optind) {
		err += process_read_list(argv[optind], &total_reads, &total_phred_count);
	}
	printf("Initial # of Reads:     %s\n", pretty_print(total_reads).c_str());
	printf("Initial # of Phred 20s: %s\n", pretty_print(total_phred_count).c_str());
	return err;
}
