#include "pretty_print.h"	/* pretty_print() */
#include "qc_read.h"	/* QC_Read, opt_print_n_quality */
#include "qc_read_lib.h"	/* qc_calc_stats(), qc_read_sequence() */
#include "read.h"	/* opt_quality_cutoff */
#include <getopt.h>	// getopt(), optind
#include <list>		/* list<> */
#include <map>		/* map<> */
#include <stdio.h>	/* EOF, fprintf(), printf(), stderr */
#include <stdlib.h>	// exit()
#include <string>	/* string */
#include <sys/types.h>	/* size_t */

/* print histogram of N-run lengths */

static void print_n_histogram(const std::map<size_t, unsigned int> &n_hist) {
	printf("Histogram of N-run Lengths\n");
	printf("--------------------------\n");
	std::map<size_t, unsigned int>::const_iterator a = n_hist.begin();
	std::map<size_t, unsigned int>::const_iterator end = n_hist.end();
	for (; a != end; ++a) {
		printf("%lu	%s\n", a->first, pretty_print(a->second).c_str());
	}
}

/* print the number of low quality bases, binned by quality */

static void print_lq_histogram(std::map<size_t, unsigned int> &lq_hist) {
	printf("Histogram of Low Quality Values\n");
	printf("-------------------------------\n");
	size_t i;
	size_t end = lq_hist.rbegin()->first;
	for (i = 0; i <= end; ++i) {
		printf("%2lu-%2lu	%s\n", i * 5, i * 5 + 4, pretty_print(lq_hist[i]).c_str());
	}
}

/* print overall stats for fasta file */

static void print_overall_stats(const std::list<QC_Read> &read_list) {
	size_t total_size = 0;
	unsigned int total_scaffolds = 0;
	unsigned int total_contigs = 0;
	unsigned int total_n1_runs = 0;
	size_t total_n1_count = 0;
	size_t total_n2_count = 0;
	size_t total_lq_bases = 0;
	std::list<QC_Read>::const_iterator a = read_list.begin();
	std::list<QC_Read>::const_iterator end = read_list.end();
	for (; a != end; ++a) {
		/*
		 * screen small contigs from overall stats - need 8k non-gap bp
		 */
		if (a->size() - a->n1_count >= 8000) {
			++total_scaffolds;
			total_size += a->size();
			total_contigs += a->contigs;
			total_n1_runs += a->n1_runs;
			total_n1_count += a->n1_count;
			total_n2_count += a->n2_count;
			total_lq_bases += a->lq_count;
		}
	}
	printf("Sequence Size:              %s\n", pretty_print(total_size).c_str());
	printf("Scaffold Number:            %s\n", pretty_print(total_scaffolds).c_str());
	if (!read_list.empty()) {
		printf("Average Scaffold Size:      %s\n", pretty_print(total_size / total_scaffolds).c_str());
	} else {
		printf("Average Scaffold Size:     -0-\n");
	}
	printf("Contig Number:              %s\n", pretty_print(total_contigs).c_str());
	if (total_contigs != 0) {
		printf("Average Contig Size:        %s\n", pretty_print((total_size - total_n1_count) / total_contigs).c_str());
	} else {
		printf("Average Contig Size:       -0-\n");
	}
	printf("Contig Gap \"N\" Runs:        %s\n", pretty_print(total_n1_runs).c_str());
	if (total_n1_runs != 0) {
		printf("Average Gap \"N\" Run Size:   %s\n", pretty_print(total_n1_count / total_n1_runs).c_str());
	} else {
		printf("Average Gap \"N\" Run Size:  -0-\n");
	}
	printf("Gap \"N\" Bases Reported:     %s\n", pretty_print(total_n1_count).c_str());
	printf("Non-Gap \"N\" Bases Reported: %s\n", pretty_print(total_n2_count).c_str());
	printf("Actual Sequence Reported:   %s\n", pretty_print(total_size - total_n1_count - total_n2_count).c_str());
	printf("Jazz Low Quality Bases:     %s\n", pretty_print(total_lq_bases).c_str());
	if (total_size != total_n1_count) {
		printf("Percentage of Bases Marked as Low Quality: %3.2f%%\n", (double)100 * total_lq_bases / (total_size - total_n1_count));
	} else {
		printf("Percentage of Bases Marked as Low Quality: -0-\n");
	}
}

/* print stats for each scaffold in the fasta file */

static void print_scaffold_stats(std::list<QC_Read> &read_list) {
	printf(" Scaffold Name  Scaffold Size  Contigs  Gap Bases  LQ Bases  QC Percentage\n");
	printf("--------------  -------------  -------  ---------  --------  -------------\n");
	std::list<QC_Read>::const_iterator a = read_list.begin();
	std::list<QC_Read>::const_iterator end = read_list.end();
	for (; a != end; ++a) {
		printf("%-14s  %13lu  %7u  %9lu  %8lu     %6.2f%%\n", a->name().c_str(), a->size(), a->contigs, a->n1_count, a->lq_count, (double)100 * a->gc_count / (a->size() - a->n1_count - a->n2_count));
	}
}

/* print usage statement and exit */

static void print_usage() {
	fprintf(stderr, "usage: qc_stats1 [options] file1 [file2] ...\n");
	fprintf(stderr, "    -h  print histogram of N run lengths\n");
	fprintf(stderr, "    -n  print quality of N bases\n");
	fprintf(stderr, "    -q  turn off all warnings\n");
	exit(1);
}

int main(int argc, char **argv) {
	bool opt_histogram = 0;
	bool opt_warnings = 1;
	opt_print_n_quality = 0;
	opt_quality_cutoff = 40;
	int c;
	while ((c = getopt(argc, argv, "hnq")) != EOF) {
		switch (c) {
		    case 'h':
			opt_histogram = 1;
			break;
		    case 'n':
			opt_print_n_quality = 1;
			break;
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
	for (; optind < argc; ++optind) {
		std::list<QC_Read> read_list;
		if (qc_read_sequence(argv[optind], read_list, opt_warnings) == -1) {
			++err;
			continue;
		}
		std::map<size_t, unsigned int> n_hist;
		std::map<size_t, unsigned int> lq_hist;
		qc_calc_stats(read_list, n_hist, lq_hist);
		if (opt_print_n_quality) {
		} else if (opt_histogram) {
			print_n_histogram(n_hist);
		} else {
			print_overall_stats(read_list);
			printf("\n");
			print_lq_histogram(lq_hist);
			printf("\n");
			print_scaffold_stats(read_list);
		}
	}
	return err;
}
