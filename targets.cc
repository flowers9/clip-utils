#include "itoa.h"	/* itoa() */
#include "range.h"	/* Range */
#include "read.h"	/* Read, opt_quality_cutoff */
#include "read_lib.h"	/* opt_readname_match, read_sequence() */
#include <getopt.h>	// getopt(), optarg, optind
#include <list>		/* list<> */
#include <map>		/* map<> */
#include <stdio.h>	/* EOF, FILE, fclose(), fopen(), fprintf(), stderr */
#include <stdlib.h>	// atoi(), exit()
#include <string>	/* string */
#include <sys/types.h>	/* size_t */
#include <unistd.h>	/* unlink() */

#define LQ_WINDOW 500

static bool opt_full_contigs;
static bool opt_qual_warning;
static size_t opt_min_lq_run;
static size_t opt_qual_length_cutoff;

/*
 * find all N's with quality 1 - these are gaps; if writing full contigs,
 * change non-gap N's to A's with quality 1
 */

static void find_gaps(Read *a, std::list<Range> &gaps) {
	size_t i;
	for (i = 0; i != a->size() && a->get_sequence(i) != 'N'; ++i) { }
	while (i != a->size()) {
		if (a->get_quality(i) == 1) {
			size_t start = i;
			for (; i != a->size() && a->get_sequence(i) == 'N'; ++i) { }
			/* back over trailing non-gap N's */
			for (--i; a->get_quality(i) != 1; --i) { }
			Range b(start, i);
			gaps.push_back(b);
		} else if (opt_full_contigs) {	/* erase non-gap N's */
			a->set_sequence(i, 'A');
			a->set_quality(i, 1);
		} else if (a->get_quality(i) >= opt_quality_cutoff) {
			fprintf(stderr, "Warning: high quality N: %s: %lu\n", a->name().c_str(), i + 1);
		}
		for (++i; i != a->size() && a->get_sequence(i) != 'N'; ++i) { }
	}
}

/*
 * make a list of low quality runs (not including gaps); if printing full
 * contigs, instead change all low qualities to 1, and return an empty list
 */

static void find_lq_runs(Read *a, const std::list<Range> &gaps, std::list<Range> &lq_runs) {
	std::list<Range>::const_iterator b = gaps.begin();
	std::list<Range>::const_iterator end_b = gaps.end();
	size_t i = a->quality_start;
	if (b != end_b && b->start == i) {
		i = b->stop + 1;
		++b;
	}
	for (;; ++b) {
		size_t end_i;
		if (b == end_b) {
			end_i = a->quality_stop;
		} else {
			end_i = b->start;
		}
		for (; i != end_i; ++i) {
			if (a->get_quality(i) >= opt_quality_cutoff) {
			} else if (opt_full_contigs) {
				a->set_quality(i, 1);
				/*
				 * avoid the check on opt_full_contigs for each
				 * position, since these tend to be in runs
				 */
				for (++i; i != end_i && a->get_quality(i) < opt_quality_cutoff; ++i) {
					a->set_quality(i, 1);
				}
				--i;
			} else {
				size_t start = i;
				for (++i; i != end_i && a->get_quality(i) < opt_quality_cutoff; ++i) { }
				--i;
				Range c(start, i);
				lq_runs.push_back(c);
			}
		}
		if (b == end_b) {
			break;
		}
		i = b->stop + 1;
	}
}

/*
 * make a list of targets based on gaps (Ns with a quality of 1) and
 * low quality runs
 *
 * The jgi program (which, among other things, made targets from jazz
 * files) grouped gaps (Ns for us) and low quality runs into windows of
 * no more than LQ_WINDOW basepairs (although possibly longer if a single
 * run or gap was larger than that).  If there were fewer than LQ_WINDOW
 * basepairs between the beginning (or end) of the contig and the adjacent
 * gap or run, the run would be extended to to edge of the contig.
 */

static void make_targets(Read *a, std::list<Range> &targets) {
	std::list<Range> gaps;
	find_gaps(a, gaps);
	std::list<Range> lq_runs;
	find_lq_runs(a, gaps, lq_runs);
	std::list<Range>::const_iterator b = gaps.begin();
	std::list<Range>::const_iterator end_b = gaps.end();
	std::list<Range>::const_iterator c = lq_runs.begin();
	std::list<Range>::const_iterator end_c = lq_runs.end();
	size_t start = a->quality_start;
	for (;; ++b) {
		size_t stop;
		if (b == end_b) {
			stop = a->quality_stop;
		} else {
			stop = b->start;
		}
		bool first = 1;
		for (; c != end_c && c->start < stop; ++c) {
			std::list<Range>::const_iterator d = c;
			for (++c; c != end_c && c->start < stop && c->stop - d->start + 1 <= LQ_WINDOW; ++c) { }
			--c;
			/*
			 * minimum size for low quality runs; don't include
			 * low quality runs that cover the entire contig
			 */
			if (c->stop - d->start + 1 >= opt_min_lq_run && (d->start > start || c->stop + 1 < stop)) {
				/*
				 * don't merge with gap - either not gap,
				 * or too far from it
				 */
				if (!first || d->start >= start + LQ_WINDOW) {
					Range x(d->start, c->stop);
					targets.push_back(x);
				/* extend to read beginning (effective gap) */
				} else if (targets.empty()) {
					Range x(start, c->stop);
					targets.push_back(x);
				} else {
					targets.back().stop = c->stop;
				}
				first = 0;
			}
		}
		if (b == end_b) {
			break;
		}
		/* merge next gap with last lq run, if close enough (if any) */
		if (!first && targets.back().stop + LQ_WINDOW >= b->start) {
			targets.back().stop = b->stop;
		} else {
			targets.push_back(*b);
		}
		start = b->stop + 1;
	}
}

/*
 * break up the existing reads by Ns and low quality runs to create new reads
 * named by the original readname and position; the order of the new reads is
 * the same as the old ones, just with more subsections
 */

static void breakup(std::list<Read> &read_list, std::list<Read> &target_read_list) {
	std::list<Read>::iterator a = read_list.begin();
	std::list<Read>::iterator end_a = read_list.end();
	for (; a != end_a; ++a) {
		std::list<Range> targets;
		make_targets(&(*a), targets);
		std::list<Range>::const_iterator c = targets.begin();
		std::list<Range>::const_iterator end_c = targets.end();
		/* start at the beginning, unless the first target is there */
		size_t start = 0;
		if (c != end_c && c->start == start) {
			start = c->stop + 1;
			++c;
		}
		for (;; ++c) {
			size_t stop;
			if (c == end_c) {
				stop = a->size();
			} else {
				stop = c->start;
			}
			target_read_list.push_back(a->subseq(start, stop));
			/*
			 * if at the end, or the last target reaches
			 * the end, then we're finished
			 */
			if (c == end_c || c->stop + 1 == a->size()) {
				break;
			}
			start = c->stop + 1;
		}
	}
}

/*
 * make an output filename from the original filename - strip leading
 * directories, strip a trailing .gz, .bz2, or .Z (if any), and add a .target
 */

static std::string make_filename(const std::string &file) {
	std::string s = file;
	size_t i = s.rfind('/');
	if (i != std::string::npos) {
		s.erase(0, i + 1);
	}
	if (s.length() > 2 && s.substr(s.length() - 2) == ".Z") {
		s.erase(s.length() - 2, s.length());
	} else if (s.length() > 3 && s.substr(s.length() - 3) == ".gz") {
		s.erase(s.length() - 3, s.length());
	} else if (s.length() > 4 && s.substr(s.length() - 4) == ".bz2") {
		s.erase(s.length() - 4, s.length());
	}
	s += ".target";
	return s;
}

/*
 * print reads in the same order they appeared in the quality file,
 * clipping as specified; print to individual files, named from read name
 */

static void print_targets(std::list<Read> &read_list) {
	std::string filename;
	FILE *fp_seq = NULL;
	FILE *fp_qual = NULL;
	std::list<Read>::const_iterator a = read_list.begin();
	std::list<Read>::const_iterator end = read_list.end();
	for (; a != end; ++a) {
		if (a->quality_stop - a->quality_start < opt_qual_length_cutoff) {
			if (opt_qual_warning) {
               	         	fprintf(stderr, "Warning: quality sequence too short, skipping %s\n", a->name().c_str());
			}
			continue;
		}
		/* strip trailing _xxx from name to get base read name */
		std::string name = a->name();
		size_t i = name.rfind('_');
		if (i != std::string::npos) {
			name.erase(i);
		}
		if (filename != name) {
			filename = name;
			if (fp_seq != NULL && fp_qual != NULL) {
				fclose(fp_seq);
				fclose(fp_qual);
			}
			fp_seq = fopen(filename.c_str(), "w");
			if (fp_seq == NULL) {
				fprintf(stderr, "Error: could not write to %s\n", filename.c_str());
				continue;
			}
			fp_qual = fopen((filename + ".qual").c_str(), "w");
			if (fp_qual == NULL) {
				fprintf(stderr, "Error: could not write to %s.qual\n", filename.c_str());
				fclose(fp_seq);
				unlink(filename.c_str());
				continue;
			}
		}
		a->print_sequence(fp_seq);
		a->print_quality(fp_qual, 96);
	}
	if (fp_seq != NULL && fp_qual != NULL) {
		fclose(fp_seq);
		fclose(fp_qual);
	}
}

/*
 * print reads in the same order they appeared in the quality file,
 * clipping as specified; print to one file, made from given given name
 */

static void print_targets(std::string file, std::list<Read> &read_list) {
	std::string filename = make_filename(file);
	FILE *fp_seq = fopen(filename.c_str(), "w");
	if (fp_seq == NULL) {
		fprintf(stderr, "Error: could not write to %s\n", filename.c_str());
		return;
	}
	filename += ".qual";
	FILE *fp_qual = fopen(filename.c_str(), "w");
	if (fp_qual == NULL) {
		fprintf(stderr, "Error: could not write to %s\n", filename.c_str());
		fclose(fp_seq);
		filename.erase(filename.length() - 5, filename.length());
		unlink(filename.c_str());
		return;
	}
	std::list<Read>::const_iterator a = read_list.begin();
	std::list<Read>::const_iterator end = read_list.end();
	for (; a != end; ++a) {
		if (a->quality_stop - a->quality_start < opt_qual_length_cutoff) {
			if (opt_qual_warning) {
               	         	fprintf(stderr, "Warning: quality sequence too short, skipping %s\n", a->name().c_str());
			}
		} else {
			a->print_sequence(fp_seq);
			a->print_quality(fp_qual, 96);
		}
	}
	fclose(fp_seq);
	fclose(fp_qual);
}

/* print usage statement and exit */

static void print_usage() {
	fprintf(stderr, "usage: targets [options] file1 [file2] ...\n");
	fprintf(stderr, "    -c ## delete sequences with less than ## basepairs\n");
	fprintf(stderr, "    -e    extract sequence without creating targets\n");
	fprintf(stderr, "    -f    include all non-gap bases, but set N's to A's and low quality to\n");
	fprintf(stderr, "          quality 1\n");
	fprintf(stderr, "    -m ## minimum length of a low quality run for targets\n");
	fprintf(stderr, "    -p    split output by read (separate file for each)\n");
	fprintf(stderr, "    -q    turn off all warnings\n");
	fprintf(stderr, "    -s XX only process read matching given string (may\n");
	fprintf(stderr, "          be specified multiple times\n");
	fprintf(stderr, "    -w    turn on short quality sequence warning\n");
	exit(1);
}

int main(int argc, char **argv) {
	/* set option defaults */
	bool opt_extract = 0;
	bool opt_split = 0;
	bool opt_warnings = 1;
	opt_full_contigs = 0;
	opt_min_lq_run = 1;
	opt_qual_length_cutoff = 0;
	opt_qual_warning = 0;
	opt_quality_cutoff = 30;
	/* read in options */
	int c;
	while ((c = getopt(argc, argv, "c:efm:pqs:w")) != EOF) {
		switch (c) {
		    case 'c':
			c = atoi(optarg);
			if (c < 0) {
				print_usage();
			}
			opt_qual_length_cutoff = c;
			break;
		    case 'e':
			opt_extract = 1;
			break;
		    case 'f':
			opt_full_contigs = 1;
			break;
		    case 'm':
			c = atoi(optarg);
			if (c < 0) {
				print_usage();
			}
			opt_min_lq_run = c;
			break;
		    case 'p':
			opt_split = 1;
			break;
		    case 'q':
			opt_warnings = 0;
			break;
		    case 's':
			opt_readname_match[optarg] = 1;
			break;
		    case 'w':
			opt_qual_warning = 1;
			break;
		    default:
			print_usage();
		}
	}
	if (optind == argc) {	/* oops, no files specified */
		print_usage();
	}
	int err = 0;
	for (; optind < argc; ++optind) {
		std::list<Read> read_list;
		if (read_sequence(argv[optind], read_list, opt_warnings) == -1) {
			++err;
			continue;
		}
		if (opt_extract) {
			if (opt_split) {
				print_targets(read_list);
			} else {
				print_targets(argv[optind], read_list);
			}
		} else {
			std::list<Read> target_read_list;
			breakup(read_list, target_read_list);
			if (opt_split) {
				print_targets(target_read_list);
			} else {
				print_targets(argv[optind], target_read_list);
			}
		}
	}
	return err;
}
