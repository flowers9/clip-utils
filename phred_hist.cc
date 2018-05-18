#include "open_compressed.h"	// close_compressed(), open_compressed(), pfgets()
#include "pretty_print.h"	// pretty_print()
#include <errno.h>	// errno
#include <exception>	// exception
#include <getopt.h>	// getopt(), optarg, optind
#include <iostream>	// cerr
#include <sstream>	// istringstream, ostringstream
#include <stdio.h>	// EOF, printf()
#include <string.h>	// strerror(), strlen()
#include <string>	// string
#include <utility>	// max(), min()
#include <vector>	// vector<>

static int opt_illumina, opt_bp_hist;
static size_t opt_failure_cutoff;

class LocalException : public std::exception {
    private:
	const std::string error;
    public:
	LocalException(const std::string &s) : error(s) { }
	virtual ~LocalException(void) throw() { }
	virtual const char * what() const throw() {
		return error.c_str();
	}
};

class Stats {
    private:
	size_t m_count;
	size_t m_hq_count;	// bases in current read with quality >= 20
	long m_last_qual;	// last qual read
    public:
	size_t good_bases;	// bases with quality >= 20
	size_t total_bases;
	size_t total_reads;
	size_t failed_reads;	// < opt_failure_cutoff bases of quality >= 20
	size_t failed_bases;	// bases with quality >= 20 for failures
	std::vector<size_t> hist;	// reads with given high quality count
	std::vector<size_t> bp_count;	// reads with this many+ bps
	std::vector<size_t> bp_qual;	// total quality at bp
	std::vector<size_t> bp_hq_hist;	// reads with quality >= 30 at bp
	Stats(void) : m_count(0), m_hq_count(0), m_last_qual(0), good_bases(0), total_bases(0), total_reads(0), failed_reads(0), failed_bases(0) { }
	~Stats(void) { }
	void add_fastq_line(const std::string &);
	void add_qual_line(const std::string &);
	void finish_read(void);
	void merge(const Stats &a) {
		good_bases += a.good_bases;
		total_bases += a.total_bases;
		total_reads += a.total_reads;
		failed_reads += a.failed_reads;
		failed_bases += a.failed_bases;
		size_t j;
		if (hist.size() < a.hist.size()) {
			hist.resize(a.hist.size(), 0);
		}
		for (j = 0; j != a.hist.size(); ++j) {
			hist[j] += a.hist[j];
		}
		if (bp_count.size() < a.bp_count.size()) {
			bp_count.resize(a.bp_count.size(), 0);
		}
		for (j = 0; j != a.bp_count.size(); ++j) {
			bp_count[j] += a.bp_count[j];
		}
		if (bp_qual.size() < a.bp_qual.size()) {
			bp_qual.resize(a.bp_qual.size(), 0);
		}
		for (j = 0; j != a.bp_qual.size(); ++j) {
			bp_qual[j] += a.bp_qual[j];
		}
		if (bp_hq_hist.size() < a.bp_hq_hist.size()) {
			bp_hq_hist.resize(a.bp_hq_hist.size(), 0);
		}
		for (j = 0; j != a.bp_hq_hist.size(); ++j) {
			bp_hq_hist[j] += a.bp_hq_hist[j];
		}
	}
	std::string::size_type size_width(void) const {
		std::ostringstream s;
		s << (hist.size() * (opt_illumina ? 10 : 50) - 1);
		return s.str().size();
	}
	std::string::size_type hist_width(void) const {
		size_t j(0);
		for (size_t i(0); i != hist.size(); ++i) {
			if (j < hist[i]) {
				j = hist[i];
			}
		}
		std::ostringstream s;
		s << j;
		return s.str().size();
	}
	std::string::size_type bp_avg_width(void) const {
		size_t j(0);
		for (size_t i(0); i != bp_count.size(); ++i) {
			const size_t k(bp_qual[i] / bp_count[i]);
			if (j < k) {
				j = k;
			}
		}
		std::ostringstream s;
		s << j;
		return s.str().size();
	}
	void trim(void) {		// trim trailing 0's off bp histograms
		if (bp_count.empty()) {
			return;
		}
		size_t i(bp_count.size() - 1);
		for (; i != static_cast<size_t>(-1) && bp_count[i] == 0; --i) { }
		if (++i != bp_count.size()) {
			bp_count.resize(i);
			bp_qual.resize(i);
			bp_hq_hist.resize(i);
		}
	}
};

void Stats::finish_read(void) {
	// ignore trailing 0 in illumina qualities
	if (m_count > 0 && opt_illumina && m_last_qual == 0) {
		--m_count;
	}
	if (m_count == 0) {
		return;
	}
	if (bp_count.size() < m_count) {
		bp_count.resize(m_count, 0);
	}
	for (size_t i(0); i < m_count; ++i) {
		++bp_count[i];
	}
	++total_reads;
	total_bases += m_count;
	good_bases += m_hq_count;
	if (m_hq_count < opt_failure_cutoff) {
		++failed_reads;
		failed_bases += m_hq_count;
	}
	m_hq_count /= opt_illumina ? 10 : 50;
	if (hist.size() < m_hq_count + 1) {
		hist.resize(m_hq_count + 1, 0);
	}
	++hist[m_hq_count];
	m_count = m_hq_count = 0;
	m_last_qual = 0;
}

void Stats::add_qual_line(const std::string &line) {
	size_t max_count(m_count);
	if (line.size() < 1048576) {	// approx size for small lines
		max_count += (line.size() + 1) / 2;
	} else {
		std::istringstream x(line);
		long i;
		for (x >> i; !x.fail(); ++max_count, x >> i) { }
	}
	if (bp_qual.size() < max_count) {
		bp_qual.resize(max_count, 0);
		bp_hq_hist.resize(max_count, 0);
	}
	std::istringstream x(line);
	long i;
	for (x >> i; !x.fail(); m_last_qual = i, ++m_count, x >> i) {
		if (i > 0 && i != 98) {
			bp_qual[m_count] += i;
			if (i >= 20) {
				++m_hq_count;
				if (i >= 30) {
					++bp_hq_hist[m_count];
				}
			}
		}
	}
}

void Stats::add_fastq_line(const std::string &line) {
	const size_t max_count(m_count + line.size());
	if (bp_qual.size() < max_count) {
		bp_qual.resize(max_count, 0);
		bp_hq_hist.resize(max_count, 0);
	}
	const size_t end_j(line.size());
	for (size_t j(0); j != end_j; ++j, ++m_count) {
		const int i(line[j] - 33);	// convert qual to int
		m_last_qual = i;
		if (i > 0 && i != 98) {
			bp_qual[m_count] += i;
			if (i >= 20) {
				++m_hq_count;
				if (i >= 30) {
					++bp_hq_hist[m_count];
				}
			}
		}
	}
}

static int count_file(const std::string &filename, Stats &stats) {
	const int fd(open_compressed(filename));
	if (fd == -1) {
		throw LocalException("open_compressed " + filename + ": " + strerror(errno));
	}
	std::string line;
	while (pfgets(fd, line) != -1) {
		if (line.empty()) {
		} else if (line[0] == '>') {
			stats.finish_read();
		} else if (line[0] == '@') {	// fastq
			stats.finish_read();
			pfgets(fd, line);	// seq
			pfgets(fd, line);	// +
			pfgets(fd, line);	// quality
			stats.add_fastq_line(line);
		} else {
			stats.add_qual_line(line);
		}
	}
	stats.finish_read();
	close_compressed(fd);
	stats.trim();
	return 1;
}

static void print_overall_stats(const Stats &stats) {
	printf("Number of reads: %s\n", pretty_print(stats.total_reads).c_str());
	printf("Total bases: %s\n", pretty_print(stats.total_bases).c_str());
	printf("Total Phred 20 bases: %s\n", pretty_print(stats.good_bases).c_str());
	printf("Average length: ");
	if (stats.total_reads == 0) {
		printf("0\n");
	} else {
		printf("%0.1f\n", static_cast<double>(stats.total_bases) / stats.total_reads);
	}
	printf("Phred average: ");
	if (stats.total_reads == 0) {
		printf("0\n");
	} else {
		printf("%0.1f\n", static_cast<double>(stats.good_bases) / stats.total_reads);
	}
	printf("Phred average without failures: ");
	if (stats.total_reads == stats.failed_reads) {
		printf("0\n");
	} else {
		printf("%0.1f\n", static_cast<double>(stats.good_bases - stats.failed_bases) / (stats.total_reads - stats.failed_reads));
	}
	printf("Percent failed: ");
	if (stats.total_reads == 0) {
		printf("0\n");
	} else {
		printf("%0.1f%%\n", static_cast<double>(100) * stats.failed_reads / stats.total_reads);
	}
	printf("\n");
}

static void print_hist(const Stats &stats) {
	const int size_width(std::max(3, static_cast<int>(stats.size_width())));
	const int hist_width(std::max(5, static_cast<int>(stats.hist_width())));
	printf("%-*s %*s %%ofReads\n%s %s --------\n", 2 * size_width + 1, "Phred20", hist_width, "Reads", std::string(2 * size_width + 1, '-').c_str(), std::string(hist_width, '-').c_str());
	if (stats.total_reads > 0) {
		const size_t n(opt_illumina ? 10 : 50);
		for (size_t i(0); i != stats.hist.size(); ++i) {
			const double x(static_cast<double>(100) * stats.hist[i] / stats.total_reads);
			printf("%*ld-%*ld %*ld %5.1f%%\t|", size_width, i * n, size_width, i * n + n - 1, hist_width, stats.hist[i], x);
			for (int j(1); j < x; j += 2) {
				printf("X");
			}
			printf("\n");
		}
	}
	printf("\n");
}

static void print_bp_hist(const Stats &stats) {
	std::ostringstream s;
	s << stats.bp_count.size();
	int count_width(s.str().size());
	int avg_width(stats.bp_avg_width());
	if (count_width < 2) {
		count_width = 2;
	}
	if (avg_width < 9) {
		avg_width = 9;
	}
	printf("%*s %*s %%Reads>=30\n%s %s ----------\n", count_width, "BP", avg_width, "Avg Score", std::string(count_width, '-').c_str(), std::string(avg_width, '-').c_str());
	size_t i;
	for (i = 0; i != stats.bp_count.size(); ++i) {
		printf("%*ld %*ld  %5.1f%%\n", count_width, i, avg_width, stats.bp_qual[i] / stats.bp_count[i], (double)100 * stats.bp_hq_hist[i] / stats.total_reads);
	}
	printf("\n");
}

static void print_stats(const std::string &name, const Stats &stats) {
	printf("%s\n%s\n\n", name.c_str(), std::string(name.size(), '=').c_str());
	print_hist(stats);
	if (opt_bp_hist) {
		print_bp_hist(stats);
	}
	print_overall_stats(stats);
}

static void print_usage() {
	std::cerr <<
		"usage: phred_hist [opts] [file1] [file2] ...\n"
		"    -b     print basepair histogram\n"
		"    -i     files are illumina sequences (strip trailing zero, size 10 bins)\n"
		"    -p ##  minimum number of phred 20s to pass [40]\n";
}

// set defaults and read in options
static int get_opts(int argc, char **argv) {
	opt_bp_hist = 0;
	opt_failure_cutoff = 40;
	opt_illumina = 0;
	int c;
	while ((c = getopt(argc, argv, "bip:")) != EOF) {
		switch (c) {
		    case 'b':
			opt_bp_hist = 1;
			break;
		    case 'i':
			opt_illumina = 1;
			break;
		    case 'p':
			std::istringstream(optarg) >> opt_failure_cutoff;
			break;
		    default:
			print_usage();
			return -1;
		}
	}
	if (optind == argc) {	// no file(s) specified
		print_usage();
		return -1;
	}
	return 0;
}

int main(int argc, char **argv) {
	if (get_opts(argc, argv) == -1) {
		return 1;
	}
	try {
		Stats overall_stats;
		int files(0);
		for (; optind < argc; ++optind) {
			std::string file(argv[optind]);
			Stats stats;
			count_file(file, stats);
			++files;
			print_stats(file, stats);
			overall_stats.merge(stats);
		}
		if (files > 1) {
			print_stats("Overall Totals", overall_stats);
		}
	} catch (std::exception &e) {
		std::cerr << "Error: " << e.what() << '\n';
		return 1;
	}
	return 0;
}
