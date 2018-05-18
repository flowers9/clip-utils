// this program first screens out all reads that are included in other reads,
// then trims the remaining reads if they have ends that match in the middle
// of other reads

#include "open_compressed.h"	// close_compressed(), open_compressed(), pfgets()
#include <getopt.h>	// getopt(), optarg, optind
#include <iostream>	// cerr, cout
#include <list>		// list<>
#include <map>		// map<>
#include <set>		// set<>
#include <sstream>	// istringstream
#include <stdlib.h>	// exit()
#include <string>	// string
#include <utility>	// pair<>
#include <vector>	// vector<>

// list of valid reads and any needed trimming (left/right)
static std::map<std::string, std::pair<int, int> > trims;
static std::set<std::string> inside_list;
static double opt_identity;
static int opt_read_offset, opt_trim_match, opt_trim_offset, opt_min_length;

static void print_usage(void) {
	std::cerr <<
		"usage: screen_reads [opts] <blat_file>\n"
		"\t-C ##\tlength of offset from end to be consider in the middle [50]\n"
		"\t-I ##\tmatch identity [.98]\n"
		"\t-i ##\tfile with list of inside reads\n"
		"\t-L ##\tfile with list of extra blat files\n"
		"\t-l ##\tminimum length of match to be used for trimming [50]\n"
		"\t-m ##\tminimum length of trimmed read to count as good [2000]\n"
		"\t-O ##\tmaximum match offset from read end [2]\n";
	exit(0);
}

static void read_blat_files(std::list<std::string> &blat_files, const char *blat_file_list) {
	const int fd = open_compressed(blat_file_list);
	if (fd == -1) {
		std::cerr << "Error reading blat file list: open: " << blat_file_list << "\n";
		exit(1);
	}
	std::string line;
	while (pfgets(fd, line) != -1) {
		blat_files.push_back(line);
	}
	close_compressed(fd);
}

static void read_inside_list(const std::string &file) {
	const int fd = open_compressed(file);
	if (fd == -1) {
		std::cerr << "Error reading inside list: open: " << file << "\n";
		exit(1);
	}
	std::string line;
	while (pfgets(fd, line) != -1) {
		inside_list.insert(line);
	}
	close_compressed(fd);
}

static int get_opts(int argc, char **argv, std::list<std::string> &blat_files) {
	opt_identity = .98;
	opt_min_length = 2000;
	opt_read_offset = 2;
	opt_trim_match = 50;
	opt_trim_offset = 50;
	std::string inside_file;
	std::istringstream x;
	int c;
	while ((c = getopt(argc, argv, "C:I:i:L:l:m:O:")) != -1) {
		switch(c) {
		    case 'C':
			x.str(optarg);
			x >> opt_trim_match;
			if (opt_trim_match < 0) {
				std::cerr << "Error: trim match is negative: " << opt_trim_match << "\n";
				print_usage();
			}
			break;
		    case 'I':
			x.str(optarg);
			x >> opt_identity;
			if (opt_identity < 0 || 1 < opt_identity) {
				std::cerr << "Error: match identity is out of range [0,1]: " << opt_identity << "\n";
				print_usage();
			}
			break;
		    case 'i':
			read_inside_list(optarg);
			break;
		    case 'L':
			read_blat_files(blat_files, optarg);
			break;
		    case 'l':
			x.str(optarg);
			x >> opt_trim_offset;
			if (opt_trim_offset < 0) {
				std::cerr << "Error: trim offset is negative: " << opt_trim_offset << "\n";
				print_usage();
			}
			break;
		    case 'm':
			x.str(optarg);
			x >> opt_min_length;
			if (opt_min_length < 0) {
				std::cerr << "Error: minimum read length is negative: " << opt_min_length << "\n";
				print_usage();
			}
			break;
		    case 'O':
			x.str(optarg);
			x >> opt_read_offset;
			if (opt_read_offset < 0) {
				std::cerr << "Error: match offset is negative: " << opt_read_offset << "\n";
				print_usage();
			}
			break;
		    default:
			print_usage();
		}
	}
	if (optind == argc && blat_files.empty()) {
		std::cerr << "Error: no blat files given\n";
		print_usage();
	}
	opt_min_length -= 2;	// to handle double offset from trimming output
	return inside_list.empty() ? 1 : 2;
}

// make list of what reads count as inside (which will then be ignored)

static void get_insides(const std::string &blat_file) {
	const int fd = open_compressed(blat_file);
	if (fd == -1) {
		std::cerr << "Error reading dup list: open: " << blat_file << "\n";
		exit(1);
	}
	std::string line;
	pfgets(fd, line);	// skip header lines
	pfgets(fd, line);
	pfgets(fd, line);
	pfgets(fd, line);
	pfgets(fd, line);
	while (pfgets(fd, line) != -1) {
		int match, mismatch, repmatch, nmatch, q_gap_count, q_gap_bases;
		int t_gap_count, t_gap_bases, q_size, q_start, q_end, t_size;
		int t_start, t_end, block_count;
		std::string strand, q_name, t_name, block_sizes, q_starts;
		std::string t_starts;
		std::istringstream x(line);
		x >> match;
		x >> mismatch;
		x >> repmatch;
		x >> nmatch;
		x >> q_gap_count;
		x >> q_gap_bases;
		x >> t_gap_count;
		x >> t_gap_bases;
		x >> strand;
		x >> q_name;
		x >> q_size;
		x >> q_start;
		x >> q_end;
		x >> t_name;
		x >> t_size;
		x >> t_start;
		x >> t_end;
		x >> block_count;
		x >> block_sizes;
		x >> q_starts;
		x >> t_starts;
		if (x.fail()) {
			std::cerr << "Warning: bad line in " << blat_file << ": " << line << "\n";
			continue;
		}
		// skip if query name == target name
		if (q_name == t_name) {
			continue;
		}
		// this is the longest possible match between these two
		const int match_length = std::min(q_size, t_size);
		const int identity = match + repmatch;
		// skip if not similar
		if (identity < opt_identity * match_length) {
			continue;
		}
		int start_offset, end_offset;
		if (strand[0] == '+') {			// forward
			start_offset = q_start - t_start;
			end_offset = q_size - q_end - (t_size - t_end);
		} else {			// reverse
			start_offset = q_start - (t_size - t_end);
			end_offset = q_size - q_end - t_start;
		}
		const bool q_inside = start_offset <= opt_read_offset && end_offset <= opt_read_offset;
		const bool t_inside = -start_offset <= opt_read_offset && -end_offset <= opt_read_offset;
		// see if one counts as inside the other
		if (q_inside && (!t_inside || q_size < t_size)) {
			inside_list.insert(q_name);
		} else if (t_inside) {
			inside_list.insert(t_name);
		}
	}
	close_compressed(fd);
}

static void get_trims(const std::string &blat_file) {
	const int fd = open_compressed(blat_file);
	if (fd == -1) {
		std::cerr << "Error reading blat file: open: " << blat_file << "\n";
		exit(1);
	}
	std::string line;
	pfgets(fd, line);	// skip header lines
	pfgets(fd, line);
	pfgets(fd, line);
	pfgets(fd, line);
	pfgets(fd, line);
	while (pfgets(fd, line) != -1) {
		int match, mismatch, repmatch, nmatch, q_gap_count, q_gap_bases;
		int t_gap_count, t_gap_bases, q_size, q_start, q_end, t_size;
		int t_start, t_end, block_count;
		std::string strand, q_name, t_name, block_sizes, q_starts;
		std::string t_starts;
		std::istringstream x(line);
		x >> match;
		x >> mismatch;
		x >> repmatch;
		x >> nmatch;
		x >> q_gap_count;
		x >> q_gap_bases;
		x >> t_gap_count;
		x >> t_gap_bases;
		x >> strand;
		x >> q_name;
		x >> q_size;
		x >> q_start;
		x >> q_end;
		x >> t_name;
		x >> t_size;
		x >> t_start;
		x >> t_end;
		x >> block_count;
		x >> block_sizes;
		x >> q_starts;
		x >> t_starts;
		if (x.fail()) {
			std::cerr << "Warning: bad line in " << blat_file << ": " << line << "\n";
			continue;
		}
		// skip if query name == target name
		if (q_name == t_name) {
			continue;
		}
		const int length = match + repmatch;
		// skip if not long enough
		if (length < opt_trim_match) {
			continue;
		}
		const int length2 = length + mismatch + nmatch;
		// skip if gaps are too large
		if (length < (length2 + q_gap_bases) * opt_identity || length < (length2 + t_gap_bases) * opt_identity) {
			continue;
		}
		const int q_stop = q_size - q_end;
		const int t_stop = t_size - t_end;
		const bool q_at_end = q_start <= opt_read_offset || q_stop <= opt_read_offset;
		const bool t_at_end = t_start <= opt_read_offset || t_stop <= opt_read_offset;
		const bool q_in_middle = q_start >= opt_trim_offset && q_stop >= opt_trim_offset;
		const bool t_in_middle = t_start >= opt_trim_offset && t_stop >= opt_trim_offset;
		// figure out what (if any) trimming is needed
		if (q_at_end && t_in_middle) {
			std::pair<int, int> &a = trims[q_name];
			if (a.first == 0) {	// initialize, if needed
				a.first = 1;
				a.second = q_size - 1;
			}
			if (q_start <= opt_read_offset) {	// front end
				if (a.first < q_end) {
					a.first = q_end;
				}
			}
			if (q_stop <= opt_read_offset) {	// back end
				if (a.second > q_start) {
					a.second = q_start;
				}
			}
		}
		if (t_at_end && q_in_middle) {
			std::pair<int, int> &a = trims[t_name];
			if (a.first == 0) {	// initialize, if needed
				a.first = 1;
				a.second = t_size - 1;
			}
			if (t_start <= opt_read_offset) {	// front end
				if (a.first < t_end) {
					a.first = t_end;
				}
			}
			if (t_stop <= opt_read_offset) {	// back end
				if (a.second > t_start) {
					a.second = t_start;
				}
			}
		}
	}
	close_compressed(fd);
}

int main(int argc, char **argv) {
	std::list<std::string> blat_files;
	const int pass = get_opts(argc, argv, blat_files);
	for (int n = optind; n != argc; ++n) {
		blat_files.push_back(argv[n]);
	}
	std::list<std::string>::const_iterator a = blat_files.begin();
	const std::list<std::string>::const_iterator end_a = blat_files.end();
	if (pass == 1) {
		for (; a != end_a; ++a) {
			get_insides(*a);
		}
		std::set<std::string>::const_iterator b = inside_list.begin();
		const std::set<std::string>::const_iterator end_b = inside_list.end();
		for (; b != end_b; ++b) {
			std::cout << *b << "\n";
		}
	} else {
		for (; a != end_a; ++a) {
			get_trims(*a);
		}
		std::map<std::string, std::pair<int, int> >::const_iterator b = trims.begin();
		const std::map<std::string, std::pair<int, int> >::const_iterator end_b = trims.end();
		for (; b != end_b; ++b) {
			std::cout << b->first;
			if (b->second.second - b->second.first >= opt_min_length) {
				std::cout << " " << (b->second.first - 1) << " " << (b->second.second + 1);
			}
			std::cout << "\n";
		}
	}
	return 0;
}
