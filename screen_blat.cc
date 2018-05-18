#include "breakup_line.h"	/* breakup_line() */
#include <getopt.h>	// getopt(), optarg, strtod()
#include <iostream>	/* cin, cout, endl */
#include <map>		/* map<> */
#include <stdio.h>	/* EOF, fprintf(), stderr */
#include <stdlib.h>	// exit()
#include <streambuf>	/* getline() */
#include <string>	/* string */
#include <vector>	/* vector<> */

static bool opt_forward_only, opt_print_matched_read;
static double opt_read_identity;

static void print_usage() {
	fprintf(stderr, "usage: screen_blat [opts]\n");
	fprintf(stderr, "\t-I ##\tmatch identity [.98]\n");
	fprintf(stderr, "\t-m\tprint reads matched against\n");
	fprintf(stderr, "\t-r\tmatch reverse as well as forward directions\n");
	exit(1);
}

static void get_opts(const int argc, char * const * argv) {
	opt_forward_only = 1;
	opt_print_matched_read = 0;
	opt_read_identity = 0.98;
	int c;
	while ((c = getopt(argc, argv, "I:mr")) != EOF) {
		switch (c) {
		    case 'I':
			opt_read_identity = strtod(optarg, NULL);
			if (opt_read_identity < 0 || 1 < opt_read_identity) {
				fprintf(stderr, "Error: read identity is out of range [0,1]: %lf\n", opt_read_identity);
				print_usage();
			}
			break;
		    case 'm':
			opt_print_matched_read = 1;
			break;
		    case 'r':
			opt_forward_only = 0;
			break;
		    default:
			print_usage();
		}
	}
}

static void parse_output_matched(std::map<std::string, std::map<std::string, char> > &reads) {
	std::string line;
	getline(std::cin, line);	// skip header lines
	getline(std::cin, line);
	getline(std::cin, line);
	getline(std::cin, line);
	getline(std::cin, line);
	while (getline(std::cin, line)) {
		std::vector<std::string> fields;
		breakup_line(line, fields);
		if (fields.size() < 17) {
			fprintf(stderr, "Warning: short line: %lu\n", fields.size());
			continue;
		}
		// skip if reversed and forward_only is set
		bool is_forward = fields[8] == "+";
		if (opt_forward_only && !is_forward) {
			continue;
		}
		// skip if not very similar
		int query_size = atoi(fields[10].c_str());
		if (atoi(fields[0].c_str()) + atoi(fields[2].c_str()) < query_size * opt_read_identity) {
			continue;
		}
		// skip if query doesn't start and end within target
		if (is_forward && (atoi(fields[15].c_str()) < atoi(fields[11].c_str()) || atoi(fields[14].c_str()) - atoi(fields[16].c_str()) < query_size - atoi(fields[12].c_str()))) {
			continue;
		}
		if (!is_forward && (atoi(fields[15].c_str()) < query_size - atoi(fields[12].c_str()) || atoi(fields[14].c_str()) - atoi(fields[16].c_str()) < atoi(fields[11].c_str()))) {
			continue;
		}
		reads[fields[9]][fields[13]] = 0;
	}
}

static void print_output_matched(std::map<std::string, std::map<std::string, char> > &reads) {
	std::map<std::string, std::map<std::string, char> >::const_iterator a = reads.begin();
	std::map<std::string, std::map<std::string, char> >::const_iterator a_end = reads.end();
	for (; a != a_end; ++a) {
		std::cout << a->first;
		std::map<std::string, char>::const_iterator b = a->second.begin();
		std::map<std::string, char>::const_iterator b_end = a->second.end();
		for (; b != b_end; ++b) {
			std::cout << ' ' << b->first;
		}
		std::cout << std::endl;
	}
}

static void parse_output(std::map<std::string, char> &reads) {
	std::string line;
	getline(std::cin, line);	// skip header lines
	getline(std::cin, line);
	getline(std::cin, line);
	getline(std::cin, line);
	getline(std::cin, line);
	while (getline(std::cin, line)) {
		std::vector<std::string> fields;
		breakup_line(line, fields);
		if (fields.size() < 17) {
			fprintf(stderr, "Warning: short line: %lu\n", fields.size());
			continue;
		}
		// skip if reversed and forward_only is set
		bool is_forward = fields[8] == "+";
		if (opt_forward_only && !is_forward) {
			continue;
		}
		// skip if not very similar
		int query_size = atoi(fields[10].c_str());
		if (atoi(fields[0].c_str()) + atoi(fields[2].c_str()) < query_size * opt_read_identity) {
			continue;
		}
		// skip if query doesn't start and end within target
		if (is_forward && (atoi(fields[15].c_str()) < atoi(fields[11].c_str()) || atoi(fields[14].c_str()) - atoi(fields[16].c_str()) < query_size - atoi(fields[12].c_str()))) {
			continue;
		}
		if (!is_forward && (atoi(fields[15].c_str()) < query_size - atoi(fields[12].c_str()) || atoi(fields[14].c_str()) - atoi(fields[16].c_str()) < atoi(fields[11].c_str()))) {
			continue;
		}
		reads[fields[9]] = 0;
	}
}

static void print_output(std::map<std::string, char> &reads) {
	std::map<std::string, char>::const_iterator a = reads.begin();
	std::map<std::string, char>::const_iterator a_end = reads.end();
	for (; a != a_end; ++a) {
		std::cout << a->first << std::endl;
	}
}

int main(const int argc, char * const * argv) {
	get_opts(argc, argv);
	if (opt_print_matched_read) {
		std::map<std::string, std::map<std::string, char> > reads;
		parse_output_matched(reads);
		print_output_matched(reads);
	} else {
		std::map<std::string, char> reads;
		parse_output(reads);
		print_output(reads);
	}
	return 0;
}
