#include "hashn.h"	/* hashn */
#include "hist_lib_hashn.h"	/* init_mer_constants(), reverse_key() */
#include "open_compressed.h"	/* close_compressed(), open_compressed(), pfgets() */
#include "read_lib.h"	/* make_read_name(), opt_strip_tracename */
#include "version.h"	/* VERSION */
#include <errno.h>	/* errno */
#include <getopt.h>	// getopt(), optarg, optind
#include <map>		/* map<> */
#include <stdio.h>	/* EOF, FILE, fprintf(), stderr, stdout */
#include <stdlib.h>	// atoi(), atol(), exit()
#include <string.h>	/* strerror() */
#include <string>	/* string */
#include <sys/types.h>	/* size_t */

static size_t opt_mer_length;
static size_t opt_nmers;
static std::string opt_output;

/*
 * return the number represented by s, which may be suffixed by a k, m, or g
 * which act as multipliers to the base amount
 */

static size_t get_value(std::string s) {
	/* validate string - digits optionally followed by k, m, or g */
	size_t i = s.find_first_not_of("0123456789");
	if (i == std::string::npos) {
		return atol(s.c_str());
	} else if (i + 1 == s.length()) {
		size_t x = atol(s.c_str());
		switch (s[i]) {
		    case 'g':
			x *= 1024;
		    case 'm':
			x *= 1024;
		    case 'k':
			x *= 1024;
			return x;
		    default:
			return 0;
		}
	} else {	/* bad value */
		return 0;
	}
}

/* print usage statement and exit */

static void print_usage() {
	fprintf(stderr, "usage: screen_pairs [options] <fasta_file>\n");
	fprintf(stderr, "\t-h\tprint this information\n");
	fprintf(stderr, "\t-m ##\tbasepair match length [48]\n");
	fprintf(stderr, "\t-o ##\tprint output to file instead of stdout\n");
	fprintf(stderr, "\t-t\tstrip first part of trace id\n");
	fprintf(stderr, "\t-V\tprint version\n");
	fprintf(stderr, "\t-z ##\tnumber of possible n-mers to allocate memory for [200m]\n");
	fprintf(stderr, "\t\t(k, m, or g may be suffixed)\n");
	exit(1);
}

static void get_opts(int argc, char **argv) {
	opt_mer_length = 48;
	opt_nmers = 200 * 1024 * 1024;
	opt_strip_tracename = 0;
	int i, c;
	while ((c = getopt(argc, argv, "hm:o:tVz:")) != EOF) {
		switch (c) {
		    case 'h':
			print_usage();
			break;
		    case 'm':
			i = atoi(optarg);
			if (i < 1) {
				fprintf(stderr, "Error: invalid mer length: %d < 1\n", i);
				print_usage();
			}
			opt_mer_length = i;
			break;
		    case 'o':
			opt_output = optarg;
			break;
		    case 't':
			opt_strip_tracename = 1;
			break;
		    case 'V':
			fprintf(stderr, "histogram_hashn version %s%s\n", VERSION,
#ifdef COMPRESS_READS
" (read compression)"
#else
""
#endif
);
			exit(0);
		    case 'z':
			opt_nmers = get_value(optarg);
			if (opt_nmers == 0) {
				fprintf(stderr, "Error: bad n-mer count %s\n", optarg);
				print_usage();
			}
			break;
		    default:
			fprintf(stderr, "Error: unknown option %c\n", c);
			print_usage();
		}
	}
	if (optind == argc) {
		fprintf(stderr, "Error: no files to process\n");
		print_usage();
	}
	if (opt_output.empty()) {
		opt_output = argv[optind];
	}
}

static int add_to_key(const std::string &s, hashn::key_type &key) {
	size_t i;
	for (i = 0; i != opt_mer_length; ++i) {
		switch (s[i]) {
		    case 'A':
		    case 'a':
			key.push_back(0);
			break;
		    case 'C':
		    case 'c':
			key.push_back(1);
			break;
		    case 'G':
		    case 'g':
			key.push_back(2);
			break;
		    case 'T':
		    case 't':
			key.push_back(3);
			break;
		    default:
			return -1;
		}
	}
	return 0;
}

static bool add_read(std::string &name, const std::string &data, std::map<std::string, std::string> &read_lookup, hashn &mer_list) {
	size_t n = name.size();
	if (n < 3) {
		return 1;
	}
	--n;
	if ((name[n - 1] != '/' && (name[n - 1] != 'R' || name[n - 2] != '-')) || (name[n] != '1' && name[n] != '2')) {
		return 1;
	}
	std::map<std::string, std::string>::iterator a = read_lookup.find(name);
	// pair not found yet, so store until it is
	if (a == read_lookup.end()) {
		name[n] = name[n] == '1' ? '2' : '1';
		// save space if read is too short by not saving the sequence;
		// need to save name, or its pair will hang out forever
		read_lookup[name] = data.size() < opt_mer_length ? "" : data.substr(0, opt_mer_length);
		return 1;
	}
	// one or both reads insufficiently long, so skip
	if (data.size() < opt_mer_length || a->second.size() < opt_mer_length) {
		read_lookup.erase(a);
		return 1;
	}
	hashn::key_type key(mer_list);
	if (add_to_key(name[n] == '1' ? data : a->second, key) == -1) {
		read_lookup.erase(a);
		return 1;
	}
	if (add_to_key(name[n] == '2' ? data : a->second, key) == -1) {
		read_lookup.erase(a);
		return 1;
	}
	read_lookup.erase(a);
	hashn::key_type comp_key(mer_list);
	reverse_key(key, comp_key);
	return mer_list.increment(key < comp_key ? key : comp_key);
}

static int add_reads(const char *filename, hashn &mer_list) {
	int fd = open_compressed(filename);
	if (fd == -1) {
		return -1;
	}
	std::map<std::string, std::string> read_lookup;
	std::string line, name, data;
	while (pfgets(fd, line) != -1) {
		if (line[0] == '>') {
			if (!add_read(name, data, read_lookup, mer_list)) {
				fprintf(stderr, "Error: n-mer list incomplete - give a larger -z value\n");
				return -1;
			}
			name = make_read_name(line);
			data.clear();
		} else {
			data += line;
		}
	}
	close_compressed(fd);
	if (!add_read(name, data, read_lookup, mer_list)) {
		fprintf(stderr, "Error: n-mer list incomplete - give a larger -z value\n");
		return -1;
	}
	return 0;
}

static void print_read(std::string &name, const std::string &data, std::map<std::string, std::string> &read_lookup, hashn &mer_list, FILE *fp_unique, FILE *fp_dup, FILE *fp_bad) {
	size_t n = name.size();
	if (n == 0) {
		return;
	} else if (n < 3) {
		fprintf(fp_bad, "%s\n", name.c_str());
		return;
	}
	--n;
	if ((name[n - 1] != '/' && (name[n - 1] != 'R' || name[n - 2] != '-')) || (name[n] != '1' && name[n] != '2')) {
		fprintf(fp_bad, "%s\n", name.c_str());
		return;
	}
	name[n] = name[n] == '1' ? '2' : '1';
	std::map<std::string, std::string>::iterator a = read_lookup.find(name);
	name[n] = name[n] == '1' ? '2' : '1';
	// pair not found yet, so store until it is
	if (a == read_lookup.end()) {
		// save space if read is too short by not saving the sequence;
		// need to save name, or its pair will hang out forever
		read_lookup[name] = data.size() < opt_mer_length ? "" : data.substr(0, opt_mer_length);
		return;
	}
	// one or both reads insufficiently long, so skip
	if (data.size() < opt_mer_length || a->second.size() < opt_mer_length) {
		fprintf(fp_bad, "%s\n%s\n", a->first.c_str(), name.c_str());
		read_lookup.erase(a);
		return;
	}
	hashn::key_type key(mer_list);
	if (add_to_key(name[n] == '1' ? data : a->second, key) == -1) {
		fprintf(fp_bad, "%s\n%s\n", a->first.c_str(), name.c_str());
		read_lookup.erase(a);
		return;
	}
	if (add_to_key(name[n] == '2' ? data : a->second, key) == -1) {
		fprintf(fp_bad, "%s\n%s\n", a->first.c_str(), name.c_str());
		read_lookup.erase(a);
		return;
	}
	hashn::key_type comp_key(mer_list);
	reverse_key(key, comp_key);
	const hashn::key_type &k = key < comp_key ? key : comp_key;
	const hashn::value_type x = mer_list.value(k);
	if (x == 0) {
		fprintf(fp_dup, "%s\n%s\n", a->first.c_str(), name.c_str());
	} else {
		fprintf(fp_unique, "%s\n%s\n", a->first.c_str(), name.c_str());
		if (x != 1) {
			mer_list.assign(k, 0);
		}
	}
	read_lookup.erase(a);
}

// second pass of file - print unique reads
static int print_reads(const std::string &filename, hashn &mer_list) {
	int fd = open_compressed(filename);
	if (fd == -1) {
		return -1;
	}
	FILE *fp_unique = fopen((opt_output + ".unique").c_str(), "w");
	if (fp_unique == NULL) {
		fprintf(stderr, "Error: could not write to unique file: %s\n", strerror(errno));
		exit(1);
	}
	FILE *fp_dup = fopen((opt_output + ".dup").c_str(), "w");
	if (fp_dup == NULL) {
		fprintf(stderr, "Error: could not write to dup file: %s\n", strerror(errno));
		exit(1);
	}
	FILE *fp_bad = fopen((opt_output + ".bad").c_str(), "w");
	if (fp_bad == NULL) {
		fprintf(stderr, "Error: could not write to bad file: %s\n", strerror(errno));
		exit(1);
	}
	std::map<std::string, std::string> read_lookup;
	std::string line, name, data;
	while (pfgets(fd, line) != -1) {
		if (line[0] == '>') {
			print_read(name, data, read_lookup, mer_list, fp_unique, fp_dup, fp_bad);
			name = make_read_name(line);
			data.clear();
		} else {
			data += line;
		}
	}
	close_compressed(fd);
	print_read(name, data, read_lookup, mer_list, fp_unique, fp_dup, fp_bad);
	std::map<std::string, std::string>::const_iterator a = read_lookup.begin();
	std::map<std::string, std::string>::const_iterator end_a = read_lookup.end();
	for (; a != end_a; ++a) {
		fprintf(fp_bad, "%s\n", a->first.c_str());
	}
	fclose(fp_unique);
	fclose(fp_dup);
	fclose(fp_bad);
	return 0;
}

int main(int argc, char **argv) {
	get_opts(argc, argv);
	init_mer_constants(2 * opt_mer_length);
	hashn mer_list(opt_nmers, 4 * opt_mer_length);
	add_reads(argv[optind], mer_list);
	print_reads(argv[optind], mer_list);
	return 0;
}
