#include "get_name.h"	/* get_name() */
#include "open_compressed.h"	/* close_compressed(), open_compressed(), pfgets() */
#include "read.h"	/* Read, opt_quality_cutoff */
#include "read_lib.h"
#include <list>		/* list<> */
#include <map>		/* map<> */
#include <string>	/* string */
#include <sys/types.h>	/* size_t */

bool opt_strip_tracename = 0;
std::map<std::string, bool> opt_readname_match;

std::string make_read_name(std::string &header) {
	if (opt_strip_tracename) {
		size_t i = header.find(' ', 1);
		if (i != std::string::npos) {
			header = ">" + header.substr(i + 1);
		}
	}
	return get_name(header);
}

/*
 * make the quality file name from the sequence file name - mainly just
 * tacks on .qual, but also handles .gz and .Z endings (since the .qual
 * needs to be added before the .gz or .Z)
 */

std::string make_qual_filename(const char *filename, bool strip_fasta) {
	std::string qual_name = filename;
	size_t i;
	if (strip_fasta && (i = qual_name.rfind(".fasta")) != std::string::npos) {
		qual_name.replace(i, 6, ".qual");
	} else {
		i = qual_name.rfind('.');
		if (i != std::string::npos && (qual_name.substr(i) == ".bz2" || qual_name.substr(i) == ".gz" || qual_name.substr(i) == ".Z")) {
			qual_name.insert(i, ".qual");
		} else {
			qual_name += ".qual";
		}
	}
	return qual_name;
}

/*
 * goes through read list and masks any basepairs with a quality less than
 * the cutoff
 */

void mask_by_phred(std::list<Read> &read_list, unsigned int opt_phred_mask_cutoff) {
	std::list<Read>::iterator a = read_list.begin();
	std::list<Read>::iterator end = read_list.end();
	for (; a != end; ++a) {
		a->mask_by_phred(opt_phred_mask_cutoff);
	}
}

/*
 * return a pointer into read_list for the current read;
 * add header to the read if it doesn't already exist
 */

static int add_read(std::string &header, std::list<Read> &read_list, std::map<std::string, Read *> &read_lookup) {
	std::string name = make_read_name(header);
	/* read name has to match string, if present */
	if (!opt_readname_match.empty() && opt_readname_match.find(name) == opt_readname_match.end()) {
		read_lookup[name] = NULL;
		return 0;
	}
	if (read_lookup.find(name) != read_lookup.end()) {
		fprintf(stderr, "Warning: duplicate read sequence: %s\n", name.c_str());
		return 0;
	} else {
		read_list.push_back(Read(header));
		read_lookup[name] = &(*read_list.rbegin());
		return 1;
	}
}

/* add quality vector of equal length to sequence, of value x */

static void set_default_quals(std::list<Read> &read_list, unsigned char x) {
	std::list<Read>::iterator a = read_list.begin();
	std::list<Read>::iterator end = read_list.end();
	for (; a != end; ++a) {
		a->set_quality(x);
	}
}

/* read in contig sequence */

int read_sequence(const char *filename, std::list<Read> &read_list, bool opt_warnings) {
	int fd = open_compressed(filename);
	if (fd == -1) {
		return -1;
	}
	std::map<std::string, Read *> read_lookup;
	int adding_read = 0;
	std::string line, data;
	while (pfgets(fd, line) != -1) {
		if (line[0] == '>') {
			if (adding_read) {
				read_list.back().add_sequence(data);
				data.clear();
			}
			adding_read = add_read(line, read_list, read_lookup);
		} else if (adding_read) {
			data += line;
		}
	}
	close_compressed(fd);
	if (adding_read) {
		read_list.back().add_sequence(data);
	}
	fd = open_compressed(make_qual_filename(filename, 0));
	if (fd == -1) {
		fd = open_compressed(make_qual_filename(filename, 1));
	}
	if (fd == -1) {
		fprintf(stderr, "Warning: %s: qual file missing, defaulting qual's to %d\n", filename, opt_quality_cutoff);
		set_default_quals(read_list, opt_quality_cutoff);
		return 1;
	}
	static Read duplicate_read;
	Read * const duplicate_read_ptr = &duplicate_read;
	Read *current_read = NULL;
	while (pfgets(fd, line) != -1) {
		if (line[0] == '>') {
			if (current_read != NULL) {
				current_read->add_quality(data, opt_warnings);
			}
			data.clear();
			std::string name = make_read_name(line);
			std::map<std::string, Read *>::iterator a = read_lookup.find(name);
			if (a == read_lookup.end()) {
				fprintf(stderr, "Warning: no sequence for quality: %s\n", name.c_str());
				current_read = NULL;
			} else if (a->second == duplicate_read_ptr) {
				fprintf(stderr, "Warning: duplicate read quality: %s\n", name.c_str());
				current_read = NULL;
			} else {
				current_read = a->second;
				a->second = duplicate_read_ptr;
			}
		} else {
			data += line + ' ';
		}
	}
	close_compressed(fd);
	if (current_read != NULL) {
		current_read->add_quality(data, opt_warnings);
	}
	return 0;
}
