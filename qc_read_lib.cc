#include "get_name.h"	/* get_name() */
#include "open_compressed.h"	/* close_compressed(), open_compressed(), pfgets() */
#include "qc_read.h"	/* QC_Read */
#include "qc_read_lib.h"
#include "read_lib.h"	/* make_qual_filename(), make_read_name(), opt_readname_match, opt_strip_tracename */
#include <list>		/* list<> */
#include <map>		/* map<> */
#include <string>	/* string */
#include <sys/types.h>	/* size_t */

/*
 * return a pointer into read_list for the current read;
 * add header to the read if it doesn't already exist
 */

static bool add_read(std::string &header, std::list<QC_Read> &read_list, std::map<std::string, QC_Read *> &read_lookup) {
	std::string name = make_read_name(header);
	/* read name has to match string, if present */
	if (!opt_readname_match.empty() && opt_readname_match.find(name) == opt_readname_match.end()) {
		return 0;
	}
	if (read_lookup.find(name) != read_lookup.end()) {
		fprintf(stderr, "Warning: duplicate read sequence: %s\n", name.c_str());
		return 0;
	} else {
		read_list.push_back(QC_Read(header));
		read_lookup[name] = &(*read_list.rbegin());
		return 1;
	}
}

/* calculates the various derived stats for qc_reads */

void qc_calc_stats(std::list<QC_Read> &read_list, std::map<size_t, unsigned int> &n_hist, std::map<size_t, unsigned int> &lq_hist) {
	std::list<QC_Read>::iterator a = read_list.begin();
	std::list<QC_Read>::iterator end = read_list.end();
	for (; a != end; ++a) {
		a->calc_stats(n_hist, lq_hist);
	}
}

/* read in contig sequence */

int qc_read_sequence(const char *filename, std::list<QC_Read> &read_list, bool opt_warnings) {
	int fd = open_compressed(filename);
	if (fd == -1) {
		return -1;
	}
	std::map<std::string, QC_Read *> read_lookup;
	bool adding_read = 0;
	std::string line, data;
	while (pfgets(fd, line) != -1) {
		if (line[0] == '>') {
			if (adding_read) {
				read_list.back().add_sequence(data);
			}
			data.clear();
			adding_read = add_read(line, read_list, read_lookup);
		} else {
			data += line;
		}
	}
	close_compressed(fd);
	if (adding_read) {
		read_list.back().add_sequence(data);
	}
	fd = open_compressed(make_qual_filename(filename));
	if (fd == -1 && opt_strip_tracename) {
		fd = open_compressed(make_qual_filename(filename, 1));
	}
	if (fd == -1) {
		return -1;
	}
	QC_Read *current_read = NULL;
	while (pfgets(fd, line) != -1) {
		if (line[0] == '>') {
			if (current_read != NULL) {
				current_read->add_quality(data, opt_warnings);
			}
			data.clear();
			std::string name = make_read_name(line);
			std::map<std::string, QC_Read *>::iterator a = read_lookup.find(name);
			if (a == read_lookup.end()) {
				fprintf(stderr, "Warning: no sequence for quality: %s\n", name.c_str());
				current_read = NULL;
			} else if (a->second == NULL) {
				fprintf(stderr, "Warning: duplicate read quality: %s\n", name.c_str());
				current_read = NULL;
			} else {
				current_read = a->second;
				a->second = NULL;
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
