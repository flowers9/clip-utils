#include "open_compressed.h"	// find_suffix(), pfgets(), pfpeek()
#include "read.h"	// Read, opt_quality_cutoff
#include "read_file.h"
#include <list>		// list<>
#include <map>		// map<>
#include <string>	// string
#include <sys/types.h>	// size_t
#include <unistd.h>	// pathconf(), _PC_PATH_MAX

bool opt_strip_tracename(0);
std::map<std::string, bool> opt_readname_match;

std::string ReadFile::make_read_name(std::string &header) {
	if (opt_strip_tracename) {
		const size_t i(header.find(' ', 1));
		if (i != std::string::npos) {
			header = header[0] + header.substr(i + 1);
		}
	}
	return get_name(header);
}

// goes through read list and masks any basepairs with a quality less than
// the cutoff

void ReadFile::mask_by_phred(const unsigned int opt_phred_mask_cutoff) {
	std::list<Read>::iterator a(read_list.begin());
	const std::list<Read>::const_iterator end_a(read_list.end());
	for (; a != end_a; ++a) {
		a->mask_by_phred(opt_phred_mask_cutoff);
	}
}

// return 1 if we find a qual file, 0 if not

static int check_file(const std::string &name, std::string &qual_file) {
	qual_file = name + ".qual";
	std::string qual_suffix;
	if (find_suffix(qual_file, qual_suffix) == 0) {
		return 1;
	}
	if (name.size() > 4 && name.compare(name.size() - 4, 4, ".fna") == 0) {
		qual_file = name.substr(0, name.size() - 3) + "qual";
		if (find_suffix(qual_file, qual_suffix) == 0) {
			return 1;
		}
	}
	if (name.size() > 6 && name.compare(name.size() - 6, 6, ".fasta") == 0) {
		qual_file = name.substr(0, name.size() - 5) + "qual";
		if (find_suffix(qual_file, qual_suffix) == 0) {
			return 1;
		}
	}
	if (name.size() > 1 && name[0] == 'f' && name.find_first_not_of("0123456789", 1) == std::string::npos) {
		qual_file = "q" + name.substr(1);
		if (find_suffix(qual_file, qual_suffix) == 0) {
			return 1;
		}
	}
	qual_file.clear();
	return 0;
}

void ReadFile::check_fastq(void) {
	fd_seq = open_compressed(seq_file);
	if (fd_seq == -1) {
		return;
	}
	char c;
	if (pfpeek(fd_seq, &c, 1) == 1 && c == '@') {
		fastq_file = 1;
	}
}

// make the quality file name from the sequence file name - mainly just
// tacks on .qual, but also handles Z/gz/bz2 endings (since the .qual
// needs to be added before the compression suffix); also opens seq_file

int ReadFile::find_qual() {
	if (seq_file.empty() || seq_file == "-") {
		check_fastq();
		return 0;
	}
	// first find actual seq name (and suffix)
	std::string suffix;
	if (find_suffix(seq_file, suffix) == -1) {
		return -1;
	}
	check_fastq();
	if (fastq_file) {
		return 0;
	}
	std::string name(seq_file, 0, seq_file.size() - suffix.size());
	if (check_file(name, qual_file)) {
		return 1;
	}
	std::string::size_type k(name.rfind("/fasta/"));
	if (k != std::string::npos) {
		name.replace(k + 1, 5, "qual");
		if (check_file(name, qual_file)) {
			return 1;
		}
	}
	const size_t buf_size(pathconf(".", _PC_PATH_MAX));
	char buf[buf_size + 1];
	const ssize_t j(readlink(seq_file.c_str(), buf, buf_size));
	if (j != -1) {
		buf[j] = 0;
		name = buf;
		get_suffix(name, suffix);
		if (check_file(name, qual_file)) {
			return 1;
		}
	}
	return 0;
}

// print errors for unmatched parts

void ReadFile::consistency_check() const {
	std::list<Read>::const_iterator a(tmp_read_list.begin());
	const std::list<Read>::const_iterator end_a(tmp_read_list.end());
	for (; a != end_a; ++a) {
		fprintf(stderr, "Warning: no quality for sequence: %s\n", a->name().c_str());
	}
	std::map<std::string, std::string>::const_iterator b(spare_quals.begin());
	const std::map<std::string, std::string>::const_iterator end_b(spare_quals.end());
	for (; b != end_b; ++b) {
		fprintf(stderr, "Warning: no sequence for quality: %s\n", b->first.c_str());
	}
}

// add quality vector of equal length to sequence, of value x

void ReadFile::set_default_quals(const unsigned char x) {
	std::list<Read>::iterator a(tmp_read_list.begin());
	const std::list<Read>::const_iterator end_a(tmp_read_list.end());
	for (; a != end_a; ++a) {
		a->set_quality(x);
	}
}

// transfer reads with quality from tmp_read_list to read_list

void ReadFile::transfer_reads(const bool opt_warnings) {
	std::list<Read>::iterator a(tmp_read_list.begin());
	const std::list<Read>::const_iterator end_a(tmp_read_list.end());
	std::list<Read>::iterator b(a);
	while (a != end_a) {
		for (; a != end_a && a->has_quality(); ++a) { }
		read_list.splice(read_list.end(), tmp_read_list, b, a);
		while (a != end_a && !a->has_quality()) {
			std::map<std::string, std::string>::iterator c(spare_quals.find(a->name()));
			if (c != spare_quals.end()) {
				a->add_quality(c->second, opt_warnings);
				spare_quals.erase(c);
				read_list.splice(read_list.end(), tmp_read_list, a++);
			} else {
				++a;
			}
		}
		b = a;
	}
	// clear out reads that got transferred
	if (!track_dups) {
		std::map<std::string, Read *>::iterator c(read_lookup.begin());
		const std::map<std::string, Read *>::const_iterator end_c(read_lookup.end());
		while (c != end_c) {
			if (c->second == duplicate_read_ptr) {
				read_lookup.erase(c++);
			} else {
				++c;
			}
		}
	}
}

void ReadFile::add_read(const std::string &data) {
	if (sheader.empty()) {
		return;
	}
	const std::string name(make_read_name(sheader));
	// read name has to match string, if present
	if (!opt_readname_match.empty() && opt_readname_match.find(name) == opt_readname_match.end()) {
		read_lookup[name] = 0;
	} else if (read_lookup.find(name) != read_lookup.end()) {
		fprintf(stderr, "Warning: duplicate read sequence: %s\n", name.c_str());
	} else {
		tmp_read_list.push_back(Read(sheader, data));
		read_lookup[name] = &tmp_read_list.back();
	}
}

void ReadFile::add_quality(const std::string &data, const bool opt_warnings) {
	if (qheader.empty()) {
		return;
	}
	std::string name(make_read_name(qheader));
	std::map<std::string, Read *>::iterator a(read_lookup.find(name));
	if (a == read_lookup.end()) {
		spare_quals[name] = data;
	} else if (a->second == duplicate_read_ptr) {
		fprintf(stderr, "Warning: duplicate read quality: %s\n", name.c_str());
	} else if (a->second == 0) {
		a->second = duplicate_read_ptr;
	} else {
		a->second->add_quality(data, opt_warnings);
		a->second = duplicate_read_ptr;
	}
}

// read in contig sequence

int ReadFile::read_all(const bool opt_warnings) {
	if (fd_seq == -1) {
		return -1;
	}
	if (fastq_file) {
		return read_all_fastq(opt_warnings);
	}
	std::string line, data;
	while (pfgets(fd_seq, line) != -1) {
		if (line[0] == '>') {
			add_read(data);
			data.clear();
			sheader = line;
		} else {
			data += line;
		}
	}
	add_read(data);
	data.clear();
	sheader.clear();
	if (fd_qual == -1) {
		set_default_quals(opt_quality_cutoff);
		read_list.splice(read_list.end(), tmp_read_list);
		this->close();
		return 1;
	}
	while (pfgets(fd_qual, line) != -1) {
		if (line[0] == '>') {
			add_quality(data, opt_warnings);
			data.clear();
			qheader = line;
		} else {
			data += line + ' ';
		}
	}
	add_quality(data, opt_warnings);
	qheader.clear();
	transfer_reads(opt_warnings);
	this->close();
	return 0;
}

int ReadFile::read_batch(const bool opt_warnings) {
	read_list.clear();
	if (batch_size == 0) {
		return read_all(opt_warnings);
	}
	if (fd_seq == -1) {
		return -1;
	}
	if (fastq_file) {
		return read_batch_fastq(opt_warnings);
	}
	size_t i(0);
	std::string line, data;
	while (pfgets(fd_seq, line) != -1) {
		if (line[0] == '>') {
			add_read(data);
			data.clear();
			sheader = line;
			if (++i == batch_size) {
				break;
			}
		} else {
			data += line;
		}
	}
	if (i != batch_size) {
		add_read(data);
		sheader.clear();
		return read_all(opt_warnings);	// read qual, close file
	}
	if (fd_qual == -1) {
		set_default_quals(opt_quality_cutoff);
		read_list.splice(read_list.end(), tmp_read_list);
		return 1;
	}
	i = 0;
	while (pfgets(fd_qual, line) != -1) {
		if (line[0] == '>') {
			add_quality(data, opt_warnings);
			data.clear();
			qheader = line;
			if (++i == batch_size) {
				break;
			}
		} else {
			data += line + ' ';
		}
	}
	transfer_reads(opt_warnings);
	return 0;
}

void ReadFile::add_read(const std::string &seq, const std::string &qual, const bool opt_warnings) {
	if (sheader.empty()) {
		return;
	}
	sheader[0] = '>';
	const std::string name(make_read_name(sheader));
	// read name has to match string, if present
	if (opt_readname_match.empty() || opt_readname_match.find(name) != opt_readname_match.end()) {
		read_list.push_back(Read(sheader, seq, qual, opt_warnings));
	}
}

int ReadFile::read_all_fastq(const bool opt_warnings) {
	while (pfgets(fd_seq, sheader) != -1) {
		if (sheader[0] == '@') {
			std::string seq, line, qual;
			if (pfgets(fd_seq, seq) == -1) {
				break;
			}
			if (pfgets(fd_seq, line) == -1) { // + line - discard
				break;
			}
			if (pfgets(fd_seq, qual) == -1) {
				break;
			}
			add_read(seq, qual, opt_warnings);
		}
	}
	this->close();
	return 0;
}

int ReadFile::read_batch_fastq(const bool opt_warnings) {
	size_t i(0);
	while (pfgets(fd_seq, sheader) != -1) {
		if (sheader[0] == '@') {
			std::string seq, line, qual;
			if (pfgets(fd_seq, seq) == -1) {
				break;
			}
			if (pfgets(fd_seq, line) == -1) { // + line - discard
				break;
			}
			if (pfgets(fd_seq, qual) == -1) {
				break;
			}
			add_read(seq, qual, opt_warnings);
			if (++i == batch_size) {
				return 0;
			}
		}
	}
	this->close();
	return 0;
}
