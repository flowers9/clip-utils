#ifndef _READ_FILE_H
#define _READ_FILE_H

#include "open_compressed.h"	// close_compressed(), open_compressed()
#include "read.h"	// Read, opt_quality_cutoff
#include <list>		// list<>
#include <map>		// map<>
#include <stdio.h>	// fprintf(), stderr
#include <string>	// string

extern bool opt_strip_tracename;
extern std::map<std::string, bool> opt_readname_match;

class ReadFile {
    private:
	int find_qual(void);
	void add_read(const std::string &);
	void add_read(const std::string &, const std::string &, bool);
	void add_quality(const std::string &, bool);
	void set_default_quals(unsigned char);
	void transfer_reads(bool);
	void consistency_check(void) const;
	std::string make_read_name(std::string &);
	void check_fastq(void);
	int read_all_fastq(bool);
	int read_batch_fastq(bool);
    private:
	int fd_seq, fd_qual;
	int track_dups, fastq_file;
	size_t batch_size;
	std::string sheader, qheader;
	std::map<std::string, Read *> read_lookup;
	std::list<Read> tmp_read_list;
	std::map<std::string, std::string> spare_quals;
	Read duplicate_read;
	Read * const duplicate_read_ptr;
    public:
	std::list<Read> read_list;
	std::string seq_file, qual_file;
	ReadFile(void) : fd_seq(-1), fd_qual(-1), track_dups(1), fastq_file(0), batch_size(0), duplicate_read_ptr(&duplicate_read) { }
	explicit ReadFile(const char * const s, const size_t b = 0, const bool t = 1) : fd_seq(-1), fd_qual(-1), track_dups(t ? 1 : 0), fastq_file(0), batch_size(b), duplicate_read_ptr(&duplicate_read), seq_file(s) {
		this->open();
	}
	explicit ReadFile(const std::string &s, const size_t b = 0, const bool t = 1) : fd_seq(-1), fd_qual(-1), track_dups(t ? 1 : 0), fastq_file(0), batch_size(b), duplicate_read_ptr(&duplicate_read), seq_file(s) {
		this->open();
	}
	~ReadFile(void) {
		this->close();
	}
	void open(void) {
		switch(find_qual()) {	// sets fd_seq
		    case 0:		// no qual file
			if (!fastq_file) {
				fprintf(stderr, "Warning: %s: qual file missing, defaulting qual's to %d\n", seq_file.c_str(), opt_quality_cutoff);
			}
		    case -1:		// no qual or seq file
			break;
		    default:
			if (fd_seq != -1) {
				fd_qual = open_compressed(qual_file);
			}
		}
	}
	void close(void) {
		if (fd_seq != -1) {
			close_compressed(fd_seq);
			fd_seq = -1;
		}
		if (fd_qual != -1) {
			close_compressed(fd_qual);
			fd_qual = -1;
		}
		consistency_check();
		read_lookup.clear();
		tmp_read_list.clear();
		spare_quals.clear();
		sheader.clear();	// in case of close before file end
		qheader.clear();
	}
	void reset(void) {
		this->close();
		this->open();
	}
	int read_all(bool);
	int read_batch(bool);
	void mask_by_phred(unsigned int);
	void set_batch_size(const size_t i) {
		batch_size = i;
		track_dups = i == 0 ? 1 : 0;
	}
	void set_track_dups(const bool i) {
		track_dups = i ? 1 : 0;
	}
};

#endif // !_READ_FILE_H
