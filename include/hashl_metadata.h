#ifndef _HASHL_METADATA_H
#define _HASHL_METADATA_H

#include "hashl.h"	// hashl
#include <map>		// map<>
#include <string>	// string
#include <utility>	// pair<>
#include <vector>	// vector<>

class hashl_metadata {
    public:
	typedef typename std::vector<hashl::base_type>::size_type size_type;
	struct position {				// for doing position lookups
		size_type file, read, read_start;
	};
    private:						// convenience variables for read_data()
	std::vector<hashl::base_type> data;
	size_type byte_offset;
	int bit_offset;
    private:
	std::vector<std::string> files;
	std::vector<std::vector<std::string> > reads;	// reads associated with each file
	// inclusive start, exclusive end		// ranges associated with each read (non-acgt basepairs create breaks in the read)
	std::vector<std::vector<std::vector<std::pair<size_type, size_type> > > > read_ranges;
    public:
	explicit hashl_metadata(void) { }
	~hashl_metadata(void) { }
	void add_filename(const std::string &file_name);	// add new file
	void add_readname(const std::string &read_name);	// add new read for current file
	void add_read_range(size_type, size_type);	// add new read range for current read
	void finalize_file(void);			// remove last adds if empty
	void read_data(std::vector<hashl::base_type> &data_out, bool feedback = 0);
	void pack(std::vector<char> &) const;		// create blob of our data
	void unpack(const std::vector<char> &);		// fill our data from blob
	void print(void) const;
	std::pair<size_type, size_type> total_reads(void) const;	// reads & read ranges
	size_type max_kmers(size_type mer_length) const;
	size_type sequence_length(void) const;
	std::vector<size_type> read_ends(void) const;
	void add(hashl_metadata &, size_type padding = 0);
	void create_lookup_map(std::map<size_type, position> &) const;
	void update_ranges(const std::vector<std::pair<size_type, size_type> > &);
	size_type file_count(void) const {
		return files.size();
	}
	size_type read_count(const size_type i) const {
		return reads[i].size();
	}
	const std::string &file(const size_type i) const {
		return files[i];
	}
	const std::string &read(const size_type i, const size_type j) const {
		return reads[i][j];
	}
    private:
	void read_file(size_type);
	void get_subreads(const std::string &, const std::vector<std::pair<size_type, size_type> > &);
};

#endif // !_HASHL_METADATA_H
