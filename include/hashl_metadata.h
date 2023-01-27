#ifndef _HASHL_METADATA_H
#define _HASHL_METADATA_H

#include "hashl.h"	// hashl
#include <stdint.h>	// uint64_t
#include <string>	// string
#include <utility>	// pair<>
#include <vector>	// vector<>

class hashl_metadata {
    private:						// convenience variables for read_data()
	std::vector<hashl::base_type> data;
	size_t byte_offset;
	int bit_offset;
    private:
	std::vector<std::string> files;
	std::vector<std::vector<std::string> > reads;
	// inclusive start, exclusive end
	std::vector<std::vector<std::vector<std::pair<uint64_t, uint64_t> > > > read_ranges;
    public:
	explicit hashl_metadata(void) { }
	~hashl_metadata(void) { }
	void add_file(const std::string &file_name);	// add new file
	void add_read(const std::string &read_name);	// add new read at current file
	void add_read_range(uint64_t, uint64_t);	// add new read range at current read
	void finalize(void);				// remove last adds if empty
	void read_data(std::vector<hashl::base_type> &data_out, bool feedback = 0);
	void pack(std::vector<char> &) const;		// create blob of our data
	void unpack(const std::vector<char> &);		// fill our data from blob
	void print(void) const;
	std::pair<size_t, size_t> total_reads(void) const;	// reads & read ranges
	size_t max_kmers(size_t mer_length) const;
	size_t sequence_length(void) const;
	std::vector<size_t> read_ends(void) const;
	void add(hashl_metadata &, size_t padding = 0);
    private:
	void read_file(size_t);
	void get_subreads(const std::string &, const std::vector<std::pair<uint64_t, uint64_t> > &);
};

#endif // !_HASHL_METADATA_H
