#include "hashl.h"	// hashl
#include "hashl_metadata.h"
#include "open_compressed.h"	// close_compressed(), open_compressed(), pfgets()
#include <iostream>	// cerr, cout
#include <stdint.h>	// uint64_t
#include <stdlib.h>	// exit()
#include <string.h>	// memcpy()
#include <string>	// string
#include <utility>	// make_pair(), pair<>
#include <vector>	// vector<>

// add metadata from another hash to ours
void hashl_metadata::add(hashl_metadata &a, const size_t padding) {
	if (padding) {				// synchronize with hashl's data padding
		add_read("__padding__");
		add_read_range(0, padding);
	}
	files.reserve(files.size() + a.files.size());
	std::move(a.files.begin(), a.files.end(), std::back_inserter(files));
	reads.reserve(reads.size() + a.reads.size());
	std::move(a.reads.begin(), a.reads.end(), std::back_inserter(reads));
	read_ranges.reserve(read_ranges.size() + a.read_ranges.size());
	std::move(a.read_ranges.begin(), a.read_ranges.end(), std::back_inserter(read_ranges));
}

// return 0-3, exit for bad sequence
static hashl::base_type convert_char(const char c) {
	switch (c) {
	    case 'A':
	    case 'a':
		return 0;
	    case 'C':
	    case 'c':
		return 1;
	    case 'G':
	    case 'g':
		return 2;
	    case 'T':
	    case 't':
		return 3;
	    default:
		std::cerr << "Error: non-ACGT basepair: " << static_cast<char>(c) << "\n";
		exit(1);
	}
}

void hashl_metadata::get_subreads(const std::string &seq, const std::vector<std::pair<uint64_t, uint64_t> > &ranges) {
	for (size_t i(0); i < ranges.size(); ++i) {
		for (size_t j(ranges[i].first); j < ranges[i].second; ++j) {
			if (bit_offset) {
				bit_offset -= 2;
				data[byte_offset] |= convert_char(seq[j]) << bit_offset;
			} else {
				bit_offset = sizeof(hashl::base_type) * 8 - 2;
				data[++byte_offset] = convert_char(seq[j]) << bit_offset;
			}
		}
	}
}

void hashl_metadata::read_file(const size_t i) {
	const int fd(open_compressed(files[i]));
	if (fd == -1) {
		std::cerr << "Error: open: " << files[i] << "\n";
		exit(1);
	}
	size_t j(0);	// which read we're on (but only incremented for reads we use)
	std::string line, seq;
	if (pfgets(fd, line) == -1) {		// empty file
		std::cerr << "Error: File is empty: " << files[i] << "\n";
		exit(1);
	} else if (line[0] == '>') {		// fasta file
		std::string header;
		do {
			const std::string &read_name(reads[i][j]);
			if (!read_name.compare(0, read_name.size(), line, 1, read_name.size()) && (line.size() == read_name.size() + 1 || isspace(line[read_name.size() + 1]))) {
				seq.clear();
				while (pfgets(fd, line) != -1 && line[0] != '>') {
					seq += line;
				}
				get_subreads(seq, read_ranges[i][j]);
				++j;
			} else {
				while (pfgets(fd, line) != -1 && line[0] != '>') { }
			}
		} while (j < reads[i].size() && line[0] == '>');
	} else if (line[0] == '@') {		// fastq file
		do {
			if (pfgets(fd, seq) == -1) {		// sequence
				std::cerr << "Error: truncated fastq file: " << files[i] << "\n";
				exit(1);
			}
			const std::string &read_name(reads[i][j]);
			if (!read_name.compare(0, read_name.size(), line, 1, read_name.size()) && (line.size() == read_name.size() + 1 || isspace(line[read_name.size() + 1]))) {
				get_subreads(seq, read_ranges[i][j]);
				++j;
			}
			// skip quality header and quality
			// (use seq because it'll be the same length as quality)
			if (pfgets(fd, line) == -1 || pfgets(fd, seq) == -1) {
				std::cerr << "Error: truncated fastq file: " << files[i] << "\n";
				exit(1);
			}
		} while (j < reads[i].size() && pfgets(fd, line) != -1);
	} else {
		std::cerr << "Error: unknown file format: " << files[i] << "\n";
		exit(1);
	}
	if (j < reads[i].size()) {
		std::cerr << "Error: File is shorter than before: " << files[i] << "\n";
		exit(1);
	}
	close_compressed(fd);
}

// mostly const - only the convenience variables get modified

void hashl_metadata::read_data(std::vector<hashl::base_type> &data_out, const bool feedback) {
	// convert size from basepairs to length of hashl::base_type array
	const size_t data_size((2 * sequence_length() + sizeof(hashl::base_type) * 8 - 1) / (sizeof(hashl::base_type) * 8));
	data.assign(data_size, 0);
	byte_offset = 0;
	bit_offset = sizeof(hashl::base_type) * 8;
	for (size_t i(0); i < files.size(); ++i) {
		if (feedback) {
			std::cerr << time(0) << ": Reading in " << files[i] << "\n";
		}
		read_file(i);
	}
	data_out.swap(data);
	data.clear();		// just in case data_out wasn't empty
}

std::pair<size_t, size_t> hashl_metadata::total_reads() const {
	size_t read_count(0), subread_count(0);
	for (size_t i(0); i < read_ranges.size(); ++i) {
		read_count += read_ranges[i].size();
		for (size_t j(0); j < read_ranges[i].size(); ++j) {
			subread_count += read_ranges[i][j].size();
		}
	}
	return std::make_pair(read_count, subread_count);
}

size_t hashl_metadata::max_kmers(const size_t mer_length) const {
	size_t x(0);
	for (size_t i(0); i < read_ranges.size(); ++i) {
		for (size_t j(0); j < read_ranges[i].size(); ++j) {
			for (size_t k(0); k < read_ranges[i][j].size(); ++k) {
				x += read_ranges[i][j][k].second - read_ranges[i][j][k].first - mer_length + 1;
			}
		}
	}
	return x;
}

size_t hashl_metadata::sequence_length() const {
	size_t x(0);
	for (size_t i(0); i < read_ranges.size(); ++i) {
		for (size_t j(0); j < read_ranges[i].size(); ++j) {
			for (size_t k(0); k < read_ranges[i][j].size(); ++k) {
				x += read_ranges[i][j][k].second - read_ranges[i][j][k].first;
			}
		}
	}
	return x;
}

std::vector<size_t> hashl_metadata::read_ends() const {
	std::vector<size_t> list;
	size_t x(0);
	for (size_t i(0); i < read_ranges.size(); ++i) {
		for (size_t j(0); j < read_ranges[i].size(); ++j) {
			for (size_t k(0); k < read_ranges[i][j].size(); ++k) {
				x += read_ranges[i][j][k].second - read_ranges[i][j][k].first;
				list.push_back(x);
			}
		}
	}
	return list;
}

void hashl_metadata::add_file(const std::string &file_name) {
	finalize();
	files.push_back(file_name);
	reads.push_back(std::vector<std::string>());
	read_ranges.push_back(std::vector<std::vector<std::pair<uint64_t, uint64_t> > >());
}

void hashl_metadata::add_read(const std::string &read_name) {
	reads.back().push_back(read_name);
	read_ranges.back().push_back(std::vector<std::pair<uint64_t, uint64_t> >());
}

void hashl_metadata::add_read_range(const uint64_t start, const uint64_t end) {
	read_ranges.back().back().push_back(std::make_pair(start, end));
}

void hashl_metadata::finalize() {
	if (!files.empty()) {
		// check if last read of last file had any ranges
		if (read_ranges.back().back().empty()) {
			read_ranges.back().pop_back();
			reads.back().pop_back();
		}
		// check if last file had any reads
		if (read_ranges.back().empty()) {
			read_ranges.pop_back();
			reads.pop_back();
			files.pop_back();
		}
	}
}

void hashl_metadata::pack(std::vector<char> &d) const {
	// count space needed
	size_t metadata_size(sizeof(uint64_t));				// number of files
	for (size_t i(0); i < files.size(); ++i) {
		metadata_size += files[i].size() + 1;			// file name and null
		metadata_size += sizeof(uint64_t);			// number of reads
		for (size_t j(0); j < reads[i].size(); ++j) {
			metadata_size += reads[i][j].size() + 1;	// read name and null
			metadata_size += sizeof(uint64_t);		// number of ranges
			metadata_size += read_ranges[i][j].size() * sizeof(uint64_t) * 2;
		}
	}
	// allocate space
	d.assign(metadata_size, 0);
	// fill space with metadata
	size_t offset(0);
	uint64_t tmp;
	memcpy(&d[offset], &(tmp = files.size()), sizeof(tmp));
	offset += sizeof(tmp);
	for (size_t i(0); i < files.size(); ++i) {
		memcpy(&d[offset], files[i].c_str(), files[i].size() + 1);
		offset += files[i].size() + 1;
		memcpy(&d[offset], &(tmp = reads[i].size()), sizeof(tmp));
		offset += sizeof(tmp);
		for (size_t j(0); j < reads[i].size(); ++j) {
			memcpy(&d[offset], reads[i][j].c_str(), reads[i][j].size() + 1);
			offset += reads[i][j].size() + 1;
			memcpy(&d[offset], &(tmp = read_ranges[i][j].size()), sizeof(tmp));
			offset += sizeof(tmp);
			for (size_t k(0); k < read_ranges[i][j].size(); ++k) {
				memcpy(&d[offset], &read_ranges[i][j][k].first, sizeof(uint64_t));
				offset += sizeof(uint64_t);
				memcpy(&d[offset], &read_ranges[i][j][k].second, sizeof(uint64_t));
				offset += sizeof(uint64_t);
			}
		}
	}
}

void hashl_metadata::unpack(const std::vector<char> &d) {
	size_t offset(0);
	uint64_t file_count;
	memcpy(&file_count, &d[offset], sizeof(file_count));
	offset += sizeof(file_count);
	files.clear();
	files.reserve(file_count);
	reads.assign(file_count, std::vector<std::string>());
	read_ranges.assign(file_count, std::vector<std::vector<std::pair<uint64_t, uint64_t> > >());
	for (size_t i(0); i < file_count; ++i) {
		files.push_back(&d[offset]);			// null terminated c-string
		offset += files.back().size() + 1;
		uint64_t read_count;
		memcpy(&read_count, &d[offset], sizeof(read_count));
		offset += sizeof(read_count);
		reads[i].reserve(read_count);
		read_ranges[i].assign(read_count, std::vector<std::pair<uint64_t, uint64_t> >());
		for (size_t j(0); j < read_count; ++j) {
			reads[i].push_back(&d[offset]);		// null terminated c-string
			offset += reads[i].back().size() + 1;
			uint64_t read_range_count;
			memcpy(&read_range_count, &d[offset], sizeof(read_range_count));
			offset += sizeof(read_range_count);
			read_ranges[i][j].reserve(read_range_count);
			for (size_t k(0); k < read_range_count; ++k) {
				uint64_t start, end;
				memcpy(&start, &d[offset], sizeof(start));
				offset += sizeof(start);
				memcpy(&end, &d[offset], sizeof(end));
				offset += sizeof(end);
				read_ranges[i][j].push_back(std::make_pair(start, end));
			}
		}
	}
	if (d.size() != offset) {
		std::cerr << "Error: metadata size mismatch: " << d.size() << " != " << offset << "\n";
		exit(1);
	}
}

void hashl_metadata::print() const {
	for (size_t i(0); i < files.size(); ++i) {
		std::cout << files[i].c_str() << "\n";
		for (size_t j(0); j < reads[i].size(); ++j) {
			std::cout << "\t" << reads[i][j].c_str() << "\n";
			for (size_t k(0); k < read_ranges[i][j].size(); ++k) {
				std::cout << "\t\t" << read_ranges[i][j][k].first << ' ' << read_ranges[i][j][k].second << "\n";
			}
		}
	}
}
