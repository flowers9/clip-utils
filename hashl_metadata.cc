#include "hashl.h"	// hashl
#include "hashl_metadata.h"
#include "open_compressed.h"	// close_compressed(), open_compressed(), pfgets()
#include <iostream>	// cerr, cout
#include <map>		// map<>
#include <stdlib.h>	// exit()
#include <string.h>	// memcpy()
#include <string>	// string
#include <utility>	// as_const(), make_pair(), pair<>
#include <vector>	// vector<>

// add metadata from another hash to ours
void hashl_metadata::add(hashl_metadata &a, const size_type padding) {
	if (padding) {				// synchronize with hashl's data padding
		add_readname("__padding__");
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

void hashl_metadata::get_subreads(const std::string &seq, const std::vector<std::pair<size_type, size_type> > &ranges) {
	for (const auto &a : ranges) {
		for (size_type i = a.first; i < a.second; ++i) {
			if (bit_offset) {
				bit_offset -= 2;
				data[byte_offset] |= convert_char(seq[i]) << bit_offset;
			} else {
				bit_offset = sizeof(hashl::base_type) * 8 - 2;
				data[++byte_offset] = convert_char(seq[i]) << bit_offset;
			}
		}
	}
}

void hashl_metadata::read_file(const size_type i) {
	const int fd = open_compressed(files[i]);
	if (fd == -1) {
		std::cerr << "Error: open: " << files[i] << "\n";
		exit(1);
	}
	size_type j = 0;	// which read we're on (but only incremented for reads we use)
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
	const size_type data_size = (2 * sequence_length() + sizeof(hashl::base_type) * 8 - 1) / (sizeof(hashl::base_type) * 8);
	data.assign(data_size, 0);
	byte_offset = 0;
	bit_offset = sizeof(hashl::base_type) * 8;
	for (size_type i = 0; i < files.size(); ++i) {
		if (feedback) {
			std::cerr << time(0) << ": Reading in " << files[i] << "\n";
		}
		read_file(i);
	}
	data_out.swap(data);
	data = std::vector<hashl::base_type>();		// just in case data_out wasn't empty
}

std::pair<hashl_metadata::size_type, hashl_metadata::size_type> hashl_metadata::total_reads() const {
	size_type read_count = 0, subread_count = 0;
	for (const auto &a : read_ranges) {
		read_count += a.size();
		for (const auto &b : a) {
			subread_count += b.size();
		}
	}
	return std::make_pair(read_count, subread_count);
}

hashl_metadata::size_type hashl_metadata::max_kmers(const size_type mer_length) const {
	size_type x = 0;
	for (const auto &a : read_ranges) {
		for (const auto &b : a) {
			for (const auto &c : b) {
				x += c.second - c.first - mer_length + 1;
			}
		}
	}
	return x;
}

hashl_metadata::size_type hashl_metadata::sequence_length() const {
	size_type x = 0;
	for (const auto &a : read_ranges) {
		for (const auto &b : a) {
			for (const auto &c : b) {
				x += c.second - c.first;
			}
		}
	}
	return x;
}

std::vector<hashl_metadata::size_type> hashl_metadata::read_ends() const {
	std::vector<size_type> list;
	size_type x = 0;
	for (const auto &a : read_ranges) {
		for (const auto &b : a) {
			for (const auto &c : b) {
				x += c.second - c.first;
				list.push_back(x);
			}
		}
	}
	return list;
}

void hashl_metadata::add_filename(const std::string &file_name) {
	files.push_back(file_name);
	reads.push_back(std::vector<std::string>());
	read_ranges.push_back(std::vector<std::vector<std::pair<size_type, size_type> > >());
}

void hashl_metadata::add_readname(const std::string &read_name) {
	reads.back().push_back(read_name);
	read_ranges.back().push_back(std::vector<std::pair<size_type, size_type> >());
}

void hashl_metadata::add_read_range(const size_type start, const size_type end) {
	read_ranges.back().back().push_back(std::make_pair(start, end));
}

void hashl_metadata::finalize_file() {
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
	size_type metadata_size = sizeof(size_type);			// number of files
	for (size_type i = 0; i < files.size(); ++i) {
		metadata_size += files[i].size() + 1;			// file name and null
		metadata_size += sizeof(size_type);			// number of reads
		for (size_type j = 0; j < reads[i].size(); ++j) {
			metadata_size += reads[i][j].size() + 1;	// read name and null
			metadata_size += sizeof(size_type);		// number of ranges
			metadata_size += read_ranges[i][j].size() * sizeof(size_type) * 2;
		}
	}
	// allocate space
	d.assign(metadata_size, 0);
	// fill space with metadata
	size_type offset = 0, tmp;
	memcpy(&d[offset], &(tmp = files.size()), sizeof(tmp));
	offset += sizeof(tmp);
	for (size_type i = 0; i < files.size(); ++i) {
		memcpy(&d[offset], files[i].c_str(), files[i].size() + 1);
		offset += files[i].size() + 1;
		memcpy(&d[offset], &(tmp = reads[i].size()), sizeof(tmp));
		offset += sizeof(tmp);
		for (size_type j = 0; j < reads[i].size(); ++j) {
			memcpy(&d[offset], reads[i][j].c_str(), reads[i][j].size() + 1);
			offset += reads[i][j].size() + 1;
			memcpy(&d[offset], &(tmp = read_ranges[i][j].size()), sizeof(tmp));
			offset += sizeof(tmp);
			for (size_type k = 0; k < read_ranges[i][j].size(); ++k) {
				memcpy(&d[offset], &read_ranges[i][j][k].first, sizeof(size_type));
				offset += sizeof(size_type);
				memcpy(&d[offset], &read_ranges[i][j][k].second, sizeof(size_type));
				offset += sizeof(size_type);
			}
		}
	}
}

void hashl_metadata::unpack(const std::vector<char> &d) {
	size_type offset = 0, file_count;
	memcpy(&file_count, &d[offset], sizeof(file_count));
	offset += sizeof(file_count);
	files.clear();
	files.reserve(file_count);
	reads.assign(file_count, std::vector<std::string>());
	read_ranges.assign(file_count, std::vector<std::vector<std::pair<size_type, size_type> > >());
	for (size_type i = 0; i < file_count; ++i) {
		files.push_back(&d[offset]);			// null terminated c-string
		offset += files.back().size() + 1;
		size_type read_count;
		memcpy(&read_count, &d[offset], sizeof(read_count));
		offset += sizeof(read_count);
		reads[i].reserve(read_count);
		read_ranges[i].assign(read_count, std::vector<std::pair<size_type, size_type> >());
		for (size_type j = 0; j < read_count; ++j) {
			reads[i].push_back(&d[offset]);		// null terminated c-string
			offset += reads[i].back().size() + 1;
			size_type read_range_count;
			memcpy(&read_range_count, &d[offset], sizeof(read_range_count));
			offset += sizeof(read_range_count);
			read_ranges[i][j].reserve(read_range_count);
			for (size_type k = 0; k < read_range_count; ++k) {
				size_type start, end;
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
	for (size_type i = 0; i < files.size(); ++i) {
		std::cout << files[i].c_str() << "\n";
		for (size_type j = 0; j < reads[i].size(); ++j) {
			std::cout << "\t" << reads[i][j].c_str() << "\n";
			for (size_type k = 0; k < read_ranges[i][j].size(); ++k) {
				std::cout << "\t\t" << read_ranges[i][j][k].first << ' ' << read_ranges[i][j][k].second << "\n";
			}
		}
	}
}

// create a map allowing translation of a data position to a file/read/read_start triplet
void hashl_metadata::create_lookup_map(std::map<size_type, hashl_metadata::position> &lookup) const {
	lookup.clear();
	size_type offset = 0;
	position x;
	for (x.file = 0; x.file < read_ranges.size(); ++x.file) {
		const std::vector<std::vector<std::pair<size_type, size_type> > > &read_list = read_ranges[x.file];
		for (x.read = 0; x.read < read_list.size(); ++x.read) {
			const std::vector<std::pair<size_type, size_type> > &range_list = read_list[x.read];
			for (size_type k = 0; k < range_list.size(); ++k) {
				x.read_start = range_list[k].first;
				lookup[offset] = x;
				offset += range_list[k].second - range_list[k].first;
			}
		}
	}
}

// update read_ranges to just the subset included in kept_offsets
// reads without ranges are removed, as are files without reads

void hashl_metadata::update_ranges(const std::vector<std::pair<size_type, size_type> > &kept_offsets) {
	std::vector<std::string> new_files;
	std::vector<std::vector<std::string> > new_reads;
	std::vector<std::vector<std::vector<std::pair<size_type, size_type> > > > new_read_ranges;
	size_type current_offset = 0;		// used to convert read ranges into data offsets
	auto kept = std::as_const(kept_offsets).begin();
	for (size_type i = 0; kept != kept_offsets.end() && i < read_ranges.size(); ++i) {
		new_reads.push_back(std::vector<std::string>());
		new_read_ranges.push_back(std::vector<std::vector<std::pair<size_type, size_type> > >());
		const auto &file_read_ranges = read_ranges[i];
		for (size_type j = 0; kept != kept_offsets.end() && j < file_read_ranges.size(); ++j) {
			new_read_ranges.back().push_back(std::vector<std::pair<size_type, size_type> >());
			for (const auto &a : file_read_ranges[j]) {
				// find any kept_range that is in the range's offset range
				const size_type stop = current_offset + a.second - a.first;
				for (; kept != kept_offsets.end() && kept->first < stop; ++kept) {
					// convert kept range to a read range and save
					const size_type start = a.first + kept->first - current_offset;
					new_read_ranges.back().back().push_back(std::make_pair(start, start + kept->second - kept->first));
				}
				current_offset = stop;
			}
			if (new_read_ranges.back().back().empty()) {
				new_read_ranges.back().pop_back();
			} else {
				new_reads.back().push_back(reads[i][j]);
			}
		}
		if (new_reads.back().empty()) {
			new_read_ranges.pop_back();
			new_reads.pop_back();
		} else {
			new_files.push_back(files[i]);
		}
	}
	files.swap(new_files);
	reads.swap(new_reads);
	read_ranges.swap(new_read_ranges);
}
