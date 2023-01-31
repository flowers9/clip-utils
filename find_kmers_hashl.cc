#include "hashl.h"	// hashl
#include "hashl_metadata.h"	// hashl_metadata
#include "open_compressed.h"	// close_compressed(), get_suffix(), open_compressed()
#include "version.h"	// VERSION
#include "write_fork.h"	// close_fork(), write_fork()
#include <fstream>	// ofstream
#include <getopt.h>	// getopt(), optarg, optind
#include <iostream>	// cerr, cout, ostream
#include <map>		// map<>
#include <stdlib.h>	// exit()
#include <string>	// string
#include <utility>	// make_pair(), pair<>
#include <vector>	// vector<>

// given a set of reference hashes and a hash of kmers to search for,
// create a fasta file of the matched sections, with markers for sections
// that may be duplicated elsewhere

static std::string opt_reference_files;

static void print_usage() {
	std::cerr << "usage: find_kmers_hashl <kmer_list_hash> <reference_hash1> [reference_hash2 [...] ]\n"
		"    -h    print this help\n"
		"    -o ## output file for base reference file names [stderr]\n"
		"    -V    print version\n";
	exit(1);
}

static void get_opts(const int argc, char * const * const argv) {
	int c;
	while ((c = getopt(argc, argv, "ho:V")) != EOF) {
		switch (c) {
		    case 'h':
			print_usage();
			break;
		    case 'o':
			opt_reference_files = optarg;
			break;
		    case 'V':
			std::cerr << "find_kmers_hashl version " << VERSION << "\n";
			exit(0);
		    default:
			std::cerr << "Error: unknown option " << char(c) << "\n";
			print_usage();
		}
	}
	if (optind + 2 > argc) {	// need at least two arguements
		print_usage();
	}
}

struct hit_data {
	uint64_t end;
	hashl::data_offset_type offset;
	hashl::small_value_type frequency;
};

// add range, either as itself or extending an existing one
// (possibly merging two)

static void add_range(const std::map<hashl::data_offset_type, hashl_metadata::position> &lookup_map, const std::pair<hashl::data_offset_type, hashl::small_value_type> &x, std::vector<std::vector<std::map<uint64_t, hit_data> > > &hits, const int mer_length) {
	// convert data offset into file/read/read_start
	// (sadly, lower_bound() doesn't return <= position, but >=)
	std::map<hashl::data_offset_type, hashl_metadata::position>::const_iterator pos(lookup_map.upper_bound(x.first));
	--pos;
	// add range, or extend overlapping one
	std::map<uint64_t, hit_data> &ranges(hits[pos->second.file][pos->second.read]);
	const uint64_t start(pos->second.read_start + x.first - pos->first);
	const uint64_t end(start + mer_length);
	// see if we should extend an existing range
	if (!ranges.empty()) {
		std::map<uint64_t, hit_data>::iterator a(ranges.upper_bound(start));
		if (a != ranges.begin()) {	// check for front overlap
			std::map<uint64_t, hit_data>::iterator b(a);
			--b;
			if (b->second.end >= start) {	// overlap
				// not complete overlap - extend end
				if (b->second.end < end) {
					b->second.end = end;
					// merge with following entry?
					if (a != ranges.end() && b->second.end >= a->first) {
						// can't completely overlap
						b->second.end = a->second.end;
						if (b->second.frequency < a->second.frequency) {
							b->second.frequency = a->second.frequency;
						}
						ranges.erase(a);
					}
				}
				if (b->second.frequency < x.second) {
					b->second.frequency = x.second;
				}
				return;
			}
		} else if (a != ranges.end() && end >= a->first) {
			// back overlap
			ranges[start] = hit_data({a->second.end, x.first, x.second > a->second.frequency ? x.second : a->second.frequency});
			ranges.erase(a);
		}
	} else {
		ranges[start] = hit_data({end, x.first, x.second});
	}
}

// output format is >F#/read_name/start_end frequency(if not unique)
// (F# is 0-offset, end is exclusive)

static void print_hits(const std::vector<std::vector<std::map<uint64_t, hit_data> > > &hits, const hashl_metadata &md, const hashl &reference, std::vector<std::string> &file_list) {
	const size_t file_offset(file_list.size());
	std::string seq;
	for (size_t i(0); i < hits.size(); ++i) {
		file_list.push_back(md.file(i));
		const std::vector<std::map<uint64_t, hit_data> > &reads(hits[i]);
		for (size_t j(0); j < reads.size(); ++j) {
			std::map<uint64_t, hit_data>::const_iterator a(reads[j].begin());
			const std::map<uint64_t, hit_data>::const_iterator end_a(reads[j].end());
			for (; a != end_a; ++a) {
				std::cout << ">F" << file_offset + i << '/' << md.read(i, j) << '/' << a->first << '_' << a->second.end;
				if (a->second.frequency != 1) {
					std::cout << ' ' << a->second.frequency;
				}
				std::cout << '\n';
				reference.get_sequence(a->second.offset, a->second.end - a->first, seq);
				std::cout << seq << '\n';
			}
		}
	}
}

// go through all kmers in lookup and find positions in reference,
// then map and combine them to form a list of ranges over reads in the
// reference file(s), then print out those ranges as a fasta file

static void check_reference(const hashl &lookup, const hashl &reference, std::vector<std::string> &file_list) {
	const int mer_length(reference.bits() / 2);
	// collect positions of all matches
	// to convert data positions into file/read/range_start
	hashl_metadata md;
	md.unpack(reference.get_metadata());
	std::map<hashl::data_offset_type, hashl_metadata::position> lookup_map;
	md.create_lookup_map(lookup_map);
	// [file][read][range_start] = (range_end, if kmer is non-unique)
	std::vector<std::vector<std::map<uint64_t, hit_data> > > hits;
	hashl::const_iterator a(lookup.begin());
	const hashl::const_iterator end_a(lookup.end());
	hashl::key_type key(lookup);
	for (; a != end_a; ++a) {
		if (a.value() && a.value() != hashl::invalid_value) {
			a.get_key(key);
			// .value() returns 0 if key not found
			const std::pair<hashl::data_offset_type, hashl::small_value_type> x(reference.entry(key));
			if (x.second) {
				add_range(lookup_map, x, hits, mer_length);
			}
		}
	}
	print_hits(hits, md, reference, file_list);
}

static void print_reference_files(const std::vector<std::string> &files) {
	std::ofstream fp_out;
	if (!opt_reference_files.empty()) {
		fp_out.open(opt_reference_files);
		if (!fp_out.is_open()) {
			std::cerr << "Error: could not write to " << opt_reference_files << "\n";
			exit(1);
		}
	}
	std::ostream &io_out(opt_reference_files.empty() ? std::cerr : fp_out);
	std::vector<std::string>::const_iterator a(files.begin());
	const std::vector<std::string>::const_iterator end_a(files.end());
	for (; a != end_a; ++a) {
		io_out << *a << '\n';
	}
}

int main(const int argc, char * const * const argv) {
	get_opts(argc, argv);
	// load lookup_hash
	int fd(open_compressed(argv[optind]));
	if (fd == -1) {
		std::cerr << "Error: could not read lokup hash: " << argv[optind] << "\n";
		return 1;
	}
	hashl lookup_hash;
	lookup_hash.init_from_file(fd);
	close_compressed(fd);
	// loop through reference_hashes to generate fasta
	std::vector<std::string> file_list;
	hashl reference_hash;
	for (int i(optind + 1); i < argc; ++i) {
		fd = open_compressed(argv[optind]);
		if (fd == -1) {
			std::cerr << "Error: could not read lokup hash: " << argv[optind] << "\n";
			return 1;
		}
		reference_hash.init_from_file(fd);
		close_compressed(fd);
		check_reference(lookup_hash, reference_hash, file_list);
	}
	print_reference_files(file_list);
	return 0;
}
