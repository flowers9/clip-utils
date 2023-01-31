#include "breakup_line.h"	// breakup_line()
#include "hashp.h"	// hashp
#include "open_compressed.h"	// close_compressed(), open_compressed(), pfgets()
#include <algorithm>	// move(), sort(), swap()
#include <condition_variable>	// condition_variable
#include <exception>	// exception
#include <getopt.h>	// getopt(), optarg, optind
#include <iostream>	// cerr, cout
#include <iterator>	// back_inserter()
#include <mutex>	// mutex, lock_guard<>
#include <set>		// set<>
#include <sstream>	// istringstream, ostringstream
#include <stdio.h>	// EOF
#include <string>	// string
#include <sys/types.h>	// size_t
#include <thread>	// thread
#include <vector>	// vector<>

#ifdef PROFILE
#include <atomic>
#include <time_used.h>
#endif

// given a list of parent loci files and a read file, will find any exact
// sequence matches (comped or not) between loci sequence and read sequence
// that contain the entire locus sequence

// parent file format: id1 id2 sequence
// read file must be in fastq format
// output format: sequence loci_list read_id_list
//
// single parent locus output format is id1_id2; multi-parent locus output
// format is id1_id2_parents; loci_list may contain multiple loci separated
// by a semicolon
// 
// read_id_list is read names separated by a semicolon

// parent file sequence must contain no non-ACGT basepairs, but those are
// fine (and skipped) in the read file

class LocalException : public std::exception {
    private:
	const std::string error_;
	const int show_usage_;
    public:
	explicit LocalException(const std::string &s) : error_(s), show_usage_(0) { }
	LocalException(const std::string &s, int i) : error_(s), show_usage_(i) { }
	virtual ~LocalException(void) throw() { }
	virtual const char * what() const throw() {
		return error_.c_str();
	}
	const int &show_usage(void) const throw() {
		return show_usage_;
	}
};

static std::vector<std::string> read_names;
static std::vector<std::vector<size_t> > matches;	// [locus] = @(read_ids)

class ThreadOutput {
    public:
	ThreadOutput(void) { }
	~ThreadOutput(void) { }
	size_t add_read(const std::string &read_name) {
		read_names_.push_back(read_name);
		// creates empty matches list for new read
		matches_.resize(read_names_.size());
		return read_names_.size() - 1;
	}
	void add_match(const size_t locus, const size_t read_id) {
		matches_[read_id].insert(locus);
	}
	void move_to_global(void) {
		const size_t offset(read_names.size());
		std::move(read_names_.begin(), read_names_.end(), std::back_inserter(read_names));
		read_names_.clear();
		for (size_t i(0); i != matches_.size(); ++i) {
			std::set<size_t>::const_iterator a(matches_[i].begin());
			const std::set<size_t>::const_iterator end_a(matches_[i].end());
			for (; a != end_a; ++a) {
				matches[*a].push_back(i + offset);
			}
		}
		matches_.clear();
	}
    private:
	std::vector<std::string> read_names_;
	// for tracking read-locus matches (using set to force uniqueness)
	std::vector<std::set<size_t> > matches_;	// [read_id] = @loci
};

static hashp lookup_list;
static hashp::key_type mer_mask;
static size_t mer_length;	// min(opt_mer_length, 32)
static size_t opt_mer_length;	// passed in through options
static size_t mer_offset;	// opt_mer_length - mer_length
static char basepair_comp[256];
static hashp::key_type basepair_lookup[256];
static hashp::key_type bp_comp[4];
// reserving space for this up front (when given a size hint by the user)
// seems to use way more cycles, for some reason; so don't do that ;)
static std::vector<std::string> loci;

// all the threading globals (mostly buffers and buffer passing locks/waits)
static enum { RUNNING, FINISH_INPUT, FINISH_OUTPUT } run_state;
static std::condition_variable input_empty_wait, input_filled_wait;
static std::condition_variable output_empty_wait, output_filled_wait;
static std::mutex input_empty_mutex, input_filled_mutex;
static std::mutex output_empty_mutex, output_filled_mutex;
static std::vector<size_t> input_empty_buffers, input_filled_buffers;
static std::vector<size_t> output_empty_buffers, output_filled_buffers;
static std::vector<std::vector<std::string> > input_buffers;
static std::vector<ThreadOutput> output_buffers;

#ifdef PROFILE
static std::atomic<uint64_t> iec, ifc, oec, ofc;
static std::atomic<uint64_t> iet, ift, oet, oft;
static std::atomic<uint64_t> input_read;
#endif

static void print_usage(void) {
	std::cerr <<
		"usage: find_kmers [opts] <fastq_file> <parent_file1> [parent_file2] ...\n"
		"    (multiple parent files are indexed 0-9A-Za-z)\n"
		"    -b ## per-thread input buffer size (in lines) [4k]\n"
		"    -j ## threads [4]\n"
		"    -m ## mer length [32]\n"
	;
}

static void init_mer(void) {
	const size_t max_mer_length(sizeof(hashp::key_type) * 4);
	mer_length = opt_mer_length < max_mer_length ? opt_mer_length : max_mer_length;
	if (mer_length == max_mer_length) {
		// have to special case this one to prevent an overflow
		mer_mask = static_cast<hashp::key_type>(-1);
	} else {
		mer_mask = (static_cast<hashp::key_type>(1) << (2 * mer_length)) - 1;
	}
	mer_offset = opt_mer_length - mer_length;
	for (int i(0); i != 256; ++i) {
		basepair_lookup[i] = -1;
	}
	basepair_lookup['A'] = basepair_lookup['a'] = 0;
	basepair_lookup['C'] = basepair_lookup['c'] = 1;
	basepair_lookup['G'] = basepair_lookup['g'] = 2;
	basepair_lookup['T'] = basepair_lookup['t'] = 3;
	for (int i(0); i != 256; ++i) {
		basepair_comp[i] = 0;
	}
	basepair_comp['A'] = basepair_comp['a'] = 'T';
	basepair_comp['C'] = basepair_comp['c'] = 'G';
	basepair_comp['G'] = basepair_comp['g'] = 'C';
	basepair_comp['T'] = basepair_comp['t'] = 'A';
	bp_comp[0] = static_cast<hashp::key_type>(3) << (2 * (mer_length - 1));
	bp_comp[1] = static_cast<hashp::key_type>(2) << (2 * (mer_length - 1));
	bp_comp[2] = static_cast<hashp::key_type>(1) << (2 * (mer_length - 1));
	bp_comp[3] = 0;
}

// return the number represented by s, which may be suffixed by a k, m, or g
// which act as multipliers to the base amount

static size_t get_value(const std::string s) {
	// validate string - digits optionally followed by k, m, or g
	const size_t i(s.find_first_not_of("0123456789"));
	if (i == std::string::npos) {
		return atol(s.c_str());
	} else if (i + 1 == s.length()) {
		size_t x(atol(s.c_str()));
		switch (s[i]) {
		    case 'g':
			return x << 30;
		    case 'm':
			return x << 20;
		    case 'k':
			return x << 10;
		    default:
			return 0;
		}
	} else {	// bad value
		return 0;
	}
}

static void get_opts(int argc, char **argv, size_t &input_buffer_size, size_t &n_threads) {
	input_buffer_size = 4 * 1024;	// lines
	n_threads = 4;
	opt_mer_length = 32;
	int c;
	while ((c = getopt(argc, argv, "b:j:m:")) != EOF) {
		switch (c) {
		    case 'b':
			input_buffer_size = get_value(optarg);
			break;
		    case 'j':
			std::istringstream(optarg) >> n_threads;
			break;
		    case 'm':
			std::istringstream(optarg) >> opt_mer_length;
			break;
		    default:
			throw LocalException("bad option: " + static_cast<char>(c), 1);
		}
	}
	if ((input_buffer_size & 1) == 1) {	// needs to be an even number
		++input_buffer_size;
	}
	if (optind + 2 > argc) {
		throw LocalException("too few files specified", 1);
	}
}

static inline void generate_key(const std::string &s, hashp::key_type &key) {
	for (size_t i(0); i != mer_length; ++i) {
		key = ((key << 2) & mer_mask) | basepair_lookup[static_cast<int>(s[i])];
	}
}

static void read_parent(const char c, const char * const filename) {
	const size_t start_size(loci.size());
	std::string parent_suffix;
	if (c != 0) {
		parent_suffix.push_back('_');
		parent_suffix.push_back(c);
	}
	const int fd(open_compressed(filename));
	if (fd == -1) {
		throw LocalException("couldn't open parent file");
	}
	std::string line;
	while (pfgets(fd, line) != EOF) {
		// first two fields are parent id info
		std::vector<std::string> list;
		breakup_line(line, list);
		if (list.size() != 3) {
			throw LocalException("bad line: " + line);
		}
		if (list[2].compare("NONE") == 0) {	// skip these
			continue;
		}
		if (list[2].size() != opt_mer_length) {
			throw LocalException("sequence wrong length: " + line);
		}
		loci.push_back(list[2] + list[0] + "_" + list[1] + parent_suffix);
	}
	close_compressed(fd);
	if (start_size == loci.size()) {
		throw LocalException("empty parent file: " + std::string(filename));
	}
}

// just concatenate loci with a ;
static inline void squash_loci_one_parent(size_t i, const size_t n) {
	std::string &write_to(loci[i]);
	for (++i; i != n; ++i) {
		write_to += ";" + loci[i].substr(opt_mer_length);
	}
}

// take loci with matching parents and reduce them to one with a longer
// parent prefix; concatenate different loci with a ;

static void squash_loci_multi_parent(size_t i, const size_t n) {
	std::string &write_to(loci[i]);
	// have to be careful writing out first set, as we're writing into
	// the first entry
	int not_first(0);
	for (;;) {
		const size_t start(i);
		const size_t length(loci[i].size());
		// -2 to avoid the parent suffix (_#)
		const size_t locus_match_length(length - opt_mer_length - 2);
		std::string parent_suffix(1, loci[i][length - 1]);
		// check lengths because if loci are equal, size is equal,
		// and we can short-circuit a compare that way
		for (++i; i != n && loci[i].size() == length && loci[start].compare(opt_mer_length, locus_match_length, loci[i], opt_mer_length, locus_match_length) == 0; ++i) {
			parent_suffix.push_back(loci[i][length - 1]);
		}
		if (not_first) {	// simply append
			write_to += ";" + loci[start].substr(opt_mer_length, locus_match_length) + "_" + parent_suffix;
		} else {		// only rewrite parent index
			not_first = 1;
			// replace because parent_suffix includes write_to's
			write_to.replace(length - 1, 1, parent_suffix);
		}
		if (i == n) {
			break;
		}
	}
}

static void order_loci(const int multi_parent) {
	std::sort(loci.begin(), loci.end());
	// now squash duplicate sequences - put all loci with the same index
	// into one entry (and further squash loci matching in all but parent
	// into one smaller locus entry)
	size_t read_index(0), write_index(0);
	const size_t end_n(loci.size());
	const size_t end_reads(end_n - 1);
	for (;;) {
		// find next duplicate
		for (; read_index != end_reads && loci[read_index].compare(0, opt_mer_length, loci[read_index + 1], 0, opt_mer_length) != 0; ++read_index, ++write_index) {
			if (write_index != read_index) {
				// swap should be faster than assignment
				// (can't use move, as we might use cell later
				// and move would leave it in a bad state)
				std::swap(loci[write_index], loci[read_index]);
			}
		}
		if (read_index == end_reads) {
			break;
		}
		// copy parent info of duplicate(s)
		size_t n(read_index + 2);
		for (; n != end_n && loci[read_index].compare(0, opt_mer_length, loci[n], 0, opt_mer_length) == 0; ++n) { }
		if (multi_parent) {
			squash_loci_multi_parent(read_index, n);
		} else {
			squash_loci_one_parent(read_index, n);
		}
		if (n == end_n) {
			break;
		}
		if (write_index != read_index) {
			std::swap(loci[write_index], loci[read_index]);
		}
		read_index = n;
		++write_index;
	}
	// last one moved outside loop to make loop easier to read
	if (write_index != read_index) {
		std::swap(loci[write_index], loci[read_index]);
	}
	loci.resize(write_index + 1);
}

static void hash_loci(void) {
	hashp::key_type key(0);		// make sure to initialize this
	size_t i(0);
	const size_t end_i(loci.size());
	while (i != end_i) {
		const std::string x(loci[i].substr(0, mer_length));
		generate_key(x, key);
		const size_t start(i);
		// find all loci with matching start sequence
		do {
			loci[i].erase(0, mer_length);
			++i;
		} while (i != end_i && loci[i].compare(0, mer_length, x) == 0);
		lookup_list.add(key, start, i);
	}
}

// find first span of at least n valid basepairs, starting at i

static int find_span(const std::string &seq, size_t &i, const size_t n) {
	size_t end_j(i + n);
	if (end_j > seq.size()) {
		return 0;
	}
	for (size_t j(i); j != end_j; ++j) {
		if (basepair_lookup[static_cast<int>(seq[j])] == static_cast<hashp::key_type>(-1)) {
			// start over
			i = j + 1;
			end_j = i + n;
			if (end_j > seq.size()) {
				return 0;
			}
		}
	}
	return 1;
}

// given the sequence, create the key for the first mer length - 1 proper
// (i.e., ACGT) base pairs, returning the current position in the sequence;
// the incoming key is presumed to have no set bits greater than the mer_mask

static int preload_key(const std::string &seq, size_t &i, hashp::key_type &key, hashp::key_type &comp_key) {
	if (!find_span(seq, i, opt_mer_length)) {
		return 0;
	}
	const size_t end_i(i + mer_length - 1);
	for (size_t j(i + mer_offset); i != end_i; ++i, ++j) {
		key = ((key << 2) & mer_mask) | basepair_lookup[static_cast<int>(seq[i])];
		comp_key = (comp_key >> 2) | bp_comp[basepair_lookup[static_cast<int>(seq[j])]];
	}
	return 1;
}

// note: always hold lock through notify() to prevent possible races

static inline void mark_input_buffer_empty(const size_t i) {
	std::lock_guard<std::mutex> lock(input_empty_mutex);
	input_empty_buffers.push_back(i);
	input_empty_wait.notify_one();
}

static inline void mark_input_buffer_filled(const size_t i) {
	std::lock_guard<std::mutex> lock(input_filled_mutex);
	input_filled_buffers.push_back(i);
	input_filled_wait.notify_one();
}

static inline size_t get_empty_input_buffer(void) {
	std::unique_lock<std::mutex> lock(input_empty_mutex);
#ifdef PROFILE
	++iec;
	iet += input_empty_buffers.size();
#endif
	while (input_empty_buffers.empty()) {
		input_empty_wait.wait(lock);
	}
	const size_t i(input_empty_buffers.back());
	input_empty_buffers.pop_back();
	return i;
}

static inline size_t get_filled_input_buffer(void) {
	std::unique_lock<std::mutex> lock(input_filled_mutex);
#ifdef PROFILE
	++ifc;
	ift += input_filled_buffers.size();
#endif
	while (input_filled_buffers.empty() && run_state != FINISH_INPUT) {
		input_filled_wait.wait(lock);
	}
	if (input_filled_buffers.empty()) {
		return static_cast<size_t>(-1);
	}
	const size_t i(input_filled_buffers.back());
	input_filled_buffers.pop_back();
	return i;
}

static inline void mark_output_buffer_empty(const size_t i) {
	std::lock_guard<std::mutex> lock(output_empty_mutex);
	output_empty_buffers.push_back(i);
	output_empty_wait.notify_one();
}

static inline void mark_output_buffer_filled(const size_t i) {
	std::lock_guard<std::mutex> lock(output_filled_mutex);
	output_filled_buffers.push_back(i);
	output_filled_wait.notify_one();
}

static inline size_t get_empty_output_buffer(void) {
	std::unique_lock<std::mutex> lock(output_empty_mutex);
#ifdef PROFILE
	++oec;
	oet += output_empty_buffers.size();
#endif
	while (output_empty_buffers.empty()) {
		output_empty_wait.wait(lock);
	}
	const size_t i(output_empty_buffers.back());
	output_empty_buffers.pop_back();
	return i;
}

static inline size_t get_filled_output_buffer(void) {
	std::unique_lock<std::mutex> lock(output_filled_mutex);
#ifdef PROFILE
	++ofc;
	oft += output_filled_buffers.size();
#endif
	while (output_filled_buffers.empty() && run_state != FINISH_OUTPUT) {
		output_filled_wait.wait(lock);
	}
	if (output_filled_buffers.empty()) {
		return static_cast<size_t>(-1);
	}
	const size_t i(output_filled_buffers.back());
	output_filled_buffers.pop_back();
	return i;
}

static std::string reverse_complement(const std::string &s) {
	std::string t;
	t.reserve(s.size());
	for (size_t i(s.size() - 1); i != std::string::npos; --i) {
		t.push_back(basepair_comp[static_cast<int>(s[i])]);
	}
	return t;
}

static void find_matches(ThreadOutput &buffer, const std::string &name, size_t &read_id, const size_t i, const hashp::key_type key, const std::string &seq) {
	hashp::value_type j, k;
	if (lookup_list.has_key(key, j, k)) {
		// check remainder against each parent starting
		// with the same mer_length basepairs
		for (; j != k; ++j) {
			if (seq.compare(i, mer_offset, loci[j], 0, mer_offset) == 0) {
				if (read_id == static_cast<size_t>(-1)) {
					read_id = buffer.add_read(name.substr(1));
				}
				buffer.add_match(j, read_id);
			}
		}
	}
}

static void grep_seq(ThreadOutput &buffer, const std::string &name, const std::string &seq) {
	size_t i(0);
	hashp::key_type key(0), comp_key(0);
	if (!preload_key(seq, i, key, comp_key)) {
		return;
	}
	size_t read_id(-1);
	const size_t rc_base(seq.size() + mer_length - mer_offset);
	// if we skipped any starting sequence, we don't need it here
	const std::string rc(reverse_complement(seq.substr(i + 1 - mer_length)));
	for (;;) {
		key = ((key << 2) & mer_mask) | basepair_lookup[static_cast<int>(seq[i])];
		comp_key = (comp_key >> 2) | bp_comp[basepair_lookup[static_cast<int>(seq[i + mer_offset])]];
		++i;
		find_matches(buffer, name, read_id, i, key, seq);
		find_matches(buffer, name, read_id, rc_base - i, comp_key, rc);
		if (i + mer_offset == seq.size()) {
			return;
		}
		if (basepair_lookup[static_cast<int>(seq[i + mer_offset])] == static_cast<hashp::key_type>(-1)) {
			// hit a bad spot, move past
			i += mer_offset + 1;
			if (!preload_key(seq, i, key, comp_key)) {
				return;
			}
		}
	}
}

static void process_input_buffer(void) {
	for (;;) {
		const size_t input_buffer_id(get_filled_input_buffer());
		if (input_buffer_id == static_cast<size_t>(-1)) {
			break;
		}
		const size_t output_buffer_id(get_empty_output_buffer());
		const std::vector<std::string> &input_buffer(input_buffers[input_buffer_id]);
		ThreadOutput &output_buffer(output_buffers[output_buffer_id]);
		for (size_t i(0); i != input_buffer.size(); i += 2) {
			grep_seq(output_buffer, input_buffer[i], input_buffer[i + 1]);
		}
		mark_input_buffer_empty(input_buffer_id);
		mark_output_buffer_filled(output_buffer_id);
	}
}

// moving the output buffer from thread temporary to global permanent
// storage (rolling this into the processing threads uses a touch less
// memory (the buffers), a touch fewer cycles, and a lot more time)

static void store_thread_buffer(void) {
	for (;;) {
		const size_t i(get_filled_output_buffer());
		if (i == static_cast<size_t>(-1)) {
			return;
		}
		output_buffers[i].move_to_global();
		mark_output_buffer_empty(i);
	}
}

// returns 0 if EOF encountered

static void fill_input_buffer(int fd, std::vector<std::string> &buffer) {
	size_t i(0);
	std::string line;
	for (;;) {
		std::string &name(buffer[i]);
		if (pfgets(fd, name) == EOF) {
			buffer.resize(i);
			return;
		}
#ifdef PROFILE
		input_read += name.size();
#endif
		if (name[0] != '@') {
			throw LocalException("bad read line" + line);
		}
		size_t j(name.find(' ', 1));
		if (j != std::string::npos) {
			name.resize(j);
		}
		std::string &seq(buffer[++i]);
		if (pfgets(fd, seq) == EOF) {
			throw LocalException("premature end of read file");
		}
#ifdef PROFILE
		input_read += seq.size();
#endif
		// switching to skip_next_line() didn't actually help
		if (pfgets(fd, line) == EOF) {
			throw LocalException("premature end of read file");
		} else if (line.compare("+") != 0) {
			throw LocalException("bad quality read name line");
		}
#ifdef PROFILE
		input_read += line.size();
#endif
		if (pfgets(fd, line) == EOF) {
			throw LocalException("premature end of read file");
		}
#ifdef PROFILE
		input_read += line.size();
#endif
		if (++i == buffer.size()) {
			return;
		}
	}
}

// converts key to sequence

static std::string convert_key(hashp::key_type key) {
	static const char values[4] = { 'A', 'C', 'G', 'T' };
	std::string sequence(mer_length, 0);
	for (size_t i(mer_length - 1); i != static_cast<size_t>(-1); --i, key >>= 2) {
		sequence[i] = values[key & 3];
	}
	return sequence;
}

static void print_output(void) {
	// iterate over lookup hash so we can reverse engineer the
	// first mer_length basepairs of the sequence from the hash key
	hashp::const_iterator a(lookup_list.begin());
	const hashp::const_iterator end_a(lookup_list.end());
	for (; a != end_a; ++a) {
		const std::string prefix(convert_key(a.key));
		hashp::value_type i(a.v1_out);
		const hashp::value_type end_i(a.v2_out);
		for (; i != end_i; ++i) {
			// no point printing out non-matching loci
			if (matches[i].empty()) {
				continue;
			}
			// print sequence, related loci, read matches
			std::cout << prefix << loci[i].substr(0, mer_offset) << "\t" << loci[i].substr(mer_offset) << "\t";
			std::vector<size_t>::const_iterator b(matches[i].begin());
			const std::vector<size_t>::const_iterator end_b(matches[i].end());
			std::cout << read_names[*b];
			for (++b; b != end_b; ++b) {
				std::cout << ";" << read_names[*b];
			}
			std::cout << "\n";
		}
	}
}

// read fastq file, checking to see if any kmers are in lookup_list

static void grep_file(const char * const filename, const size_t input_buffer_size, const size_t n_threads) {
	run_state = RUNNING;
	std::thread threads[n_threads];
	input_buffers.resize(n_threads + 1);
	output_buffers.resize(n_threads + 1);
	input_empty_buffers.reserve(n_threads + 1);
	input_filled_buffers.reserve(n_threads + 1);
	output_empty_buffers.reserve(n_threads + 1);
	output_filled_buffers.reserve(n_threads + 1);
	for (size_t i(0); i != n_threads + 1; ++i) {
		input_buffers[i].resize(input_buffer_size);
		input_empty_buffers.push_back(i);
		output_empty_buffers.push_back(i);
	}
#ifdef PROFILE
	iet = iec = ift = ifc = oet = oec = oft = ofc = 0;
	input_read = 0;
	start_time();
#endif
	std::thread storage_thread(store_thread_buffer);
	for (size_t i(0); i != n_threads; ++i) {
		threads[i] = std::thread(process_input_buffer);
	}
	const int fd(open_compressed(filename));
	if (fd == -1) {
		throw LocalException("couldn't open read file");
	}
	// cycle buffers
	for (;;) {
#ifdef PROFILE
		const double x(elapsed_time());
		if (x >= 10) {
			std::cerr <<
				(double)iet / iec << " " <<
				(double)ift / ifc << " " <<
				(double)oet / oec << " " <<
				(double)oft / ofc << " " <<
				input_read / 1048576 / x << " (";
			for (size_t i(0); i < output_filled_buffers.size(); ++i) {
				std::cerr << " " << output_filled_buffers[i];
			}
			std::cerr << ")\n";
			iet = iec = ift = ifc = oet = oec = oft = ofc = 0;
			input_read = 0;
			start_time();
		}
#endif
		const size_t i(get_empty_input_buffer());
		fill_input_buffer(fd, input_buffers[i]);
		if (input_buffers[i].empty()) {
			break;
		}
		mark_input_buffer_filled(i);
	}
	close_compressed(fd);
	{	// lock to prevent very unlikely race condition
		std::lock_guard<std::mutex> lock(input_filled_mutex);
		run_state = FINISH_INPUT; // finish processing input buffers
	}
	input_filled_wait.notify_all();
	for (size_t i(0); i != n_threads; ++i) {
		threads[i].join();
	}
	{	// lock to prevent very unlikely race condition
		std::lock_guard<std::mutex> lock(output_filled_mutex);
		run_state = FINISH_OUTPUT; // finish printing output buffers
	}
	output_filled_wait.notify_one();
	storage_thread.join();
	print_output();
}

int main(int argc, char **argv) {
	int had_error(0);
	try {
		size_t input_buffer_size, n_threads;
		get_opts(argc, argv, input_buffer_size, n_threads);
		init_mer();
		const char * const fastq_file(argv[optind]);
		++optind;
		int multi_parent(optind + 1 != argc);
		if (multi_parent) {
			const std::string parent_index("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
			if (static_cast<size_t>(argc - optind) > parent_index.size()) {
				throw LocalException("Too many parent files");
			}
			for (size_t i(0); optind != argc; ++optind, ++i) {
				read_parent(parent_index[i], argv[optind]);
			}
		} else {
			read_parent(0, argv[optind]);
		}
		order_loci(multi_parent);
		matches.resize(loci.size());
		// *2 for hash lookup efficiency - it seems a good
		// tradeoff of memory versus speed
		lookup_list.init(loci.size() * 2);
		hash_loci();
		// read read file, print matches
		grep_file(fastq_file, input_buffer_size, n_threads);
	} catch (std::exception &e) {
		std::cerr << "Error: " << e.what() << "\n";
		LocalException *x(dynamic_cast<LocalException *>(&e));
		if (x != NULL && x->show_usage()) {
			print_usage();
		}
		had_error = 1;
	}
	return had_error;
}
