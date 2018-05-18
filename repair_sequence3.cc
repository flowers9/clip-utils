#include "breakup_line.h"	/* breakup_line() */
#include "itoa.h"	/* itoa() */
#include "open_compressed.h"	/* close_compressed(), open_compressed(), pfgets() */
#include "strtostr.h"	/* strtostr() */
#include "version.h"	/* VERSION */
#include "write_fork.h"	/* close_fork(), pfputs(), write_fork() */
#include <algorithm>	/* max(), min() */
#include <ctype.h>	/* toupper() */
#include <dirent.h>	/* DIR, closedir(), opendir(), readdir() */
#include <errno.h>	/* errno */
#include <getopt.h>	// getopt(), optarg, optind
#include <list>		/* list<> */
#include <map>		/* map<> */
#include <math.h>	/* floor(), sqrt() */
#include <stdio.h>	/* EOF, fprintf(), stderr */
#include <stdlib.h>	// atoi(), exit()
#include <string.h>	/* strerror() */
#include <string>	/* string */
#include <sys/stat.h>	/* stat(), struct stat */
#include <unistd.h>	/* unlink() */
#include <vector>	/* vector<> */

#define FASTA_WIDTH 50

static bool opt_cleanup, opt_pair_match, opt_indel, opt_report;
static double opt_overlap_cutoff;
static int fd_log, opt_abort, opt_change, opt_confirm, opt_length_cutoff;
static std::map<std::string, bool> opt_exclude;
static std::string opt_strip_trace;

static std::string find_qual_file(const std::string &, bool = 0);

class Match {
    public:
	long length, score, start, stop, s_start, s_stop;
	std::string id, qs, ss, tag;
	explicit Match(void) { }
	Match(const Match &__a) : length(__a.length), score(__a.score), start(__a.start), stop(__a.stop), s_start(__a.s_start), s_stop(__a.s_stop), id(__a.id), qs(__a.qs), ss(__a.ss), tag(__a.tag) { }
	~Match(void) { }
	Match operator=(const Match &__a) {
		length = __a.length;
		score = __a.score;
		start = __a.start;
		stop = __a.stop;
		s_start = __a.s_start;
		s_stop = __a.s_stop;
		id = __a.id;
		qs = __a.qs;
		ss = __a.ss;
		tag = __a.tag;
		return *this;
	}
	bool operator<(const Match &__a) const {
		if (start != __a.start) {
			return start < __a.start;
		} else if (stop != __a.stop) {
			return stop < __a.stop;
		} else {
			return id < __a.id;
		}
	}
	bool trim_alignment(void);
	bool identity_check(long, long) const;
	void make_tag(void);
	void slide(void);
	void scrub(size_t);
	bool post_scrub_trim(void);
	void add_pads(const std::vector<int> &, const std::vector<int> &);
	void update_tag(const std::map<long, char> &);
	void print(void) const {
		fprintf(stderr, ">%s %ld\n%ld %ld %ld %ld %ld\n", id.c_str(), length, score, start, stop, s_start, s_stop);
		fprintf(stderr, "%s\n", qs.c_str());
		//fprintf(stderr, "%s\n", tag.c_str());
		fprintf(stderr, "%s\n", ss.c_str());
	}
};

class CompareMatchScore {
    public:
	bool operator()(const Match &__a, const Match &__b) {
		if (__a.score != __b.score) {
			return __a.score < __b.score;
		} else if (__a.stop - __a.start != __b.stop - __b.start) {
			return __a.stop - __a.start < __b.stop - __b.start;
		} else {
			return __a < __b;
		}
	}
};

class BlatFile {
    private:
	bool _is_current;
	int _fd;
	std::string _filename, _current_id, _line;
    public:
	explicit BlatFile(void) : _fd(-1) { }
	explicit BlatFile(const std::string &__s) {
		open_blat(__s);
	}
	~BlatFile(void) {
		close_blat();
	}
	void open_blat(const std::string &__s) {
		_is_current = 0;
		if ((_fd = open_compressed(_filename = __s)) == -1) {
			fprintf(stderr, "Error: open_compressed: %s: %s\n", __s.c_str(), strerror(errno));
			exit(1);
		}
	}
	void close_blat(void) {
		_is_current = 0;
		if (_fd != -1) {
			close_compressed(_fd);
			_fd = -1;
		}
	}
	void unlink_blat(void) {
		close_blat();
		if (!_filename.empty()) {
			unlink(_filename.c_str());
			_filename.clear();
		}
	}
	bool is_current(void) const {
		return _is_current;
	}
	bool find_next(const std::string &);
	bool read_next(size_t, Match &);
	void finish_current(void) {
		_current_id.clear();
	}
};

class FastaFile {
    private:
	bool _unlink;
	int _fd_fasta, _fd_qual;
	std::string _fasta_file, _qual_file;
	std::string _header, _id, _seq, _line, _line_qual, _seq_orig;
	std::vector<int> _qual;
	void read_next_qual(void);
    public:
	explicit FastaFile(void) : _fd_fasta(-1), _fd_qual(-1) { }
	explicit FastaFile(const std::string &__s, bool __b, const std::string &__t) : _fd_qual(-1), _fasta_file(__s) {
		if ((_fd_fasta = open_compressed(_fasta_file)) == -1) {
			fprintf(stderr, "Error: open_compressed: %s: %s\n", _fasta_file.c_str(), strerror(errno));
			exit(1);
		}
		if (__b) {
			_qual_file = find_qual_file(_fasta_file);
			if (!_qual_file.empty()) {
				if ((_fd_qual = open_compressed(_qual_file)) == -1) {
					fprintf(stderr, "Error: open_compressed: %s: %s\n", _qual_file.c_str(), strerror(errno));
					exit(1);
				}
			}
		}
		_unlink = _fasta_file.find(__t) != std::string::npos;
	}
	~FastaFile(void) {
		fasta_close();
	}
	const std::string &fasta_file(void) const {
		return _fasta_file;
	}
	const std::string &header(void) const {
		return _header;
	}
	const std::string &id(void) const {
		return _id;
	}
	const std::string &seq(void) const {
		return _seq;
	}
	bool is_open(void) const {
		return _fd_fasta != -1;
	}
	bool has_qual(void) const {
		return _fd_qual != -1;
	}
	void fasta_close(void) {
		if (_fd_fasta != -1) {
			close_compressed(_fd_fasta);
			_fd_fasta = -1;
		}
		if (_fd_qual != -1) {
			close_compressed(_fd_qual);
			_fd_qual = -1;
		}
	}
	void fasta_unlink(void) {
		fasta_close();
		if (_unlink) {
			if (!_fasta_file.empty()) {
				unlink(_fasta_file.c_str());
				_fasta_file.clear();
			}
			if (!_qual_file.empty()) {
				unlink(_qual_file.c_str());
				_qual_file.clear();
			}
		}
	}
	void add_trace(const std::string &__s) {
		_header.insert(1, __s);
	}
	bool read_next(void);
	void add_pads(const std::vector<int> &, size_t);
	void print(int, int = -1) const;
	void update_seq(std::string::size_type __i, char __c) {
		_seq[__i] = toupper(__c);
	}
	void update_quality(const std::string &);
	void remove_pads(void);
	void revert(void) {
		_seq = _seq_orig;
		_seq_orig.clear();
	}
};

/*
 * trim homopolymer from ends (longest of query and search sequences);
 * discard matches with mismatches in the starting or ending homopolymers,
 * as I haven't found a good way to assure the rest of the match is acceptable
 * (^aaaggg vs ^aagggg, when it should really be aaggg- vs aagggg, for example)
 */

// returns if match should be used
bool Match::trim_alignment() {
	// find lengths of starting and ending homopolymers
	std::string::size_type i = qs.find_first_not_of(qs[0], 1);
	if (i == std::string::npos || i != ss.find_first_not_of(ss[0], 1)) {
		return 0;
	}
	std::string::size_type j = qs.find_last_not_of(*qs.rbegin());
	if (j != ss.find_last_not_of(*ss.rbegin())) {
		return 0;
	}
	++j;
	start += i;
	stop -= qs.size() - j;
	qs.resize(j);
	qs.erase(0, i);
	ss.resize(j);
	ss.erase(0, i);
	tag.resize(j);
	tag.erase(0, i);
	return 1;
}

// return if potential match has good identity
bool Match::identity_check(long identity, long seq_length) const {
	long c, d;
	long b = seq_length - stop;
	if (s_start < s_stop) {
		c = s_start;
		d = length - s_stop;
	} else {
		c = length - s_start;
		d = s_stop;
	}
	// check for going off ends, and don't count that in the length
	long x = length;	// potential match length
	if (start < c) {
		x -= c - start;
	}
	if (b < d) {
		x -= d - b;
	}
	return identity >= opt_overlap_cutoff * x;
}

void Match::make_tag() {
	tag.assign(qs.size(), '|');
	if (ss.empty()) {	// ss isn't given if it's identical to qs
		ss = qs;
		return;
	}
	// mismatches to non-standard basepairs are not counted as mismatches
	const std::string good("-acgt");
	std::string::size_type i;
	for (i = 0; i != qs.size(); ++i) {
		if (qs[i] != ss[i] && good.find(ss[i]) != std::string::npos) {
			tag[i] = ' ';
		}
	}
}

// replace existing contents of string from pos to pos + n - 1 with c
static inline void replace(std::string &s, std::string::size_type pos, std::string::size_type n, char c) {
	for (n += pos; pos != n; ++pos) {
		s[pos] = c;
	}
}

// rearrange string from i to j so that c comes first, - after;
// return number of c's
static size_t shift(std::string &s, std::string::size_type s_begin, std::string::size_type s_end, char c) {
	size_t n = 0;
	std::string::size_type i;
	for (i = s_begin; i != s_end; ++i) {
		if (s[i] == c) {
			s[i] = '-';
			++n;
		}
	}
	replace(s, s_begin, n, c);
	return n;
}

// for any homopolymer indels, shift them to the end of the repeat
// (blat seems to randomly place the indel at the beginning or end)
void Match::slide() {
	// check over insertions
	std::string::size_type i = static_cast<std::string::size_type>(-1);
	while ((i = qs.find('-', i + 1)) != std::string::npos) {
		char c = ss[i];
		std::string s = "-" + std::string(1, c);
		std::string::size_type j = qs.find_first_not_of(s, i + 1);
		if (j == std::string::npos) {
			j = qs.size();
		}
		std::string::size_type k = ss.find_first_not_of(s, i + 1);
		if (k == std::string::npos) {
			k = ss.size();
		}
		if (j > k) {
			j = k;
		}
		if (j == i + 1) {
			continue;
		}
		size_t x = shift(qs, i, j, c);
		size_t y = shift(ss, i, j, c);
		if (x > y) {
			k = x - y;
			x = y;
		} else {
			k = y - x;
		}
		i = j - 1;
	}
	// check over deletions
	i = static_cast<std::string::size_type>(-1);
	while ((i = ss.find('-', i + 1)) != std::string::npos) {
		char c = qs[i];
		std::string s = "-" + std::string(1, c);
		std::string::size_type j = qs.find_first_not_of(s, i + 1);
		if (j == std::string::npos) {
			j = qs.size();
		}
		std::string::size_type k = ss.find_first_not_of(s, i + 1);
		if (k == std::string::npos) {
			k = ss.size();
		}
		if (j > k) {
			j = k;
		}
		if (j == i + 1) {
			continue;
		}
		size_t x = shift(qs, i, j, c);
		size_t y = shift(ss, i, j, c);
		if (x > y) {
			k = x - y;
			x = y;
		} else {
			k = y - x;
		}
		i = j - 1;
	}
}

// turn homopolymers surrounding k into n's
void Match::scrub(size_t k) {
	// find end of adjacent leading homopolymer
	std::string::size_type i, j;
	if (k == 0) {
		i = 0;
	} else {
		i = qs.find_last_not_of(std::string(1, '-') + qs[k], k - 1);
		if (i == std::string::npos) {
			i = 0;
		} else if (i != 0) {
			i = qs.find_last_not_of(std::string(1, '-') + qs[i], i - 1);
			if (i == std::string::npos) {
				i = 0;
			} else {
				std::string::size_type m = ss.find_last_not_of(std::string(1, '-') + ss[k], k - 1);
				if (m == std::string::npos || m == 0) {
					i = 0;
				} else {
					m = ss.find_last_not_of(std::string(1, '-') + ss[m], m - 1);
					if (m == std::string::npos) {
						i = 0;
					} else {
						if (i > m) {
							i = m;
						}
						++i;
					}
				}
			}
		}
	}
	// find end of adjacent trailing homopolymer
	std::string::size_type n = qs.size() - 1;
	if (k == n) {
		j = n;
	} else {
		j = qs.find_first_not_of(std::string(1, '-') + qs[k], k + 1);
		if (j == std::string::npos) {
			j = n;
		} else if (j != n) {
			j = qs.find_first_not_of(std::string(1, '-') + qs[j], j + 1);
			if (j == std::string::npos) {
				j = n;
			} else {
				std::string::size_type m = ss.find_first_not_of(std::string(1, '-') + ss[k], k + 1);
				if (m == std::string::npos || m == n) {
					j = n;
				} else {
					m = ss.find_first_not_of(std::string(1, '-') + ss[m], m + 1);
					if (m == std::string::npos) {
						j = n;
					} else {
						if (j < m) {
							j = m;
						}
						--j;
					}
				}
			}
		}
	}
	j -= i - 1;
	replace(ss, i, j, 'n');
	replace(tag, i, j, '|');
}

// remove any leading/trailing n's to speed up later routines
// returns if entire match is now n's
bool Match::post_scrub_trim() {
	// remove trailing n's (do first, so string is smaller)
	std::string::size_type i = ss.find_last_not_of('n');
	if (i == std::string::npos) {		// whoops, all n's
		return 1;
	}
	if (i != ss.size() - 1) {		// trailing n's found
		++i;
		stop -= ss.size() - i;
		qs.erase(i);
		ss.erase(i);
		tag.erase(i);
	}
	// remove leading n's
	i = ss.find_first_not_of('n');
	if (i != 0) {				// leading n's found
		start += i;
		qs.erase(0, i);
		ss.erase(0, i);
		tag.erase(0, i);
	}
	return 0;
}

// only insert pads that aren't already there
void Match::add_pads(const std::vector<int> &pads, const std::vector<int> &all_pads) {
	std::string::size_type j = qs.find_first_not_of('-');
	if (j == std::string::npos) {
		return;
	}
	long i;
	for (i = start + 1, ++j; i != stop; ++i) {
		int k = pads[i];
		if (k != 0) {
			std::string::size_type m = qs.find_first_not_of('-', j);
			if (m == std::string::npos) {
				break;
			}
			if (j + k > m) {
				int n = j + k - m;
				qs.insert(j, n, '-');
				ss.insert(j, n, '-');
				tag.insert(j, n, '|');
				j = m + 1 + n;
			} else {
				j = m + 1;
			}
		} else {
			if (qs[j] == '-') {	// that one didn't count
				--i;
			}
			++j;
		}
	}
	start += all_pads[start];
	stop += all_pads[stop - 1];
}

// change tag: | -> ., ' ' -> X, or Y for matching changes
void Match::update_tag(const std::map<long, char> &changes) {
	std::string::size_type i;
	for (i = 0; i != tag.size(); ++i) {
		if (tag[i] == '|') {
			tag[i] = '.';
		} else {
			std::map<long, char>::const_iterator a = changes.find(start + i);
			if (a != changes.end() && a->second == ss[i]) {
				tag[i] = 'Y';
			} else {
				tag[i] = 'X';
			}
		}
	}
}

bool BlatFile::find_next(const std::string &id) {
	if (_fd == -1) {
		return 0;
	}
	if (_current_id.empty()) {
		if (_line.empty()) {
			if (pfgets(_fd, _line) == -1) {
				close_blat();
				return 0;
			}
		}
		while (_line[0] != '=' && pfgets(_fd, _line) != -1) { }
		if (_line.empty()) {
			close_blat();
			return 0;
		}
		_current_id = _line.substr(1);
		pfgets(_fd, _line);		// preload next line
	}
	return _is_current = id.size() >= _current_id.size() && id.compare(0, _current_id.size(), _current_id) == 0 && (id.size() == _current_id.size() || id[_current_id.size()] == ' ');
}

// return if a good current match was found;
bool BlatFile::read_next(size_t seq_length, Match &b) {
	for (;;) {		// loop over all current matches
		std::string::size_type i;
		for (;;) {	// loop until next good match
			if (_line[0] == '>') {
				b.id = strtostr(_line, &(i = 1));
				b.length = atoi(strtostr(_line, &i).c_str());
				if (opt_exclude.find(b.id) == opt_exclude.end() && b.length >= opt_length_cutoff) {
					break;
				}
			} else if (_line[0] == '=') {
				return 0;
			}
			if (pfgets(_fd, _line) == -1) {
				return 0;
			}
		}
		pfgets(_fd, _line);		// read info line
		int identity = atoi(strtostr(_line, &(i = 0)).c_str());
		b.score = atoi(strtostr(_line, &i).c_str());
		b.start = atoi(strtostr(_line, &i).c_str());
		b.stop = atoi(strtostr(_line, &i).c_str());
		b.s_start = atoi(strtostr(_line, &i).c_str());
		b.s_stop = atoi(strtostr(_line, &i).c_str());
		pfgets(_fd, _line);
		if (b.identity_check(identity, seq_length)) {
			b.qs = strtostr(_line, &(i = 0));
			b.ss = strtostr(_line, &i);
			b.slide();
			b.make_tag();
			--b.start;		// change to zero offset
			pfgets(_fd, _line);	// preload next line
			return 1;
		}
		if (pfgets(_fd, _line) == -1) {
			return 0;
		}
	}
}

void FastaFile::read_next_qual() {
	if (_line_qual.empty()) {
		pfgets(_fd_qual, _line_qual);
	}
	if (_line_qual != _header) {
		fprintf(stderr, "Error: read name mismatch between read and qual file: %s != %s\n", _header.c_str(), _line_qual.c_str());
		exit(1);
	}
	if (pfgets(_fd_qual, _line_qual) == -1) {
		return;
	}
	while (_line_qual[0] != '>') {
		std::vector<std::string> list;
		breakup_line(_line_qual, list);
		std::vector<std::string>::const_iterator a = list.begin();
		std::vector<std::string>::const_iterator end_a = list.end();
		for (; a != end_a; ++a) {
			if (a->compare("") != 0) {
				_qual.push_back(atoi(a->c_str()));
			}
		}
		if (pfgets(_fd_qual, _line_qual) == -1) {
			break;
		}
	}
	if (_qual.size() != 0 && _qual.size() != _seq.size()) {
		fprintf(stderr, "Warning: length mismatch between sequence and qual: %lu != %lu: %s\n", _seq.size(), _qual.size(), _id.c_str());
	}
}

bool FastaFile::read_next() {
	if (_line.empty()) {
		pfgets(_fd_fasta, _header);
	} else {
		_header = _line;
		_line.clear();
	}
	_id.clear();
	_seq.clear();
	_qual.clear();
	if (!_header.empty()) {
		std::string::size_type i;
		_id = strtostr(_header, &(i = 1));	// skip leading >
		while (pfgets(_fd_fasta, _line) != -1 && _line[0] != '>') {
			_seq += _line;
		}
		if (_fd_qual != -1) {
			read_next_qual();
		}
		_header += "\n";
	}
	return !_header.empty();
}

// pad sequence and quality
void FastaFile::add_pads(const std::vector<int> &pads, size_t total_pads) {
	_seq.resize(_seq.size() + total_pads);
	std::string::size_type i, j;
	for (i = pads.size() - 1, j = _seq.size() - 1; i != j; --i, --j) {
		_seq[j] = _seq[i];
		if (pads[i] != 0) {
			replace(_seq, j -= pads[i], pads[i], '-');
		}
	}
	if (opt_abort != -1) {	// make backup in case we have to back out
		_seq_orig = _seq;
	}
	if (!_qual.empty()) {
		_qual.resize(_qual.size() + total_pads);
		for (i = pads.size() - 1, j = _qual.size() - 1; i != j; --i) {
			_qual[j--] = _qual[i];
			if (pads[i] != 0) {
				size_t k = j - pads[i];
				for (; j != k; --j) {
					_qual[j] = 0;
				}
			}
		}
	}
}

void FastaFile::print(int fd_fasta_out, int fd_qual_out) const {
	if (pfputs(fd_fasta_out, _header) == -1) {
		exit(1);
	}
	size_t i;
	for (i = 0; i + FASTA_WIDTH < _seq.size(); i += FASTA_WIDTH) {
		if (pfputs(fd_fasta_out, _seq.substr(i, FASTA_WIDTH) + "\n") == -1) {
			exit(1);
		}
	}
	if (!_seq.empty() && pfputs(fd_fasta_out, _seq.substr(i) + "\n") == -1) {
		exit(1);
	}
	if (fd_qual_out != -1) {
		if (pfputs(fd_qual_out, _header) == -1) {
			exit(1);
		}
		if (!_qual.empty()) {
			for (i = 0; i + FASTA_WIDTH < _qual.size();) {
				size_t n = i + FASTA_WIDTH - 1;
				for (; i != n; ++i) {
					if (pfputs(fd_qual_out, itoa(_qual[i]) + " ") == -1) {
						exit(1);
					}
				}
				if (pfputs(fd_qual_out, itoa(_qual[i++]) + "\n") == -1) {
					exit(1);
				}
			}
			for (; i < _qual.size() - 1; ++i) {
				if (pfputs(fd_qual_out, itoa(_qual[i]) + " ") == -1) {
					exit(1);
				}
			}
			if (pfputs(fd_qual_out, itoa(_qual[i]) + "\n") == -1) {
				exit(1);
			}
		}
	}
}

// upgrade quality of confirmed basepairs to 40
void FastaFile::update_quality(const std::string &confirms) {
	if (_qual.empty()) {
		return;
	}
	std::string::size_type i = static_cast<std::string::size_type>(-1);
	while ((i = confirms.find('1', i + 1)) != std::string::npos) {
		_qual[i] = 40;
	}
}

void FastaFile::remove_pads() {
	std::string::size_type i = _seq.find('-');
	if (i == std::string::npos) {
		return;
	}
	std::string::size_type j = _seq.find_first_not_of('-', i + 1);
	while (j != std::string::npos) {
		_seq[i] = _seq[j];
		if (!_qual.empty()) {
			_qual[i] = _qual[j];
		}
		++i;
		j = _seq.find_first_not_of('-', j + 1);
	}
	_seq.resize(i);
	if (!_qual.empty()) {
		_qual.resize(i);
	}
}

static std::string find_qual_file(const std::string &z, bool skip_check) {
	std::string filename(z);	// so I can modify it
	std::string suffix;
	const char *suffix_list[3] = { ".gz", ".bz2", ".Z" };
	const size_t suffix_size[3] = { 3, 4, 2 };
	int i;
	for (i = 0; i != 3; ++i) {
		const size_t j = suffix_size[i];
		if (filename.size() > j && filename.compare(filename.size() - j, j, suffix_list[i]) == 0) {
			filename.erase(filename.size() - j);
			suffix = suffix_list[i];
			break;
		}
	}
	std::string s;
	struct stat buf;
	s = filename + ".qual";
	if (stat(s.c_str(), &buf) == 0) {
		return s;
	}
	std::string t = s + suffix;
	if (stat(t.c_str(), &buf) == 0) {
		return t;
	}
	for (i = 0; i != 3; ++i) {
		t = s + suffix_list[i];
		if (stat(t.c_str(), &buf) == 0) {
			return t;
		}
	}
	const char *fasta_list[2] = { ".fna", ".fasta" };
	const size_t fasta_size[2] = { 4, 6 };
	for (i = 0; i != 2; ++i) {
		const size_t j = fasta_size[i];
		if (filename.size() > j && filename.compare(filename.size() - j, j, fasta_list[i]) == 0) {
			filename.erase(filename.size() - j);
			filename += ".qual";
			break;
		}
	}
	if (i == 2) {	// no matching file, but this is as close as it gets
		return skip_check ? s + suffix : "";
	}
	if (stat(filename.c_str(), &buf) == 0) {
		return filename;
	}
	s = filename + suffix;
	if (stat(s.c_str(), &buf) == 0) {
		return s;
	}
	for (i = 0; i != 3; ++i) {
		t = filename + suffix_list[i];
		if (stat(t.c_str(), &buf) == 0) {
			return t;
		}
	}
	return skip_check ? s : "";
}

static void read_excludes(const char *filename) {
	int fd = open_compressed(filename);
	if (fd == -1) {
		fprintf(stderr, "Error: open_compressed: %s: %s\n", filename, strerror(errno));
		exit(1);
	}
	std::string line;
	while (pfgets(fd, line) != -1) {
		opt_exclude[line] = 1;
	}
	close_compressed(fd);
}

static void print_usage() {
	fprintf(stderr, "usage: repair_sequence.remote [opts] <-m ##> <tmp_dir> <do_qual> <fasta> <index>\n");
	fprintf(stderr, "\t-A ##\tpercentage of sequence changed that causes an abort [off]\n");
	fprintf(stderr, "\t-c\tdelete data files once finished processing\n");
	fprintf(stderr, "\t-I\tonly make indel changes; if given twice, n's can also be changed\n");
	fprintf(stderr, "\t-k ##\tminimum percent of match overlap [90]\n");
	fprintf(stderr, "\t-l ##\tminimum length of matching read [50]\n");
	fprintf(stderr, "\t-m ##\tminimum number of confirming sequences from db\n");
	fprintf(stderr, "\t-n ##\tminimum number of non-confirming sequences from db\n");
	fprintf(stderr, "\t\t[same as confirming sequences]\n");
	fprintf(stderr, "\t-r\treport matching sequence instead of processing it\n");
	fprintf(stderr, "\t-S\tonly use db reads where both pairs match a given query sequence\n");
	fprintf(stderr, "\t-t ##\tadd trace marker back to read ids\n");
	fprintf(stderr, "\t-x ##\tfile with list of database reads to exclude from matching\n");
	exit(1);
}

static void get_opts(int argc, char **argv) {
	opt_abort = -1;
	opt_change = -1;
	opt_cleanup = 0;
	opt_confirm = -1;
	opt_pair_match = 0;
	opt_indel = 0;
	opt_length_cutoff = 50;
	opt_overlap_cutoff = 90;
	opt_report = 0;
	int c;
	while ((c = getopt(argc, argv, "A:cIk:l:m:n:rSt:Vx:")) != EOF) {
		switch (c) {
		    case 'A':
			opt_abort = atoi(optarg);
			if (opt_abort <= 0) {
				fprintf(stderr, "Error: -C option non-positive: %d\n", opt_abort);
				print_usage();
			}
			break;
		    case 'c':
			opt_cleanup = 1;
			break;
		    case 'I':
			++opt_indel;
			break;
		    case 'k':
			opt_overlap_cutoff = atoi(optarg);
			break;
		    case 'l':
			opt_length_cutoff = atoi(optarg);
			break;
		    case 'm':
			opt_confirm = atoi(optarg);
			if (opt_confirm <= 0) {
				fprintf(stderr, "Error: -m option non-positive: %d\n", opt_confirm);
				print_usage();
			}
			break;
		    case 'n':
			opt_change = atoi(optarg);
			if (opt_change <= 0) {
				fprintf(stderr, "Error: -n option non-positive: %d\n", opt_confirm);
				print_usage();
			}
			break;
		    case 'r':
			opt_report = 1;
			break;
		    case 'S':
			opt_pair_match = 1;
			break;
		    case 't':
			opt_strip_trace = optarg;
			break;
		    case 'V':
			fprintf(stderr, "repair_sequence version %s\n", VERSION);
			exit(0);
		    case 'x':
			read_excludes(optarg);
			break;
		    default:
			print_usage();
		}
	}
	if (opt_confirm == -1 || optind + 4 != argc) {
		print_usage();
	}
	opt_overlap_cutoff /= 100;	// convert to a fraction
	if (opt_change == -1) {
		opt_change = opt_confirm;
	}
}

static void open_blats(const std::string &tmp_dir, const std::string &index, std::list<BlatFile> &blats) {
	DIR *dir = opendir(tmp_dir.c_str());
	if (dir == NULL) {
		fprintf(stderr, "Error: opendir (%s): %s\n", tmp_dir.c_str(), strerror(errno));
		exit(1);
	}
	struct dirent *a;
	std::string::size_type n = index.size();
	errno = 0;
	while ((a = readdir(dir)) != NULL) {
		std::string s(a->d_name);
		// pattern match for ^b$n\.\d+(\.gz|\.bz2)?$, where $n == index
		if (s.size() >= n + 3 && s[0] == 'b' && s[n + 1] == '.' && s.compare(1, n, index) == 0) {
			std::string::size_type i = s.find_first_not_of("0123456789", n + 2);
			if (i == std::string::npos || (i == s.size() - 3 && s.compare(s.size() - 3, 3, ".gz") == 0) || (i == s.size() - 4 && s.compare(s.size() - 4, 4, ".bz2") == 0)) {
				blats.resize(blats.size() + 1);
				blats.back().open_blat(tmp_dir + "/" + s);
			}
		}
		errno = 0;
	}
	if (errno != 0) {
		fprintf(stderr, "Error: readdir (%s): %s\n", tmp_dir.c_str(), strerror(errno));
		exit(1);
	}
	closedir(dir);
}

static bool find_next_blats(std::list<BlatFile> &blats, const std::string &id) {
	bool found_one = 0;
	std::list<BlatFile>::iterator a = blats.begin();
	std::list<BlatFile>::iterator end_a = blats.end();
	for (; a != end_a; ++a) {
		found_one |= a->find_next(id);
	}
	return found_one;
}

static void get_matches1(std::list<BlatFile> &blats, size_t seq_length, std::list<Match> &matches, std::list<Match> &all_matches) {
	std::list<BlatFile>::iterator a = blats.begin();
	std::list<BlatFile>::iterator end_a = blats.end();
	for (; a != end_a; ++a) {
		if (!a->is_current()) {
			continue;
		}
		Match b;
		while (a->read_next(seq_length, b)) {
			all_matches.push_back(b);
			if (b.trim_alignment()) {
				matches.push_back(b);
			}
		}
		a->finish_current();
	}
}

// as above, except only adds matches for which both forward and
// reverse parts of the read match

static void get_matches2(std::list<BlatFile> &blats, size_t seq_length, std::list<Match> &matches, std::list<Match> &all_matches) {
	all_matches.push_back(Match());
	std::list<Match>::const_iterator b = all_matches.end();
	const std::list<Match>::const_iterator end_b = b--;
	std::map<std::string, int> pairs;
	std::list<BlatFile>::iterator a = blats.begin();
	std::list<BlatFile>::iterator end_a = blats.end();
	for (; a != end_a; ++a) {
		if (!a->is_current()) {
			continue;
		}
		while (a->read_next(seq_length, all_matches.back())) {
			const std::string &s = all_matches.back().id;
			const size_t n = s.size();
			if (n < 3 || s[n - 2] != '/') {
			} else if (s[n - 1] == '1') {
				pairs[s.substr(0, n - 2)] |= 1;
			} else if (s[n - 1] == '2') {
				pairs[s.substr(0, n - 2)] |= 2;
			}
			all_matches.push_back(Match());
		}
		a->finish_current();
	}
	if (++b == end_b) {	// check to make sure at least one was added
		all_matches.pop_back();
		return;
	}
	--b;
	all_matches.pop_back();
	matches.push_back(Match());
	for (; b != end_b; ++b) {
		const std::string &s = b->id;
		const size_t n = s.size();
		if (n > 2) {
			const std::map<std::string, int>::const_iterator d = pairs.find(s.substr(0, n - 2));
			if (d != pairs.end() && d->second == 3) {
				matches.back() = *b;
				if (matches.back().trim_alignment()) {
					matches.push_back(Match());
				}
			}
		}
	}
	matches.pop_back();
}

static void add_coverage(std::vector<int> &coverage, const std::list<Match> &matches) {
	std::list<Match>::const_iterator a = matches.begin();
	std::list<Match>::const_iterator end_a = matches.end();
	for (; a != end_a; ++a) {
		long i;
		for (i = a->start; i != a->stop; ++i) {
			++coverage[i];
		}
	}
}

// if coverage is sufficiently uneven, toss matches to even coverage
// (lowest score first, in excessively covered region)
static void fix_coverage(std::list<Match> &matches, size_t seq_length) {
	// calculate basepair coverage levels
	std::vector<int> coverage(seq_length);
	add_coverage(coverage, matches);
	int min = std::max(opt_confirm, opt_change);
	min *= min;
	size_t i;
	for (i = 0; i != seq_length && coverage[i] < min; ++i) { }
	if (i == seq_length) {		// coverage is too low to matter
		return;
	}
	// find lowest (at least min) and highest coverage levels
	int low = coverage[i];
	int high = low;
	long highest = i;
	for (++i; i != seq_length; ++i) {
		if (coverage[i] < min) {
		} else if (low > coverage[i]) {
			low = coverage[i];
		} else if (high < coverage[i]) {
			high = coverage[i];
			highest = i;
		}
	}
	// max is level to reduce coverage to
	int max = std::max(2 * low, (int)floor(sqrt(high)));
	if (high <= max) {		// not enough variation to matter
		return;
	}
	matches.sort(CompareMatchScore());
	while (high > max) {
		// do in chunks for time efficiency
		int x = (high - max + 1) / 2;
		std::list<Match>::iterator a = matches.begin();
		for (; x != 0; --x) {	// remove lowest scores in region
			// skip matches covering basepair
			for (; highest < a->start || a->stop < highest; ++a) { }
			long j;
			for (j = a->start; j != a->stop; ++j) {
				--coverage[j];
			}
			a = matches.erase(a);
		}
		// find new highest coverage
		for (i = 0; coverage[i] < min; ++i) { }
		high = coverage[i];
		highest = i;
		for (++i; i != seq_length; ++i) {
			if (high < coverage[i]) {
				high = coverage[i];
				highest = i;
			}
		}
	}
}

// scrub matches to replace error regions with n's,
// then sort matches into good and bad ones
static void sort_matches(std::list<Match> &matches, std::list<Match> &good_matches, std::list<Match> &bad_matches) {
	// find polymorphism coverage rate to evaluate basepairs;
	// get rid of non-indel matches if option is set
	std::map<int, std::map<char, int> > polys;	// pos -> bp -> count
	std::list<Match>::iterator a = matches.begin();
	std::list<Match>::iterator end_a = matches.end();
	while (a != end_a) {
		std::string::size_type i = a->tag.find(' ');
		if (i == std::string::npos) {
			good_matches.splice(good_matches.end(), matches, a++);
		} else if (opt_indel && a->qs[i] != '-' && a->ss[i] != '-' && (opt_indel == 1 || a->qs[i] != 'n')) {
			a = matches.erase(a);
		} else {
			std::list<std::string::size_type> list(1, i);
			for (i = a->tag.find(' ', i + 1); i != std::string::npos; i = a->tag.find(' ', i + 1)) {
				if (opt_indel && a->qs[i] != '-' && a->ss[i] != '-' && (opt_indel == 1 || a->qs[i] != 'n')) {
					a = matches.erase(a);
					break;
				}
				list.push_back(i);
			}
			if (i == std::string::npos) {
				std::list<std::string::size_type>::const_iterator b = list.begin();
				std::list<std::string::size_type>::const_iterator end_b = list.end();
				for (; b != end_b; ++b) {
					++polys[a->start + *b][a->ss[*b]];
				}
				++a;
			}
		}
	}
	// remove high coverage polymorphisms from list
	int min_coverage = (std::min(opt_confirm, opt_change) + 1) / 2;
	std::map<int, std::map<char, int> >::iterator b = polys.begin();
	std::map<int, std::map<char, int> >::iterator end_b = polys.end();
	while (b != end_b) {
		std::map<char, int>::iterator c = b->second.begin();
		std::map<char, int>::iterator end_c = b->second.end();
		while (c != end_c) {
			if (c->second >= min_coverage) {
				b->second.erase(c++);
			} else {
				++c;
			}
		}
		if (b->second.empty()) {
			polys.erase(b++);
		} else {
			++b;
		}
	}
	// take low support polymorphisms, and change to n's
	a = matches.begin();
	while (a != end_a) {
		bool scrubbed = 0;
		std::string::size_type i = static_cast<std::string::size_type>(-1);
		while ((i = a->tag.find(' ', i + 1)) != std::string::npos) {
			b = polys.find(a->start + i);
			if (b != polys.end() && b->second.find(a->ss[i]) != b->second.end()) {
				a->scrub(i);
				scrubbed = 1;
			}
		}
		if (scrubbed && a->post_scrub_trim()) {
			a = matches.erase(a);
		} else {
			++a;
		}
	}
	// sort resulting matches
	a = matches.begin();
	while (a != end_a) {
		if (a->tag.find(' ') == std::string::npos) {
			good_matches.splice(good_matches.end(), matches, a++);
		} else {
			bad_matches.splice(bad_matches.end(), matches, a++);
		}
	}
}

// counts number of pads before a given location
static void count_pads(const std::list<Match> &matches, std::vector<int> &pads) {
	std::list<Match>::const_iterator a = matches.begin();
	std::list<Match>::const_iterator end_a = matches.end();
	for (; a != end_a; ++a) {
		long k = a->start;
		std::string::size_type i = 0;
		while ((i = a->qs.find('-', i)) != std::string::npos) {
			std::string::size_type j = a->qs.find_first_not_of('-', i + 1);
			if (j == std::string::npos) {
				j = a->qs.size();
			}
			int n = j - i;
			k -= n;
			if (pads[k + j] < n) {
				pads[k + j] = n;
			}
			i = j + 1;
		}
	}
}

// create a padded -> unpadded transformation vector; also counts total pads
// before each location; returns number of pads
static size_t make_unpadded(const std::vector<int> &pads, std::vector<int> &unpad, std::vector<int> &all_pads) {
	all_pads.reserve(pads.size());
	unpad.reserve(pads.size());	// underestimate
	size_t i, j, total;
	for (i = j = total = 0; i != pads.size(); ++i) {
		int n = pads[i];
		total += n;
		all_pads.push_back(total);
		int k;
		for (k = -1; k != n; ++k, ++j) {
			unpad.push_back(i + 1);		// change to 1 offset
		}
	}
	return unpad.size() - pads.size();
}

// pad matches
static void add_pads_match(const std::vector<int> &pads, const std::vector<int> &all_pads, std::list<Match> &matches) {
	std::list<Match>::iterator a = matches.begin();
	std::list<Match>::iterator end_a = matches.end();
	for (; a != end_a; ++a) {
		a->add_pads(pads, all_pads);
	}
}

static void initialize_confirms(std::string &confirms, const std::vector<int> &confirm) {
	size_t i;
	for (i = 0; i != confirm.size(); ++i) {
		if (confirm[i] >= opt_confirm) {
			confirms[i] = '1';
		}
	}
}

static bool all_confirmed(const std::string &confirms, const Match &match) {
	long i;
	for (i = match.start; i != match.stop; ++i) {
		if (confirms[i] == '0') {
			return 0;
		}
	}
	return 1;
}

// remove all good matches covering all confirmed sequence
static void clean_good_matches(std::list<Match> &matches, const std::string &confirms) {
	std::list<Match>::iterator a = matches.begin();
	std::list<Match>::iterator end_a = matches.end();
	while (a != end_a) {
		if (all_confirmed(confirms, *a)) {
			a = matches.erase(a);
		} else {
			++a;
		}
	}
}

// remove bad matches that disagree with a confirmed basepair;
// make list of possibly supported changes
static void make_change_list(std::list<Match> &matches, const std::string &confirms, std::map<long, std::map<char, int> > &changes) {
	std::list<Match>::iterator a = matches.begin();
	std::list<Match>::iterator end_a = matches.end();
	while (a != end_a) {
		std::list<long> list;
		std::string::size_type i = static_cast<std::string::size_type>(-1);
		while ((i = a->tag.find(' ', i + 1)) != std::string::npos) {
			if (confirms[a->start + i] == '1') {
				a = matches.erase(a);
				break;
			}
			list.push_back(i);
		}
		if (i == std::string::npos) {
			std::list<long>::const_iterator b = list.begin();
			std::list<long>::const_iterator end_b = list.end();
			for (; b != end_b; ++b) {
				++changes[a->start + *b][a->ss[*b]];
			}
			++a;
		}
	}
	// remove possible changes with insufficient support
	std::map<long, std::map<char, int> >::iterator b = changes.begin();
	std::map<long, std::map<char, int> >::iterator end_b = changes.end();
	while (b != end_b) {
		std::map<char, int>::iterator c = b->second.begin();
		std::map<char, int>::iterator end_c = b->second.end();
		while (c != end_c) {
			if (c->second < opt_change) {
				b->second.erase(c++);
			} else {
				++c;
			}
		}
		if (b->second.empty()) {
			changes.erase(b++);
		} else {
			++b;
		}
	}
}

// subtract match mis-matches from changes
static void remove_changes(std::map<long, std::map<char, int> > &changes, const Match &match) {
	std::string::size_type i = static_cast<std::string::size_type>(-1);
	while ((i = match.tag.find(' ', i + 1)) != std::string::npos) {
		std::map<long, std::map<char, int> >::iterator a = changes.find(match.start + i);
		if (a != changes.end()) {
			std::map<char, int>::iterator b = a->second.find(match.qs[i]);
			if (b != a->second.end() && --b->second < opt_change) {
				a->second.erase(b);
				if (a->second.empty()) {
					changes.erase(a);
				}
			}
		}
	}
}

// find coverage of good reads (not including n's)
static void add_confirm(std::vector<int> &confirm, std::list<Match>::const_iterator a, std::list<Match>::const_iterator end_a) {
	const std::string good("-acgt");
	for (; a != end_a; ++a) {
		long i;
		std::string::size_type j;
		for (i = a->start, j = 0; i != a->stop; ++i, ++j) {
			if (good.find(a->ss[j]) != std::string::npos) {
				++confirm[i];
			}
		}
	}
}

static void update_matches(std::list<Match> &good_matches, std::list<Match> &bad_matches, std::vector<int> &confirm, std::string &confirms, std::map<long, std::map<char, int> > &changes, long position, char c) {
	const std::string good("-acgt");
	std::string non_match(good);
	non_match.erase(non_match.find(c), 1);
	// remove good matches that no longer match the changed basepair;
	// also, remove new perfect matches
	std::list<Match>::iterator a = good_matches.begin();
	std::list<Match>::iterator end_a = good_matches.end();
	while (a != end_a) {
		if (a->start <= position && position < a->stop) {
			if (non_match.find(a->ss[position - a->start]) != std::string::npos) {
				long i;
				std::string::size_type j;
				for (i = a->start, j = 0; i != a->stop; ++i, ++j) {
					if (good.find(a->ss[j]) != std::string::npos) {
						--confirm[i];
					}
				}
				a = good_matches.erase(a);
			} else if (all_confirmed(confirms, *a)) {
				a = good_matches.erase(a);
			} else {
				++a;
			}
		} else {
			++a;
		}
	}
	// update bad matches for the change, move to good if they now agree
	a = bad_matches.begin();
	end_a = bad_matches.end();
	long check_start = position;
	long check_stop = position;
	while (a != end_a) {
		if (a->start <= position && position < a->stop) {
			if (non_match.find(a->ss[position - a->start]) != std::string::npos) {
				// whoops, will never match now
				remove_changes(changes, *a);
				a = bad_matches.erase(a);
			} else {
				a->tag[position - a->start] = '|';
				if (a->tag.find(' ') == std::string::npos) {
					// it's a good match now
					if (check_start > a->start) {
						check_start = a->start;
					}
					if (check_stop < a->stop) {
						check_stop = a->stop;
					}
					std::list<Match>::iterator b = a++;
					add_confirm(confirm, b, a);
					if (all_confirmed(confirms, *b)) {
						bad_matches.erase(b);
					} else {
						good_matches.splice(good_matches.end(), bad_matches, b);
					}
				} else {
					++a;
				}
			}
		} else {
			++a;
		}
	}
	// find new confirmations that could knock out mismatches
	std::map<long, bool> new_changes;
	long revise_start = -1;
	long revise_stop = 0;			// to avoid compiler warning
	long i;
	for (i = check_start; i != check_stop; ++i) {
		if (confirm[i] >= opt_confirm && confirms[i] == '0') {
			confirms[i] = '1';
			std::map<long, std::map<char, int> >::iterator b = changes.find(i);
			if (b != changes.end()) {
				changes.erase(b);
				new_changes[i] = 1;
			}
			if (revise_start == -1) {
				revise_start = i;
			}
			revise_stop = i;
		}
	}
	if (revise_start == -1) {
		return;
	}
	// check for new perfect matches
	a = good_matches.begin();
	end_a = good_matches.end();
	while (a != end_a) {
		if (a->start <= revise_stop && revise_start < a->stop && all_confirmed(confirms, *a)) {
			a = good_matches.erase(a);
		} else {
			++a;
		}
	}
	// check for bad matches with new confirmed mismatches
	if (new_changes.empty()) {
		return;
	}
	a = bad_matches.begin();
	end_a = bad_matches.end();
	while (a != end_a) {
		if (a->start <= revise_stop && revise_start < a->stop) {
			std::string::size_type j = static_cast<std::string::size_type>(-1);
			while ((j = a->tag.find(' ', j + 1)) != std::string::npos && new_changes.find(a->start + j) == new_changes.end()) { }
			if (j != std::string::npos) {
				remove_changes(changes, *a);
				a = bad_matches.erase(a);
			} else {
				++a;
			}
		} else {
			++a;
		}
	}
}

static void make_changes(std::list<Match> &good_matches, std::list<Match> &bad_matches, std::vector<int> &confirm, std::string &confirms, std::map<long, std::map<char, int> > &changes, FastaFile &fasta, std::map<long, char> &changes_made, const std::vector<int> &unpad, const size_t change_cutoff) {
	std::string change_log;
	while (change_cutoff != changes_made.size()) {
		// find best supported change
		char c = 0;			// to avoid compiler warning
		long position = -1;		// to avoid compiler warning
		int score = opt_change - 1;
		std::map<long, std::map<char, int> >::iterator a = changes.begin();
		std::map<long, std::map<char, int> >::iterator end_a = changes.end();
		for (; a != end_a; ++a) {
			std::map<char, int>::iterator b = a->second.begin();
			std::map<char, int>::iterator end_b = a->second.end();
			for (; b != end_b; ++b) {
				if (score < b->second) {
					score = b->second;
					position = a->first;
					c = b->first;
				}
			}
		}
		if (score == opt_change - 1) {
			if (!change_log.empty() && pfputs(fd_log, itoa(unpad[position]) + ": " + fasta.seq()[position] + " -> " + (char)toupper(c) + "\n") == -1) {
				exit(1);
			}
			return;
		}
		// make change
		changes.erase(position);
		changes_made[position] = c;
		if (fd_log != -1) {
			if (change_cutoff == static_cast<size_t>(-1)) {
				if (pfputs(fd_log, itoa(unpad[position]) + ": " + fasta.seq()[position] + " -> " + (char)toupper(c) + "\n") == -1) {
					exit(1);
				}
			} else {
				change_log += itoa(unpad[position]) + ": " + fasta.seq()[position] + " -> " + (char)toupper(c) + "\n";
			}
		}
		fasta.update_seq(position, c);
		// force good to prevent cycles
		confirms[position] = '1';
		update_matches(good_matches, bad_matches, confirm, confirms, changes, position, c);
	}
	fprintf(stderr, "Warning: %s: hit cutoff limit, reverting changes\n", fasta.id().c_str());
	fasta.revert();
	changes_made.clear();
	confirms.assign(confirms.size(), '0');	// to prevent quality changes
}

static void print_report(int fd, std::list<Match> &matches, const std::map<long, char> &changes) {
	if (matches.empty()) {
		return;
	}
	matches.sort();
	long offset = matches.front().start;
	std::list<Match>::iterator a = matches.begin();
	std::list<Match>::iterator end_a = matches.end();
	for (; a != end_a; ++a) {
		long i = itoa(std::max(a->start + 1, a->s_start)).size() + 1;
		long j = a->start - offset;
		if (j < i) {
			offset += j - i;
			j = i;
		}
		a->update_tag(changes);
		if (pfputs(fd, std::string(j - i, ' ') + ">" + a->id + "\n") == -1) {
			exit(1);
		}
		std::string s = itoa(a->start + 1);
		if (pfputs(fd, std::string(j - s.size() - 1, ' ') + s + " " + a->qs + " " + itoa(a->stop) + "\n") == -1) {
			exit(1);
		}
		if (pfputs(fd, std::string(j, ' ') + a->tag + "\n") == -1) {
			exit(1);
		}
		s = itoa(a->s_start);
		if (pfputs(fd, std::string(j - s.size() - 1, ' ') + s + " " + a->ss + " " + itoa(a->s_stop) + "\n") == -1) {
			exit(1);
		}
	}
	if (pfputs(fd, "\n") == -1) {
		exit(1);
	}
}

static void process_blats(FastaFile &fasta, std::list<BlatFile> &blats, int fd_fasta, int fd_qual) {
	if (fd_log != -1) {
		if (pfputs(fd_log, ">" + fasta.id() + "\n") == -1) {
			exit(1);
		}
	}
	if (!find_next_blats(blats, fasta.id())) {
		fasta.add_trace(opt_strip_trace);
		if (opt_report) {
			if (pfputs(fd_fasta, fasta.header()) == -1) {
				exit(1);
			}
		} else {
			fasta.print(fd_fasta, fd_qual);
		}
		return;
	}
	fasta.add_trace(opt_strip_trace);
	std::list<Match> matches, all_matches;
	if (opt_pair_match) {
		get_matches2(blats, fasta.seq().size(), matches, all_matches);
	} else {
		get_matches1(blats, fasta.seq().size(), matches, all_matches);
	}
	fix_coverage(matches, fasta.seq().size());
	std::vector<int> pads(fasta.seq().size());
	count_pads(matches, pads);
	std::vector<int> unpad, all_pads;
	size_t total_pads = make_unpadded(pads, unpad, all_pads);
	size_t change_cutoff = opt_abort == -1 ? static_cast<size_t>(-1) : fasta.seq().size() * opt_abort / 100;
	fasta.add_pads(pads, total_pads);
	add_pads_match(pads, all_pads, matches);
	std::list<Match> good_matches, bad_matches;
	sort_matches(matches, good_matches, bad_matches);
	std::vector<int> confirm(fasta.seq().size());
	add_confirm(confirm, good_matches.begin(), good_matches.end());
	// confirms is sticky - once set to 1, will never change back to 0
	std::string confirms(confirm.size(), '0');
	initialize_confirms(confirms, confirm);
	clean_good_matches(good_matches, confirms);
	std::map<long, std::map<char, int> > changes;
	make_change_list(bad_matches, confirms, changes);
	std::map<long, char> changes_made;
	make_changes(good_matches, bad_matches, confirm, confirms, changes, fasta, changes_made, unpad, change_cutoff);
	if (opt_report) {
		add_pads_match(pads, all_pads, all_matches);
		if (pfputs(fd_fasta, fasta.header()) == -1) {
			exit(1);
		}
		print_report(fd_fasta, all_matches, changes_made);
	} else {
		fasta.update_quality(confirms);
		fasta.remove_pads();
		fasta.print(fd_fasta, fd_qual);
	}
}

static void close_blats(std::list<BlatFile> &blats) {
	std::list<BlatFile>::iterator a = blats.begin();
	std::list<BlatFile>::iterator end_a = blats.end();
	for (; a != end_a; ++a) {
		a->close_blat();
	}
}

static void unlink_blats(std::list<BlatFile> &blats) {
	std::list<BlatFile>::iterator a = blats.begin();
	std::list<BlatFile>::iterator end_a = blats.end();
	for (; a != end_a; ++a) {
		a->unlink_blat();
	}
}

int main(int argc, char **argv) {
	get_opts(argc, argv);
	std::string tmp_dir(argv[optind++]);
	bool do_qual = atoi(argv[optind++]) && !opt_report;
	// open input files
	FastaFile fasta(argv[optind++], do_qual, tmp_dir);
	if (!fasta.is_open()) {
		fprintf(stderr, "Error: could not open fasta file %s: %s\n", fasta.fasta_file().c_str(), strerror(errno));
		return 1;
	}
	if (do_qual && !fasta.has_qual()) {
		fprintf(stderr, "Warning: no qual file found for %s\n", fasta.fasta_file().c_str());
	}
	std::string index(argv[optind]);
	std::list<BlatFile> blats;
	open_blats(tmp_dir, index, blats);
	// open output files
	std::list<std::string> fork_args(1, "bzip2");
	fork_args.push_back("-c");
	int fd_fasta_out = write_fork(fork_args, tmp_dir + "/f" + index + ".bz2");
	if (fd_fasta_out == -1) {
		fprintf(stderr, "Error: could not write fasta output file\n");
		return 1;
	}
	int fd_qual_out = -1;
	if (fasta.has_qual()) {
		fd_qual_out = write_fork(fork_args, tmp_dir + "/q" + index + ".bz2");
		if (fd_fasta_out == -1) {
			fprintf(stderr, "Error: could not write qual output file\n");
			return 1;
		}
	}
	fd_log = write_fork(fork_args, tmp_dir + "/l" + index + ".bz2");
	if (fd_log == -1) {
		fprintf(stderr, "Warning: could not write log file\n");
	}
	// process files
	while (fasta.read_next()) {
		process_blats(fasta, blats, fd_fasta_out, fd_qual_out);
	}
	// close files
	close_blats(blats);
	fasta.fasta_close();
	close_fork(fd_fasta_out);
	if (fd_qual_out != -1) {
		close_fork(fd_qual_out);
	}
	if (fd_log != -1) {
		close_fork(fd_log);
	}
	if (opt_cleanup) {	// delete temporary files, if asked to
		fasta.fasta_unlink();
		unlink_blats(blats);
	}
	return 0;
}
