#include "itoa.h"	// itoa()
#include "open_compressed.h"	// close_compressed(), open_compressed(), pfgets()
#include "strtostr.h"	// strtostr()
#include "write_fork.h"	// close_fork(), pfputs(), write_fork()
#include <algorithm>	// swap()
#include <dirent.h>	// DIR, closedir(), opendir(), readdir()
#include <errno.h>	// errno
#include <exception>	// exception
#include <fstream>	// ofstream
#include <getopt.h>	// getopt(), optarg, optind
#include <glob.h>	// GLOB_NOSORT, glob(), globfree(), glob_t
#include <iostream>	// cerr
#include <list>		// list<>
#include <map>		// map<>, multimap<>
#include <set>		// set<>
#include <stdlib.h>	// abs()
#include <string.h>	// strerror()
#include <string>	// string
#include <sys/stat.h>	// S_IFDIR, mkdir(), stat(), struct stat
#include <unistd.h>	// getcwd(), pathconf(), unlink(), _PC_PATH_MAX
#include <utility>	// make_pair()
#include <vector>	// vector<>

#define MAX_DB_SIZE 745378110UL
#define SEARCH_MAX (1UL << 50)
#define INSERT_LENGTH 48

#define BZIP2 "/usr/bin/bzip2"
#define GZIP "/usr/bin/gzip"
#define SSD_DIR "/mnt/ssd/tmp"

class LocalException : public std::exception {
    private:
	const std::string error_;
	const int show_usage_;
    public:
	explicit LocalException(const std::string &s, int i = 0) : error_(s), show_usage_(i) { }
	virtual ~LocalException(void) throw() { }
	virtual const char * what() const throw() {
		return error_.c_str();
	}
	const int &show_usage(void) const throw() {
		return show_usage_;
	}
};

std::string PWD;

static void set_pwd() {
	errno = 0;	// allow for testing of no limit given
	long size(pathconf(".", _PC_PATH_MAX));
	if (size != -1) {
		++size;			// add room for trailing null
	} else if (errno == 0) {	// no limit given
		size = 1 << 16;
	} else {
		throw LocalException("pathconf: " + std::string(strerror(errno)));
	}
	char * const buf(new char[size]);
	if (getcwd(buf, size) == 0) {
		throw LocalException("getcwd: " + std::string(strerror(errno)));
	}
	PWD = buf;
	delete[] buf;
}

static void print_usage() {
	std::cerr <<
		"usage: copy_dbs [opts] <run_dir> <min_dbs>\n"
		"\t-d ##\tfiles to blat\n"
		"\t\t(may be specified multiple times; will expand globs)\n"
		"\t-l ##\tlinker (strips this and any sequence past it)\n"
		"\t-M ##\tminimum length of read\n"
		"\t-x ##\tfile containing a list of reads to not blat\n"
		"\t\t(may be specified multiple times)\n";
}

// returns zero normally, 1 if help was asked for
static int get_opts(const int argc, char * const argv[], std::list<std::string> &db_list, std::vector<std::string> &exclude_names, std::string &run_dir, size_t &opt_length_cutoff, std::string &opt_linker, int &min_dbs) {
	opt_length_cutoff = 1;
	int c;
	while ((c = getopt(argc, argv, "d:hl:M:x:")) != EOF) {
		switch(c) {
		    case 'd':
			db_list.push_back(optarg);
			if (db_list.back().empty()) {
				std::cerr << "Warning: empty database name\n";
				db_list.pop_back();
			}
			break;
		    case 'h':
			print_usage();
			return 1;
		    case 'l':
			opt_linker = optarg;
			break;
		    case 'M':
			c = atoi(optarg);
			if (c < 1) {
				throw LocalException("bad length cutoff: " + std::string(optarg));
			}
			opt_length_cutoff = c;
			break;
		    case 'x':
			exclude_names.push_back(optarg);
			break;
		    default:
			throw LocalException("bad option: " + char(c), 1);
		}
	}
	if (db_list.empty()) {
		throw LocalException("no -d options given", 1);
	}
	if (optind + 2 != argc) {
		throw LocalException("incorrect number of arguments", 1);
	}
	run_dir = argv[optind];
	if (run_dir.empty()) {
		throw LocalException("blank run_dir", 1);
	}
	min_dbs = atoi(argv[optind + 1]);
	if (min_dbs < 1) {
		throw LocalException("bad minimum number of dbs: " + std::string(argv[optind + 1]));
	}
	return 0;
}

static void read_excludes(const std::vector<std::string> &exclude_names, std::map<std::string, int> &exclude_list) {
	std::vector<std::string>::const_iterator a(exclude_names.begin());
	const std::vector<std::string>::const_iterator end_a(exclude_names.end());
	for (int i(-1); a != end_a; ++a, --i) {
		const int fd(open_compressed(*a));
		if (fd != -1) {
			std::string line;
			while (pfgets(fd, line) != -1) {
				exclude_list[line] = i;
			}
			close_compressed(fd);
		} else {
			throw LocalException("failed to open " + *a);
		}
	}
}

// clean up a directory name
static void cleanup_dir(std::string &dir) {
	if (dir.empty()) {
		return;
	}
	// remove duplicate /'s
	size_t j(0);
	const size_t end_i(dir.size());
	for (size_t i(1); i != end_i; ++i) {
		if (dir[i] != '/' || dir[j] != '/') {
			if (i != ++j) {
				dir[j] = dir[i];
			}
		}
	}
	// don't keep trailing /, if any
	dir.resize(dir[j] == '/' ? j : j + 1);
}

// clean up run_dir and convert it to an ssd dir for tmp_dir
static std::string get_tmp_dir(std::string &run_dir) {
	cleanup_dir(run_dir);
	if (run_dir.empty()) {
		throw LocalException("bad run_dir: /");
	}
	const std::string::size_type j(run_dir.rfind('/'));
	return std::string(SSD_DIR) + "/" + run_dir.substr(j == std::string::npos ? 0 : j + 1);
}

// do glob expansion on database file names; also, if directory, replace
// with all non-dot files in the directory, and check if each file is a
// list of files, rather than a fasta file (and expand those names, if so)

static void expand_included_files(std::list<std::string> &files) {
	std::list<std::string>::iterator a(files.begin());
	const std::list<std::string>::const_iterator end_a(files.end());
	while (a != end_a) {
		// convert to absolute path name
		if ((*a)[0] != '/') {
			a->insert(0, PWD + "/");
		}
		// shell expansion
		glob_t b;
		if (glob(a->c_str(), GLOB_NOSORT, 0, &b) != 0) {
			throw LocalException("glob: " + *a);
		}
		if (b.gl_pathc != 0) {
			// retain pointer to start of list
			const std::list<std::string>::iterator c(a);
			++a;
			for (size_t i(0); i != b.gl_pathc; ++i) {
				files.insert(a, b.gl_pathv[i]);
			}
			// remove old entry and continue processing
			// starting with the newly expanded values
			a = files.erase(c);
		}
		globfree(&b);
		struct stat buf;
		if (stat(a->c_str(), &buf) == -1) {
			throw LocalException("stat: " + *a + ": " + strerror(errno));
		} else if (buf.st_mode & S_IFDIR) {	// directory expansion
			DIR * const dir(opendir(a->c_str()));
			if (dir == 0) {
				throw LocalException("opendir: " + *a + ": " + strerror(errno));
			}
			// retain pointer to start of list
			const std::list<std::string>::iterator c(a);
			for (++a;;) {
				errno = 0;
				struct dirent * const d(readdir(dir));
				if (d) {
					if (d->d_name && d->d_name[0] && d->d_name[0] != '.') {
						files.insert(a, *c + "/" + d->d_name);
					}
				} else if (errno) {
					throw LocalException("readdir: " + *c + ": " + strerror(errno));
				} else {
					break;
				}
			}
			closedir(dir);
			// remove old entry and continue processing
			// starting with the newly expanded values
			a = files.erase(c);
		} else {		// file expansion
			const int fd(open_compressed(*a));
			if (fd == -1) {
				throw LocalException("open_compressed: " + *a);
			}
			std::string line;
			while (pfgets(fd, line) != -1 && line.empty()) { }
			if (!line.empty() && line[0] != '>') {
				// retain pointer to start of list
				const std::list<std::string>::iterator c(a);
				++a;
				// get directory of expanded file
				const size_t i(c->rfind('/'));
				c->resize(i == std::string::npos ? 0 : i);
				cleanup_dir(*c);
				*c += '/';
				do {
					if (!line.empty()) {
						files.insert(a, line[0] == '/' ? line : *c + line);
					}
				} while (pfgets(fd, line) != -1);
				a = files.erase(c);
			} else {
				++a;
			}
			close_compressed(fd);
		}
	}
}

// returns 1 if it got something, 0 if out of file
static int read_fasta_next(const int fd, std::string &id, std::string &seq, std::string &line) {
	id = line;
	while (id.empty() && pfgets(fd, id) != -1) { }
	if (id.empty()) {
		return 0;
	}
	std::string::size_type i(1);	// skip leading >
	id = strtostr(id, &i);		// get id from header
	seq.clear();
	while (pfgets(fd, line) != -1 && (line.empty() || line[0] != '>')) {
		seq += line;
	}
	return 1;
}

static void complement(std::string &s) {
	std::string::size_type i(0);
	std::string::size_type j(s.size());
	for (; i != j; ++i) {
		switch (s[i]) {
		     case 'A':
			s[i] = 'T';
			break;
		     case 'C':
			s[i] = 'G';
			break;
		     case 'G':
			s[i] = 'C';
			break;
		     case 'T':
			s[i] = 'A';
			break;
		     case 'a':
			s[i] = 't';
			break;
		     case 'c':
			s[i] = 'g';
			break;
		     case 'g':
			s[i] = 'c';
			break;
		     case 't':
			s[i] = 'a';
			break;
		}
	}
	for (i = 0, --j; i < j; ++i, --j) {
		std::swap(s[i], s[j]);
	}
}

static void split_db(const std::string &run_dir, const std::string &tmp_dir, const std::list<std::string> &target_files, std::map<std::string, int> &exclude_list, const size_t opt_length_cutoff, const std::string &opt_linker, std::multimap<size_t, std::string> &q_list) {
	const size_t max_size(SEARCH_MAX / MAX_DB_SIZE);
	const std::string insert_seq(INSERT_LENGTH, 'N');
	std::string filename;
	size_t filenum(0);
	size_t printed(0);
	std::set<std::string> skipped;
						  // reads looking for partners
	std::map<std::string, std::string> reads; // id, seq
	std::list<std::string> args(1, BZIP2);
	args.push_back("-c");
	const int fd_reads(write_fork(args, run_dir + "/read_names.bz2"));
	if (fd_reads == -1) {
		throw LocalException("could not write read_names.bz2");
	}
	args.front() = GZIP;
	int fd_out(-1);
	std::list<std::string>::const_iterator a(target_files.begin());
	const std::list<std::string>::const_iterator end_a(target_files.end());
	for (; a != end_a; ++a) {
		const int fd(open_compressed(*a));
		if (fd == -1) {
			throw LocalException("could not open " + *a);
		}
		std::string id, seq, line;
		while (read_fasta_next(fd, id, seq, line)) {
			const std::map<std::string, int>::iterator b(exclude_list.find(id));
			if (b != exclude_list.end()) {
				b->second = abs(b->second);
				continue;
			}
			// trim linker and any trailing sequence
			if (!opt_linker.empty()) {
				const std::string::size_type i(seq.find(opt_linker));
				if (i != std::string::npos) {
					seq.resize(i);
				}
			}
			if (seq.size() < opt_length_cutoff) {
				skipped.insert(id);
				continue;
			}
			const char ending(*id.rbegin());
			if (ending != '1' && ending != '2') {
				skipped.insert(id);
				continue;
			}
			std::string x(id);
			*x.rbegin() = ending == '1' ? '2' : '1';
			const std::map<std::string, std::string>::iterator c(reads.find(x));
			// other read pair not present yet
			if (c == reads.end()) {
				reads[id] = seq;
				continue;
			}
			if (ending == '1') {
				id.replace(id.size() - 2, 2, "-" + itoa(seq.size()));
				complement(seq);
				seq += insert_seq + c->second;
			} else {
				id.replace(id.size() - 2, 2, "-" + itoa(c->second.size()));
				complement(c->second);
				seq = c->second + insert_seq + seq;
			}
			reads.erase(c);
			const size_t size(id.size() + seq.size() + 3);
			if (filename.empty() || printed + size > max_size) {
				if (!filename.empty()) {
					//out.close();
					close_fork(fd_out);
					q_list.insert(make_pair(printed, filename));
				}
				filename = tmp_dir + "/q" + itoa(filenum) + ".gz";
				++filenum;
				fd_out = write_fork(args, filename);
				//out.open(filename.c_str());
				printed = 0;
			}
			id += "\n";
			pfputs(fd_reads, id);
			pfputs(fd_out, ">" + id + seq + "\n");
			//out << ">" << id << seq << "\n";
			printed += size;
		}
		close_compressed(fd);
	}
	close_fork(fd_reads);
	if (!filename.empty()) {
		//out.close();
		close_fork(fd_out);
		q_list.insert(make_pair(printed, filename));
	}
	std::ofstream out;
	out.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	out.open((run_dir + "/no_match").c_str());
	if (!reads.empty()) {
		std::map<std::string, std::string>::const_iterator b(reads.begin());
		const std::map<std::string, std::string>::const_iterator end_b(reads.end());
		for (; b != end_b; ++b) {
			out << b->first << "\n";
		}
	}
	if (!skipped.empty()) {
		std::set<std::string>::const_iterator b(skipped.begin());
		const std::set<std::string>::const_iterator end_b(skipped.end());
		for (; b != end_b; ++b) {
			out << *b << "\n";
		}
	}
	out.close();
}

static void bin_files(const size_t files, const std::multimap<size_t, std::string> &sizes, const size_t max, std::multimap<size_t, std::list<std::string> > &bins) {
	std::multimap<size_t, std::string>::const_reverse_iterator a(sizes.rbegin());
	for (size_t i(0); i != files; ++i, ++a) {
		bins.insert(make_pair(a->first, std::list<std::string>(1, a->second)));
	}
	const std::multimap<size_t, std::string>::const_reverse_iterator end_a(sizes.rend());
	for (; a != end_a; ++a) {
		const std::multimap<size_t, std::list<std::string> >::iterator b(bins.begin());
		// bad binning; try again with more files
		if (b->first + a->first > max) {
			bins.clear();
			bin_files(files + 1, sizes, max, bins);
			return;
		}
		// add new entry with increased size and extra file
		b->second.push_back(a->second);
		bins.insert(make_pair(b->first + a->first, b->second));
		// delete old entry
		bins.erase(b);
	}
}

static double score_bin(const std::multimap<size_t, std::list<std::string> > &bins, unsigned long long total) {
	double avg_score(total / bins.size());
	double score(0);
	std::multimap<size_t, std::list<std::string> >::const_iterator a(bins.begin());
	const std::multimap<size_t, std::list<std::string> >::const_iterator end_a(bins.end());
	for (; a != end_a; ++a) {
		const double x(avg_score - a->first);
		score += x * x;
	}
	return score;
}

// group dbs into sets that are as even as possible, and less than MAX_DB_SIZE
static size_t combine_dbs(const std::string &tmp_dir, const std::multimap<size_t, std::string> &q_list, const int min_dbs) {
	unsigned long long total_size(0);
	std::multimap<size_t, std::string>::const_iterator a(q_list.begin());
	const std::multimap<size_t, std::string>::const_iterator end_a(q_list.end());
	for (; a != end_a; ++a) {
		total_size += a->first;
	}
	// smallest number of files possible
	size_t files(size_t(total_size / MAX_DB_SIZE));
	// keep files a multiple of min_dbs
	if (files % min_dbs) {
		files += min_dbs - (files % min_dbs);
	}
	size_t size;	// largest size that gives the right number of files
	if (files != 1) {
		size = size_t(total_size / (files - 1) - 1);
	} else {
		size = total_size;
	}
	std::multimap<size_t, std::list<std::string> > best_output_set;
	bin_files(files, q_list, size, best_output_set);
	double best_score(score_bin(best_output_set, total_size));
	// smallest size that gives the right number of files;
	const size_t diff((size - total_size / (files + 1) + 1) / 10);
	int i(0);
	for (size -= diff; i != 10; ++i, size -= diff) {
		std::multimap<size_t, std::list<std::string> > output_set;
		bin_files(files, q_list, size, output_set);
		const double score(score_bin(output_set, total_size));
		if (best_score > score && (output_set.size() % min_dbs == 0 || best_output_set.size() % min_dbs != 0)) {
			best_score = score;
			best_output_set = output_set;
		}
	}
	size_t filenum(0);
	std::ofstream out;
	out.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	std::multimap<size_t, std::list<std::string> >::const_iterator b(best_output_set.begin());
	const std::multimap<size_t, std::list<std::string> >::const_iterator end_b(best_output_set.end());
	for (; b != end_b; ++b) {
		const std::string filename(tmp_dir + "/db" + itoa(filenum));
		++filenum;
		out.open(filename.c_str());
		std::list<std::string>::const_iterator c(b->second.begin());
		const std::list<std::string>::const_iterator end_c(b->second.end());
		for (; c != end_c; ++c) {
			out << *c << "\n";
		}
		out.close();
	}
	return filenum;
}

static void print_exclude_log(const std::string &file, const std::map<std::string, int> &exclude_list, const std::vector<std::string> &exclude_names) {
	const size_t end_i(exclude_names.size());
	std::vector<size_t> count(0, end_i);
	std::map<std::string, int>::const_iterator a(exclude_list.begin());
	const std::map<std::string, int>::const_iterator end_a(exclude_list.end());
	for (; a != end_a; ++a) {
		if (a->second > 0) {
			++count[a->second - 1];
		}
	}
	std::ofstream out(file.c_str());
	for (size_t i(0); i != end_i; ++i) {
		out << exclude_names[i] << " " << count[i] << "\n";
	}
	out.close();
}

static void copy_dbs(const std::string &run_dir, const std::string &tmp_dir, std::list<std::string> &files, const std::vector<std::string> &exclude_names, const size_t opt_length_cutoff, const std::string &opt_linker, const int min_dbs) {
	if (mkdir(tmp_dir.c_str(), 0777) == -1) {
		throw LocalException("could not create temporary directory " + tmp_dir + ": " + strerror(errno));
	}
	std::map<std::string, int> exclude_list;
	read_excludes(exclude_names, exclude_list);
	expand_included_files(files);
	std::multimap<size_t, std::string> q_list;
	split_db(run_dir, tmp_dir, files, exclude_list, opt_length_cutoff, opt_linker, q_list);
	const size_t db_files(combine_dbs(tmp_dir, q_list, min_dbs));
	std::ofstream out;
	out.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	out.open((run_dir + "/db_list").c_str());
	out << "q" << q_list.size() << "\n";
	out << "db" << db_files << "\n";
	out.close();
	print_exclude_log(run_dir + "/exclude_count", exclude_list, exclude_names);
}

int main(const int argc, char * const argv[]) {
	int had_error(0);
	try {
		set_pwd();
		int min_dbs;
		size_t opt_length_cutoff;
		std::string run_dir, opt_linker;
		std::list<std::string> db_list;
		std::vector<std::string> exclude_names;
		if (get_opts(argc, argv, db_list, exclude_names, run_dir, opt_length_cutoff, opt_linker, min_dbs)) {
			return 0;
		}
		const std::string tmp_dir(get_tmp_dir(run_dir));
		copy_dbs(run_dir, tmp_dir, db_list, exclude_names, opt_length_cutoff, opt_linker, min_dbs);
	} catch (std::exception &e) {
		std::cerr << "Error: " << e.what() << "\n";
		const LocalException * const x(dynamic_cast<LocalException *>(&e));
		if (x != 0 && x->show_usage()) {
			print_usage();
		}
		had_error = 1;
	}
	return had_error;
}
