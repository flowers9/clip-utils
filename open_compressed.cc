#include "breakup_line.h"	// breakup_line()
#include "open_compressed.h"
#include "refcount_array.h"	// refcount_array<>
#include <cassert>	// assert()
#include <errno.h>	// EINVAL, EISDIR, ENFILE, ENOENT, errno
#include <fcntl.h>	// O_RDONLY, open()
#include <iostream>	// cerr
#include <list>		// list<>
#include <map>		// map<>
#include <stdlib.h>	// getenv()
#include <string.h>	// memcpy(), memmove(), strerror()
#include <string>	// string
#include <sys/stat.h>	// S_IFDIR, stat(), struct stat
#include <sys/types.h>	// pid_t, size_t, ssize_t
#include <sys/wait.h>	// WNOHANG, waitpid()
#include <unistd.h>	// _SC_OPEN_MAX, _exit(), close(), dup2(), execlp(), fork(), pipe(), read(), sysconf()
#include <vector>	// vector<>

// XXX - consider switching to popen instead of fork/exec

#define BUFSIZE 32768
#define BZIP2_COMMAND_DEFAULT "bzip2"
#define GZIP_COMMAND_DEFAULT "gzip"
#define XZ_COMMAND_DEFAULT "xz"

class OpenCompressedLocalData {
    private:
	bool m_already_closed_stdin;	// prevent reuse of closed stream
	// map of open files to decompression process id's
	std::map<int, pid_t> m_open_processes;
	// list of closed processes that need to be waited on
	std::list<pid_t> m_closed_processes;
	std::vector<std::string> m_gzip, m_bzip2, m_xz;
    private:
	void finish_nohang(void) {
		std::list<pid_t>::iterator a(m_closed_processes.begin());
		const std::list<pid_t>::const_iterator end_a(m_closed_processes.end());
		while (a != end_a) {
			if (waitpid(*a, 0, WNOHANG) > 0) {
				a = m_closed_processes.erase(a);
			} else {
				++a;
			}
		}
	}
    public:
	const ssize_t open_max;
	// internal file buffers - has to be a raw array, sadly
	refcount_array<char> * const buffers;
	ssize_t * const buffer_start;
	ssize_t * const buffer_length;
    public:
	OpenCompressedLocalData(void) :
		m_already_closed_stdin(0),
		m_open_processes(),
		m_closed_processes(),
		m_gzip(), m_bzip2(), m_xz(),
		open_max(sysconf(_SC_OPEN_MAX)),
		buffers(new refcount_array<char>[open_max]),
		buffer_start(new ssize_t[open_max]),
		buffer_length(new ssize_t[open_max])
	{
		assert(open_max);
		assert(buffers);
		assert(buffer_start);
		assert(buffer_length);
		// zero may get reused as something else, but always starts as stdin
		buffer_start[0] = buffer_length[0] = 0;
	}
	~OpenCompressedLocalData(void) {
		std::map<int, pid_t>::const_iterator a(m_open_processes.begin());
		const std::map<int, pid_t>::const_iterator end_a(m_open_processes.end());
		for (; a != end_a; ++a) {
			close(a->first);
		}
		delete[] buffers;
		delete[] buffer_start;
		delete[] buffer_length;
		std::list<pid_t>::const_iterator b(m_closed_processes.begin());
		const std::list<pid_t>::const_iterator end_b(m_closed_processes.end());
		for (; b != end_b; ++b) {
			waitpid(*b, 0, 0);
		}
		for (a = m_open_processes.begin(); a != end_a; ++a) {
			waitpid(a->second, 0, 0);
		}
	}
	void add_open(const int i, const pid_t j) {
		m_open_processes[i] = j;
	}
	void close_process(const int i) {
		assert(-1 < i && i < open_max);
		close(i);
		if (i == 0) {
			m_already_closed_stdin = 1;
		}
		buffers[i].resize(0);
		const std::map<int, pid_t>::iterator a(m_open_processes.find(i));
		if (a != m_open_processes.end()) {
			m_closed_processes.push_back(a->second);
			m_open_processes.erase(a);
		}
		finish_nohang();
	}
	void close_process_wait(const int i) {
		assert(-2 < i && i < open_max);
		if (i == -1) {		// wait for all closed processes
			std::list<pid_t>::const_iterator b(m_closed_processes.begin());
			const std::list<pid_t>::const_iterator end_b(m_closed_processes.end());
			for (; b != end_b; ++b) {
				waitpid(*b, 0, 0);
			}
			m_closed_processes.clear();
		} else {
			close(i);
			if (i == 0) {
				m_already_closed_stdin = 1;
			}
			buffers[i].resize(0);
			const std::map<int, pid_t>::iterator a(m_open_processes.find(i));
			if (a != m_open_processes.end()) {
				waitpid(a->second, 0, 0);
				m_open_processes.erase(a);
			}
			finish_nohang();
		}
	}
	const bool &already_closed_stdin(void) const {
		return m_already_closed_stdin;
	}
	// return array suitable to passing to execvp()
	char **gzip_args(const std::string &);
	char **bzip2_args(const std::string &);
	char **xz_args(const std::string &);
    private:
	OpenCompressedLocalData(const OpenCompressedLocalData &);
	OpenCompressedLocalData &operator=(const OpenCompressedLocalData &);
};

char **OpenCompressedLocalData::gzip_args(const std::string &file) {
	if (m_gzip.empty()) {
		char * const env_cmd(getenv("GZIP_COMMAND"));
		breakup_line(env_cmd ? env_cmd : GZIP_COMMAND_DEFAULT, m_gzip);
		m_gzip.push_back("-d");		// an extra -d or -c doesn't hurt
		m_gzip.push_back("-c");
	}
	m_gzip.push_back(file);
	size_t n(0);
	for (size_t i(0); i != m_gzip.size(); ++i) {
		n += m_gzip[i].size() + 1;	// +1 for trailing null
	}
	// this would be a memory leak, but we only ever allocate
	// it immediately before calling execvp()
	char *m_gzip_string(new char[n]);
	char **m_gzip_args(new char*[m_gzip.size() + 1]);
	char *s(m_gzip_string);
	for (size_t i(0); i != m_gzip.size(); ++i) {
		memcpy(s, m_gzip[i].c_str(), m_gzip[i].size() + 1);
		m_gzip_args[i] = s;
		s += m_gzip[i].size() + 1;
	}
	m_gzip_args[m_gzip.size()] = 0;	// trailing null
	m_gzip.pop_back();		// remove file
	return m_gzip_args;
}

char **OpenCompressedLocalData::bzip2_args(const std::string &file) {
	if (m_bzip2.empty()) {
		char * const env_cmd(getenv("BZIP2_COMMAND"));
		breakup_line(env_cmd ? env_cmd : BZIP2_COMMAND_DEFAULT, m_bzip2);
		m_bzip2.push_back("-d");		// an extra -d or -c doesn't hurt
		m_bzip2.push_back("-c");
	}
	m_bzip2.push_back(file);
	size_t n(0);
	for (size_t i(0); i != m_bzip2.size(); ++i) {
		n += m_bzip2[i].size() + 1;	// +1 for trailing null
	}
	// this would be a memory leak, but we only ever allocate
	// it immediately before calling execvp()
	char *m_bzip2_string(new char[n]);
	char **m_bzip2_args(new char*[m_bzip2.size() + 1]);
	char *s(m_bzip2_string);
	for (size_t i(0); i != m_bzip2.size(); ++i) {
		memcpy(s, m_bzip2[i].c_str(), m_bzip2[i].size() + 1);
		m_bzip2_args[i] = s;
		s += m_bzip2[i].size() + 1;
	}
	m_bzip2_args[m_bzip2.size()] = 0;	// trailing null
	m_bzip2.pop_back();		// remove file
	return m_bzip2_args;
}

char **OpenCompressedLocalData::xz_args(const std::string &file) {
	if (m_xz.empty()) {
		char * const env_cmd(getenv("XZ_COMMAND"));
		breakup_line(env_cmd ? env_cmd : XZ_COMMAND_DEFAULT, m_xz);
		m_xz.push_back("-d");		// an extra -d or -c doesn't hurt
		m_xz.push_back("-c");
	}
	m_xz.push_back(file);
	size_t n(0);
	for (size_t i(0); i != m_xz.size(); ++i) {
		n += m_xz[i].size() + 1;	// +1 for trailing null
	}
	// this would be a memory leak, but we only ever allocate
	// it immediately before calling execvp()
	char *m_xz_string(new char[n]);
	char **m_xz_args(new char*[m_xz.size() + 1]);
	char *s(m_xz_string);
	for (size_t i(0); i != m_xz.size(); ++i) {
		memcpy(s, m_xz[i].c_str(), m_xz[i].size() + 1);
		m_xz_args[i] = s;
		s += m_xz[i].size() + 1;
	}
	m_xz_args[m_xz.size()] = 0;	// trailing null
	m_xz.pop_back();		// remove file
	return m_xz_args;
}

static OpenCompressedLocalData local;

// simply return the suffix of the file name, if it matches one of the list

void get_suffix(const std::string &filename, std::string &suffix) {
	const char *suffix_list[4] = { ".gz", ".bz2", ".xz", ".Z" };
	const size_t suffix_size[4] = { 3, 4, 3, 2 };
	suffix.clear();
	for (int i(0); i != 4; ++i) {
		const size_t &j(suffix_size[i]);
		if (filename.size() > j && filename.compare(filename.size() - j, j, suffix_list[i]) == 0) {
			suffix = suffix_list[i];
			break;
		}
	}
}

// returns empty, .Z, .gz, .xz, or .bz2; checks to see if the filename
// ends in any given suffix, if not checks to see if a file with the given
// suffix exists

int find_suffix(std::string &filename, std::string &suffix) {
	const char *suffix_list[4] = { ".gz", ".bz2", ".xz", ".Z" };
	const size_t suffix_size[4] = { 3, 4, 3, 2 };
	suffix.clear();
	for (int i(0); i != 4; ++i) {
		const size_t &j(suffix_size[i]);
		if (filename.size() > j && filename.compare(filename.size() - j, j, suffix_list[i]) == 0) {
			suffix = suffix_list[i];
			break;
		}
	}
	struct stat buf;
	if (stat(filename.c_str(), &buf) == 0) {
		// only open regular files
		if ((buf.st_mode & S_IFDIR) != 0) {
			errno = EISDIR;
			return -1;
		} else {
			return 0;
		}
	} else if (errno != ENOENT) {
		std::cerr << "Error: stat: " << filename << ": " << strerror(errno) << '\n';
		return -1;
	} else if (!suffix.empty()) {
		return -1;
	}
	for (int i(0); i != 4; ++i) {
		const std::string s(filename + suffix_list[i]);
		if (stat(s.c_str(), &buf) == -1) {
			if (errno != ENOENT) {
				std::cerr << "Error: stat: " << filename << ": " << strerror(errno) << '\n';
				return -1;
			}
		} else if ((buf.st_mode & S_IFDIR) == 0) {
			suffix = suffix_list[i];
			filename += suffix;
			return 0;
		}
	}
	errno = ENOENT;
	return -1;
}

// open a file, with gzip/bzip/xz if the filename ends in .gz, .bz2, .xz, or .Z;
// if file is not found, .gz, .bz2, .xz, and .Z are added to end, to see
// if file is compressed

int open_compressed(const std::string &filename, bool force_uncompressed) {
	int fd(-1);
	std::string s(filename);
	std::string suffix;
	// see if file exists
	if (!s.empty() && s.compare("-") != 0 && (force_uncompressed || find_suffix(s, suffix) == -1)) {
		return -1;
	}
	if (!suffix.empty()) {
		pid_t pid;
		int pipefd[2];
		if (pipe(pipefd) == -1) {
			std::cerr << "Error: pipe: " << strerror(errno) << '\n';
			return -1;
		} else if (pipefd[0] >= local.open_max) { // too many open files
			close(pipefd[0]);
			close(pipefd[1]);
			errno = ENFILE;
			std::cerr << "Error: open: " << strerror(errno) << '\n';
			return -1;
		} else if ((pid = fork()) == -1) {
			std::cerr << "Error: fork: " << strerror(errno) << '\n';
			close(pipefd[0]);
			close(pipefd[1]);
			return -1;
		} else if (pid == 0) {	// child
			close(0);
			close(pipefd[0]);
			if (dup2(pipefd[1], 1) == -1) {
				std::cerr << "Error: dup2: " << strerror(errno) << '\n';
				_exit(1);
			}
			close(pipefd[1]);
			// close everything except stdin, stdout, and stderr
			for (int i(3); i < local.open_max; ++i) {
				close(i);
			}
			char **cmd;
			if (suffix == ".bz2") {
				cmd = local.bzip2_args(s);
			} else if (suffix == ".xz") {
				cmd = local.xz_args(s);
			} else {
				cmd = local.gzip_args(s);
			}
			if (execvp(cmd[0], cmd) == -1) {
				std::cerr << "Error: execvp " << cmd[0] << ": " << strerror(errno) << '\n';
			}
			_exit(1);
		} else {		// parent
			close(pipefd[1]);
			fd = pipefd[0];
			local.add_open(fd, pid);
		}
	} else if (s.empty() || s.compare("-") == 0) {
		if (local.already_closed_stdin()) {
			return -1;
		}
		fd = 0;
		local.buffers[fd].resize(BUFSIZE);
		// don't reset start/length, to allow rereading of
		// buffered parts of the stream; resize won't change
		// buffer if old size == new size
		return fd;
	} else if ((fd = open(s.c_str(), O_RDONLY)) == -1) {
		std::cerr << "Error: open: " << strerror(errno) << '\n';
		return -1;
	} else if (fd >= local.open_max) {	// too many open files
		close(fd);
		errno = ENFILE;
		std::cerr << "Error: open: " << strerror(errno) << '\n';
		return -1;
	}
	local.buffers[fd].resize(BUFSIZE);
	local.buffer_start[fd] = 0;
	local.buffer_length[fd] = 0;
	return fd;
}

// close the file and wait on the gzip process, if any

void close_compressed(const int fd) {
	local.close_process(fd);
}

// wait for all remaining writes to finish

void close_compressed_wait(const int fd) {
	local.close_process_wait(fd);
}

// read input from a file descriptor - returns data up to end of line or
// end of file (and strips the end of line character); returns -1 on error,
// otherwise number of characters read (minus end of line, if any)

ssize_t pfgets(const int fd, std::string &line, const char delim) {
	if (fd < 0 || local.open_max <= fd) {
		std::cerr << "Error: pfgets: fd out of range: " << fd << '\n';
		return -1;
	}
	if (local.buffers[fd].empty()) {
		std::cerr << "Error: pfgets: buffer unallocated\n";
		return -1;
	}
	line.clear();
	char * const buf(local.buffers[fd].array());
	ssize_t &i(local.buffer_start[fd]);
	ssize_t &j(local.buffer_length[fd]);
	for (;;) {
		ssize_t k;
		for (k = i; k != j && buf[k] != delim; ++k) { }
		line.append(buf + i, k - i);
		if (k != j) {
			i = k + 1;
			return static_cast<ssize_t>(line.size());
		}
		i = 0;
		j = read(fd, buf, BUFSIZE);
		if (j <= 0) {
			if (j == -1) {
				std::cerr << "Error: read(" << fd << "): " << strerror(errno) << '\n';
			}
			j = 0;
			// only return -1 if we don't return anything else
			return line.empty() ? -1 : static_cast<ssize_t>(line.size());
		}
	}
}

ssize_t skip_next_line(const int fd, const char delim) {
	if (fd < 0 || local.open_max <= fd) {
		std::cerr << "Error: pfgets: fd out of range: " << fd << '\n';
		return -1;
	}
	if (local.buffers[fd].empty()) {
		std::cerr << "Error: pfgets: buffer unallocated\n";
		return -1;
	}
	char * const buf(local.buffers[fd].array());
	ssize_t &i(local.buffer_start[fd]);
	ssize_t &j(local.buffer_length[fd]);
	ssize_t amount_read(0);
	for (;;) {
		ssize_t k;
		for (k = i; k != j && buf[k] != delim; ++k) { }
		amount_read += k - i;
		if (k != j) {
			i = k + 1;
			return amount_read;
		}
		i = 0;
		j = read(fd, buf, BUFSIZE);
		if (j <= 0) {
			if (j == -1) {
				std::cerr << "Error: read(" << fd << "): " << strerror(errno) << '\n';
			}
			j = 0;
			// only return -1 if we don't return anything else
			return amount_read == 0 ? -1 : amount_read;
		}
	}
}

// read and discard next size chars from fd

ssize_t skip_next_chars(const int fd, const size_t size) {
	if (fd < 0 || local.open_max <= fd) {
		std::cerr << "Error: skip_next_chars: fd out of range: " << fd << '\n';
		return -1;
	}
	if (local.buffers[fd].empty()) {
		std::cerr << "Error: skip_next_chars: buffer unallocated\n";
		return -1;
	}
	size_t k(size);
	char * const buf(local.buffers[fd].array());
	ssize_t &i(local.buffer_start[fd]);
	ssize_t &j(local.buffer_length[fd]);
	for (;;) {
		const size_t n(j - i);
		if (n >= k) {
			i += k;
			return size;
		}
		k -= n;
		i = 0;
		j = read(fd, buf, BUFSIZE);
		if (j <= 0) {
			if (j == -1) {
				std::cerr << "Error: read(" << fd << "): " << strerror(errno) << '\n';
			}
			j = 0;
			// only return -1 if we don't skip anything
			return k == size ? -1 : static_cast<ssize_t>(size - k);
		}
	}
}

// read up to size bytes from fd and put them into ptr

ssize_t pfread(const int fd, void * const ptr, const size_t size) {
	if (fd < 0 || local.open_max <= fd) {
		std::cerr << "Error: pfread: fd out of range: " << fd << '\n';
		return -1;
	}
	if (local.buffers[fd].empty()) {
		std::cerr << "Error: pfread: buffer unallocated\n";
		return -1;
	}
	char *s(static_cast<char *>(ptr));
	size_t k(size);
	char * const buf(local.buffers[fd].array());
	ssize_t &i(local.buffer_start[fd]);
	ssize_t &j(local.buffer_length[fd]);
	const size_t n(j - i);
	if (n >= k) {		// have enough in buffer
		memcpy(s, buf + i, k);
		i += k;
		return size;
	}
	// start by copying buffer out
	memcpy(s, buf + i, n);
	s += n;
	k -= n;
	i = j = 0;
	do {	// now just read directly into ptr
		const ssize_t m(read(fd, s, k));
		if (m <= 0) {
			if (m == -1) {
				std::cerr << "Error: read(" << fd << "): " << strerror(errno) << '\n';
			}
			// only return -1 if we don't return anything else
			return k == size ? -1 : static_cast<ssize_t>(size - k);
		}
		s += m;
		k -= m;
	} while (k);
	return size;
}

// like pfread, but don't remove anything unread from the buffer
// (however, you can add to it, and move stuff around a bit)

ssize_t pfpeek(const int fd, void * const ptr, size_t size) {
	if (fd < 0 || local.open_max <= fd) {
		std::cerr << "Error: pfpeek: fd out of range: " << fd << '\n';
		return -1;
	}
	if (local.buffers[fd].empty()) {
		std::cerr << "Error: pfpeek: buffer unallocated\n";
		return -1;
	}
	if (size > BUFSIZE) {
		std::cerr << "Error: pfpeek: request for " << size << " bytes, buffer is only " << BUFSIZE << " long\n";
		return -1;
	}
	char * const s(static_cast<char *>(ptr));
	char * const buf(local.buffers[fd].array());
	ssize_t &i(local.buffer_start[fd]);
	ssize_t &j(local.buffer_length[fd]);
	if (size > BUFSIZE - static_cast<size_t>(i)) {	// need to make space
		// move unread section to the front
		j -= i;
		memmove(buf, buf + i, j);
		i = 0;
	}
	for (;;) {
		const size_t n(j - i);
		if (n >= size) {
			memcpy(s, buf + i, size);
			return size;
		}
		// we only need (size - n), but fill the buffer anyway
		const ssize_t k(read(fd, buf + j, BUFSIZE - j));
		if (k <= 0) {
			if (k == -1) {
				std::cerr << "Error: read(" << fd << "): " << strerror(errno) << '\n';
			}
			// only return -1 if we don't return anything else
			if (n == 0) {
				return -1;
			} else {
				memcpy(s, buf + i, n);
				return n;
			}
		}
		j += k;
	}
}
