#include "open_compressed.h"	// get_suffix()
#include "write_fork.h"
#include <cassert>	// assert()
#include <errno.h>	// errno
#include <fcntl.h>	// O_CREATE, O_TRUNC, O_WRONLY, open()
#include <iostream>	// cerr, cout
#include <list>		// list<>
#include <map>		// map<>
#include <new>		// new
#include <string.h>	// strerror()
#include <string>	// string
#include <sys/types.h>	// mode_t, pid_t, size_t, ssize_t
#include <sys/wait.h>	// WNOHANG, waitpid()
#include <unistd.h>	// _SC_OPEN_MAX, _exit(), STDOUT_FILENO, close(), dup2(), execvp(), fork(), pipe(), sysconf, write()

class WriteForkLocalData {
    private:
	// map of open files to forked process id's
	std::map<int, pid_t> m_open_processes;
	// list of closed processes that need to be waited on
	std::list<pid_t> m_closed_processes;
	void finish_nohang(void) {
		std::list<pid_t>::iterator a(m_closed_processes.begin());
		const std::list<pid_t>::iterator end_a(m_closed_processes.end());
		while (a != end_a) {
			if (waitpid(*a, 0, WNOHANG) > 0) {
				a = m_closed_processes.erase(a);
			} else {
				++a;
			}
		}
	}
    public:
	const long open_max;
	WriteForkLocalData(void) : m_open_processes(), m_closed_processes(), open_max(sysconf(_SC_OPEN_MAX)) {
		assert(open_max != 0);
	}
	void add_open(const int i, const pid_t j) {
		m_open_processes[i] = j;
	}
	void close_process(const int i) {
		assert(-1 < i && i < open_max);
		close(i);
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
			const std::map<int, pid_t>::iterator a(m_open_processes.find(i));
			if (a != m_open_processes.end()) {
				waitpid(a->second, 0, 0);
				m_open_processes.erase(a);
			}
			finish_nohang();
		}
	}
	~WriteForkLocalData(void) {
		std::map<int, pid_t>::const_iterator a(m_open_processes.begin());
		const std::map<int, pid_t>::const_iterator end_a(m_open_processes.end());
		for (; a != end_a; ++a) {
			close(a->first);
			m_closed_processes.push_back(a->second);
		}
		std::list<pid_t>::const_iterator b(m_closed_processes.begin());
		const std::list<pid_t>::const_iterator end_b(m_closed_processes.end());
		for (; b != end_b; ++b) {
			waitpid(*b, 0, 0);
		}
	}
};

static WriteForkLocalData local;

// pipe output through another program into a file

int write_fork(const std::list<std::string> &args, const std::string &filename, const mode_t mode) {
	int fd(-1);
	if (args.empty()) {			// direct write, no pipe
		if (filename.empty() || filename == "-") { // it's just stdout
			return STDOUT_FILENO;
		}
		fd = open(filename.c_str(), O_WRONLY | O_CREAT | O_TRUNC, mode);
		if (fd == -1) {
			std::cerr << "Error: open: " << strerror(errno) << '\n';
		}
		return fd;
	}
	pid_t pid;
	int pipefd[2];
	if (pipe(pipefd) == -1) {
		std::cerr << "Error: pipe: " << strerror(errno) << '\n';
		return -1;
	} else if (pipefd[0] >= local.open_max) {	// too many open files
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
		close(pipefd[1]);
		if (dup2(pipefd[0], 0) == -1) {		// set up stdin (pipe)
			std::cerr << "Error: dup2: " << strerror(errno) << '\n';
			_exit(1);
		}
		close(pipefd[0]);
		for (int j = 3; j < local.open_max; ++j) {
			close(j);
		}
		// convert vector to array of char *
		char **argv(new char *[args.size() + 1]);
		size_t i(0);
		std::list<std::string>::const_iterator a(args.begin());
		const std::list<std::string>::const_iterator end_a(args.end());
		for (; a != end_a; ++a, ++i) {
			argv[i] = const_cast<char *>(a->c_str());
		}
		argv[i] = 0;
		if (filename.empty() || filename == "-") {
			fd = 1;		// stdout is already what we want
		} else {		// set up stdout (file)
			close(1);
			fd = open(filename.c_str(), O_WRONLY | O_CREAT | O_TRUNC, mode);
			if (fd == -1) {
				std::cerr << "Error: open: " << strerror(errno) << '\n';
			} else if (fd != 1) {
				if (dup2(fd, 1) == -1) {
					std::cerr << "Error: dup2: " << strerror(errno) << '\n';
					close(fd);
					fd = -1;
				} else {
					close(fd);
					fd = 1;
				}
			}
		}
		if (fd == 1 && execvp(argv[0], argv) == -1) {
			std::cerr << "Error: execvp: " << strerror(errno) << '\n';
		}
		_exit(1);
	} else {		// parent
		close(pipefd[0]);
		fd = pipefd[1];
		local.add_open(fd, pid);
	}
	return fd;
}

// try and guess if we should pipe this through a compressor
int write_fork(const std::string &filename, const mode_t mode) {
	std::string suffix;
	get_suffix(filename, suffix);
	std::list<std::string> args;
	if (suffix == ".gz") {
		args.push_back("gzip");
		args.push_back("-c");
	} else if (suffix == ".bz2") {
		args.push_back("bzip2");
		args.push_back("-c");
	} else if (suffix == ".Z") {
		args.push_back("compress");
		args.push_back("-c");
	}
	return write_fork(args, filename, mode);
}

// close the file and wait on the forked process

void close_fork(const int fd) {
	local.close_process(fd);
}

// wait for all remaining writes to finish

void close_fork_wait(const int fd) {
	local.close_process_wait(fd);
}

// write one character to a file descriptor; return -1 on error, 1 on success
ssize_t pfputc(const int fd, const char c) {
	assert(-1 < fd && fd < local.open_max);
	if (write(fd, &c, 1) == -1) {
		std::cerr << "Error: write(" << fd << ' ' << c << "): " << strerror(errno) << '\n';
		return -1;
	}
	return 1;
}

// write output to a file descriptor; returns -1 on error, otherwise
// number of characters written

ssize_t pfputs(const int fd, const std::string &line) {
	assert(-1 < fd && fd < local.open_max);
	const char *buf(line.c_str());
	ssize_t left_to_write(line.size());
	while (left_to_write != 0) {
		const ssize_t written(write(fd, buf, left_to_write));
		if (written == -1) {
			std::cerr << "Error: write(" << fd << ' ' << line << "): " << strerror(errno) << '\n';
			return -1;
		}
		left_to_write -= written;
		buf += written;
	}
	return line.size();
}

ssize_t pfwrite(const int fd, const void * const ptr, const size_t size) {
	assert(-1 < fd && fd < local.open_max);
	const char *buf(static_cast<const char *>(ptr));
	ssize_t i(size);
	while (i != 0) {
		const ssize_t j(write(fd, buf, i));
		if (j == -1) {
			std::cerr << "Error: write (" << fd << ' ' << ptr << "): " << strerror(errno) << '\n';
			return -1;
		}
		i -= j;
		buf += j;
	}
	return size;
}
