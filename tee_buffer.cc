// like tee, except it fills a buffer first, starts up processes it'll
// write to, and then writes from the buffer while also continuing to
// accept input

#include "breakup_line.h"	// breakup_line(), breakup_line_quoted()
#include <algorithm>	// min()
#include <errno.h>	// errno
#include <exception>	// exception
#include <fcntl.h>      // O_CREATE, O_TRUNC, O_WRONLY, open()
#include <getopt.h>	// getopt(), optarg, optind
#include <iostream>	// cerr
#include <signal.h>	// SIGUSR1, kill(), sigset_t, sigaction(), sigaddset(), sigemptyset(), sigwait()
#include <sstream>	// istringstream
#include <stdio.h>
#include <stdio.h>	// EOF
#include <string.h>	// strerror()
#include <string>	// string
#include <sys/stat.h>	// S_IRUSR, S_IWUSR, S_IRGRP, S_IWGRP, S_IROTH, S_IWOTH
#include <sys/time.h>	// FD_SET(), fd_set, pselect()
#include <sys/types.h>	// pid_t, ssize_t
#include <sys/wait.h>	// wait()
#include <unistd.h>	// STDIN_FILENO, STDOUT_FILENO, _exit(), close(), dup2(), execvp(), fork(), pipe(), read(), write()
#include <vector>	// vector<>

class LocalException : public std::exception {
    private:
	const std::string error_;
	const int show_usage_;
    public:
	explicit LocalException(const std::string &s) : error_(s), show_usage_(0) { }
	explicit LocalException(const std::string &s, const int i) : error_(s), show_usage_(i) { }
	virtual ~LocalException(void) throw() { }
	virtual const char * what() const throw() {
		return error_.c_str();
	}
	const int &show_usage(void) const throw() {
		return show_usage_;
	}
};

class Buffer {
    private:
	std::vector<char> buffer_;
	size_t cycle_size_;
	size_t read_offset_, write_offset_;	// write to buffer, read from buffer
	int was_filled_;
    public:
	explicit Buffer(const size_t i, const size_t j) : buffer_(i), cycle_size_(j), read_offset_(0), write_offset_(0), was_filled_(0) {
		// fill buffer on initialization
		size_t n(buffer_.size());
		while (n != 0) {
			const ssize_t k(read(STDIN_FILENO, &buffer_[write_offset_], n));
			if (k <= 0) {
				if (k == -1) {
					throw LocalException("read(stdin): " + std::string(strerror(errno)));
				}
				return;
			}
			write_offset_ += k;
			n -= k;
		}
		write_offset_ = 0;		// handle wrap
		was_filled_ = 1;
	}
	~Buffer() { }
	int was_filled() const {
		return was_filled_;
	}
    private:
	// send exactly n bytes to fd_list (need this to keep output in sync)
	void write_exactly(const std::vector<int> &fd_list, const size_t n) {
		const char * const buf_start(&buffer_[read_offset_]);
		std::vector<int>::const_iterator a(fd_list.begin());
		const std::vector<int>::const_iterator end_a(fd_list.end());
		for (; a != end_a; ++a) {
			const char *buf(buf_start);
			for (size_t j(n); j != 0;) {
				const ssize_t i(write(*a, buf, j));
				if (i == -1) {
					throw LocalException("write: " + std::string(strerror(errno)));
				}
				j -= i;
				buf += i;
			}
		}
		if ((read_offset_ += n) == buffer_.size()) {
			read_offset_ = 0;
		}
	}
    public:
	void loop(const std::vector<int> &fd_list) {
		// use this for checking for read blocking
		const struct timespec timeout({0, 0});
		fd_set stdin_fd;
		FD_SET(STDIN_FILENO, &stdin_fd);
		// since we start filled, start with a write, then read
		for (;;) {
			// limit write sizes to prevent read pipe from filling
			// (need <= to handle full buffer condition)
			const size_t n(std::min((write_offset_ <= read_offset_ ? buffer_.size() : write_offset_) - read_offset_, cycle_size_));
			write_exactly(fd_list, n);
			// if buffer isn't empty, check for blocking
			if (read_offset_ != write_offset_ && pselect(1, &stdin_fd, 0, 0, &timeout, 0) != 1) {
				// need to reset, as it was cleared
				FD_SET(STDIN_FILENO, &stdin_fd);
				// nothing available, so keep going through buffer
				continue;
			}
			const size_t buffer_left((write_offset_ < read_offset_ ? read_offset_ : buffer_.size()) - write_offset_);
			const ssize_t i(read(STDIN_FILENO, &buffer_[write_offset_], buffer_left));
			if (i <= 0) {
				if (i == -1) {
					throw LocalException("read(stdin): " + std::string(strerror(errno)));
				}
				return;
			}
			if ((write_offset_ += i) == buffer_.size()) {
				write_offset_ = 0;
			}
		}
	}
	void empty(const std::vector<int> &fd_list) {
		if (read_offset_ == write_offset_) {		// empty
			return;
		}
		// first, empty to end of buffer
		const size_t n((write_offset_ < read_offset_ ? buffer_.size() : write_offset_) - read_offset_);
		write_exactly(fd_list, n);
		// next, empty anything at beginning of buffer
		if (read_offset_ != write_offset_) {
			write_exactly(fd_list, write_offset_);
		}
	}
};

static void print_usage() {
	std::cerr <<
		"usage: tee [opts] <file1> [<file2> ...]\n"
		"    -b ##  buffer size [1mb]\n"
		"    -c ##  buffer cycle size [32kb]\n"
		"    -n     don't write to stdout\n";
}

static void get_opts(int argc, char **argv, size_t &buffer_size, size_t &buffer_cycle_size, std::vector<int> &fd_list) {
	size_t x;
	int write_stdout(1);
	int c;
	while ((c = getopt(argc, argv, "b:c:hn")) != EOF) {
		switch (c) {
		    case 'b':
			std::istringstream(optarg) >> x;
			if (x != 0) {			// if zero, use default
				buffer_size = x;
			}
			break;
		    case 'c':
			std::istringstream(optarg) >> x;
			if (x != 0) {			// if zero, use default
				buffer_cycle_size = x;
			}
			break;
		    case 'h':
			throw LocalException("", 1);
			break;
		    case 'n':
			write_stdout = 0;
			break;
		    default:
			throw LocalException("bad option: " + static_cast<char>(c), 1);
		}
	}
	if (buffer_cycle_size > buffer_size) {
		buffer_cycle_size = buffer_size;
	}
	if (write_stdout) {
		fd_list.push_back(STDOUT_FILENO);
	} else {
		close(STDOUT_FILENO);
	}
}

// used to signal child that buffer has been filled
static int caught_signal(0);
static void handle_signal(const int signal) {
	(void)signal;		// prevent unused warning
	caught_signal = 1;
}

// list1[1] is the output file to write to (if any)
// list2 contains a list of pipeline commands
static void start_child(const int pipe_in, const std::vector<std::string> &list1, const std::vector<std::string> &list2) {
	if (dup2(pipe_in, STDIN_FILENO) == -1) {			// set up stdin (pipe)
		throw LocalException("dup2: " + std::string(strerror(errno)));
	}
	close(pipe_in);
	// convert list2 into format suitable for execvp
	// args_list hold the actrual data, argv_list points to it
	std::vector<std::vector<std::string> > args_list(list2.size());
	std::vector<std::vector<char *> > argv_list(list2.size());
	for (size_t i(0); i != list2.size(); ++i) {
		std::vector<std::string> &args(args_list[i]);
		std::vector<char *> &argv(argv_list[i]);
		breakup_line_quoted(list2[i], args);
		std::vector<std::string>::const_iterator a(args.begin());
		const std::vector<std::string>::const_iterator end_a(args.end());
		for (; a != end_a; ++a) {
			argv.push_back(const_cast<char *>(a->c_str()));
		}
		argv.push_back(0);
	}
	// wait for signal to start (sent when buffer is filled)
	sigset_t wait_signal;
	if (sigemptyset(&wait_signal) == -1) {
		throw LocalException("sigemptyset: " + std::string(strerror(errno)));
	}
	if (sigaddset(&wait_signal, SIGUSR1) == -1) {
		throw LocalException("sigaddset: " + std::string(strerror(errno)));
	}
	int sig;
	if (!caught_signal && sigwait(&wait_signal, &sig) != 0) {
		throw LocalException("sigwait");
	}
	// start up pipeline
	for (size_t i(0); i != argv_list.size() - 1; ++i) {
		pid_t pid;
		int pipefd[2];
		if (pipe(pipefd) == -1) {
			throw LocalException("pipe: " + std::string(strerror(errno)));
		} else if ((pid = fork()) == -1) {
			close(pipefd[0]);
			close(pipefd[1]);
			throw LocalException("fork: " + std::string(strerror(errno)));
		} else if (pid == 0) {					// child
			close(pipefd[1]);
			if (dup2(pipefd[0], STDIN_FILENO) == -1) {	// set up stdin
				throw LocalException("dup2: " + std::string(strerror(errno)));
			}
			close(pipefd[0]);
		} else {						// parent
			close(pipefd[0]);
			if (dup2(pipefd[1], STDOUT_FILENO) == -1) {	// set up stdout
				throw LocalException("dup2: " + std::string(strerror(errno)));
			}
			if (execvp(argv_list[i][0], argv_list[i].data()) == -1) {
				throw LocalException("execvp: " + std::string(strerror(errno)));
			}
		}
	}
	// set up stdout for final segment of pipeline
	if (list1.size() == 2 && !list1[1].empty() && list1[1] != "-") {	// not stdout
		close(STDOUT_FILENO);
		const int fd(open(list1[1].c_str(), O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH));
		if (fd == -1) {
			throw LocalException("open: " + list1[1] + ": " + std::string(strerror(errno)));
		}
		if (dup2(fd, STDOUT_FILENO) == -1) {
			close(fd);
			throw LocalException("dup2: " + std::string(strerror(errno)));
		}
	}
	if (execvp(argv_list.back()[0], argv_list.back().data()) == -1) {
		throw LocalException("execvp: " + std::string(strerror(errno)));
	}
}

// children wait for signal before starting the pipeline
static int spawn_outputs(const std::string &command, std::vector<pid_t> &children) {
	struct sigaction act;
	act.sa_handler = handle_signal;
	if (sigemptyset(&act.sa_mask) == -1) {
		throw LocalException("sigemptyset: " + std::string(strerror(errno)));
	}
	if (sigaddset(&act.sa_mask, SIGUSR1) == -1) {
		throw LocalException("sigaddset: " + std::string(strerror(errno)));
	}
	act.sa_flags = 0;
	sigaction(SIGUSR1, &act, 0);
	// segment command into separate pipes and ending file
	std::vector<std::string> list1, list2;
	breakup_line(command, list1, '>');
	if (list1.size() > 2) {
		throw LocalException("bad command: multiple > in pipeline: " + command);
	}
	breakup_line(list1.front(), list2, '|');
	if (list1.size() == 1 && list2.size() == 1) {		// simple file
		if (command.empty() || command == "-") {	// it's just stdout
			return STDOUT_FILENO;
		}
		const int fd(open(list1[0].c_str(), O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH));
		if (fd == -1) {
			throw LocalException("open: " + command + ": " + std::string(strerror(errno)));
		}
		return fd;
	}
	pid_t pid;
	int pipefd[2];
	if (pipe(pipefd) == -1) {
		throw LocalException("pipe: " + std::string(strerror(errno)));
	} else if ((pid = fork()) == -1) {
		close(pipefd[0]);
		close(pipefd[1]);
		throw LocalException("fork: " + std::string(strerror(errno)));
	} else if (pid == 0) {  // child
		try {
			close(pipefd[1]);
			start_child(pipefd[0], list1, list2);	// doesn't return
		} catch (std::exception &e) {
			if (e.what()[0] != 0) {
				std::cerr << "Error: " << e.what() << "\n";
			}
		}
		_exit(1);
	} else {		// parent
		close(pipefd[0]);
		children.push_back(pid);
		return pipefd[1];
	}
}

int main(int argc, char **argv) {
	int had_error(0);
	try {
		size_t buffer_size(1 << 24);
		size_t buffer_cycle_size(1 << 15);
		std::vector<int> fd_list;
		get_opts(argc, argv, buffer_size, buffer_cycle_size, fd_list);
		// fork off outputs before we have a large memory footprint
		std::vector<pid_t> children;
		for (; optind != argc; ++optind) {
			fd_list.push_back(spawn_outputs(argv[optind], children));
		}
		Buffer buffer(buffer_size, buffer_cycle_size);
		// now that buffer is filled, send outputs signal to start
		std::vector<pid_t>::const_iterator b(children.begin());
		const std::vector<pid_t>::const_iterator end_b(children.end());
		for (; b != end_b; ++b) {
			kill(*b, SIGUSR1);
		}
		if (buffer.was_filled()) {
			buffer.loop(fd_list);
		}
		buffer.empty(fd_list);
		// close outputs
		std::vector<int>::const_iterator a(fd_list.begin());
		const std::vector<int>::const_iterator end_a(fd_list.end());
		for (; a != end_a; ++a) {
			close(*a);
		}
		// now wait for them all to finish
		while (wait(0) != -1) { }
	} catch (std::exception &e) {
		if (e.what()[0] != 0) {
			std::cerr << "Error: " << e.what() << "\n";
		}
		LocalException * const x(dynamic_cast<LocalException *>(&e));
		if (x != NULL && x->show_usage()) {
			print_usage();
		}
		had_error = 1;
	}
	return had_error;
}
