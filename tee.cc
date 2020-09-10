// like tee, except it fills the buffer first, then starts up any
// process it'll write to (and/or open files)

#include "breakup_line.h"	// breakup_line(), breakup_line_quoted()
#include <algorithm>	// min()
#include <errno.h>	// errno
#include <exception>	// exception
#include <fcntl.h>      // O_CREATE, O_TRUNC, O_WRONLY, open()
#include <getopt.h>	// getopt(), optarg, optind
#include <iostream>	// cerr
#include <netdb.h>	// gethostbyname(), h_errno, struct hostent
#include <sstream>	// istringstream
#include <stdio.h>	// EOF
#include <string.h>	// memcpy(), memset(), strerror()
#include <string>	// string
#include <sys/socket.h>	// AF_INET, SOCK_STREAM, connect(), htons(), socket(), struct socketaddr, struct sockaddr_in
#include <sys/stat.h>	// S_IRUSR, S_IWUSR, S_IRGRP, S_IWGRP, S_IROTH, S_IWOTH
#include <sys/types.h>	// pid_t, ssize_t
#include <sys/wait.h>	// wait()
#include <unistd.h>	// STDIN_FILENO, STDOUT_FILENO, _exit(), close(), dup(), dup2(), execvp(), fork(), pipe(), read(), write()
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
	ssize_t length_;	// amount of buffer currently used
    public:
	explicit Buffer(const size_t i) : buffer_(i), length_(0) {
		// fill buffer on initialization
		size_t n(buffer_.size());
		while (n != 0) {
			const ssize_t k(read(STDIN_FILENO, &buffer_[length_], n));
			if (k <= 0) {
				if (k == -1) {
					throw LocalException("read(stdin): " + std::string(strerror(errno)));
				}
				return;
			}
			length_ += k;
			n -= k;
		}
	}
	~Buffer() { }
	void loop(const std::vector<int> &fd_list) {
		std::vector<int>::const_iterator a;
		const std::vector<int>::const_iterator end_a(fd_list.end());
		for (;;) {
			// start with a write, since we pre-filled
			for (a = fd_list.begin(); a != end_a; ++a) {
				const char *buf(&buffer_[0]);
				ssize_t j(length_);
				for (;;) {
					const ssize_t i(write(*a, buf, j));
					if (i == -1) {
						throw LocalException("write: " + std::string(strerror(errno)));
					}
					if ((j -= i) == 0) {
						break;
					}
					buf += i;
				}
			}
			// just take however much we get from here on
			length_ = read(STDIN_FILENO, &buffer_[0], buffer_.size());
			if (length_ <= 0) {
				if (length_ == -1) {
					throw LocalException("read(stdin): " + std::string(strerror(errno)));
				}
				return;
			}
		}
	}
};

static void print_usage() {
	std::cerr <<
		"usage: tee [opts] <pipeline/file1> [<pipeline/file2> ...]\n" <<
		"    -b ##         buffer size [32kb]\n" <<
		"    -f host:port  send flag once buffer is full \n" <<
		"    -h            print this help\n" <<
		"    -n            don't write to stdout\n";
}

static void get_opts(int argc, char **argv, size_t &buffer_size, std::vector<int> &fd_list, std::string &flag_host, int &flag_port) {
	int write_stdout(1);
	size_t x;
	const char *s;
	int c;
	while ((c = getopt(argc, argv, "b:f:hn")) != EOF) {
		switch (c) {
		    case 'b':
			std::istringstream(optarg) >> x;
			if (x != 0) {			// if zero, use default
				buffer_size = x;
			}
			break;
		    case 'f':
			s = strchr(optarg, ':');
			if (s == 0) {
				throw LocalException("bad -f option: " + std::string(optarg), 1);
			}
			flag_host.assign(optarg, s - optarg);
			std::istringstream(s + 1) >> flag_port;
			if (flag_port < 0) {
				throw LocalException("bad flag port", 1);
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
	if (write_stdout) {
		// stdout gets used for pipelines, so use a duplicate for our output
		const int fd(dup(STDOUT_FILENO));
		if (fd == -1) {
			throw LocalException("dup: stdout: " + std::string(strerror(errno)));
		}
		fd_list.push_back(fd);
	}
}

// convert list into format suitable for execvp
// args_list holds the actual data, argv_list points to it
static void convert_to_argv(const std::vector<std::string> &list, std::vector<std::vector<std::string> > &args_list, std::vector<std::vector<char *> > &argv_list) {
	args_list.clear();
	args_list.resize(list.size());
	for (size_t i(0); i != list.size(); ++i) {
		breakup_line_quoted(list[i], args_list[i]);
	}
	argv_list.clear();
	argv_list.resize(list.size());
	for (size_t i(0); i != list.size(); ++i) {
		const std::vector<std::string> &args(args_list[i]);
		std::vector<char *> &argv(argv_list[i]);
		std::vector<std::string>::const_iterator a(args.begin());
		const std::vector<std::string>::const_iterator end_a(args.end());
		for (; a != end_a; ++a) {
			argv.push_back(const_cast<char *>(a->c_str()));
		}
		argv.push_back(0);
	}
}

// spawn proceses back to front so they're all main process's children
// list1[1] is the output file to write to (if any)
// list2 contains a list of pipeline commands
static int start_child(const std::vector<std::string> &list1, const std::vector<std::string> &list2) {
	std::vector<std::vector<std::string> > args_list;
	std::vector<std::vector<char *> > argv_list;
	convert_to_argv(list2, args_list, argv_list);
	// set up stdout for final segment of pipeline
	if (list1.size() == 2 && !list1[1].empty() && list1[1] != "-") {	// not stdout
		const int fd(open(list1[1].c_str(), O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH));
		if (fd == -1) {
			throw LocalException("open: " + list1[1] + ": " + std::string(strerror(errno)));
		}
		if (dup2(fd, STDOUT_FILENO) == -1) {
			throw LocalException("dup2: " + std::string(strerror(errno)));
		}
		close(fd);
	}
	// start up pipeline
	for (size_t i(argv_list.size() - 1);; --i) {
		pid_t pid;
		int pipefd[2];				// we only use 1(write) -> 0(read)
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
			if (execvp(argv_list[i][0], argv_list[i].data()) == -1) {
				throw LocalException("execvp: " + std::string(strerror(errno)));
			}
		} else {						// parent
			close(pipefd[0]);
			if (i == 0) {
				return pipefd[1];
			}
			if (dup2(pipefd[1], STDOUT_FILENO) == -1) {	// set up next stdout
				throw LocalException("dup2: " + std::string(strerror(errno)));
			}
			close(pipefd[1]);
		}
	}
}

static int spawn_child(const std::string &command) {
	if (command.empty()) {
		throw LocalException("empty command");
	}
	// see if we're writing to a file at the end of a pipeline
	std::vector<std::string> list1;
	breakup_line(command, list1, '>');
	if (list1.size() > 2) {
		throw LocalException("bad command: multiple > in pipeline: " + command);
	} else if (list1.size() == 2 || command.find('|') != std::string::npos) {
		// segment command into separate pipes
		std::vector<std::string> list2;
		breakup_line(list1[0], list2, '|');
		return start_child(list1, list2);
	} else {					// not a pipeline - simple file
		const int fd(open(command.c_str(), O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH));
		if (fd == -1) {
			throw LocalException("open: " + command + ": " + std::string(strerror(errno)));
		}
		return fd;
	}
}

static void notify_flag_server(const std::string &host, const int port) {
	if (host.empty()) {
		return;
	}
	const int socket_fd(socket(AF_INET, SOCK_STREAM, 0));
	if (socket_fd == -1) {
		std::cerr << "Warning: socket: " << strerror(errno) << '\n';
		return;
	}
	const struct hostent * const server(gethostbyname(host.c_str()));
	if (server == NULL) {
		std::cerr << "Warning: gethostbyname: " << h_errno << '\n';
		return;
	}
	struct sockaddr_in server_address;
	memset(&server_address, 0, sizeof(server_address));
	server_address.sin_family = AF_INET;
	memcpy(&server_address.sin_addr, server->h_addr, server->h_length);
	server_address.sin_port = htons(port);
	if (connect(socket_fd, (struct sockaddr *)&server_address, sizeof(server_address)) == -1) {
		std::cerr << "Warning: connect: " << strerror(errno) << '\n';
		return;
	}
	if (write(socket_fd, "f", 1) == -1) {
		std::cerr << "Warning: write: " << strerror(errno) << '\n';
	}
	close(socket_fd);
}

int main(int argc, char **argv) {
	int had_error(0);
	try {
		size_t buffer_size(1 << 15);
		std::vector<int> fd_list;
		std::string flag_host;
		int flag_port;
		get_opts(argc, argv, buffer_size, fd_list, flag_host, flag_port);
		Buffer buffer(buffer_size);	// fills buffer
		// send flag to flag server, if given
		notify_flag_server(flag_host, flag_port);
		// fork off children
		for (; optind != argc; ++optind) {
			fd_list.push_back(spawn_child(argv[optind]));
		}
		// close this or children later in the pipeline won't see an EOF
		close(STDOUT_FILENO);
		buffer.loop(fd_list);
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
