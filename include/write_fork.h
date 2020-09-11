#ifndef _WRITE_FORK_H
#define _WRITE_FORK_H

#include <list>		// list<>
#include <string>	// string
#include <sys/stat.h>	// S_IRUSR, S_IWUSR, S_IRGRP, S_IWGRP, S_IROTH, S_IWOTH
#include <sys/types.h>	// mode_t, size_t, ssize_t

extern int write_fork(const std::list<std::string> &, const std::string &, mode_t = S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
// this version tries to guess the arguments from the file ending
extern int write_fork(const std::string &, mode_t = S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
extern void close_fork(int);
extern void close_fork_wait(int);
extern ssize_t pfputc(int, char);
extern ssize_t pfputs(int, const std::string &);
extern ssize_t pfwrite(int, const void *, const size_t);

#endif // !_WRITE_FORK_H
