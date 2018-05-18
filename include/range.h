#ifndef _RANGE_H
#define _RANGE_H

#include <sys/types.h>	// size_t

class Range {
    public:
	size_t start, stop;
	Range(const size_t start_in, const size_t stop_in) : start(start_in), stop(stop_in) { }
	~Range() { }
	// check to see if new range immediately follows this one -
	// if so, extend
	bool extend(const Range &a) {
		if (stop + 1 != a.start) {
			return 0;
		} else {
			stop = a.stop;
			return 1;
		}
	}
	// extend a range by a single position, if it immediately follows
	// the current end
	bool extend(size_t position) {
		if (stop + 1 != position) {
			return 0;
		} else {
			++stop;
			return 1;
		}
	}
	bool operator<(const Range &a) const {
		if (start != a.start) {
			return start < a.start;
		} else {
			return stop < a.stop;
		}
	}
};

#endif // !_RANGE_H
