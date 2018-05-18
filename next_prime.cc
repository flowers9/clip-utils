#include "next_prime.h"
#include <list>		// list<>
#include <math.h>	// sqrt()
#include <sys/types.h>	// size_t

// return the smallest prime greater than or equal to x;
// we keep a list of primes up to sqrt(x)

size_t next_prime(size_t x) {
	if (x <= 2) {
		return 2;
	}
	if ((x & 1) == 0) {			// avoid multiples of two
		++x;
	}
	size_t n(static_cast<size_t>(sqrt(x)));
	size_t n2((n + 1) * (n + 1));
	// maintain list of known primes to speed up subsequent calls
	static std::list<size_t> primes(1, 3);
	// used for growing list of primes
	static size_t y(3);	// biggest prime found so far
	static size_t m(2);	// m is sqrt(y), more or less
	static size_t m2(9);	// (m + 1) ^ 2
	for (;;) {
		// grow primes, if needed, as we need primes up to sqrt(x)
		if (n >= y) {
			for (;;) {
				y += 2;
				if (y >= m2) {
					++m;
					m2 += 2 * m + 1;
				}
				std::list<size_t>::const_iterator a(primes.begin());
				// only need to search to sqrt(y) (i.e., m)
				for (; *a <= m && y % *a != 0; ++a) { }
				if (*a > m) {
					primes.push_back(y);
					if (y > n) {	// we've found enough
						break;
					}
				}
			}
		}
		for (;;) {	// search for next prime, starting with x
			std::list<size_t>::const_iterator a(primes.begin());
			for (; *a <= n && x % *a != 0; ++a) { }
			if (*a > n) {
				return x;
			}
			x += 2;
			if (x >= n2) {
				++n;
				n2 += 2 * n + 1;
				if (n == primes.back()) { // need to grow list
					break;
				}
			}
		}
	}
}
