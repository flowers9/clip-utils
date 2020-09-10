#include "hash_read_hits.h"	// hash_read_hits::read_type
#include "kmer_lookup_info.h"
#include "open_compressed.h"	// pfread()
#include "write_fork.h"	// pfwrite()
#include <new>		// new[]

void KmerLookupInfo::save(const int fd) const {
	kmer_hash.save(fd);
	pfwrite(fd, &count, sizeof(count));
	pfwrite(fd, &data_size, sizeof(data_size));
	pfwrite(fd, list, count * sizeof(hash_read_hits::read_type));
	pfwrite(fd, data, data_size);
}

void KmerLookupInfo::restore(const int fd) {
	kmer_hash.restore(fd);
	pfread(fd, &count, sizeof(count));
	pfread(fd, &data_size, sizeof(data_size));
	list = new hash_read_hits::read_type[count];
	data = new char[data_size];
	pfread(fd, list, count * sizeof(hash_read_hits::read_type));
	pfread(fd, data, data_size);
}
