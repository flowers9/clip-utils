#include "breakup_line.h"	/* breakup_line_exact() */
#include <list>		/* list<> */
#include <stdio.h>	/* fprint(), stderr */
#include <stdlib.h>	/* getenv() */
#include <string.h>	/* strncmp() */
#include <string>	/* string */
#include <sys/stat.h>	/* stat(), struct stat */
#include <vector>	/* vector<> */

#define RAN_DIR "/home/raid2/MB/ranblocks/"

int main() {
	char *s = getenv("REMOTE_ADDR");
	if (s == NULL || (strncmp(s, "172.26.2", 8) != 0 && strncmp(s, "127.0.0.1", 9) != 0)) {
		if (s == NULL) {
			fprintf(stderr, "check_barcodes: no address\n");
		} else {
			fprintf(stderr, "check_barcodes: incorrect address: %s\n", s);
		}
		return 1;
	}
	s = getenv("QUERY_STRING");
	if (s == NULL || strncmp(s, "barcodes=", 9) != 0) {
		if (s == NULL) {
			fprintf(stderr, "check_barcodes: no query\n");
		} else {
			fprintf(stderr, "check_barcodes: incorrect query: %s\n", s);
		}
		return 2;
	}
	std::list<std::string> used_list;
	std::vector<std::string> list;
	breakup_line_exact(s + 9, ",", list);
	struct stat buf;
	std::vector<std::string>::const_iterator a = list.begin();
	std::vector<std::string>::const_iterator end_a = list.end();
	for (; a != end_a; ++a) {
		if (stat(std::string(RAN_DIR + *a).c_str(), &buf) != -1) {
			used_list.push_back(*a);
		}
	}
	printf("Content-Type: text/plain\n\n");
	int first = 0;
	std::list<std::string>::const_iterator b = used_list.begin();
	std::list<std::string>::const_iterator end_b = used_list.end();
	for (; b != end_b; ++b) {
		if (first) {
			printf(" ");
		} else {
			first = 1;
		}
		printf("%s", b->c_str());
	}
	printf("\n");
	return 0;
}
