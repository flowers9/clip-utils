// Sometimes we get a qual file that has been hard reformated, with lines
// exactly a certain length, but, in doing so, qualities originally at the
// end of a line get mashed up against the quality at the beginning of the
// next line (and some qualities get split down the middle).  This program
// fixes the lines so numbers don't get wrapped in the middle, and tries
// its best to split up numbers that got mashed together.  It makes the
// assumptions that the only valid quality scores are 0-56, 98, and 99,
// and that the original file had consistent line lengths for each read.
//
// Currently, it also does a fixup pass to check for #9 8 splits where the
// second bp is an N, and turns it into a # 98.

#include <fstream>	// ifstream
#include <iostream>	// cerr, cout
#include <list>		// list<>
#include <map>		// map<>
#include <sstream>	// istringstream
#include <string>	// getline(), string

static const int g_line_width(24);

class Qual {
    public:
	int qual;
	bool split;
	Qual(void) : qual(-1), split(0) { }
	explicit Qual(int x) : qual(x), split(0) { }
	Qual(int x, bool i) : qual(x), split(i) { }
	explicit Qual(const std::string &s) : split(0) {
		std::istringstream(s) >> qual;
	}
	Qual(const std::string &s, bool i) : split(i) {
		std::istringstream(s) >> qual;
	}
	~Qual(void) { }
};

class Read {
    private:
	std::list<Qual> qual;
	void choose(std::list<Qual>::iterator);
	int find_mode(void) const;
    public:
	std::string header, seq, qual_data;
	Read(void) { }
	~Read(void) { }
	void reset(const std::string &s) {
		header = s;
		seq.clear();
		qual_data.clear();
		qual.clear();
	}
	bool empty(void) const {
		return header.empty();
	}
	bool get_seq(std::ifstream &, std::string &);
	void convert_qual(void);
	void print_qual(void) const;
	void print_frequency(void) const;
};

bool Read::get_seq(std::ifstream &f_seq, std::string &seq_line) {
	if (seq_line.empty() && !getline(f_seq, seq_line)) {
		std::cerr << "Error: could not read seq\n";
		return 1;
	} else if (seq_line != header) {
		std::cerr << "Error: " << seq_line << " != " << header << "\n";
		return 1;
	}
	while (getline(f_seq, seq_line)) {
		if (seq_line.empty()) {
		} else if (seq_line[0] == '>') {
			break;
		} else {
			seq += seq_line;
		}
	}
	return 0;
}

// for non-forced three digit splits, examine nearby values to choose
// best split (### -> # ## or ## #)

void Read::choose(std::list<Qual>::iterator a) {
	const std::list<Qual>::iterator end_a(qual.end());
	std::list<Qual>::iterator start_b(a), end_b(a), b;
	// find range to nearest splits (but limit it to g_line_width)
	// go to half that distance for the chi^2 calculation
	int i(-1), j(-1);
	for (; start_b != end_a && !start_b->split && i != g_line_width; --start_b, ++i) { }
	for (; end_b != end_a && !end_b->split && j != g_line_width; ++end_b, ++j) { }
	// if ranges are unequal, use the lesser of the two
	if (start_b != end_a && end_b != end_a && i != j) {
		if (i > j) {
			i = j;
		} else {
			j = i;
		}
	}
	for (i /= 2, start_b = a; i != 0; --start_b, --i) { }
	for (j /= 2, end_b = a; j != 0; ++end_b, --j) { }
	++end_b;
	double x1(0), x2;	// averages
	// start i at 2 to account for split of a
	for (i = 2, b = start_b; b != a; ++i, ++b) {
		x1 += b->qual;
	}
	for (++b; b != end_b; ++i, ++b) {
		x1 += b->qual;
	}
	x2 = (x1 + a->qual / 100 + a->qual % 100) / i;
	x1 = (x1 + a->qual / 10 + a->qual % 10) / i;
	double y1(0), y2(0);	// diff^2
	for (b = start_b; b != a; ++b) {
		y1 += (x1 - b->qual) * (x1 - b->qual);
		y2 += (x2 - b->qual) * (x2 - b->qual);
	}
	for (++b; b != end_b; ++b) {
		y1 += (x1 - b->qual) * (x1 - b->qual);
		y2 += (x2 - b->qual) * (x2 - b->qual);
	}
	y1 += (x1 - a->qual / 10) * (x1 - a->qual / 10) + (x1 - a->qual % 10) * (x1 - a->qual % 10);
	y2 += (x2 - a->qual / 100) * (x2 - a->qual / 100) + (x2 - a->qual % 100) * (x2 - a->qual % 100);
	if (y1 < y2) {
		qual.insert(a, Qual(a->qual / 10, 1));
		a->qual %= 10;
	} else {
		qual.insert(a, Qual(a->qual / 100, 1));
		a->qual %= 100;
	}
}

// return frequency of splits, which should correspond to the line width
// of the original file (and thus be consistent at least within each read)

int Read::find_mode(void) const {
	std::map<int, int> count;
	std::list<Qual>::const_iterator a(qual.begin());
	const std::list<Qual>::const_iterator end_a(qual.end());
	int i(0);
	for (; a != end_a; ++a) {
		if (a->split) {
			++count[i];
			i = 0;
		} else {
			++i;
		}
	}
	if (count.empty()) {
		return g_line_width;
	}
	std::map<int, int>::const_iterator b(count.begin()), best(b);
	const std::map<int, int>::const_iterator end_b(count.end());
	for (++b; b != end_b; ++b) {
		if (best->second < b->second) {
			best = b;
		}
	}
	return best->first;
}

void Read::convert_qual(void) {
	if (qual_data.empty()) {
		return;
	}
	if (*qual_data.rbegin() != ' ') { // make sure it ends with a space
		qual_data += ' ';
	}
	// be careful not to ignore leading zeros during conversion
	std::string::size_type i(0), j(qual_data.find(' '));
	do {
		if (j != i + 1 && qual_data[i] == '0') {
			qual.push_back(Qual(0, 1));
			qual.push_back(Qual(qual_data.substr(i + 1, j - i - 1)));
		} else {
			qual.push_back(Qual(qual_data.substr(i, j - i)));
		}
		j = qual_data.find(' ', i = j + 1);
	} while (j != std::string::npos);
	// first, cover impossible cases with no choices
	std::list<std::list<Qual>::iterator> list;	// ones with choices
	std::list<Qual>::iterator a(qual.begin());
	const std::list<Qual>::iterator end_a(qual.end());
	for (; a != end_a; ++a) {
		if (a->qual > 999) {			// four digits
			qual.insert(a, Qual(a->qual / 100, 1));
			a->qual %= 100;
		} else if (a->qual > 99) {		// three digits
			const int x(a->qual / 10);
			const int y(a->qual % 100);
			if (56 < x && x < 98) {
				qual.insert(a, Qual(a->qual / 100, 1));
				a->qual = y;
			} else if ((56 < y && y < 98) || x % 10 == 0) {
				qual.insert(a, Qual(x, 1));
				a->qual %= 10;
			} else {			// not forced
				list.push_back(a);
			}
		} else if (56 < a->qual && a->qual < 98) {	// two digits
			qual.insert(a, Qual(a->qual / 10, 1));
			a->qual %= 10;
		}
	}
	// next, choose for impossible cases with a choice
	std::list<std::list<Qual>::iterator>::const_iterator b(list.begin());
	const std::list<std::list<Qual>::iterator>::const_iterator end_b(list.end());
	for (; b != end_b; ++b) {
		choose(*b);
	}
	// finally, go through remainder of quals
	// doing ## -> # # splits at the former line breaks
	int k(0);
	for (a = qual.begin(); a != end_a; ++a) {
		if (a->split) {
			k = 0;
		} else if (k == g_line_width) {
			std::list<Qual>::iterator c(a);
			if (++c != end_a || qual.size() != seq.size()) {
				if (a->qual < 10) {
					std::cerr << "Error: bad split: " << header << ": " << a->qual << "\n";
				}
				qual.insert(a, Qual(a->qual / 10, 1));
				a->qual %= 10;
				k = 1;
			}
		} else {
			++k;
		}
	}
	// now go back and check for incorrect #98 -> #9 8 splits
	// (I'm assuming 98 is only used for N bps)
	for (i = 0, a = qual.begin(); a != end_a; ++a, ++i) {
		if (a->split && a->qual > 9 && a->qual % 10 == 9 && seq[i + 1] == 'N') {
			std::list<Qual>::iterator c(a);
			++c;
			if (c->qual == 8) {
				a->qual /= 10;
				c->qual = 98;
			}
		}
	}
}

void Read::print_qual(void) const {
	std::cout << header << "\n";
	std::list<Qual>::const_iterator a(qual.begin());
	const std::list<Qual>::const_iterator end_a(qual.end());
	if (a == end_a) {
		return;
	}
	std::cout << a->qual;
	for (++a; a != end_a; ++a) {
		std::cout << " " << a->qual;
	}
	std::cout << "\n";
}

// used to find the original line width

void Read::print_frequency(void) const {
	std::cout << header << "\n";
	std::list<Qual>::const_iterator a(qual.begin());
	const std::list<Qual>::const_iterator end_a(qual.end());
	if (a == end_a) {
		return;
	}
	int i(0);
	for (; a != end_a; ++a) {
		if (a->split) {
			std::cout << i << (i != g_line_width ? "!\n" : "\n");
			i = 0;
		} else {
			++i;
		}
	}
	if (i != 0) {
		std::cout << i << (i > g_line_width ? "!\n" : "\n");
	}
}

int main(int argc, char *argv[]) {
	if (argc != 3) {
		std::cerr << "usage: " << argv[0] << " <fasta_file> <qual_file>\n";
		return 1;
	}
	std::ifstream f_seq(argv[1]);
	std::ifstream f_qual(argv[2]);
	Read read;
	std::string line, seq_line;
	while (getline(f_qual, line)) {
		if (line.empty()) {
		} else if (line[0] == '>') {
			if (!read.empty()) {
				if (read.get_seq(f_seq, seq_line)) {
					return 1;
				}
				read.convert_qual();
				read.print_qual();
			}
			read.reset(line);
		} else {
			read.qual_data += line;
		}
	}
	if (!read.empty()) {
		if (read.get_seq(f_seq, seq_line)) {
			return 1;
		}
		read.convert_qual();
		read.print_qual();
	}
	return 0;
}
