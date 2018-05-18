//#include "ReadMetaInfo.h"	// ReadMetaInfo
//#include "String.h"		// String
#include "parse_readnames.h"
//#include "system/Assert.h"	// AssertEq()
//#include "system/System.h"	// FatalErr()
#include <ctype.h>	// isalnum(), isdigit(), islower(), isupper()
#include <sstream>	// ostringstream
#include <string>	// string

#define AssertEq(x,y) { }
#define FatalErr(x) { }

// to handle new read name patterns, make a new Parser sub-class
// with a parse() method to turn the header line into the trace name,
// template id, plate name, well name, and direction, then add a
// check for it in pick_parser(); the parse() method should return 1
// on success, 0 on failure

// plate name is currently being set to "unknown", so the bits that
// set the plate name are commented out

// convert alphanumeric well to number
static void convert_3well(std::string &well) {
	int x(0);
	std::string::const_iterator a(well.begin());
	const std::string::const_iterator end_a(well.end());
	for (; a != end_a; ++a) {
		x *= 36;
		if (isdigit(*a)) {
			x += *a - '0';
		} else if (isupper(*a)) {
			x += *a - 'A' + 10;
		} else if (islower(*a)) {
			x += *a - 'a' + 10;
		} else {
			FatalErr("non-alphanumeric well: " + well);
		}
	}
	std::ostringstream s;
	s << x;
	well = s.str();
}

// set lib using base filename up to (but not including) first _ or .

void ReadNameParser::reset_filename(const String &fn) {
	const std::string filename(fn.c_str());	// rather not mess with Strings
	const std::string::size_type k(filename.find_last_of('/'));
	const std::string::size_type i(k == std::string::npos ? 0 : k + 1);
	const std::string::size_type j(filename.find_first_of("_.", i));
	if (j == std::string::npos) {
		lib_ = filename.substr(i) + "_";
	} else if (j != i) {
		lib_ = filename.substr(i, j - i) + "_";
	} else {
		lib_.clear();
	}
}

const ReadMetaInfo &ReadNameParser::rmi() {
	rmi_.setName(trace_);
	rmi_.setCenter("GSC");
	rmi_.setPlate("unknown");
	rmi_.setWell(well_);
	rmi_.setTemplateId(id_);
	rmi_.setTiNumber(11394);
	rmi_.setDirection(direction_);
	return rmi_;
}

// >(([[:alnum:]]{11})([[:alnum:]]*([[:alnum:]])))
// id/trace is lib + _ + $1
// plate is lib + _ + $2
// well is $3
// direction is $4 (either R or, if not R, then F)

int ReadNameParser_454::parse(const std::string &line) {
	if (line.size() < 12) {
		return 0;
	}
	std::string::size_type i(0);
	for (; i != line.size() && isalnum(line[i]); ++i) { }
	if (i < 12) {
		return 0;
	}
	trace_ = id_ = lib_ + line.substr(0, i - 1);
	//plate_ = lib_ + line.substr(0, 11);
	well_ = line.substr(11, i - 11);
	direction_ = line[i - 1] == 'R' ? 'R' : 'F';
	return 1;
}

void ReadNameParser_454::extractNameFromBuffer(char * const buffer, String &name) const {
	name.clear();
	if (buffer == 0 || *buffer == 0) {
		return;
	}
	AssertEq(*buffer, '>');
	char * const name_start(buffer + 1);
	char *a(name_start);
	for (; *a != 0 && isalnum(*a); ++a) { }
	if (a - name_start > 11) {
		*a = 0;
		name = lib_ + name_start;
	}
}

// well is converted to base 10 from base 36

int ReadNameParser_454_3well::parse(const std::string &line) {
	if (ReadNameParser_454::parse(line) == 0) {
		return 0;
	}
	convert_3well(well_);
	return 1;
}

// >((([[:alnum:]]{11})([[:alnum:]]+))\.([RF]))
// trace is lib + _ + $1
// id is lib + _ + $2
// plate is $3
// well is $4
// direction is $5

int ReadNameParser_454FR::parse(const std::string &line) {
	if (line.size() < 14) {
		return 0;
	}
	std::string::size_type i(0);
	for (; i != line.size() && isalnum(line[i]); ++i) { }
	if (i < 12 || line[i] != '.') {
		return 0;
	}
	++i;
	if (line[i] != 'R' && line[i] != 'F') {
		return 0;
	}
	++i;
	trace_ = lib_ + line.substr(0, i - 1);
	id_ = trace_.substr(0, trace_.size() - 2);
	//plate_ = lib_ + line.substr(0, 11);
	well_ = line.substr(11, i - 13);
	direction_ = trace_[trace_.size() - 1] == 'R' ? 'R' : 'F';
	return 1;
}

void ReadNameParser_454FR::extractNameFromBuffer(char * const buffer, String &name) const {
	name.clear();
	if (buffer == 0 || *buffer == 0) {
		return;
	}
	AssertEq(*buffer, '>');
	char * const name_start(buffer + 1);
	char *a(name_start);
	for (; *a != 0 && isalnum(*a); ++a) { }
	if (a - name_start < 12 || *a != '.') {
		return;
	}
	++a;
	if (*a == 'R' || *a == 'F') {
		++a;
		*a = 0;
		name = lib_ + name_start;
	}
}

// well is converted to base 10 from base 36

int ReadNameParser_454FR_3well::parse(const std::string &line) {
	if (ReadNameParser_454FR::parse(line) == 0) {
		return 0;
	}
	convert_3well(well_);
	return 1;
}

// >(.*_(([^_]*)_[^_]*)(?:-R|/)([12]))
// trace is lib_$1-R$4
// id is lib_$1
// plate is $2
// well is $3
// direction is F if $4 == 1, else R

int ReadNameParser_Ill::parse(const std::string &line) {
	if (line.size() < 7) {
		return 0;
	}
	// remove any cruft from end of line, as it may contain
	// characters we search on; also remove leading >
	const std::string::size_type l(line.find(' '));
	const std::string s(line.substr(0, l == std::string::npos ? l : l - 1));
	std::string::size_type j(s.rfind("-R"));
	if (j == std::string::npos || j + 2 == s.size() || (s[j + 2] != '1' && s[j + 2] != '2')) {
		j = s.rfind('/');
		if (j == std::string::npos || j + 1 == s.size() || (s[j + 1] != '1' && s[j + 1] != '2')) {
			return 0;
		}
	}
	const std::string::size_type k(s.rfind('_', j));
	if (k == std::string::npos || k == 0) {
		return 0;
	}
	std::string::size_type i(s.rfind('_', k - 1));
	if (i == std::string::npos) {
		return 0;
	}
	++i;
	id_ = lib_ + s.substr(0, j);
	const char c(s[j] == '-' ? s[j + 2] : s[j + 1]);
	trace_ = id_ + "-R";
	trace_ += c;
	well_ = s.substr(i, k - i);
	// plate_ = s.substr(i, j - i);
	direction_ = c == '1' ? 'F' : 'R';
	return 1;
}

void ReadNameParser_Ill::extractNameFromBuffer(char * const buffer, String &name) const {
	name.clear();
	if (buffer == 0 || *buffer == 0) {
		return;
	}
	AssertEq(*buffer, '>');
	char * const name_start(buffer + 1);
	char *a(name_start);
	// remove any cruft from end of line, as it may contain
	// characters we search on
	for (; *a != 0 && *a != ' '; ++a) { }
	if (a == name_start) {
		return;
	}
	*a = 0;
	for (a -= 2; a >= name_start && *a != '-' && *a != '/'; --a) { }
	if (a < name_start) {
		return;
	} else if (*a == '-') {
		if (a[1] != 'R' || (a[2] != '1' && a[2] != '2')) {
			return;
		}
	} else if (*a == '/') {
		if (a[1] != '1' && a[1] != '2') {
			return;
		}
	}
	const char c(*a == '-' ? a[2] : a[1]);
	*a = 0;
	name = lib_ + name_start + "-R";
	name += c;
}

// >(.*-([^-]*))
// trace is $1
// id is $1
// plate is unknown
// well is $2
// direction is F

int ReadNameParser_mol::parse(const std::string &line) {
	if (line.size() < 3) {
		return 0;
	}
	// remove any cruft from end of line, as it may contain
	// characters we search on; also remove leading >
	const std::string::size_type i(line.find(' '));
	const std::string s(line.substr(0, i == std::string::npos ? i : i - 1));
	const std::string::size_type j(s.rfind('-'));
	if (j == std::string::npos || j + 1 == s.size()) {
		return 0;
	}
	trace_ = id_ = lib_ + s;
	well_ = s.substr(j + 1);
	direction_ = 'F';
	return 1;
}

void ReadNameParser_mol::extractNameFromBuffer(char * const buffer, String &name) const {
	name.clear();
	if (buffer == 0 || *buffer == 0) {
		return;
	}
	AssertEq(*buffer, '>');
	char * const name_start(buffer + 1);
	char *a(name_start);
	// remove any cruft from end of line, as it may contain
	// characters we search on
	for (; *a != 0 && *a != ' '; ++a) { }
	if (a == name_start) {
		return;
	}
	*a = 0;
	// now finish check for name correctness
	for (--a; a >= name_start && *a != '-'; --a) { }
	if (a >= name_start) {
		name = lib_ + name_start;
	}
}

ReadNameParser * pick_readname_parser(const std::string &read, const bool opt_454_3well) {
	// order matters if multiple parsers can parse the same line;
	// for example, 454FR has to be tested for before 454, as 454 will
	// match 454FR lines, but won't set the direction properly
	ReadNameParser * parser = new ReadNameParser_454FR;
	if (parser->parse(read)) {
		if (opt_454_3well) {
			delete parser;
			return new ReadNameParser_454FR_3well;
		} else {
			return parser;
		}
	}
	delete parser;
	parser = new ReadNameParser_454;
	if (parser->parse(read)) {
		if (opt_454_3well) {
			delete parser;
			return new ReadNameParser_454_3well;
		} else {
			return parser;
		}
	}
	delete parser;
	parser = new ReadNameParser_Ill;
	if (parser->parse(read)) {
		return parser;
	}
	delete parser;
	parser = new ReadNameParser_mol;
	if (parser->parse(read)) {
		return parser;
	}
	delete parser;
	return 0;
}
