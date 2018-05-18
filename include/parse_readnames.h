#ifndef _PARSE_READNAME_H
#define _PARSE_READNAME_H

//#include "FastaNameParser.h"	// FastaNameParser
//#include "ReadMetaInfo.h"	// ReadMetaInfo
#include <string>	// string

class FastaNameParser { };
class ReadMetaInfo {
    public:
	void setName(const std::string &) { }
	void setCenter(const std::string &) { }
	void setPlate(const std::string &) { }
	void setWell(const std::string &) { }
	void setTemplateId(const std::string &) { }
	void setTiNumber(int) { }
	void setDirection(char) { }
};
typedef std::string String;

// to handle new read name patterns, make a new Parser sub-class
// with a parse() method to turn the header line into the trace name,
// template id, plate name, well name, and direction, then add a
// check for it in pick_parser(); the parse() method should return 1
// on success, 0 on failure

// plate name is currently being set to "unknown"

class ReadNameParser : public FastaNameParser {
    protected:
	std::string lib_;
	std::string trace_, id_, plate_, well_;
	char direction_;
    private:	// only used to return an rmi
	ReadMetaInfo rmi_;
    public:
	virtual ~ReadNameParser(void) { }
	void reset_filename(const String &);
	const ReadMetaInfo &rmi(void);
	const std::string &trace(void) const {
		return trace_;
	}
	const std::string &id(void) const {
		return id_;
	}
	const std::string &plate(void) const {
		return plate_;
	}
	const std::string &well(void) const {
		return well_;
	}
	const char &direction(void) const {
		return direction_;
	}
    public:	// these get set by derived classes
	// parse takes a plain readname (no leading >)
	virtual int parse(const std::string &read) = 0;
	// extract takes a line from a file
	virtual void extractNameFromBuffer(char *line, String &name) const = 0;
};

class ReadNameParser_454 : public ReadNameParser {
    public:
	virtual int parse(const std::string &read);
	void extractNameFromBuffer(char *line, String &name) const;
};

class ReadNameParser_454_3well : public ReadNameParser_454 {
    public:
	int parse(const std::string &read);
};

class ReadNameParser_454FR : public ReadNameParser {
    public:
	virtual int parse(const std::string &read);
	void extractNameFromBuffer(char *line, String &name) const;
};

class ReadNameParser_454FR_3well : public ReadNameParser_454FR {
    public:
	int parse(const std::string &read);
};

class ReadNameParser_Ill : public ReadNameParser {
    public:
	int parse(const std::string &read);
	void extractNameFromBuffer(char *line, String &name) const;
};

class ReadNameParser_mol : public ReadNameParser {
    public:
	int parse(const std::string &read);
	void extractNameFromBuffer(char *line, String &name) const;
};

ReadNameParser * pick_readname_parser(const std::string &read, bool opt_454_3well = false);

#endif // !_PARSE_READNAME_H
