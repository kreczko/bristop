/*
 * Logger.hh
 *
 *  Created on: May 12, 2010
 *      Author: lkreczko
 */
#include <string>
#include <iostream>
#include <cstdlib>
#include <sstream>

#ifndef LOGGER_HH_
#define LOGGER_HH_
//template<typename Char, typename Traits = char_traits<Char> >
//struct LogStream{
//    typedef std::basic_ostream<Char, Traits> ostream_type;
//    typedef ostream_type& (*manip_type)(ostream_type&);
//    LogStream(ostream_type& os):os(os){}
//    LogStream &operator<<(manip_type pfn) {
//        if(pfn == static_cast<manip_type>(std::endl)) {
//            time_t t = time(0);
//            os << " --- " << ctime(&t) << pfn;
//        } else
//            os << pfn;
//        return *this;
//    }
//    template<typename T>
//    LogStream &operator<<(T const& t) {
//        os << t;
//        return *this;
//    }
//private:
//    ostream_type & os;
//};

class Logger {//TODO: change it to be a struct since it can be stored with 'new' only

public:
	enum LogLevel {
		LogLevel_DEBUG, LogLevel_NORMAL
	};

	enum OutputDevice {
		OutputDevice_SCREEN, OutputDevice_FILE
	};
	Logger();
	Logger(LogLevel level, OutputDevice device);
	~Logger();

	void log(std::string message);
	void log(std::string message, LogLevel level);
	void logError(std::string errorMessage);

	void printOnScreen(std::string message);
	void printErrorOnScreen(std::string errorMessage);
	void writeIntoFile(std::string message);

	void setLogLevel(LogLevel level);
	void setOutPutFile(std::string outputFilename);
	void setOutputDevice(OutputDevice device);

	void recreateLogfile();
	template<typename T>
	Logger &operator<<(T const& message) {
		this->outputStream << message;
		return *this;
	}

	Logger &operator=(Logger const &){
		return *this;
	}

//	Logger &operator<<(manip_type pfn) {
//		if (pfn == static_cast<manip_type> (std::endl)) {
//
//		} else
//			this->outputStream << pfn;
//		return *this;
//	}


private:
	std::stringstream outputStream;
	LogLevel logLevel;
	OutputDevice outputDevice;
	std::string outputFilename;

};

#endif /* LOGGER_HH_ */
