/*
 * Logger.cc
 *
 *  Created on: May 12, 2010
 *      Author: lkreczko
 */

#include "Logger.hh"
#include <iostream>
#include <fstream>

Logger::Logger() {
	Logger(Logger::LogLevel_NORMAL, Logger::OutputDevice_SCREEN);
}

Logger::~Logger() {

}

Logger::Logger(LogLevel level, OutputDevice device) {
	setOutputDevice(device);
	setLogLevel(level);
}

void Logger::setOutputDevice(OutputDevice device) {
	this->outputDevice = device;
}

void Logger::setLogLevel(LogLevel level) {
	this->logLevel = level;
}

void Logger::setOutPutFile(std::string outputFilename) {
	if (outputFilename == "")
		std::cerr << "ERROR: output filename for Logger is empty!" << std::endl;
	else {
		this->outputFilename = outputFilename;
		recreateLogfile();
	}
}

void Logger::log(std::string message) {
	if (this->outputDevice == OutputDevice_SCREEN)
		printOnScreen(message);
	else if (this->outputDevice == OutputDevice_FILE)
		writeIntoFile(message);
}

void Logger::logError(std::string errorMessage) {
	if (this->outputDevice == OutputDevice_SCREEN)
		printErrorOnScreen(errorMessage);
	else if (this->outputDevice == OutputDevice_FILE)
		writeIntoFile(errorMessage);
}
void Logger::printOnScreen(std::string message) {
	std::cout << message << std::endl;
}

void Logger::printErrorOnScreen(std::string errorMessage) {
	std::cerr << errorMessage << std::endl;
}

void Logger::writeIntoFile(std::string message) {
	if (this->outputFilename != "") {
		std::ofstream logfile;
		logfile.open(this->outputFilename.c_str(), std::ios::app);
		logfile << message + "\n";
		logfile.close();
	} else {
		printErrorOnScreen("ERROR: Could not write into file. File name is not set.");
	}
}

void Logger::recreateLogfile(){
	std::ofstream logfile;
	logfile.open(this->outputFilename.c_str(), std::ios::trunc);
	logfile << "Starting log file \n";
	logfile.close();
}
