/*
 * Analysis.cpp
 *
 *  Created on: 12 Jul 2010
 *      Author: kreczko
 */

#include "Analysis.h"
#include <iostream>
using namespace BAT;
using namespace std;
Analysis::Analysis(): eventReader(new BAT::NTupleEventReader()), eventFilter(Filter::makeStandardFilter()) {

}

Analysis::~Analysis() {
	// TODO Auto-generated destructor stub
}

void Analysis::addInputFile(const char* fileName){
	eventReader->addInputFile(fileName);
}

void Analysis::analyze(){
	unsigned long numberOfEvents = eventReader->getNumberOfEvents();
	cout << "total number of events to analyse: " << numberOfEvents << endl;
	for(unsigned long eventIndex = 0; eventIndex < numberOfEvents; eventIndex++){
		cout << "analysing event No. " << eventIndex << endl;
		Event* event = eventReader->getNextEvent();
		cout << event->getOtherElectrons()->at(0).energy() << endl;
	}
}
