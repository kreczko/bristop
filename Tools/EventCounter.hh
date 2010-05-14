/*
 * EventCounter.hh
 *
 *  Created on: May 9, 2010
 *      Author: lkreczko
 */

#ifndef EVENTCOUNTER_HH_
#define EVENTCOUNTER_HH_

//template<class T>
//typedef std::vector<std::vector<std::vector<T> > > 3DArray;

template<class T>
struct Arrays {
typedef std::vector<std::vector<std::vector<T> > > ThreeDimensional;
};
class EventCounter {
public:
	EventCounter();
	EventCounter(unsigned short sizeOfFirstDimension, unsigned short sizeOfSecondDimension, unsigned short sizeOfThirdDimension);
	~EventCounter();
	void printEventTable();
private:
	Arrays::ThreeDimensional<double> counter;
	unsigned short sizeOfFirstDimension;
	unsigned short sizeOfSecondDimension;
	unsigned short sizeOfThirdDimension;
};
#endif /* EVENTCOUNTER_HH_ */
