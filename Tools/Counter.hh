/*
 * EventCounter.hh
 *
 *  Created on: May 9, 2010
 *      Author: lkreczko
 */

#ifndef EVENTCOUNTER_HH_
#define EVENTCOUNTER_HH_

template<class T>
struct Arrays {
	typedef std::vector<std::vector<T> > TwoDimensional;
	typedef std::vector<std::vector<std::vector<T> > > ThreeDimensional;
};
class Counter {
public:
	Counter();
	Counter(unsigned short sizeOfFirstDimension,
			unsigned short sizeOfSecondDimension);
	~Counter();
	void printEventTable();//TODO: move this to EventCounter
private:
	Arrays<double>::TwoDimensional counter;
	Arrays<double>::TwoDimensional weighted_counter;
	unsigned short sizeOfFirstDimension;
	unsigned short sizeOfSecondDimension;
};
#endif /* EVENTCOUNTER_HH_ */
