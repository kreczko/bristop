/*
 * Event.hh
 *
 *  Created on: May 9, 2010
 *      Author: lkreczko
 */

#ifndef EVENT_HH_
#define EVENT_HH_

class Event{
private:
	std::vector<TLorentzVector> goodJets;
	std::vector<TLorentzVector> goodElectrons;
	std::vector<TLorentzVector> goodBarrelElectrons;
	std::vector<TLorentzVector> goodEndcapElectrons;
	std::vector<TLorentzVector> goodIsolatedElectrons;
	std::vector<TLorentzVector> goodIsolatedBarrelElectrons;
	std::vector<TLorentzVector> goodIsolatedEndcapElectrons;
	std::vector<TLorentzVector> missingET;
};

#endif /* EVENT_HH_ */
