/*
 * Event.hh
 *
 *  Created on: May 9, 2010
 *      Author: lkreczko
 */
#include <vector>
#include "TLorentzVector.h"
#include "RecoObjects.hh"
#include "TopPairEventCandidate.hh"

class NTupleReader;
class EventFilter;

#ifndef EVENTFILTER_HH_
#include "EventFilter.hh"
#endif

#ifndef NTUPLEREADER_H_
#include "../NTupleReader.hh"
#endif

#ifndef EVENT_HH_
#define EVENT_HH_

typedef unsigned short ushort;

class Event {
public:
	//big number since it will be minimized
	static const double initial_chi2 = 9999999.;
	static const double initial_highValue = 9999999.;
	static const short initial_id = -1;

	enum BtagTypes {
		btag_type_none,
		btag_type_fake,
		btag_type_TrackCount_highEff,
		btag_type_TrackCount_highPur,
		btag_type_JetBProb,
		btag_type_JetProb,
		btag_type_secondaryVertex,
		btag_softEle,
		btag_type_softMuon,
		btag_type_softMuonNoIP,
		NUMBER_OF_BTAG_TYPES
	};
	Event();
	Event(NTupleReader*, EventFilter*);
	~Event();
	ushort getNumberOfGoodJets();

	ushort getNumberOfGoodElectrons();
	ushort getNumberOfGoodBarrelElectrons();
	ushort getNumberOfGoodEndcapElectrons();

	ushort getNumberOfIsolatedElectrons();
	ushort getNumberOfIsolatedBarrelElectrons();
	ushort getNumberOfIsolatedEndcapElectrons();

	Jet getGoodJet(ushort index);

	Electron getGoodElectrons(ushort index);
	Electron getGoodBarrelElectrons(ushort index);
	Electron getGoodEndcapElectrons(ushort index);

	Electron getIsolatedElectrons(ushort index);
	Electron getIsolatedBarrelElectrons(ushort index);
	Electron getIsolatedEndcapElectrons(ushort index);

	Event& operator=(Event const &) {
		return *this;
	}

	void reconstructTopPairInvariantMass();//move to TTbarEventCandidate

private:
	TopPairEventCandidate ttbarCandidate;
	NTupleReader* reader;
	EventFilter* filter;

	JetCollection goodJets;

	ElectronCollection goodElectrons;
	ElectronCollection goodBarrelElectrons;
	ElectronCollection goodEndcapElectrons;

	ElectronCollection goodIsolatedElectrons;
	ElectronCollection goodIsolatedBarrelElectrons;
	ElectronCollection goodIsolatedEndcapElectrons;

	ElectronCollection goodNonIsolatedElectrons;
	ElectronCollection goodNonIsolatedBarrelElectrons;
	ElectronCollection goodNonIsolatedEndcapElectrons;

	std::vector<TLorentzVector> missingET;

	double total_chi2, leptonic_chi2, hadronic_chi2, chi2Global;
	short leptonic_bjet_id, hadronic_bjet_id, lightquark_1, chosen_lightquark_2;
	ushort chosen_neutrino;

	void initilizeValues();
	void loadJets();
	void loadAntiKT5Jets();
	void loadSISCone5Jets();
	void loadParticleFlowJets();
	void loadKTJets();

	void loadElectrons();
	void addGoodElectron(Electron electron);
	double computeCorrectedElectronD0(Electron electron);

	void loadMuons();

	void loadMissingTransverseEnergies();

	double getCombinedRelativeIsolation(Electron);
	double getOldCombinedRelativeIsolation(Electron);

};

#endif /* EVENT_HH_ */
