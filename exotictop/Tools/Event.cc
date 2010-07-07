/*
 * Event.cc
 *
 *  Created on: 14 May 2010
 *      Author: kreczko
 */

#include "Event.hh"
#include <iostream>
#include <string>

Event::Event() :
	total_chi2(Event::initial_chi2), leptonic_chi2(Event::initial_chi2), hadronic_chi2(Event::initial_chi2),
			chi2Global(Event::initial_chi2) {
}


Event::Event(NTupleReader* reader, EventFilter* filter) {
	Event();
	if (reader == NULL || filter == NULL) {
		std::cerr << "Event received  NULL pointer" << std::endl;
		exit(-1);
	} else {
		this->reader = reader;
		this->filter = filter;
		loadJets();
		loadElectrons();
		loadMuons();
		loadMissingTransverseEnergies();
	}
}

void Event::loadJets() {
	std::string jetAlgorithm = reader->GetJetAlgorithm();
	if (jetAlgorithm == "Default")
		loadAntiKT5Jets();
	else if (jetAlgorithm == "SC5")
		loadSISCone5Jets();
	else if (jetAlgorithm == "pfjet")
		loadParticleFlowJets();
	else if (jetAlgorithm == "KT4")
		loadKTJets();
}

void Event::loadAntiKT5Jets() {
	for (ushort jetIndex = 0; jetIndex < this->reader->Njets; ++jetIndex) {
		Jet currentJet;
		currentJet.fourVector = TLorentzVector(reader->jets_px->at(jetIndex), reader->jets_py->at(jetIndex),
				reader->jets_pz->at(jetIndex), reader->jets_energy->at(jetIndex));
		currentJet.electromagnetic_fraction = reader->jets_emf->at(jetIndex);
		if (filter->isGoodJet(currentJet)) {
			currentJet.btags.resize(NUMBER_OF_BTAG_TYPES);
			currentJet.btags[btag_type_TrackCount_highEff] = reader->jets_btag_TC_highEff->at(jetIndex);
			currentJet.btags[btag_type_TrackCount_highPur] = reader->jets_btag_TC_highPur->at(jetIndex);
			currentJet.btags[btag_type_JetBProb] = reader->jets_btag_jetBProb->at(jetIndex);
			currentJet.btags[btag_type_JetProb] = reader->jets_btag_jetProb->at(jetIndex);
//			currentJet.btags[btag_softEle] = reader->jets_btag_softEle->at(jetIndex);
//			currentJet.btags[btag_type_softMuon] = reader->jets_btag_softMuon->at(jetIndex);
//			currentJet.btags[btag_type_softMuonNoIP] = reader->jets_btag_softMuonNoIP->at(jetIndex);
			currentJet.btags[btag_type_secondaryVertex] = reader->jets_btag_secVertex->at(jetIndex);
			currentJet.btags[btag_type_fake] = reader->jets_parton_Id->at(jetIndex);
			this->goodJets.push_back(currentJet);
		}
	}
}

void Event::loadSISCone5Jets() {

}

void Event::loadParticleFlowJets() {

}

void Event::loadKTJets() {

}

void Event::loadElectrons() {
	for (ushort electronIndex = 0; electronIndex < reader->Nels; ++electronIndex) {
		Electron electron;
		electron.fourVector = TLorentzVector(reader->els_px->at(electronIndex), reader->els_py->at(electronIndex),
				reader->els_pz->at(electronIndex), reader->els_energy->at(electronIndex));
		electron.dr04EcalRecHitSumEt = reader->els_dr04EcalRecHitSumEt->at(electronIndex);
		electron.dr04HcalTowerSumEt = reader->els_dr04HcalTowerSumEt->at(electronIndex);
		electron.tightID = reader->els_tightId->at(electronIndex);
		electron.looseID = reader->els_looseId->at(electronIndex);
		electron.robustTightID = reader->els_robustTightId->at(electronIndex);
		electron.robustLooseID = reader->els_robustLooseId->at(electronIndex);
		electron.robustHighEnergyID = reader->els_robustHighEnergyId->at(electronIndex);
		electron.vertex.x = reader->els_vx->at(electronIndex);
		electron.vertex.y = reader->els_vy->at(electronIndex);
		electron.vertex.px = reader->els_vpx->at(electronIndex);
		electron.vertex.py = reader->els_vpy->at(electronIndex);

		electron.d0 = computeCorrectedElectronD0(electron);
		electron.combinedRelativeIsolation = getCombinedRelativeIsolation(electron);
		electron.oldCombinedRelativeIsolation = getOldCombinedRelativeIsolation(electron);

		if (filter->isGoodElectron(electron))
			addGoodElectron(electron);
	}
}

void Event::addGoodElectron(Electron goodElectron){
	this->goodElectrons.push_back(goodElectron);
	if(this->filter->isInBarrelRegion(goodElectron))
		this->goodBarrelElectrons.push_back(goodElectron);
	else if(this->filter->isInEndcapRegion(goodElectron))
		this->goodEndcapElectrons.push_back(goodElectron);
}

double Event::computeCorrectedElectronD0(Electron electron) {
	float referenceX = reader->beamSpot_x->at(0);
	float referenceY = reader->beamSpot_y->at(0);
	double correctedD0 = -(-(electron.vertex.x - referenceX) * electron.vertex.px + (electron.vertex.y - referenceY)
			* electron.vertex.py) / TMath::Hypot(electron.vertex.px, electron.vertex.py);
	return correctedD0;
}

void Event::loadMuons() {

}

double Event::getCombinedRelativeIsolation(Electron lepton) {
	return (lepton.dr04EcalRecHitSumEt + lepton.dr04HcalTowerSumEt + lepton.trackerIsolation) / lepton.fourVector.Et();
}

double Event::getOldCombinedRelativeIsolation(Electron lepton) {
	return lepton.fourVector.Et() / (lepton.fourVector.Et() + lepton.dr04EcalRecHitSumEt + lepton.dr04HcalTowerSumEt
			+ lepton.trackerIsolation);
}

void Event::loadMissingTransverseEnergies() {

}

ushort Event::getNumberOfGoodJets() {
	return this->goodJets.size();
}

ushort Event::getNumberOfGoodElectrons() {
	return this->goodElectrons.size();
}
ushort Event::getNumberOfGoodBarrelElectrons() {
	return this->goodBarrelElectrons.size();
}
ushort Event::getNumberOfGoodEndcapElectrons() {
	return this->goodEndcapElectrons.size();
}

ushort Event::getNumberOfIsolatedElectrons() {
	return this->goodIsolatedElectrons.size();
}
ushort Event::getNumberOfIsolatedBarrelElectrons() {
	return this->goodIsolatedBarrelElectrons.size();
}
ushort Event::getNumberOfIsolatedEndcapElectrons() {
	return this->goodIsolatedEndcapElectrons.size();
}

Jet Event::getGoodJet(ushort index) {
	return this->goodJets.at(index);
}

Electron Event::getGoodElectrons(ushort index) {
	return this->goodElectrons.at(index);
}
Electron Event::getGoodBarrelElectrons(ushort index) {
	return this->goodBarrelElectrons.at(index);
}
Electron Event::getGoodEndcapElectrons(ushort index) {
	return this->goodEndcapElectrons.at(index);
}

Electron Event::getIsolatedElectrons(ushort index) {
	return this->goodIsolatedElectrons.at(index);
}
Electron Event::getIsolatedBarrelElectrons(ushort index) {
	return this->goodIsolatedBarrelElectrons.at(index);
}
Electron Event::getIsolatedEndcapElectrons(ushort index) {
	return this->goodIsolatedEndcapElectrons.at(index);
}

Event::~Event() {
	this->goodJets.clear();

	this->goodElectrons.clear();
	this->goodEndcapElectrons.clear();
	this->goodBarrelElectrons.clear();

	this->goodIsolatedElectrons.clear();
	this->goodIsolatedEndcapElectrons.clear();
	this->goodIsolatedBarrelElectrons.clear();

	this->goodNonIsolatedElectrons.clear();
	this->goodNonIsolatedEndcapElectrons.clear();
	this->goodNonIsolatedBarrelElectrons.clear();

	this->missingET.clear();
}

