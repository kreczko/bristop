/*
 * EventFilter.cc
 *
 *  Created on: 14 May 2010
 *      Author: kreczko
 */

#include "EventFilter.hh"

EventFilter::EventFilter() :
	jet_pt_cut(30.), jet_eta_cut(2.1), jet_minimal_electromagnetic_fraction(0.01),
			minimalCombinedRelativeIsolation(0.1), oldMaximalCombinedRelativeIsolation(0.85), electron_Et_cut(30),
			electron_eta_cut(2.4), electron_d0_cut(0.02) {
}

EventFilter::~EventFilter() {

}

bool EventFilter::isGoodJet(Jet jet) {
	bool passPtRequirement = jet.fourVector.Pt() > jet_pt_cut;
	bool passEtaRequirement = fabs(jet.fourVector.Eta()) < jet_eta_cut;
	bool passEMFractionRequirement = jet.electromagnetic_fraction > jet_minimal_electromagnetic_fraction;
	return passPtRequirement && passEtaRequirement && passEMFractionRequirement;
}

bool EventFilter::isGoodElectron(Electron electron) {
	bool idRequirement = electron.robustTightID > 0;
	bool EtRequirement = electron.fourVector.Et() > electron_Et_cut;
	bool EtaRequirement = fabs(electron.fourVector.Eta()) < electron_eta_cut;
	bool d0Requirement = fabs(electron.d0) < electron_d0_cut;
	return idRequirement && EtRequirement && EtaRequirement && d0Requirement;
}

bool EventFilter::isInBarrelRegion(Electron recoObject) {
	return fabs(recoObject.fourVector.Eta()) < barrel_eta_region;
}

bool EventFilter::isInEndcapRegion(Electron recoObject) {
	double absolute_eta = fabs(recoObject.fourVector.Eta());
	return absolute_eta > endcap_minimal_eta && absolute_eta < endcap_maximal_eta;
}

bool EventFilter::isIsolated(Electron electron){
	return electron.combinedRelativeIsolation < minimalCombinedRelativeIsolation;
}

