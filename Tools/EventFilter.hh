/*
 * EventFiler.hh
 *
 *  Created on: May 9, 2010
 *      Author: lkreczko
 */



#ifndef EVENTFILTER_HH_
#define EVENTFILTER_HH_

#include "Event.hh"

class EventFilter {

private:
	double jet_pt_cut, jet_eta_cut, jet_minimal_electromagnetic_fraction;
	double minimalCombinedRelativeIsolation, oldMaximalCombinedRelativeIsolation;
	double electron_Et_cut, electron_eta_cut, electron_d0_cut;
	float electronID;
	bool passHLT_Trigger(const Event evt);
public:
	static const double barrel_eta_region = 1.442;
	static const double endcap_minimal_eta = 1.560;
	static const double endcap_maximal_eta = 2.5;
	enum DEFINED_STAGES {
		STAGE_0_INITIAL,
		STAGE_1_GOODRUN,
		STAGE_2_PASS_HLT_TRIGGER,
		STAGE_3_MORE_THAN_ONE_TIGHT_ELECTRON,
		STAGE_4_THAN_ONE_TIGHT_ISOLATED_ELECTRON,
		STAGE_5_EXACT_ONE_TIGHT_ISOLATED_ELECTRON,
		STAGE_6_MUON_VETO,
		STAGE_7_MORE_THAN_4_GOODJETS,
		STAGE_8_MET_CUT,
		STAGE_9_Z_VETO,
		STAGE_10_CONVERSION_VETO,
		STAGE_11_ONLY_BARREL_REGION,
		STAGE_12_AT_LEAST_ONE_BTAG,
		STAGE_13_AT_LEAST_TWO_BTAG,
		NUMBER_OF_STAGES
	};
	EventFilter();
	~EventFilter();
	bool doesEventPassStage(const Event* event, DEFINED_STAGES stage);
	bool doesEventPassStages(const Event* event, std::vector<DEFINED_STAGES> stages);
	bool isGoodJet(Jet);
	bool isGoodElectron(Electron);
	bool isInBarrelRegion(Electron recoObject);
	bool isInEndcapRegion(Electron recoObject);
	bool isIsolated(Electron electron);
	bool passMuonVeto(const Event* event);
	bool passZVeto(const Event* event);
	bool hasAtLeast4GoodJets(const Event* event);
	bool hasAtLeast2Btags(const Event* event, ushort btag_type);
};

#endif /* EVENTFILTER_HH_ */
