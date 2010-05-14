/*
 * EventFiler.hh
 *
 *  Created on: May 9, 2010
 *      Author: lkreczko
 */

#ifndef EVENTFILER_HH_
#define EVENTFILER_HH_

class EventFilter {

private:
	bool passHLT_Trigger(const Event& evt);
public:
	enum DEFINED_CUTS {
		NO_CUT, HLT_TRIGGER_CUT, MORE_THAN_ONE_TIGHT_ELECTRON_CUT, MORE_THAN_ONE_TIGHT_ISOLATED_ELECTRON_CUT,
		EXACT_ONE_TIGHT_ISOLATED_ELECTRON_CUT, MUON_VETO_CUT, MORE_THAN_4_JETS_CUT, NUMBEROFCUTS;
	};
	bool doesEventPassCut(const Event& evt, DEFINED_CUTS cut);
	bool doesEventPassCuts(const Event& evt, std::vector<DEFINED_CUTS> cuts);
};

#endif /* EVENTFILER_HH_ */
