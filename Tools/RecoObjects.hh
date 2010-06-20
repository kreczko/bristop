/*
 * Structures.hh
 *
 *  Created on: 23 May 2010
 *      Author: kreczko
 */

#ifndef STRUCTURES_HH_
#define STRUCTURES_HH_
struct Vertex {
	double x, y;
	double px, py;
};

struct RecoObject {
	TLorentzVector fourVector;
};
struct Jet: public RecoObject {
	double electromagnetic_fraction;
	std::vector<double> btags;
};

struct Electron: public RecoObject {
	double combinedRelativeIsolation;
	double oldCombinedRelativeIsolation;
	double dr04EcalRecHitSumEt;
	double dr04HcalTowerSumEt;
	double trackerIsolation;

	float tightID;
	float looseID;
	float robustTightID;
	float robustLooseID;
	float robustHighEnergyID;

	double d0;
	Vertex vertex;
	double computeCorrectedElectronD0(Vertex beamspot) {
		double correctedD0 = -(-(this->vertex.x - beamspot.x) * this->vertex.px + (this->vertex.y - beamspot.y)
				* this->vertex.py) / TMath::Hypot(this->vertex.px, this->vertex.py);
		return correctedD0;
	}
};

struct Muon: public RecoObject {
	Vertex vertex;
};

typedef std::vector<Jet> JetCollection;
typedef std::vector<Electron> ElectronCollection;
typedef std::vector<Muon> MuonCollection;

#endif /* STRUCTURES_HH_ */
