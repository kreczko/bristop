/*
 * TopPairEventCandidate.hh
 *
 *  Created on: 23 May 2010
 *      Author: kreczko
 */

#include "RecoObjects.hh"

#ifndef TOPPAIREVENTCANDIDATE_HH_
#define TOPPAIREVENTCANDIDATE_HH_
struct TopPairEventCandidate {
	static const double matched_angle = 0.945666;
	static const double matched_angle_sigma = 0.311091;
	static const double matched_leptonic_top_mass = 178.377;
	static const double matched_leptonic_top_mass_sigma = 31.050;
	static const double matched_hadronic_W_mass = 89.9153;
	static const double matched_hadronic_W_mass_sigma = 13.8711;
	static const double matched_hadronic_top_mass = 182.191;
	static const double matched_hadronic_top_mass_sigma = 22.1484;
	static const double matched_ptratio = 0.18552;
	static const double matched_ptratio_sigma = 0.401973;
	static const double matched_pt_ttbarSystem = 0.0760939;
	static const double matched_pt_ttbarSystem_sigma = 0.0700391;
	static const double matched_HTSystem = 1;
	static const double matched_HTSystem_sigma = 0.1;

	Jet leptonic_b, hadronic_b, quarkFromWdecay, quark2FromWdecay;
	Electron electron;
	RecoObject neutrino;
	double HT;

	TopPairEventCandidate& operator=(TopPairEventCandidate& a) {
		return a;
	}

	TLorentzVector leptonic_W() {
		return TLorentzVector(electron.fourVector + neutrino.fourVector);
	}

	TLorentzVector hadronic_W() {
		return TLorentzVector(quarkFromWdecay.fourVector + quark2FromWdecay.fourVector);
	}
	TLorentzVector hadronic_top() {
		return TLorentzVector(hadronic_W() + hadronic_b.fourVector);
	}

	TLorentzVector leptonic_top() {
		return TLorentzVector(leptonic_W() + leptonic_b.fourVector);
	}

	double mttbar() {
		return TLorentzVector(leptonic_top() + hadronic_top()).M();
	}

	double SumPt() {
		return leptonic_b.fourVector.Pt() + hadronic_b.fourVector.Pt() + quarkFromWdecay.fourVector.Pt()
				+ quark2FromWdecay.fourVector.Pt();
	}

	double HTSystem() {
		return SumPt() / HT;
	}

	double angleBetweenBjetAndElectron() {
		return leptonic_b.fourVector.Angle(electron.fourVector.Vect());
	}

	double PtRatio() {
		return TMath::Log(hadronic_top().Pt() / hadronic_W().Pt());
	}

	double Pt_ttbarSystem() {
		return TLorentzVector(leptonic_top() + hadronic_top()).Pt() / HT;
	}

	double chi2Global() {
		double pt_ttbarSystem_term = TMath::Power(Pt_ttbarSystem(), 2) / (matched_pt_ttbarSystem_sigma
				* matched_pt_ttbarSystem_sigma * 2);
		double HTSystem_term = TMath::Power(HTSystem() - matched_HTSystem, 2) / (2 * matched_HTSystem_sigma
				* matched_HTSystem_sigma);
		double normalization = TMath::Sqrt(2);
		return (pt_ttbarSystem_term + HTSystem_term) / normalization;
	}

	double chi2Leptonic() {
		double angle_term = TMath::Power(angleBetweenBjetAndElectron() - matched_angle, 2) / (2 * matched_angle_sigma
				* matched_angle_sigma);
		double top_mass_term = TMath::Power(leptonic_top().M() - matched_leptonic_top_mass, 2) / (2
				* matched_leptonic_top_mass_sigma * matched_leptonic_top_mass_sigma);
		double normalization = TMath::Sqrt(2);
		return (angle_term + top_mass_term) / normalization;
	}

	double chi2Hadronic() {
		double W_mass_term = TMath::Power(hadronic_W().M() - matched_hadronic_W_mass, 2) / (2
				* matched_hadronic_W_mass_sigma * matched_hadronic_W_mass_sigma);
		double top_mass_term = TMath::Power(hadronic_top().M() - matched_hadronic_top_mass, 2) / (2
				* matched_hadronic_top_mass_sigma * matched_hadronic_top_mass_sigma);
		double pt_ratio_term = TMath::Power(PtRatio() - matched_ptratio, 2) / (2 * matched_ptratio_sigma
				* matched_ptratio_sigma);
		double normalization = TMath::Sqrt(3);
		return (W_mass_term + top_mass_term + pt_ratio_term) / normalization;
	}
};

#endif /* TOPPAIREVENTCANDIDATE_HH_ */
