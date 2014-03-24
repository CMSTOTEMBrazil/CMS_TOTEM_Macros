//STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TDirectory.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <TRandom.h>

//STANDARD C++ INCLUDES
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>

//OUR OWN CLASSES TO READ THE TREE
#include "MassParticles.h"
#include "MyBaseJet.h"
#include "MyBeamSpot.h"
#include "MyCaloJet.h"
#include "MyCastorDigi.h"
#include "MyCastorJet.h"
#include "MyCastorRecHit.h"
#include "MyDiJet.h"
#include "MyElectron.h"
#include "MyEvtId.h"
#include "MyFwdGap.h"
#include "MyGenJet.h"
#include "MyGenKin.h"
#include "MyGenMet.h"
#include "MyGenPart.h"
#include "MyHLTrig.h"
#include "MyJet.h"
#include "MyL1Trig.h"
#include "MyL1TrigOld.h"
//#include "MyMITEvtSel.h"
#include "MyMet.h"
#include "MyMuon.h"
#include "MyPFCand.h"
#include "MyPFJet.h"
#include "MyPUSumInfo.h"
#include "MyPart.h"
#include "MySimVertex.h"
#include "MyTracks.h"
#include "MyVertex.h"
#include "RPRootDumpReconstructedProton.h"

#include "rp_aperture_config.h"

#define PI 3.141592653589793
using namespace std;

double etaBinsHCALBoundaries[] = {-5.205, -4.903, -4.730,
	-4.552, -4.377, -4.204, -4.027, -3.853, -3.677, -3.503, -3.327, -3.152,
	-3.000, -2.868, -2.650, -2.500,
	-2.322, -2.172, -2.043, -1.930, -1.830, -1.740, -1.653, -1.566, -1.479,
	-1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.870, -0.783,
	-0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087,
	0.000, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696,
	0.783, 0.870, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392,
	1.479, 1.566, 1.653, 1.740, 1.830, 1.930, 2.043, 2.172, 2.322,
	2.500, 2.650, 2.868, 3.000,
	3.152, 3.327, 3.503, 3.677, 3.853, 4.027, 4.204, 4.377, 4.552,
	4.730, 4.903, 5.205}; // 41 + 41 bins

class MyZeroBiasData {
	public:
		MyZeroBiasData() {}
		~MyZeroBiasData() {}

		bool   proton_rec_left_valid;
		double proton_rec_left_t;
		double proton_rec_left_xi;
		bool   proton_rec_right_valid;
		double proton_rec_right_t;
		double proton_rec_right_xi;
		bool   vtx_valid;
		double vtx_ndof;
		double vtx_x;
		double vtx_y;
		double vtx_z;
};


void PYTHIA_SDJpsi_PU_2(string const& outputFileName = "pythia_18March_pu.root", const Double_t t_proton_down_=0.0, const Double_t t_proton_up_=1.0,const Int_t Bin_mass=200,  const Int_t nevt_max = -1){

	bool verbose = false;
	string treeName = "evt";//"cms_totem";
	bool selectVertex = true;
	bool Vertex = false;
	bool selectMuons = true;
	bool selectTrack = true;
	double ptMax = 9999.0;
	bool selectRPPlusAccept = false;
	bool selectRPMinusAccept =true;
	bool sdplus = true;
	bool sdminus = false;

	const Int_t nevt_max_corr = (nevt_max >= 0) ? nevt_max : 99999999;


	// Declaration of histograms
	map<string,TH1F*> histosTH1F;
	map<string,TH2F*> histosTH2F;


	histosTH1F["vtx_xpos"] = new TH1F("vtx_xpos", "x(vtx)" , 150 , -1.5 , 1.5);
	histosTH1F["vtx_ypos"] = new TH1F("vtx_ypos", "y(vtx)" , 150 , -1.5 , 1.5);
	histosTH1F["vtx_ndof"] = new TH1F("vtx_ndof", "ndof(vtx)" , 100 , 0. , 15.);
	histosTH1F["vtx_chi2"] = new TH1F("vtx_chi2", "chi2(vtx)" , 100 , 0. , 10.);

	histosTH1F["vertex_multiplicity"] = new TH1F("vertex_multiplicity", "n vertices" , 30 , 0 , 30);
	histosTH1F["vertex_multiplicity_after_vtx_sel"] = new TH1F("vertex_multiplicity_after_vtx_sel", "n vertices after vtx sel" , 30 , 0 , 30);
	histosTH1F["prim_vtx_zpos"] = new TH1F("prim_vtx_zpos", "z(vtx)" , 150 , -30. , 30.);
	histosTH1F["prim_vtx_xpos"] = new TH1F("prim_vtx_xpos", "x(vtx)" , 150 , -1.5 , 1.5);
	histosTH1F["prim_vtx_ypos"] = new TH1F("prim_vtx_ypos", "y(vtx)" , 150 , -1.5 , 1.5);

	histosTH1F["prim_vtx_ndof"] = new TH1F("prim_vtx_ndof", "ndof(vtx)" , 100 , 0. , 15.);
	histosTH1F["prim_vtx_chi2"] = new TH1F("prim_vtx_chi2", "chi2(vtx)" , 100 , 0. , 10.);
	histosTH1F["prim_vtx_chi2n"] = new TH1F("prim_vtx_chi2n", "chi2n(vtx)" , 100 , 0. , 10.);
	histosTH1F["prim_vtx_ntracks"] = new TH1F("prim_vtx_ntracks", "n_{trk}(vtx)" , 30 , 0 , 30);
	histosTH1F["prim_vtx_sumpt"] = new TH1F("prim_vtx_sumpt", "sum(p_{T})(vtx)" , 100 , 0. , 100.);

	histosTH1F["prim_vtx_zpos_after_vtx_sel"] = new TH1F("prim_vtx_zpos_after_vtx_sel", "z(vtx)" , 150 , -30. , 30.);
	histosTH1F["prim_vtx_xpos_after_vtx_sel"] = new TH1F("prim_vtx_xpos_after_vtx_sel", "x(vtx)" , 150 , -1.5 , 1.5);
	histosTH1F["prim_vtx_ypos_after_vtx_sel"] = new TH1F("prim_vtx_ypos_after_vtx_sel", "y(vtx)" , 150 , -1.5 , 1.5);

	histosTH1F["prim_vtx_ndof_after_vtx_sel"] = new TH1F("prim_vtx_ndof_after_vtx_sel", "ndof(vtx)" , 100 , 0. , 15.);
	histosTH1F["prim_vtx_chi2_after_vtx_sel"] = new TH1F("prim_vtx_chi2_after_vtx_sel", "chi2(vtx)" , 100 , 0. , 10.);
	histosTH1F["prim_vtx_chi2n_after_vtx_sel"] = new TH1F("prim_vtx_chi2n_after_vtx_sel", "chi2n(vtx)" , 100 , 0. , 10.);
	histosTH1F["prim_vtx_ntracks_after_vtx_sel"] = new TH1F("prim_vtx_ntracks_after_vtx_sel", "n_{trk}(vtx)" , 30 , 0 , 30);
	histosTH1F["prim_vtx_sumpt_after_vtx_sel"] = new TH1F("prim_vtx_sumpt_after_vtx_sel", "sum(p_{T})(vtx)" , 100 , 0. , 100.);

	//histosTH1F["pt_gen"] = new TH1F("pt_gen" , "pt_gen;pt;nTracks" , 120 , 0 , 6);
	histosTH1F["track_pt"] = new TH1F("track_pt", "p_{T}(trk)" , 150 , 0. , 15.);
	histosTH1F["track_eta"] = new TH1F("track_eta", "#eta(trk)" , 200 , -5.2 , 5.2);
	histosTH1F["track_phi"] = new TH1F("track_phi", "#phi(trk)" , 200 , -M_PI , M_PI);
	histosTH1F["track_rapidity"] = new TH1F("track_rapidity", "#rapidity(trk)" , 100 , -5.2 , 5.2);
	histosTH1F["track_multiplicity"] = new TH1F("track_multiplicity", "n tracks" , 100 , 0 , 100);

	histosTH1F["muon_pt"] = new TH1F("muon_pt", "p_{T}(muon)" , 100 , 0. , 100.);
	histosTH1F["muon_eta"] = new TH1F("muon_eta", "#eta(muon)" , 100 , -5.2 , 5.2);
	histosTH1F["muon_phi"] = new TH1F("muon_phi", "#phi(muon)" , 100 , -1.2*M_PI , 1.2*M_PI);
	histosTH1F["muon_multiplicity"] = new TH1F("muon_multiplicity", "n muons" , 100 , 0 , 100);

	//Muon1 info
	histosTH1F["muon1_pt"] = new TH1F("muon1_pt", "p_{T}(mu1)" , 100 , 0. , 100.);
	histosTH1F["muon1_eta"] = new TH1F("muon1_eta", "#eta(mu1)" , 100 , -5.2 , 5.2);
	histosTH1F["muon1_phi"] = new TH1F("muon1_phi", "#phi(muon1)" , 100 , -1.2*M_PI ,1.2*M_PI);
	histosTH1F["muon1_rapidity"] = new TH1F("muon1_rapidity", "y(mu1)" , 100 , -15. , 15.);

	//Muon2 info
	histosTH1F["muon2_pt"] = new TH1F("muon2_pt", "p_{T}(mu2)" , 100 , 0. , 100.);
	histosTH1F["muon2_eta"] = new TH1F("muon2_eta", "#eta(mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["muon2_phi"] = new TH1F("muon2_phi", "#phi(muon2)" , 100 , -1.2*M_PI ,1.2*M_PI);
	histosTH1F["muon2_rapidity"] = new TH1F("muon2_rapidity", "y(mu2)" , 100 , -15. , 15.);

	//Deltas info
	histosTH1F["muonDeltaPt"] = new TH1F("muonDeltaPt", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta"] = new TH1F("muonDeltaEta", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);	
	histosTH1F["muonDeltaPhi"] = new TH1F("muonDeltaPhi", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);
	histosTH1F["muonDeltaY"] = new TH1F("muonDeltaY", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["lorentzdphi"] = new TH1F("lorentzDPhi", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);    
	//Dimuon
	histosTH1F["dimuon_mass"] = new TH1F("dimuon_mass", "mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["proton_right_dimuon_mass"] = new TH1F("proton_right_dimuon_mass", "mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["proton_left_dimuon_mass"] = new TH1F("proton_left_dimuon_mass", "mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["dimuon_pt2"] = new TH1F("dimuon_pt2", "p_{T}^{2}(mu1,mu2)" , 100 , 0. , 1000.);
	histosTH1F["dimuon_pt"] = new TH1F("dimuon_pt", "p_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["dimuon_eta"] = new TH1F("dimuon_eta", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["dimuon_rapidity"] = new TH1F("dimuon_rapidity", "y(mu1,mu2)" , 100 , -15. , 15.);
	histosTH1F["dimuon_multiplicity"] = new TH1F("dimuon_multiplicity", "n dimuons" , 100 , 0 , 100); 

	//jpsi mass selection
	histosTH1F["jpsi_mass"] = new TH1F("jpsi_mass", "jpsi_mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["jpsi_pt"] = new TH1F("jpsi_pt", "jpsi p_{T}(mu1,mu2)" , 100 , 0. , 20.);
	histosTH1F["jpsi_pt2"] = new TH1F("jpsi_pt2", "jpsi p_{T}^{2}(mu1,mu2)" , 100 , 0. , 1000.);
	histosTH1F["jpsi_eta"] = new TH1F("jpsi_eta", "jpsi #eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["jpsi_rapidity"] = new TH1F("jpsi_rapidity", "jpsi y(mu1,mu2)" , 100 , -15. , 15.);
	histosTH1F["jpsi_multiplicity"] = new TH1F("jpsi_multiplicity", "n jpsi(dimuons)" , 100 , 0 , 100);

	//Deltas info
	histosTH1F["muonDeltaPt_jpsi"] = new TH1F("muonDeltaPt_jpsi", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta_jpsi"] = new TH1F("muonDeltaEta_jpsi", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["muonDeltaPhi_jpsi"] = new TH1F("muonDeltaPhi_jpsi", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);
	histosTH1F["muonDeltaY_jpsi"] = new TH1F("muonDeltaY_jpsi", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["lorentzdphi_jpsi"] = new TH1F("lorentzDPhi_jpsi", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);

	histosTH1F["Eta_max"] = new TH1F("Eta_max", "#eta^{max}" , 82 , etaBinsHCALBoundaries);
	histosTH1F["Eta_min"] = new TH1F("Eta_min", "#eta^{min}" , 82 , etaBinsHCALBoundaries);
	histosTH1F["Delta_eta_maxmin"] = new TH1F("Delta_eta_maxmin", "#eta^{max} - #eta^{min}" , 50 , 0,11);
	histosTH1F["xi_plus_Reco"] = new TH1F("xi+_Reco", "#xi^{+}" , 82 , 0,4);
	histosTH1F["xi_minus_Reco"] = new TH1F("xi-_Reco", "#xi^{-}" , 82 , 0,4);
	histosTH1F["logxi_plus"] = new TH1F("logxi+_Reco", "Log #xi^{+}" , 82 , -3,0.5);
	histosTH1F["logxi_plus_gen"] = new TH1F("logxi+_gen", "Log #xi_{+}^{gen}" , 82 , -3,0.5);
	histosTH1F["logxi_minus_gen"] = new TH1F("logxi-_gen", "Log #xi_{-}^{gen}" , 82 , -3,0.5);
	histosTH1F["correction"] = new TH1F("correction", "Correction factor" , 82 , 0,2);
	histosTH1F["resolution_after"] = new TH1F("resolution_after", "Resolution" , 82 , -2,2);
	histosTH1F["resolution_before"] = new TH1F("resolution_before", "Resolution" , 82 , -2,2);

	histosTH1F["pf_xiPlus_minus_proton_left_xi"] = new TH1F("pf_xiPlus_minus_proton_left_xi", "#xi diff." , 200 , -1. , 1.);
	histosTH1F["pf_xiMinus_minus_proton_right_xi"] = new TH1F("pf_xiMinus_minus_proton_right_xi", "#xi diff." , 200 , -1. , 1.);
	///
	Float_t tbins[12] = { 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.22, 0.30, 0.40, 0.50, 0.65, 1.};
	histosTH1F["xi_proton_minus"] = new TH1F("xi_proton_minus", "xi_proton_minus" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_minus_rec"] = new TH1F("xi_proton_minus_rec", "xi_proton_minus_sel" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_plus"] = new TH1F("xi_proton_plus", "xi_proton_plus" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_plus_rec"] = new TH1F("xi_proton_plus_rec", "xi_proton_plus_sel" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_right_zb"] = new TH1F("xi_proton_right_zb", "xi_proton_right" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_left_zb"] = new TH1F("xi_proton_left_zb", "xi_proton_left" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_signal_right"] = new TH1F("xi_proton_signal_right", "xi_proton_right" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_signal_right_kint"] = new TH1F("xi_proton_signal_right_kint", "xi_proton_right" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_signal_right_kin_cut"] = new TH1F("xi_proton_signal_right_kin_cut", "xi_proton_right" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_signal_right_kin_cut_halo"] = new TH1F("xi_proton_signal_right_kin_cut_halo", "xi_proton_right" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_signal_right_kint_cut"] = new TH1F("xi_proton_signal_right_kint_cut", "xi_proton_right" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_backg_right"] = new TH1F("xi_proton_backg_right", "xi_proton_right" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_backg_right_kin"] = new TH1F("xi_proton_backg_right_kin", "xi_proton_right" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_backg_right_kin_cut"] = new TH1F("xi_proton_backg_right_kin_cut", "xi_proton_right" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_backg_right_kin_cut_halo"] = new TH1F("xi_proton_backg_right_kin_cut_halo", "xi_proton_right" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_backg_right_kint_cut"] = new TH1F("xi_proton_backg_right_kint_cut", "xi_proton_right" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_backg_right_kint"] = new TH1F("xi_proton_backg_right_kint", "xi_proton_right" , 20, -0.1, 0.3);
	histosTH1F["xi_right_zb"] = new TH1F("xi_right_zb", "xi_proton_right" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_signal_left"] = new TH1F("xi_proton_signal_left", "xi_proton_left" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_signal_left_kint"] = new TH1F("xi_proton_signal_left_kint", "xi_proton_left" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_signal_left_kin_cut"] = new TH1F("xi_proton_signal_left_kin_cut", "xi_proton_left" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_signal_left_kin_cut_halo"] = new TH1F("xi_proton_signal_left_kin_cut_halo", "xi_proton_left" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_backg_left"] = new TH1F("xi_proton_backg_left", "xi_proton_left" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_backg_left_kin"] = new TH1F("xi_proton_backg_left_kin", "xi_proton_left" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_backg_left_kin_cut"] = new TH1F("xi_proton_backg_left_kin_cut", "xi_proton_left" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_backg_left_kin_cut_halo"] = new TH1F("xi_proton_backg_left_kin_cut_halo", "xi_proton_left" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_backg_left_kint"] = new TH1F("xi_proton_backg_left_kint", "xi_proton_left" , 20, -0.1, 0.3); 
	histosTH1F["xi_proton_right_both"] = new TH1F("xi_proton_right_both", "xi_proton_right" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_both_right_kint"] = new TH1F("xi_proton_both_right_kint", "xi_proton_right" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_both_right_kint_cut"] = new TH1F("xi_proton_both_right_kint_cut", "xi_proton_right" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_both_right_kin_cut"] = new TH1F("xi_proton_both_right_kin_cut", "xi_proton_right" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_both_right_kin_cut_halo"] = new TH1F("xi_proton_both_right_kin_cut_halo", "xi_proton_right" , 20, -0.1, 0.3);

	histosTH1F["xi_proton_left_both"] = new TH1F("xi_proton_left_both", "xi_proton_left" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_both_left_kint"] = new TH1F("xi_proton_both_left_kint", "xi_proton_left" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_both_left_kint_cut"] = new TH1F("xi_proton_both_left_kint_cut", "xi_proton_left" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_both_left_kin_cut"] = new TH1F("xi_proton_both_left_kin_cut", "xi_proton_left" , 20, -0.1, 0.3);
	histosTH1F["xi_proton_both_left_kin_cut_halo"] = new TH1F("xi_proton_both_left_kin_cut_halo", "xi_proton_left" , 20, -0.1, 0.3);

	histosTH1F["t_proton_plus"] = new TH1F("t_proton_plus", "t_proton_plus" , 11, tbins);
	histosTH1F["t_proton_plus_rec"] = new TH1F("t_proton_plus_rec", "t_proton_plus" , 11, tbins);
	histosTH1F["t_proton_minus"] = new TH1F("t_proton_minus", "t_proton_minus" , 11, tbins);
	histosTH1F["t_proton_minus_rec"] = new TH1F("t_proton_minus_rec", "t_proton_minus" , 11, tbins);
	histosTH1F["t_proton_right_zb"] = new TH1F("t_proton_right_zb", "t_proton_right" , 11, tbins);
	histosTH1F["t_proton_right_zb_sel"] = new TH1F("t_proton_right_zb_sel", "t_proton_right" , 11, tbins);
	histosTH1F["t_proton_left_zb"] = new TH1F("t_proton_left_zb", "t_proton_left" , 11, tbins);
	histosTH1F["t_proton_left_zb_sel"] = new TH1F("t_proton_left_zb_sel", "t_proton_left" , 11, tbins);
	histosTH1F["t_proton_signal_right"] = new TH1F("t_proton_signal_right", "t_proton_right" , 11, tbins);
	histosTH1F["t_proton_signal_right_kin_cut"] = new TH1F("t_proton_signal_right_kin_cut", "t_proton_right" , 11, tbins);
	histosTH1F["t_proton_signal_right_kin_cut_halo"] = new TH1F("t_proton_signal_right_kin_cut_halo", "t_proton_right" , 11, tbins);
	histosTH1F["t_proton_backg_right"] = new TH1F("t_proton_backg_right", "t_proton_right" , 11, tbins);
	histosTH1F["t_proton_backg_right_sel"] = new TH1F("t_proton_backg_right_sel", "t_proton_right" , 11, tbins);
	histosTH1F["t_proton_backg_right_kin_cut"] = new TH1F("t_proton_backg_right_kin_cut", "t_proton_right" , 11, tbins);
	histosTH1F["t_proton_backg_right_kin_cut_halo"] = new TH1F("t_proton_backg_right_kin_cut_halo", "t_proton_right" , 11, tbins);
	histosTH1F["t_proton_backg_right_kin"] = new TH1F("t_proton_backg_right_kin", "t_proton_right" , 11, tbins);
	histosTH1F["t_proton_signal_left"] = new TH1F("t_proton_signal_left", "t_proton_left" , 11, tbins);
	histosTH1F["t_proton_backg_left"] = new TH1F("t_proton_backg_left", "t_proton_left" , 11, tbins);
	histosTH1F["t_proton_right_both"] = new TH1F("t_proton_right_both", "t_proton_right" , 11, tbins);
	histosTH1F["t_proton_left_both"] = new TH1F("t_proton_left_both", "t_proton_left" , 11, tbins);
	histosTH1F["t_proton_both_right_kin_cut"] = new TH1F("t_proton_backg_both_kin_cut", "t_proton_minus" , 11, tbins);
	histosTH1F["t_proton_both_right_kin_cut_halo"] = new TH1F("t_proton_backg_both_kin_cut_halo", "t_proton_minus" , 11, tbins);
	histosTH1F["t_proton_both_left_kin_cut"] = new TH1F("t_proton_left_backg_both_kin_cut", "t_proton_plus" , 11, tbins);
	histosTH1F["t_proton_both_left_kin_cut_halo"] = new TH1F("t_proton_left_backg_both_kin_cut_halo", "t_proton_plus" , 11, tbins);


	histosTH1F["t_proton_signal_left_kin_cut"] = new TH1F("t_proton_signal_left_kin_cut", "t_proton_left" , 11, tbins);
	histosTH1F["t_proton_signal_left_kin_cut_halo"] = new TH1F("t_proton_signal_left_kin_cut_halo", "t_proton_left" , 11, tbins);
	histosTH1F["t_proton_backg_left_sel"] = new TH1F("t_proton_backg_left_sel", "t_proton_left" , 11, tbins);
	histosTH1F["t_proton_backg_left_kin_cut"] = new TH1F("t_proton_backg_left_kin_cut", "t_proton_left" , 11, tbins);
	histosTH1F["t_proton_backg_left_kin_cut_halo"] = new TH1F("t_proton_backg_left_kin_cut_halo", "t_proton_left" , 11, tbins);
	histosTH1F["t_proton_backg_left_kin"] = new TH1F("t_proton_backg_left_kin", "t_proton_left" , 11, tbins);

	histosTH1F["xi_cms_totem_right_signal"] = new TH1F("xi_cms_totem_right_signal", "xi_cms_totem_right" , 20, -0.1, 0.3);
	histosTH1F["xi_cms_totem_right_signal_kin"] = new TH1F("xi_cms_totem_right_signal_kin", "xi_cms_totem_right" , 20, -0.4, 0.4);
	histosTH1F["xi_cms_totem_right_zb"] = new TH1F("xi_cms_totem_right_zb", "xi_cms_totem_right" , 20, -0.1, 0.3);
	histosTH1F["xi_cms_totem_right_zb_kin"] = new TH1F("xi_cms_totem_right_zb_kin", "xi_cms_totem_right" , 20, -0.4, 0.4);
	histosTH1F["xi_cms_totem_right_both"] = new TH1F("xi_cms_totem_right_both", "xi_cms_totem_right" , 20, -0.1, 0.3);
	histosTH1F["xi_cms_totem_right_both_kin"] = new TH1F("xi_cms_totem_right_both_kin", "xi_cms_totem_right" , 20, -0.4, 0.4);
	histosTH1F["xi_cms_totem_left_both"] = new TH1F("xi_cms_totem_left_both", "xi_cms_totem_left" , 20, -0.1, 0.3);
	histosTH1F["xi_cms_totem_left_both_kin"] = new TH1F("xi_cms_totem_left_both_kin", "xi_cms_totem_left" , 20, -0.4, 0.4);
	histosTH1F["xi_cms_totem_left_zb"] = new TH1F("xi_cms_totem_left_zb", "xi_cms_totem_left" , 20, -0.1, 0.3);
	histosTH1F["xi_cms_totem_left_signal"] = new TH1F("xi_cms_totem_left_signal", "xi_cms_totem_left" , 20, -0.1, 0.3);
	histosTH1F["xi_cms_totem_left_signal_kin"] = new TH1F("xi_cms_totem_left_signal_kin", "xi_cms_totem_left" , 20, -0.4, 0.4);
	histosTH1F["xi_cms_totem_left_zb_kin"] = new TH1F("xi_cms_totem_left_zb_kin", "xi_cms_totem_left" , 20, -0.4, 0.4);

	histosTH1F["beta_proton_minus_signal"] = new TH1F("beta_proton_minus_signal", "beta_proton_minus" , 20, 0, 1 );
	histosTH1F["beta_proton_minus_signal_kin_cut"] = new TH1F("beta_proton_minus_signal_kin_cut", "beta_proton_minus" , 20, 0, 1 );
	histosTH1F["beta_proton_minus_backg"] = new TH1F("beta_proton_minus_backg", "beta_proton_minus" , 20, 0, 1 );
	histosTH1F["beta_proton_minus_backg_kin_cut"] = new TH1F("beta_proton_minus_backg_kin_cut", "beta_proton_minus" , 20, 0, 1 );
	histosTH1F["beta_proton_plus_signal"] = new TH1F("beta_proton_plus_signal", "beta_proton_plus" , 20, 0, 1 );
	histosTH1F["beta_proton_plus_signal_kin_cut"] = new TH1F("beta_proton_plus_signal_kin_cut", "beta_proton_plus" , 20, 0, 1 );
	histosTH1F["beta_proton_plus_backg"] = new TH1F("beta_proton_plus_backg", "beta_proton_plus" , 20, 0, 1 );
	histosTH1F["beta_proton_plus_backg_kin_cut"] = new TH1F("beta_proton_plus_backg_kin_cut", "beta_proton_plus" , 20, 0, 1 );
	histosTH1F["beta_proton_both_right"] = new TH1F("beta_proton_both_right", "beta_proton_plus" , 20, 0, 1 );
	histosTH1F["beta_proton_both_right_kin_cut"] = new TH1F("beta_proton_both_right_kin_cut", "beta_proton_plus" , 20, 0, 1 );

	///
	histosTH1F["jpsi_cms_sumEHFminus"] = new TH1F("jpsi_cms_sumEHFminus","sumEHF-",500, 0.,500.);
	histosTH1F["jpsi_cms_sumEHFplus"] = new TH1F("jpsi_cms_sumEHFplus","sumEHF+",500,0.,500.);
	histosTH1F["cms_sumEHFminus"] = new TH1F("cms_sumEHFminus","sumEHF-",500, 0.,500.);
	histosTH1F["cms_sumEHFplus"] = new TH1F("cms_sumEHFplus","sumEHF+",500, 0.,500.);

	histosTH1F["pf_xiPlus"] = new TH1F("pf_xiPlus","#xi^{+}",200,-1.,1.);
	histosTH1F["pf_xiMinus"] = new TH1F("pf_xiMinus","#xi^{-}",200,-1.,1.);
	histosTH1F["xi_cms_pfplus"] =  new TH1F("xi_cms_pfplus","#xi^{+}",200,-1.,1.);
	histosTH1F["xi_cms_pfminus"] =  new TH1F("xi_cms_pfminus","#xi^{-}",200,-1.,1.);


	histosTH1F["jpsi_xi_cms_pfplus"] =  new TH1F("jpsi_xi_cms_pfplus","#xi^{+}",200,-1.,1.);
	histosTH1F["jpsi_xi_cms_pfminus"] =  new TH1F("jpsi_xi_cms_pfminus","#xi^{-}",200,-1.,1.);
	histosTH1F["jpsi_Eta_max"]= new TH1F("Jpsi_Eta_max", "#eta^{max}" , 82 , etaBinsHCALBoundaries);
	histosTH1F["jpsi_Eta_min"]= new TH1F("Jpsi_Eta_min", "#eta^{min}" , 82 , etaBinsHCALBoundaries);
	histosTH1F["jpsi_Delta_eta_maxmin"] = new TH1F("Jpsi_Delta_eta_maxmin", "#eta^{max} - #eta^{min}" , 20 , 0,11);
	//histosTH1F["xi_cms_pfminus"] =  new TH1F("xi_cms_pfminus","#xi^{-}",20,0.,1.);
	histosTH1F["pf_logXiPlus"] = new TH1F("pf_logXiPlus","log(#xi^{+})",20,-4.,0.);
	histosTH1F["pf_logXiMinus"] = new TH1F("pf_logXiMinus","log(#xi^{-})",20,-4.,0.);

	histosTH1F["thx_proton_plus"] = new TH1F("thx_proton_plus", "thx_proton_plus" , 200, -5e-4, 5e-4);
	histosTH1F["thy_proton_plus"] = new TH1F("thy_proton_plus", "thy_proton_plus" , 200, -5e-4, 5e-4);


	histosTH1F["thx_proton_minus"] = new TH1F("thx_proton_minus", "thx_proton_minus" , 200, -5e-4, 5e-4);
	histosTH1F["thy_proton_minus"] = new TH1F("thy_proton_minus", "thy_proton_minus" , 200, -5e-4, 5e-4);

	histosTH1F["xi_proton_t_range_plus"] = new TH1F("xi_proton_t_range_plus", "xi_proton_plus" , 200, -1., 1.);
	histosTH1F["t_proton_xi_range_plus"] = new TH1F("t_proton_xi_range_plus", "t_proton_plus" , 100, 0., 5.);

	histosTH1F["xi_proton_t_range_minus"] = new TH1F("xi_proton_t_range_minus", "xi_proton_minus" , 200, -1., 1.);
	histosTH1F["t_proton_xi_range_minus"] = new TH1F("t_proton_xi_range_minus", "t_proton_minus" , 100, 0., 5.);

	//FIXME
	histosTH1F["xi_proton_plus_accepted"] = new TH1F("xi_proton_plus_accepted", "xi_proton_plus" , 200, -1., 1.);
	histosTH1F["t_proton_plus_accepted"] = new TH1F("t_proton_plus_accepted", "t_proton_plus" , 100, 0., 5.);
	histosTH1F["xi_proton_minus_accepted"] = new TH1F("xi_proton_minus_accepted", "xi_proton_minus" , 200, -1., 1.);
	histosTH1F["t_proton_minus_accepted"] = new TH1F("t_proton_minus_accepted", "t_proton_minus" , 100, 0., 5.);

	histosTH1F["xi_proton_t_range_plus_accepted"] = new TH1F("xi_proton_t_range_plus_accepted", "xi_proton_plus" , 200, -1., 1.);
	histosTH1F["t_proton_xi_range_plus_accepted"] = new TH1F("t_proton_xi_range_plus_accepted", "t_proton_plus" , 100, 0., 5.);
	histosTH1F["xi_proton_t_range_minus_accepted"] = new TH1F("xi_proton_t_range_minus_accepted", "xi_proton_minus" , 200, -1., 1.);
	histosTH1F["t_proton_xi_range_minus_accepted"] = new TH1F("t_proton_xi_range_minus_accepted", "t_proton_minus" , 100, 0., 5.);

	histosTH1F["xi_proton_plus_selected"] = new TH1F("xi_proton_plus_selected", "xi_proton_plus" , 200, -1., 1.);
	histosTH1F["t_proton_plus_selected"] = new TH1F("t_proton_plus_selected", "t_proton_plus" , 100, 0., 5.);
	histosTH1F["xi_proton_minus_selected"] = new TH1F("xi_proton_minus_selected", "xi_proton_minus" , 200, -1., 1.);
	histosTH1F["t_proton_minus_selected"] = new TH1F("t_proton_minus_selected", "t_proton_minus" , 100, 0., 5.);

	histosTH1F["xi_proton_t_range_plus_selected"] = new TH1F("xi_proton_t_range_plus_selected", "xi_proton_plus" , 200, -1., 1.);
	histosTH1F["t_proton_xi_range_plus_selected"] = new TH1F("t_proton_xi_range_plus_selected", "t_proton_plus" , 100, 0., 5.);
	histosTH1F["xi_proton_t_range_minus_selected"] = new TH1F("xi_proton_t_range_minus_selected", "xi_proton_minus" , 200, -1., 1.);
	histosTH1F["t_proton_xi_range_minus_selected"] = new TH1F("t_proton_xi_range_minus_selected", "t_proton_minus" , 100, 0., 5.);

	// RP stations
	histosTH1F["xi_proton_plus_accepted_020"] = new TH1F("xi_proton_plus_accepted_020", "xi_proton_plus" , 200, -1., 1.);

	histosTH1F["t_proton_plus_accepted_020"] = new TH1F("t_proton_plus_accepted_020", "t_proton_plus" , 100, 0., 5.);
	histosTH1F["posx_proton_plus_accepted_020"] = new TH1F("posx_proton_plus_accepted_020", "posx_proton_plus" , 200, -0.05, 0.05);
	histosTH1F["posy_proton_plus_accepted_020"] = new TH1F("posy_proton_plus_accepted_020", "posy_proton_plus" , 200, -0.05, 0.05);

	histosTH1F["xi_proton_plus_accepted_021"] = new TH1F("xi_proton_plus_accepted_021", "xi_proton_plus" , 200, -1., 1.);
	histosTH1F["t_proton_plus_accepted_021"] = new TH1F("t_proton_plus_accepted_021", "t_proton_plus" , 100, 0., 5.);
	histosTH1F["xi_proton_plus_accepted_022"] = new TH1F("xi_proton_plus_accepted_022", "xi_proton_plus" , 200, -1., 1.);
	histosTH1F["t_proton_plus_accepted_022"] = new TH1F("t_proton_plus_accepted_022", "t_proton_plus" , 100, 0., 5.);
	histosTH1F["xi_proton_plus_accepted_023"] = new TH1F("xi_proton_plus_accepted_023", "xi_proton_plus" , 200, -1., 1.);
	histosTH1F["t_proton_plus_accepted_023"] = new TH1F("t_proton_plus_accepted_023", "t_proton_plus" , 100, 0., 5.);
	histosTH1F["xi_proton_plus_accepted_024"] = new TH1F("xi_proton_plus_accepted_024", "xi_proton_plus" , 200, -1., 1.);
	histosTH1F["t_proton_plus_accepted_024"] = new TH1F("t_proton_plus_accepted_024", "t_proton_plus" , 100, 0., 5.);
	histosTH1F["xi_proton_plus_accepted_025"] = new TH1F("xi_proton_plus_accepted_025", "xi_proton_plus" , 200, -1., 1.);
	histosTH1F["t_proton_plus_accepted_025"] = new TH1F("t_proton_plus_accepted_025", "t_proton_plus" , 100, 0., 5.);
	histosTH1F["xi_proton_plus_accepted_120"] = new TH1F("xi_proton_plus_accepted_120", "xi_proton_plus" , 200, -1., 1.);
	histosTH1F["t_proton_plus_accepted_120"] = new TH1F("t_proton_plus_accepted_120", "t_proton_plus" , 100, 0., 5.);

	histosTH1F["xi_proton_minus_accepted_120"] = new TH1F("xi_proton_minus_accepted_120", "xi_proton_minus" , 200, -1., 1.);
	histosTH1F["t_proton_minus_accepted_120"] = new TH1F("t_proton_minus_accepted_120", "t_proton_minus" , 100, 0., 5.);
	histosTH1F["posx_proton_minus_accepted_120"] = new TH1F("posx_proton_minus_accepted_120", "posx_proton_minus" , 200, -0.05, 0.05);
	histosTH1F["posy_proton_minus_accepted_120"] = new TH1F("posy_proton_minus_accepted_120", "posy_proton_minus" , 200, -0.05, 0.05);

	histosTH1F["xi_proton_minus_accepted_020"] = new TH1F("xi_proton_minus_accepted_020", "xi_proton_minus" , 200, -1., 1.);
	histosTH1F["t_proton_minus_accepted_020"] = new TH1F("t_proton_minus_accepted_020", "t_proton_minus" , 100, 0., 5.);

	histosTH2F["t2_track_multiplicity_vs_track_multiplicity"] = new TH2F("t2_track_multiplicity_vs_track_multiplicity","t2_track_multiplicity_vs_track_multiplicity", 100 , 0 , 100, 100 , 0 , 100);
	histosTH2F["t2_track_multiplicity_vs_leadingJet_pt"] = new TH2F("t2_track_multiplicity_vs_leadingJet_pt","t2_track_multiplicity_vs_leadingJet_pt", 150 , 0. , 150., 100 , 0 , 100);
	histosTH2F["t2_track_entryY_vs_entryX_zplus"] = new TH2F("t2_track_entryY_vs_entryX_zplus","t2_track_entryY_vs_entryX", 160 , -160. , 160., 160 , -160. , 160.);
	histosTH2F["t2_track_entryY_vs_entryX_zminus"] = new TH2F("t2_track_entryY_vs_entryX_zminus","t2_track_entryY_vs_entryX", 160 , -160. , 160., 160 , -160. , 160.);

	histosTH2F["proton_right_logXi_vs_pf_logXiPlus"] = new TH2F("proton_right_logXi_vs_pf_logXiPlus","proton_right_logXi_vs_pf_logXiPlus", 200, -5., 0., 200, -5., 0.);
	histosTH2F["proton_left_logXi_vs_pf_logXiMinus"] = new TH2F("proton_left_logXi_vs_pf_logXiMinus","proton_left_logXi_vs_pf_logXiMinus", 200, -5., 0., 200, -5., 0.);
	histosTH2F["proton_right_logXi_vs_pf_logXiMinus"] = new TH2F("proton_right_logXi_vs_pf_logXiMinus","proton_right_logXi_vs_pf_logXiMinus", 200, -5., 0., 200, -5., 0.);
	histosTH2F["proton_left_logXi_vs_pf_logXiPlus"] = new TH2F("proton_left_logXi_vs_pf_logXiPlus","proton_left_logXi_vs_pf_logXiPlus", 200, -5., 0., 200, -5., 0.);
	histosTH2F["proton_right_logXi_vs_t"] = new TH2F("proton_right_logXi_vs_t","proton_right_logXi_vs_t", 200, 0., 5., 200, -5., 0.);
	histosTH2F["proton_left_logXi_vs_t"] = new TH2F("proton_left_logXi_vs_t","proton_left_logXi_vs_t", 200, 0., 5., 250, -5., 0.);
	histosTH2F["proton_right_xi_vs_pf_xiMinus"] = new TH2F("proton_right_xi_vs_pf_xiMinus","proton_right_xi_vs_pf_xiMinus", 125, 0., 1., 125, 0., 1.);
	histosTH2F["proton_right_xi_vs_pf_xiPlus"] = new TH2F("proton_right_xi_vs_pf_xiPlus","proton_right_xi_vs_pf_xiPlus", 125, 0., 1., 125, 0., 1.);
	histosTH2F["proton_left_xi_vs_pf_xiMinus"] = new TH2F("proton_left_xi_vs_pf_xiMinus","proton_left_xi_vs_pf_xiMinus", 125, 0., 1., 125, 0., 1.);
	histosTH2F["proton_left_xi_vs_pf_xiPlus"] = new TH2F("proton_left_xi_vs_pf_xiPlus","proton_left_xi_vs_pf_xiPlus", 125, 0., 1., 125, 0., 1.);

	//histosTH2F["proton_right_t_vs_leadingJet_pt"] = new TH2F("proton_right_t_vs_leadingJet_pt","proton_right_t_vs_leadingJet_pt", 150 , 0. , 150., 200 , 0. , 5.);
	//histosTH2F["proton_left_t_vs_leadingJet_pt"] = new TH2F("proton_left_t_vs_leadingJet_pt","proton_left_t_vs_leadingJet_pt", 150 , 0. , 150., 200 , 0. , 5.);

	histosTH2F["rp_track_pos_y_vs_x_020"] = new TH2F("rp_track_pos_y_vs_x_020", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);
	histosTH2F["rp_track_pos_y_vs_x_021"] = new TH2F("rp_track_pos_y_vs_x_021", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);
	histosTH2F["rp_track_pos_y_vs_x_022"] = new TH2F("rp_track_pos_y_vs_x_022", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);
	histosTH2F["rp_track_pos_y_vs_x_023"] = new TH2F("rp_track_pos_y_vs_x_023", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);
	histosTH2F["rp_track_pos_y_vs_x_024"] = new TH2F("rp_track_pos_y_vs_x_024", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);
	histosTH2F["rp_track_pos_y_vs_x_025"] = new TH2F("rp_track_pos_y_vs_x_025", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);
	histosTH2F["rp_track_pos_y_vs_x_120"] = new TH2F("rp_track_pos_y_vs_x_120", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);
	histosTH2F["rp_track_pos_y_vs_x_121"] = new TH2F("rp_track_pos_y_vs_x_121", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);
	histosTH2F["rp_track_pos_y_vs_x_122"] = new TH2F("rp_track_pos_y_vs_x_122", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);
	histosTH2F["rp_track_pos_y_vs_x_123"] = new TH2F("rp_track_pos_y_vs_x_123", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);
	histosTH2F["rp_track_pos_y_vs_x_124"] = new TH2F("rp_track_pos_y_vs_x_124", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);
	histosTH2F["rp_track_pos_y_vs_x_125"] = new TH2F("rp_track_pos_y_vs_x_125", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);

	histosTH2F["proton_plus_xi_vs_t_accepted"] = new TH2F("proton_plus_xi_vs_t_accepted","proton_plus_xi_vs_t", 100, 0., 5., 200, -1., 1.);
	histosTH2F["proton_minus_xi_vs_t_accepted"] = new TH2F("proton_minus_xi_vs_t_accepted","proton_minus_xi_vs_t", 100, 0., 5., 200, -1., 1.);

	histosTH2F["proton_plus_xi_vs_t_selected"] = new TH2F("proton_plus_xi_vs_t_selected","proton_plus_xi_vs_t", 100, 0., 5., 200, -1., 1.);
	histosTH2F["proton_minus_xi_vs_t_selected"] = new TH2F("proton_minus_xi_vs_t_selected","proton_minus_xi_vs_t", 100, 0., 5., 200, -1., 1.);

	histosTH2F["pos_y_vs_x_proton_plus_020"] = new TH2F("pos_y_vs_x_proton_plus_020", "pos_y_vs_x_proton_plus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_plus_021"] = new TH2F("pos_y_vs_x_proton_plus_021", "pos_y_vs_x_proton_plus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_plus_022"] = new TH2F("pos_y_vs_x_proton_plus_022", "pos_y_vs_x_proton_plus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_plus_023"] = new TH2F("pos_y_vs_x_proton_plus_023", "pos_y_vs_x_proton_plus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_plus_024"] = new TH2F("pos_y_vs_x_proton_plus_024", "pos_y_vs_x_proton_plus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_plus_025"] = new TH2F("pos_y_vs_x_proton_plus_025", "pos_y_vs_x_proton_plus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_minus_120"] = new TH2F("pos_y_vs_x_proton_minus_120", "pos_y_vs_x_proton_minus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_minus_121"] = new TH2F("pos_y_vs_x_proton_minus_121", "pos_y_vs_x_proton_minus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_minus_122"] = new TH2F("pos_y_vs_x_proton_minus_122", "pos_y_vs_x_proton_minus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_minus_123"] = new TH2F("pos_y_vs_x_proton_minus_123", "pos_y_vs_x_proton_minus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_minus_124"] = new TH2F("pos_y_vs_x_proton_minus_124", "pos_y_vs_x_proton_minus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_minus_125"] = new TH2F("pos_y_vs_x_proton_minus_125", "pos_y_vs_x_proton_minus" , 400, -0.05, 0.05, 400, -0.05, 0.05);

	histosTH2F["pos_thy_vs_thx_proton_plus_020"] = new TH2F("pos_thy_vs_thx_proton_plus_020", "pos_thy_vs_thx_proton_plus" , 100, -5e-4, 5e-4, 100, -5e-4, 5e-4);
	histosTH2F["pos_thy_vs_thx_proton_plus_021"] = new TH2F("pos_thy_vs_thx_proton_plus_021", "pos_thy_vs_thx_proton_plus" , 100, -5e-4, 5e-4, 100, -5e-4, 5e-4);
	histosTH2F["pos_thy_vs_thx_proton_plus_022"] = new TH2F("pos_thy_vs_thx_proton_plus_022", "pos_thy_vs_thx_proton_plus" , 100, -5e-4, 5e-4, 100, -5e-4, 5e-4);
	histosTH2F["pos_thy_vs_thx_proton_plus_023"] = new TH2F("pos_thy_vs_thx_proton_plus_023", "pos_thy_vs_thx_proton_plus" , 100, -5e-4, 5e-4, 100, -5e-4, 5e-4);
	histosTH2F["pos_thy_vs_thx_proton_plus_024"] = new TH2F("pos_thy_vs_thx_proton_plus_024", "pos_thy_vs_thx_proton_plus" , 100, -5e-4, 5e-4, 100, -5e-4, 5e-4);
	histosTH2F["pos_thy_vs_thx_proton_plus_025"] = new TH2F("pos_thy_vs_thx_proton_plus_025", "pos_thy_vs_thx_proton_plus" , 100, -5e-4, 5e-4, 100, -5e-4, 5e-4);
	histosTH2F["pos_thy_vs_thx_proton_minus_120"] = new TH2F("pos_thy_vs_thx_proton_minus_120", "pos_thy_vs_thx_proton_minus" , 100, -5e-4, 5e-4, 100, -5e-4, 5e-4);
	histosTH2F["pos_thy_vs_thx_proton_minus_121"] = new TH2F("pos_thy_vs_thx_proton_minus_121", "pos_thy_vs_thx_proton_minus" , 100, -5e-4, 5e-4, 100, -5e-4, 5e-4);
	histosTH2F["pos_thy_vs_thx_proton_minus_122"] = new TH2F("pos_thy_vs_thx_proton_minus_122", "pos_thy_vs_thx_proton_minus" , 100, -5e-4, 5e-4, 100, -5e-4, 5e-4);
	histosTH2F["pos_thy_vs_thx_proton_minus_123"] = new TH2F("pos_thy_vs_thx_proton_minus_123", "pos_thy_vs_thx_proton_minus" , 100, -5e-4, 5e-4, 100, -5e-4, 5e-4);
	histosTH2F["pos_thy_vs_thx_proton_minus_124"] = new TH2F("pos_thy_vs_thx_proton_minus_124", "pos_thy_vs_thx_proton_minus" , 100, -5e-4, 5e-4, 100, -5e-4, 5e-4);
	histosTH2F["pos_thy_vs_thx_proton_minus_125"] = new TH2F("pos_thy_vs_thx_proton_minus_125", "pos_thy_vs_thx_proton_minus" , 100, -5e-4, 5e-4, 100, -5e-4, 5e-4);
	//========================================================

	histosTH2F["pos_y_vs_x_proton_plus_accepted_020"] = new TH2F("pos_y_vs_x_proton_plus_accepted_020", "pos_y_vs_x_proton_plus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_plus_accepted_021"] = new TH2F("pos_y_vs_x_proton_plus_accepted_021", "pos_y_vs_x_proton_plus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_plus_accepted_022"] = new TH2F("pos_y_vs_x_proton_plus_accepted_022", "pos_y_vs_x_proton_plus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_plus_accepted_023"] = new TH2F("pos_y_vs_x_proton_plus_accepted_023", "pos_y_vs_x_proton_plus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_plus_accepted_024"] = new TH2F("pos_y_vs_x_proton_plus_accepted_024", "pos_y_vs_x_proton_plus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_plus_accepted_025"] = new TH2F("pos_y_vs_x_proton_plus_accepted_025", "pos_y_vs_x_proton_plus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_minus_accepted_120"] = new TH2F("pos_y_vs_x_proton_minus_accepted_120", "pos_y_vs_x_proton_minus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_minus_accepted_121"] = new TH2F("pos_y_vs_x_proton_minus_accepted_121", "pos_y_vs_x_proton_minus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_minus_accepted_122"] = new TH2F("pos_y_vs_x_proton_minus_accepted_122", "pos_y_vs_x_proton_minus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_minus_accepted_123"] = new TH2F("pos_y_vs_x_proton_minus_accepted_123", "pos_y_vs_x_proton_minus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_minus_accepted_124"] = new TH2F("pos_y_vs_x_proton_minus_accepted_124", "pos_y_vs_x_proton_minus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_minus_accepted_125"] = new TH2F("pos_y_vs_x_proton_minus_accepted_125", "pos_y_vs_x_proton_minus" , 400, -0.05, 0.05, 400, -0.05, 0.05);

	histosTH2F["pos_thy_vs_thx_proton_plus_accepted_020"] = new TH2F("pos_thy_vs_thx_proton_plus_accepted_020", "pos_thy_vs_thx_proton_plus" ,  100, -5e-4, 5e-4, 100, -5e-4, 5e-4);
	histosTH2F["pos_thy_vs_thx_proton_plus_accepted_021"] = new TH2F("pos_thy_vs_thx_proton_plus_accepted_021", "pos_thy_vs_thx_proton_plus" ,  100, -5e-4, 5e-4, 100, -5e-4, 5e-4);
	histosTH2F["pos_thy_vs_thx_proton_plus_accepted_022"] = new TH2F("pos_thy_vs_thx_proton_plus_accepted_022", "pos_thy_vs_thx_proton_plus" , 100, -5e-4, 5e-4, 100, -5e-4, 5e-4);
	histosTH2F["pos_thy_vs_thx_proton_plus_accepted_023"] = new TH2F("pos_thy_vs_thx_proton_plus_accepted_023", "pos_thy_vs_thx_proton_plus" ,  100, -5e-4, 5e-4, 100, -5e-4, 5e-4);
	histosTH2F["pos_thy_vs_thx_proton_plus_accepted_024"] = new TH2F("pos_thy_vs_thx_proton_plus_accepted_024", "pos_thy_vs_thx_proton_plus" ,  100, -5e-4, 5e-4, 100, -5e-4, 5e-4);
	histosTH2F["pos_thy_vs_thx_proton_plus_accepted_025"] = new TH2F("pos_thy_vs_thx_proton_plus_accepted_025", "pos_thy_vs_thx_proton_plus" , 100, -5e-4, 5e-4, 100, -5e-4, 5e-4);
	histosTH2F["pos_thy_vs_thx_proton_minus_accepted_120"] = new TH2F("pos_thy_vs_thx_proton_minus_accepted_120", "pos_thy_vs_thx_proton_minus" , 100, -5e-4, 5e-4, 100, -5e-4, 5e-4);
	histosTH2F["pos_thy_vs_thx_proton_minus_accepted_121"] = new TH2F("pos_thy_vs_thx_proton_minus_accepted_121", "pos_thy_vs_thx_proton_minus" , 100, -5e-4, 5e-4, 100, -5e-4, 5e-4);
	histosTH2F["pos_thy_vs_thx_proton_minus_accepted_122"] = new TH2F("pos_thy_vs_thx_proton_minus_accepted_122", "pos_thy_vs_thx_proton_minus" ,  100, -5e-4, 5e-4, 100, -5e-4, 5e-4);
	histosTH2F["pos_thy_vs_thx_proton_minus_accepted_123"] = new TH2F("pos_thy_vs_thx_proton_minus_accepted_123", "pos_thy_vs_thx_proton_minus" ,  100, -5e-4, 5e-4, 100, -5e-4, 5e-4);
	histosTH2F["pos_thy_vs_thx_proton_minus_accepted_124"] = new TH2F("pos_thy_vs_thx_proton_minus_accepted_124", "pos_thy_vs_thx_proton_minus" ,  100, -5e-4, 5e-4, 100, -5e-4, 5e-4);
	histosTH2F["pos_thy_vs_thx_proton_minus_accepted_125"] = new TH2F("pos_thy_vs_thx_proton_minus_accepted_125", "pos_thy_vs_thx_proton_minus" ,  100, -5e-4, 5e-4, 100, -5e-4, 5e-4);
	//==========================================================


	double energyMin = -10.;
	double energyMax = 190.;
	int nBinsEnergy = 1000;
	histosTH2F["energyVsEtaAllTypes"] = new TH2F("energyVsEtaAllTypes","energyVsEtaAllTypes",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
	histosTH2F["energyVsEtaUndefined"] = new TH2F("energyVsEtaUndefined","energyVsEtaUndefined",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
	histosTH2F["energyVsEtaChargedHadron"] = new TH2F("energyVsEtaChargedHadron","energyVsEtaChargedHadron",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
	histosTH2F["energyVsEtaElectron"] = new TH2F("energyVsEtaElectron","energyVsEtaElectron",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
	histosTH2F["energyVsEtaMuon"] = new TH2F("energyVsEtaMuon","energyVsEtaMuon",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
	histosTH2F["energyVsEtaGamma"] = new TH2F("energyVsEtaGamma","energyVsEtaGamma",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
	histosTH2F["energyVsEtaNeutralHadron"] = new TH2F("energyVsEtaNeutralHadron","energyVsEtaNeutralHadron",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
	histosTH2F["energyVsEtaHadronHF"] = new TH2F("energyVsEtaHadronHF","energyVsEtaHadronHF",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
	histosTH2F["energyVsEtaHadronHFEcalEnergy"] = new TH2F("energyVsEtaHadronHFEcalEnergy","energyVsEtaHadronHFEcalEnergy",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
	histosTH2F["energyVsEtaHadronHFNoEcalEnergy"] = new TH2F("energyVsEtaHadronHFNoEcalEnergy","energyVsEtaHadronHFNoEcalEnergy",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
	histosTH2F["energyVsEtaEGammaHF"] = new TH2F("energyVsEtaEGammaHF","energyVsEtaEGammaHF",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
	histosTH2F["xi_plus_reco_gen"] = new TH2F("xi+","xi+",82,0,0.5,82,0,0.5);
	histosTH2F["xi_minus_reco_gen"] = new TH2F("xi-","xi-",82,0,0.5,82,0,0.5);
	histosTH2F["logxi_plus_reco_gen"] = new TH2F("logxi+","xi+",82,-3,0.5,82,-3,0.5);
	histosTH2F["logxi_minus_reco_gen"] = new TH2F("logxi-","xi-",82,-3,0.5,82,-3,0.5);


	histosTH2F["proton_plus_xi_vs_t"] = new TH2F("proton_plus_xi_vs_t","proton_plus_xi_vs_t", 200, 0., 5., 200, 0., 1.);
	histosTH2F["proton_minus_xi_vs_t"] = new TH2F("proton_minus_xi_vs_t","proton_minus_xi_vs_t", 200, 0., 5., 200, 0., 1.);

	for(map<string,TH1F*>::const_iterator it = histosTH1F.begin(); it != histosTH1F.end(); ++it)
		it->second->Sumw2();
	for(map<string,TH2F*>::const_iterator it = histosTH2F.begin(); it != histosTH2F.end(); ++it)
		it->second->Sumw2();
	gStyle->SetPalette(1);
	//===================
	int i_tot = 0 , nevt_tot = 0;

	//MC Files

    const char *dirname="/storage2/eliza/TOTEM/PYTHIA_Jpsi_UATree/";
	const char *ext=".root";
	vector<TString>* vfiles = new vector<TString>; 

	TSystemDirectory dir(dirname, dirname);
	TList *files = dir.GetListOfFiles();
	if (files) {
		TSystemFile *file;
		TString fname;
		TIter next(files);
		while ((file=(TSystemFile*)next())) {
			fname = file->GetName();
			if (!file->IsDirectory() && fname.EndsWith(ext)) {
				TString root_file = dirname + string(fname.Data());
				vfiles->push_back(root_file); cout<<root_file<<endl;      
			}
		}
	}

	//Declaration of tree and its branches variables
	TTree* tree = new TTree(treeName.c_str(),"");
	MyEvtId*           evtId        = NULL;
	//   MyL1TrigOld*       l1Trig       = NULL;  
	//   MyHLTrig*          hltTrig      = NULL;
	vector<MyGenPart>* genPart      = NULL;
	vector<MyTracks>*  track_coll   = NULL;
	vector<MyVertex>*  vertex_coll  = NULL;
	//vector<MyPFJet>*   pfJet_coll   = NULL;
	vector<MyMuon>*    muon_coll   = NULL;
	vector<MyPFCand>*  pFlow_coll   = NULL;
	MyGenKin*  genKin   = NULL;
	//===============================

	//ZeroBias Files
	string treeNameZB = "cms_totem";
	TChain treeZB("cms_totem");
	vector<TString>* vdirs_zb = new vector<TString>;
	//vdirs_zb->push_back("/storage1/lhuertas/Analysis/CMSTotem/data/MergedNtuples/HighBeta/ZeroBias/");
	vdirs_zb->push_back("/storage1/eliza/TOTEM/ZeroBias/total/");
	vector<TString>* vfiles_zb = new vector<TString>;
	for(vector<TString>::iterator itdirs_zb = vdirs_zb->begin(); itdirs_zb != vdirs_zb->end(); ++itdirs_zb){
		TString& dirname_zb = *itdirs_zb;
		TSystemDirectory dir(dirname_zb, dirname_zb);
		TList *files = dir.GetListOfFiles();
		if (files) {
			TSystemFile *file;
			TString fname;
			TIter next(files);
			while ((file=(TSystemFile*)next())) {
				fname = file->GetName();
				if (!file->IsDirectory() && fname.EndsWith(ext)) {
					TString root_file_zb = dirname_zb + string(fname.Data());
					vfiles_zb->push_back(root_file_zb); cout<<root_file_zb<<endl;
					treeZB.Add(root_file_zb);
				}
			}
		}
	}//cout<<"n_events_total ZB "<<treeZB.GetEntries()<<endl;

	//TTree* treeZB = new TTree(treeNameZB.c_str(),"");
	MyEvtId*           evtIdZB        = NULL;
	RPRootDumpReconstructedProton* rec_proton_left  = NULL;
	RPRootDumpReconstructedProton* rec_proton_right = NULL;
	vector<MyVertex>*  vertex_coll_ZB  = NULL;
	//  TFile* fileZB = TFile::Open(fileNameZB.c_str(),"READ");

	//  treeZB = (TTree*)fileZB->Get( treeNameZB.c_str() );
	int nevZB = int(treeZB.GetEntries());

	treeZB.SetBranchAddress("cmsEvtUA",&evtIdZB);
	treeZB.SetBranchAddress("rec_prot_left.",&rec_proton_left);
	treeZB.SetBranchAddress("rec_prot_right.",&rec_proton_right);
	treeZB.SetBranchAddress("cmsVerticesUA",&vertex_coll_ZB);

	std::vector<MyZeroBiasData> zeroBias;
	//  int nevt_max_ZB = 500000;
	for(int i_evt = 0; i_evt < nevZB && i_evt < nevt_max_corr; ++i_evt){

		if( ((i_evt+1) % 10000) == 0) cout <<int(double(i_evt+1)/1000)<<"k done"<<endl;

		//Filling the variables defined setting branches
		treeZB.GetEntry(i_evt);

		MyZeroBiasData mydata;
		mydata.proton_rec_left_valid = rec_proton_left->valid;
		mydata.proton_rec_left_t = rec_proton_left->t;
		mydata.proton_rec_left_xi = rec_proton_left->xi;
		mydata.proton_rec_right_valid = rec_proton_right->valid;
		mydata.proton_rec_right_t = rec_proton_right->t;
		mydata.proton_rec_right_xi = rec_proton_right->xi;

		MyVertex& primaryVertex = vertex_coll_ZB->at(0);
		/*for(vector<MyVertex>::iterator it_vtx_zb = vertex_coll_ZB->begin() ; it_vtx_zb != vertex_coll_ZB->end() ; ++it_vtx_zb){
		  mydata.vtx_ndof = it_vtx_zb->ndof;    
		  mydata.vtx_valid = it_vtx_zb->validity;
		  }*/
		mydata.vtx_ndof = primaryVertex.ndof;
		mydata.vtx_valid = primaryVertex.validity;
		mydata.vtx_x = primaryVertex.x;
		mydata.vtx_y = primaryVertex.y;
		mydata.vtx_z = primaryVertex.z;

		zeroBias.push_back(mydata);
	}
	//fileZB->Close();

	cout << zeroBias.size() << " events analyzed" << endl;

	// for(vector<MyTOTEMData>::const_iterator it_zb = zeroBias.begin(); it_zb != zeroBias.end(); it_zb++){
	//    cout<<it_zb->proton_rec_right_t<<endl;  
	// }

	//================================================


	rp_aperture_config();

	double nevents_muons = 0; 
	//double nevents_jets = 0; 
	double nevents_pf = 0; 
	double nevents_gen = 0; 
	double nevents_total = 0; 
	//double events_jets = 0;
	double events_muons = 0;
	double events_pf = 0;
	double events_gen = 0;

	double nweight_total = 0; 
	double weight_total_leadingMuon = 0; 
	double weight_total_secondMuon = 0; 
	double weight_total_Muon_selected = 0; 
	double weight_total_PF_selected = 0; 


	int n_vertices_selected = 0;
	int n_select_Vertex_After_vtx_cut =0;
	int n_tracks_selected = 0;
	int n_evt_selectTrack = 0;
	int n_muon_selected = 0;
	int n_dimuons_selected = 0;
	int n_jpsi_selected = 0;
	int n_jpsi_mass = 0;
	int n_dimuons_rp_selected_plus_side = 0;
	int n_dimuons_rp_selected_tsel_plus_side =0;
	int n_dimuons_rp_selected_minus_side = 0;
	int n_dimuons_rp_selected_tsel_minus_side =0;
	int n_evt_PF = 0; 
	int n_jpsi_selected_rp_accept_plus_side = 0;
	int n_jpsi_selected_rp_accept_tsel_plus_side = 0;
	int n_jpsi_selected_rp_accept_minus_side = 0;
	int n_jpsi_selected_rp_accept_tsel_minus_side = 0;
	int nprotons = 0;
	int nevtxisignalleft = 0;
	int nevtxisignalright = 0;
	int nevtxisignalleftjpsi = 0;
	int nevtxisignalrightjpsi = 0; 


	int nBothsLeft=0;
	int nBckgRegionBothsLeft=0;
	int nSignalRegionBothsLeft=0;

	int nBothsRight=0;
	int nBckgRegionBothsRight=0;
	int nSignalRegionBothsRight=0;

	int nBckgRegionZBLeft=0;
	int nSignalRegionZBLeft=0;
	int nZBLeft=0;

	int nBckgRegionZBRight=0;
	int nSignalRegionZBRight=0;
	int nZBRight=0;

	int nBckgRegionMCLeft=0;
	int nSignalRegionMCLeft=0;
	int nMCLeft=0;

	int nBckgRegionMCRight=0;
	int nSignalRegionMCRight=0;
	int nMCRight=0;

	int nZBconditions=0;

	bool PF_eta_max = false;
	bool PF_eta_min = false;
	bool mu1_selected = false;
	bool mu2_selected = false;
	bool jpsi_mass=false;
	//bool chmuons = false;
	bool select_Track = false;	
	//starting Loop over files, stops at end of list of files or when reached nevt_max
	for(vector<TString>::iterator itfiles = vfiles->begin(); itfiles != vfiles->end() && i_tot < nevt_max_corr; ++itfiles){

		TFile* file = TFile::Open(*itfiles,"READ");

		//getting the tree form the current file
		tree = (TTree*) file->Get( treeName.c_str() );

		//Getting number of events
		int nev = int(tree->GetEntriesFast());
		nevt_tot += nev;
		cout <<"The current file has " << nev << " entries : "<< endl << *itfiles << endl;

		//adding branches to the tree ----------------------------------------------------------------------
		tree->SetBranchAddress("evtId",&evtId);
		tree->SetBranchAddress("generalTracks",&track_coll); 
		tree->SetBranchAddress("offlinePrimaryVertices",&vertex_coll);
		tree->SetBranchAddress("muons",&muon_coll);
		tree->SetBranchAddress("particleFlow",&pFlow_coll);
		tree->SetBranchAddress("genKin",&genKin);
		tree->SetBranchAddress("genPart",&genPart);

		//starting loop over events, stops when reached end of file or nevt_max
		for(int i_evt = 0; i_evt < nev && i_tot < nevt_max_corr; ++i_evt , ++i_tot){

			//printing the % of events done every 10k evts
			if( ((i_tot+1) % 10000) == 0) cout <<int(double(i_tot+1)/1000)<<"k done"<<endl;

			//Filling the variables defined setting branches
			tree->GetEntry(i_evt);

			bool passedHLT = false;
			bool passedvtx = false;
			//bool mu1_selected = false;
			//bool mu2_selected = false;
			bool pz_proton_max = false;
			bool PF_eta_max = false;
			bool PF_eta_min = false;
			bool xi_negat_gen = false;
			bool xi_posit_gen = false;

			//AT THIS POINT ON, CAN START USING THE VARIABLES LIKE TRACKS, VERTEX ETC !
			//double event_weight = genKin->genWeight; 
			double event_weight = 1.0;
			nweight_total += event_weight; 
			++nevents_total;


			//-------------------------------------------------------------------------------------------------


			// Vertices
			//double xpos; double ypos; double zpos;
			if(verbose)cout<<"MyVertex"<<endl;
			// Find number of good vertices
			int nGoodVertices = 0;
			// Vertices
			for(vector<MyVertex>::iterator it_vtx = vertex_coll->begin() ; it_vtx != vertex_coll->end() ; ++it_vtx){
				if (it_vtx!=vertex_coll->begin()) continue;
				if( it_vtx->ndof>4 ) passedvtx = true;
				++nGoodVertices;
				if(verbose)cout<<" nGoodVertices : "<<nGoodVertices<<endl;
				//histosTH1F["vtx_zpos"]->Fill( it_vtx->z, event_weight );
				//histosTH1F["vtx_xpos"]->Fill( it_vtx->x, event_weight );
				//histosTH1F["vtx_ypos"]->Fill( it_vtx->y, event_weight );}

		}	
		//histosTH1F["vertex_multiplicity"]->Fill(nGoodVertices, event_weight );
		if(verbose)cout<<"Before nGoodVertices: "<<nGoodVertices<<endl;
		bool select_OneVertex =(nGoodVertices > 0 && nGoodVertices <= 1);
		if (selectVertex && !select_OneVertex) continue;
		if(verbose)cout<<"After nGoodVertices= 1 : "<<nGoodVertices<<endl;
		++n_vertices_selected;	
		if(!passedvtx) continue;


		MyVertex& primaryVertex = vertex_coll->at(0);
		//histosTH1F["prim_vtx_zpos"]->Fill( primaryVertex.z, event_weight );
		//histosTH1F["prim_vtx_xpos"]->Fill( primaryVertex.x, event_weight );
		//histosTH1F["prim_vtx_ypos"]->Fill( primaryVertex.y, event_weight );
		double prim_vtx_r = sqrt( primaryVertex.x*primaryVertex.x + primaryVertex.y*primaryVertex.y );
		bool select_Vertex = ( !primaryVertex.fake && primaryVertex.validity && primaryVertex.ndof > 4 && fabs( primaryVertex.z ) < 15.0 && prim_vtx_r < 2.0);

		if(selectVertex && !select_Vertex) continue;
		//passedvtx = true;  
		++n_select_Vertex_After_vtx_cut;


		/*histosTH1F["vertex_multiplicity_after_vtx_sel"]->Fill( nGoodVertices, event_weight );  
		  histosTH1F["prim_vtx_zpos_after_vtx_sel"]->Fill( primaryVertex.z, event_weight );
		  histosTH1F["prim_vtx_xpos_after_vtx_sel"]->Fill( primaryVertex.x, event_weight );
		  histosTH1F["prim_vtx_ypos_after_vtx_sel"]->Fill( primaryVertex.y, event_weight );

		  histosTH1F["prim_vtx_ndof_after_vtx_sel"]->Fill( primaryVertex.ndof, event_weight );
		  histosTH1F["prim_vtx_chi2_after_vtx_sel"]->Fill( primaryVertex.chi2, event_weight );
		  histosTH1F["prim_vtx_chi2n_after_vtx_sel"]->Fill( primaryVertex.chi2n(), event_weight );
		  histosTH1F["prim_vtx_ntracks_after_vtx_sel"]->Fill( primaryVertex.ntracks, event_weight );
		  histosTH1F["prim_vtx_sumpt_after_vtx_sel"]->Fill( primaryVertex.SumPtTracks, event_weight );*/

		//histosTH1F["EventSelection"]->Fill( "Vertex", event_weight );

		int prim_vtx_id = primaryVertex.id;
		if(verbose)cout<<"final vertex : "<<endl;
		//-------------------------------------------------------------------------------------------------     

		// Tracks

		int n_tracks_selected = 0;
		for(vector<MyTracks>::iterator it_trk = track_coll->begin() ; it_trk != track_coll->end() ; ++it_trk){
			histosTH1F["track_pt"]->Fill( it_trk->Pt(), event_weight );
			histosTH1F["track_eta"]->Fill( it_trk->Eta(), event_weight );
			histosTH1F["track_phi"]->Fill( it_trk->Phi(), event_weight );
			histosTH1F["track_rapidity"]->Fill( it_trk->Rapidity(), event_weight );
			if( it_trk->Pt() < 0.5 ) continue;
			if( fabs( it_trk->Eta() ) > 2.5 ) continue;
			if( ( it_trk->dz / it_trk->edz ) > 5. ) continue;
			if( ( it_trk->d0 / it_trk->ed0 ) > 5. ) continue;
			if( !it_trk->quality[2] ) continue;
			select_Track = true;

			++n_tracks_selected;
		}
		histosTH1F["track_multiplicity"]->Fill( n_tracks_selected, event_weight );
		if(selectTrack && !select_Track) continue;
		++n_evt_selectTrack;
		//-------------------------------------------------------------------------------------------------
		// Muons variables 
		int n_muons = 0; int n_dimuons = 0.;
		double chmu1 = 0.; double chmu2 = 0.;
		double phimu1 = 0.; double ptmu1 = 0.;double etamu1 = 0.;double ymu1 = 0.;
		double phimu2 = 0.; double ptmu2 = 0.;double etamu2 = 0.;double ymu2 = 0.;
		double deltaphi = 0.; double deltaeta = 0.; double deltapt = 0.; double deltay = 0.; double Dphi = 0.;
		double dimuon_mass = 0.; double dimuon_pt=0.; double dimuon_pt2=0.; double dimuon_eta =0.;
		double dimuon_rapidity = 0.; double jpsi_mass = 0.; double jpsi_pt=0.;
		double jpsi_pt2=0.; double jpsi_eta =0.; double jpsi_rapidity = 0.; double dphijpsi = 0.;
		vector<MyMuon> muons_selected;
		for(vector<MyMuon>::iterator it_muon = muon_coll->begin() ; it_muon != muon_coll->end() ; ++it_muon){
			if( !(it_muon->IsTrackerMuon || it_muon->IsGlobalMuon) ) continue;
			MyTracks const& muon_innerTrack = it_muon->innerTrack;
			bool muon_id = it_muon->TMOneStationAngTight &&
				muon_innerTrack.chi2n < 1.8 &&
				muon_innerTrack.nValidPixelHits > 0 &&
				muon_innerTrack.vtxdxy[prim_vtx_id] < 3. &&
				muon_innerTrack.vtxdz[prim_vtx_id] < 30.;

			if( !muon_id ) continue;
			++n_muons; 
			if(verbose) cout <<"nr de muons = " <<n_muons << endl;
			++n_muon_selected;
			if(verbose) cout <<"nr de muons selected = " <<n_muon_selected<< endl;

			histosTH1F["muon_pt"]->Fill( it_muon->Pt(), event_weight );
			histosTH1F["muon_eta"]->Fill( it_muon->Eta(), event_weight );
			histosTH1F["muon_phi"]->Fill( it_muon->Phi(), event_weight );

			muons_selected.push_back( *it_muon );

			//histosTH1F["muon_multiplicity"]->Fill( n_muons, event_weight );
		}
		bool select_Muons = ( muons_selected.size() >= 2 );
		if(selectMuons && !select_Muons) continue;
		//histosTH1F["EventSelection"]->Fill( "Muons", event_weight );
		if(verbose)cout<<"muons selected size (>=2) :"<<muons_selected.size()<<endl;
		//----------------------------------------------------------------------------------------
		// Muon pairs
		for(vector<MyMuon>::iterator it_mu1 = muons_selected.begin() ;
				it_mu1 != muons_selected.end() ; ++it_mu1){
			histosTH1F["muon1_pt"]->Fill( it_mu1->Pt(), event_weight );
			histosTH1F["muon1_eta"]->Fill( it_mu1->Eta(), event_weight );
			histosTH1F["muon1_rapidity"]->Fill( it_mu1->Rapidity(), event_weight );
			histosTH1F["muon1_phi"]->Fill(it_mu1->Phi(), event_weight );
			phimu1 = it_mu1->Phi();  ptmu1 = it_mu1->Pt(); 
			etamu1 = it_mu1->Eta(); ymu1 = it_mu1->Rapidity();
			chmu1 = it_mu1->charge;
			mu1_selected = true;

			for(vector<MyMuon>::iterator it_mu2 = muons_selected.begin() ;
					it_mu2 != muons_selected.end() ; ++it_mu2){
				bool os_muons = ( it_mu1->charge*it_mu2->charge < 0. );
				if( !os_muons ) continue;
				++n_dimuons_selected;
				++n_dimuons;
				mu2_selected = true;

				histosTH1F["muon2_pt"]->Fill( it_mu2->Pt(), event_weight );
				histosTH1F["muon2_eta"]->Fill( it_mu2->Eta(), event_weight );
				histosTH1F["muon2_rapidity"]->Fill( it_mu2->Rapidity(), event_weight );
				histosTH1F["muon2_phi"]->Fill(it_mu2->Phi(), event_weight );
				phimu2 = it_mu2->Phi(); ptmu2 = it_mu2->Pt();
				etamu2 = it_mu2->Eta(); ymu2 = it_mu2->Rapidity();
				chmu2 = it_mu2->charge;

				//...
				TLorentzVector& muon1_lorentz = *it_mu1;
				TLorentzVector& muon2_lorentz = *it_mu2;
				TLorentzVector dimuon_lorentz(0.,0.,0.,0.);
				dimuon_lorentz += muon1_lorentz;
				dimuon_lorentz += muon2_lorentz;
				//histosTH1F["dimuon_multiplicity"]->Fill(n_dimuons_selected, event_weight );
				histosTH1F["dimuon_mass"]->Fill( dimuon_lorentz.M(), event_weight );
				histosTH1F["dimuon_pt"]->Fill( dimuon_lorentz.Pt(), event_weight );
				histosTH1F["dimuon_pt2"]->Fill( dimuon_lorentz.Pt()*dimuon_lorentz.Pt(), event_weight );
				histosTH1F["dimuon_eta"]->Fill( dimuon_lorentz.Eta(), event_weight );
				histosTH1F["dimuon_rapidity"]->Fill( dimuon_lorentz.Rapidity(), event_weight );
				dimuon_mass = dimuon_lorentz.M();  
				dimuon_pt=dimuon_lorentz.Pt(); 
				dimuon_pt2=dimuon_lorentz.Pt()*dimuon_lorentz.Pt(); 
				dimuon_eta =dimuon_lorentz.Eta();
				dimuon_rapidity = dimuon_lorentz.Rapidity();
				dphijpsi = dimuon_lorentz.Phi();
				//cout<<"Delta Definitions :"<<endl;		
				deltapt = fabs(muon1_lorentz.Pt() - muon2_lorentz.Pt());
				deltaeta = fabs(muon1_lorentz.Eta() - muon2_lorentz.Eta());
				deltaphi = fabs(phimu1 - phimu2);
				if(deltaphi > M_PI)deltaphi = (2*M_PI - deltaphi);
				deltay = fabs(muon1_lorentz.Rapidity() - muon2_lorentz.Rapidity());
				//Dphi = std::fabs(deltaPhi(phimu1 ,phimu2));  
				//cout<<"Delta Phi :"<<deltaphi<<endl;
				histosTH1F["muonDeltaPt"]->Fill(deltapt, event_weight );
				histosTH1F["muonDeltaEta"]->Fill(deltaeta, event_weight );	
				histosTH1F["muonDeltaPhi"]->Fill(deltaphi, event_weight );
				histosTH1F["muonDeltaY"]->Fill(deltay, event_weight ); 
				//histosTH2F["DeltaPhi_vs_dimuon_pt"]->Fill(dimuon_pt,deltaphi, event_weight ); 
				//histosTH1F["muonDphi"]->Fill(Dphi, event_weight );
				//histosTH1F["lorentzdphi"]->Fill(dphijpsi, event_weight );

			}
		}
		if(!mu2_selected)continue;if(verbose)cout<<mu2_selected<<endl;

		if(((dimuon_mass > 3.05) && (dimuon_mass < 3.15))){
			++n_jpsi_selected;
			jpsi_mass=true;
			//double Dphi_jpsi = std::fabs(deltaPhi(phimu1 ,phimu2));

			histosTH1F["jpsi_mass"]->Fill( dimuon_mass, event_weight );
			histosTH1F["jpsi_pt"]->Fill( dimuon_pt, event_weight );
			histosTH1F["jpsi_pt2"]->Fill( (dimuon_pt2), event_weight );
			histosTH1F["jpsi_eta"]->Fill( dimuon_eta, event_weight );
			histosTH1F["jpsi_rapidity"]->Fill( dimuon_rapidity, event_weight );
			histosTH1F["muonDeltaPt_jpsi"]->Fill(deltapt, event_weight );
			histosTH1F["muonDeltaEta_jpsi"]->Fill(deltaeta, event_weight );	
			histosTH1F["muonDeltaPhi_jpsi"]->Fill(deltaphi, event_weight );
			//histosTH1F["muonDphi_jpsi"]->Fill(Dphi, event_weight );
			//histosTH1F["lorentzdphi_jpsi"]->Fill(dphijpsi, event_weight );
			histosTH1F["muonDeltaY_jpsi"]->Fill(deltay, event_weight ); 
			//histosTH2F["jpsi_DeltaPhi_vs_dimuon_pt"]->Fill(dimuon_pt,deltaphi, event_weight );
			/*jpsi_mass = dimuon_lorentz.M();  
			  jpsi_pt=dimuon_lorentz.Pt(); 
			  jpsi_pt2=dimuon_lorentz.Pt()*dimuon_lorentz.Pt(); 
			  jpsi_eta =dimuon_lorentz.Eta();
			  jpsi_rapidity = dimuon_lorentz.Rapidity();*/
			if(verbose)cout<<"finalizando jpsi mass cut"<<endl; 
		}   
		if(!jpsi_mass)continue;if(verbose)cout<<dimuon_mass<<endl;
		++n_jpsi_mass;
		//Dimuons*/
		//GenPart
		double genEPlusPz = 0;
		double genEMinusPz = 0;
		double cm = 8000;
		double proton_pi = 4000;
		double proton_pz_plus=-999;
		double proton_px_plus = -999.;
		double proton_py_plus = -999.;
		double proton_energy_plus = 0.;
		double proton_pz_minus=999;
		double proton_px_minus = 999.;
		double proton_py_minus = 999.;
		double proton_energy_minus = 0.;
		double px_gen;
		double py_gen;
		double pz_gen;
		double energy_gen;
		double proton_pf;

		std::vector<double> eta_gen_vec;

		for(vector<MyGenPart>::iterator it_genpart = genPart->begin(); it_genpart != genPart->end(); ++it_genpart){

			double eta_gen = it_genpart->Eta();
			int status = it_genpart->status;
			int id = it_genpart->pdgId;

			if (status != 1) continue; //if(verbose)cout<<"final state for the particles"<<endl;
			if (id != 2212) continue;
			//	 if (eta_gen<4.9 && eta_gen>-4.9){ xi_posit_gen = true;
			energy_gen = it_genpart->Energy();
			px_gen = it_genpart->Px();
			py_gen = it_genpart->Py();
			pz_gen = it_genpart->Pz();
			proton_pf = sqrt(px_gen*px_gen+py_gen*py_gen+pz_gen*pz_gen);  
			double pz_cut = 0.7*proton_pi;
			if (fabs(pz_gen) < pz_cut) continue;

			genEPlusPz += (energy_gen + pz_gen);
			genEMinusPz += (energy_gen - pz_gen);

			if (pz_gen > proton_pz_plus) {
				proton_pz_plus = pz_gen; proton_energy_plus = energy_gen;
				proton_px_plus = px_gen; proton_py_plus = py_gen;   
				if(verbose)cout<<"(pz_gen > proton_pz_plus)"<<endl;    
			}
			if (pz_gen < proton_pz_minus) {
				proton_pz_minus = pz_gen; proton_energy_minus = energy_gen;
				proton_px_minus = px_gen; proton_py_minus = py_gen;
				if(verbose)cout<<"(pz_gen < proton_pz_minus)"<<endl;  
			}
			//eta_gen_vec.push_back( eta_gen);
		}


		double xi_plus_gen = genEPlusPz/cm;
		double xi_minus_gen = genEMinusPz/cm;
		double xi_proton_plus = -1.;
		double xi_proton_minus = -1.;
		double t_proton_plus = 0.;
		double t_proton_minus = 0.;
		double thx_proton_plus = 0.;
		double thy_proton_plus = 0.;
		double thx_proton_minus = 0.;
		double thy_proton_minus = 0.;

		++nevents_gen ;
		if(verbose)cout<<"(nevents_gen) ="<<nevents_gen<<endl;  
		bool proton_minus_rp_accept_120 = false;
		bool proton_minus_rp_accept_121 = false;
		bool proton_minus_rp_accept_122 = false;
		bool proton_minus_rp_accept_123 = false;
		bool proton_minus_rp_accept_124 = false;
		bool proton_minus_rp_accept_125 = false;
		bool proton_minus_rp_accept_020 = false;

		bool proton_plus_rp_accept_020 = false;
		bool proton_plus_rp_accept_021 = false;
		bool proton_plus_rp_accept_022 = false;
		bool proton_plus_rp_accept_023 = false;
		bool proton_plus_rp_accept_024 = false;
		bool proton_plus_rp_accept_025 = false;
		bool proton_plus_rp_accept_120 = false;
		bool fiducial_plus_020024 = false;
		bool fiducial_plus_021025 = false;
		bool fiducial_minus_120124 = false;
		bool fiducial_minus_121125 = false;

		std::map<int,std::vector<double> > proton_plus_pars;
		std::map<int,std::vector<double> > proton_minus_pars;

		if(proton_pz_plus > 0.){
			xi_proton_plus =  ( 1 - (proton_pz_plus/proton_pi) );
			TLorentzVector vec_pi(0.,0.,proton_pi,proton_pi);
			TLorentzVector vec_pf(proton_px_plus,proton_py_plus,proton_pz_plus,proton_energy_plus);
			TLorentzVector vec_t = (vec_pf - vec_pi);
			t_proton_plus = vec_t.Mag2();
			thx_proton_plus = atan(-proton_px_plus/proton_pi);
			thy_proton_plus = atan(proton_py_plus/proton_pi);
			if(verbose)cout<<"(proton_pz_plus > 0.)"<<endl;  
			//FIXME
			double out_x, out_thx, out_y, out_thy, out_xi;
			proton_plus_rp_accept_020 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 20, out_x, out_thx, out_y, out_thy, out_xi);
			proton_plus_pars[20] = std::vector<double>(5,0.);
			proton_plus_pars[20][0] = out_x; proton_plus_pars[20][1] = out_y;
			proton_plus_pars[20][2] = out_thx; proton_plus_pars[20][3] = out_thy;
			proton_plus_pars[20][4] = out_xi;
			histosTH2F["pos_y_vs_x_proton_plus_020"]->Fill( proton_plus_pars[20][0], proton_plus_pars[20][1] , event_weight );
			histosTH2F["pos_thy_vs_thx_proton_plus_020"]->Fill( proton_plus_pars[20][2], proton_plus_pars[20][3] , event_weight );

			proton_plus_rp_accept_021 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 21, out_x, out_thx, out_y, out_thy, out_xi);
			proton_plus_pars[21] = std::vector<double>(5,0.);
			proton_plus_pars[21][0] = out_x; proton_plus_pars[21][1] = out_y;
			proton_plus_pars[21][2] = out_thx; proton_plus_pars[21][3] = out_thy;
			proton_plus_pars[21][4] = out_xi;
			//cout<< "rp_021_x ; rp_021_y:"<<proton_plus_pars[21][0]<<" ; "<< proton_plus_pars[21][1]<<endl;
			histosTH2F["pos_y_vs_x_proton_plus_021"]->Fill( proton_plus_pars[21][0], proton_plus_pars[21][1] , event_weight );
			histosTH2F["pos_thy_vs_thx_proton_plus_021"]->Fill( proton_plus_pars[21][2], proton_plus_pars[21][3] , event_weight );
			//proton_plus_rp_accept_022 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 22);
			proton_plus_rp_accept_022 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 22, out_x, out_thx, out_y, out_thy, out_xi);
			proton_plus_pars[22] = std::vector<double>(5,0.);
			proton_plus_pars[22][0] = out_x; proton_plus_pars[22][1] = out_y;
			proton_plus_pars[22][2] = out_thx; proton_plus_pars[22][3] = out_thy;
			proton_plus_pars[22][4] = out_xi;
			histosTH2F["pos_y_vs_x_proton_plus_022"]->Fill( proton_plus_pars[22][0], proton_plus_pars[22][1] , event_weight );
			histosTH2F["pos_thy_vs_thx_proton_plus_022"]->Fill( proton_plus_pars[22][2], proton_plus_pars[22][3] , event_weight );

			proton_plus_rp_accept_023 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 23, out_x, out_thx, out_y, out_thy, out_xi);
			proton_plus_pars[23] = std::vector<double>(5,0.);
			proton_plus_pars[23][0] = out_x; proton_plus_pars[23][1] = out_y;
			proton_plus_pars[23][2] = out_thx; proton_plus_pars[23][3] = out_thy;
			proton_plus_pars[23][4] = out_xi;
			histosTH2F["pos_y_vs_x_proton_plus_023"]->Fill( proton_plus_pars[23][0], proton_plus_pars[23][1] , event_weight );
			histosTH2F["pos_thy_vs_thx_proton_plus_023"]->Fill( proton_plus_pars[23][2], proton_plus_pars[23][3] , event_weight );

			proton_plus_rp_accept_024 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 24, out_x, out_thx, out_y, out_thy, out_xi);
			proton_plus_pars[24] = std::vector<double>(5,0.);
			proton_plus_pars[24][0] = out_x; proton_plus_pars[24][1] = out_y;
			proton_plus_pars[24][2] = out_thx; proton_plus_pars[24][3] = out_thy;
			proton_plus_pars[24][4] = out_xi;
			//cout<< "rp_024_x ; rp_024_y:"<<proton_plus_pars[24][0]<<" ; "<< proton_plus_pars[24][1]<<endl;
			histosTH2F["pos_y_vs_x_proton_plus_024"]->Fill( proton_plus_pars[24][0], proton_plus_pars[24][1] , event_weight );
			histosTH2F["pos_thy_vs_thx_proton_plus_024"]->Fill( proton_plus_pars[24][2], proton_plus_pars[24][3] , event_weight );

			proton_plus_rp_accept_025 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 25, out_x, out_thx, out_y, out_thy, out_xi);
			proton_plus_pars[25] = std::vector<double>(5,0.);
			proton_plus_pars[25][0] = out_x; proton_plus_pars[25][1] = out_y;
			proton_plus_pars[25][2] = out_thx; proton_plus_pars[25][3] = out_thy;
			proton_plus_pars[25][4] = out_xi;
			//cout<< "rp_025_x ; rp_025_y:"<<proton_plus_pars[25][0]<<" ; "<< proton_plus_pars[25][1]<<endl;
			histosTH2F["pos_y_vs_x_proton_plus_025"]->Fill( proton_plus_pars[25][0], proton_plus_pars[25][1] , event_weight );
			histosTH2F["pos_thy_vs_thx_proton_plus_025"]->Fill( proton_plus_pars[25][2], proton_plus_pars[25][3] , event_weight );

			proton_plus_rp_accept_120 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 120);
			if(verbose)cout<<" RP020-025"<<endl;

		}


		if(proton_pz_minus < 0.){
			xi_proton_minus = (proton_pz_minus < 0.) ? ( 1 + (proton_pz_minus/proton_pi) ) : -1.;
			TLorentzVector vec_pi(0.,0.,-proton_pi,proton_pi);
			TLorentzVector vec_pf(proton_px_minus,proton_py_minus,proton_pz_minus,proton_energy_minus);
			TLorentzVector vec_t = (vec_pf - vec_pi);
			t_proton_minus = vec_t.Mag2();
			thx_proton_minus = atan(-proton_px_minus/proton_pi);
			thy_proton_minus = atan(proton_py_minus/proton_pi);
			if(verbose)cout<<"(proton_pz_minus > 0.)"<<endl;  
			double out_x, out_thx, out_y, out_thy, out_xi;
			proton_minus_rp_accept_120 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 120, out_x, out_thx, out_y, out_thy, out_xi);
			proton_minus_pars[120] = std::vector<double>(5,0.);
			proton_minus_pars[120][0] = out_x; proton_minus_pars[120][1] = out_y;
			proton_minus_pars[120][2] = out_thx; proton_minus_pars[120][3] = out_thy;
			proton_minus_pars[120][4] = out_xi;


			histosTH2F["pos_y_vs_x_proton_minus_120"]->Fill( proton_minus_pars[120][0], proton_minus_pars[120][1], event_weight );
			cout<< "rp_120_x ; rp_120_y:"<<proton_minus_pars[120][0]*1000.0<<" ; "<< proton_minus_pars[120][1]*1000.0<<endl;
			cout<< "rp_120_thx ; rp_120_thy:"<<proton_minus_pars[120][2]<<" ; "<< proton_minus_pars[120][3]<<endl;
			histosTH2F["pos_thy_vs_thx_proton_minus_120"]->Fill( proton_minus_pars[120][2], proton_minus_pars[120][3] , event_weight );
			proton_minus_rp_accept_121 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 121, out_x, out_thx, out_y, out_thy, out_xi);
			//if(Vertex)proton_minus_rp_accept_121 = protonRPDetected(pxpos_random, thx_proton_minus, pypos_random, thy_proton_minus, -xi_proton_minus, 121, out_x, out_thx, out_y, out_thy, out_xi);
			proton_minus_pars[121] = std::vector<double>(5,0.);
			proton_minus_pars[121][0] = out_x; proton_minus_pars[121][1] = out_y;
			proton_minus_pars[121][2] = out_thx; proton_minus_pars[121][3] = out_thy;
			proton_minus_pars[121][4] = out_xi;
			histosTH2F["pos_y_vs_x_proton_minus_121"]->Fill( proton_minus_pars[121][0], proton_minus_pars[121][1] , event_weight );
			histosTH2F["pos_thy_vs_thx_proton_minus_121"]->Fill( proton_minus_pars[121][2], proton_minus_pars[121][3] , event_weight );

			proton_minus_rp_accept_122 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 122, out_x, out_thx, out_y, out_thy, out_xi);
			//if(Vertex)proton_minus_rp_accept_122 = protonRPDetected(pxpos_random, thx_proton_minus, pypos_random, thy_proton_minus, -xi_proton_minus, 122, out_x, out_thx, out_y, out_thy, out_xi);
			proton_minus_pars[122] = std::vector<double>(5,0.);
			proton_minus_pars[122][0] = out_x; proton_minus_pars[122][1] = out_y;
			proton_minus_pars[122][2] = out_thx; proton_minus_pars[122][3] = out_thy;
			proton_minus_pars[122][4] = out_xi;
			histosTH2F["pos_y_vs_x_proton_minus_122"]->Fill( proton_minus_pars[122][0], proton_minus_pars[122][1] , event_weight );
			histosTH2F["pos_thy_vs_thx_proton_minus_122"]->Fill( proton_minus_pars[122][2], proton_minus_pars[122][3] , event_weight );

			proton_minus_rp_accept_123 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 123, out_x, out_thx, out_y, out_thy, out_xi);
			//if(Vertex)proton_minus_rp_accept_123 = protonRPDetected(pxpos_random, thx_proton_minus, pypos_random, thy_proton_minus, -xi_proton_minus, 123, out_x, out_thx, out_y, out_thy, out_xi);
			proton_minus_pars[123] = std::vector<double>(5,0.);
			proton_minus_pars[123][0] = out_x; proton_minus_pars[123][1] = out_y;
			proton_minus_pars[123][2] = out_thx; proton_minus_pars[123][3] = out_thy;
			proton_minus_pars[123][4] = out_xi;
			histosTH2F["pos_y_vs_x_proton_minus_123"]->Fill( proton_minus_pars[123][0], proton_minus_pars[123][1] , event_weight );
			histosTH2F["pos_thy_vs_thx_proton_minus_123"]->Fill( proton_minus_pars[123][2], proton_minus_pars[123][3] , event_weight );

			proton_minus_rp_accept_124 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 124, out_x, out_thx, out_y, out_thy, out_xi);
			//if(Vertex)proton_minus_rp_accept_124 = protonRPDetected(pxpos_random, thx_proton_minus, pypos_random, thy_proton_minus, -xi_proton_minus, 124, out_x, out_thx, out_y, out_thy, out_xi);
			proton_minus_pars[124] = std::vector<double>(5,0.);
			proton_minus_pars[124][0] = out_x; proton_minus_pars[124][1] = out_y;
			proton_minus_pars[124][2] = out_thx; proton_minus_pars[124][3] = out_thy;
			proton_minus_pars[124][4] = out_xi;
			histosTH2F["pos_y_vs_x_proton_minus_124"]->Fill( proton_minus_pars[124][0], proton_minus_pars[124][1] , event_weight );

			histosTH2F["pos_thy_vs_thx_proton_minus_124"]->Fill( proton_minus_pars[124][2], proton_minus_pars[124][3] , event_weight );


			proton_minus_rp_accept_125 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 125, out_x, out_thx, out_y, out_thy, out_xi);
			//if(Vertex)proton_minus_rp_accept_125 = protonRPDetected(pxpos_random, thx_proton_minus, pypos_random, thy_proton_minus, -xi_proton_minus, 125, out_x, out_thx, out_y, out_thy, out_xi);
			proton_minus_pars[125] = std::vector<double>(5,0.);
			proton_minus_pars[125][0] = out_x; proton_minus_pars[125][1] = out_y;
			proton_minus_pars[125][2] = out_thx; proton_minus_pars[125][3] = out_thy;
			proton_minus_pars[125][4] = out_xi;
			histosTH2F["pos_y_vs_x_proton_minus_125"]->Fill( proton_minus_pars[125][0], proton_minus_pars[125][1] , event_weight );  

			histosTH2F["pos_thy_vs_thx_proton_minus_125"]->Fill( proton_minus_pars[125][2], proton_minus_pars[125][3] , event_weight );



			proton_minus_rp_accept_020 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 20);

			if(verbose)cout<<" RP120-125"<<endl;

		}

		if(verbose)cout<<" starting RP combinations"<<endl;

		if( proton_plus_rp_accept_020 && proton_plus_rp_accept_024){
			cout<< "rp_020_x ; rp_020_y:"<<proton_plus_pars[20][0]*1000.0<<" ; "<< proton_plus_pars[20][1]*1000.0<<endl;
			cout<< "rp_024_x ; rp_024_y:"<<proton_plus_pars[24][0]*1000.0<<" ; "<< proton_plus_pars[24][1]*1000.0<<endl;
			fiducial_plus_020024 = ((proton_plus_pars[20][0] > -1.5*1000.0 && proton_plus_pars[20][0]*1000.0 < 6.5) &&
					(proton_plus_pars[20][1]*1000.0 > 7.0 && proton_plus_pars[20][1]*1000.0 < 29.0) &&
					(proton_plus_pars[24][0]*1000.0 > -1.5 && proton_plus_pars[24][0]*1000.0 < 6.5) &&
					(proton_plus_pars[24][1]*1000.0 > 7.0 && proton_plus_pars[24][1]*1000.0 < 29.0));}
		else if(( proton_plus_rp_accept_021 && proton_plus_rp_accept_025)){
			//cout<< "rp_021_x ; rp_021_y:"<<proton_plus_pars[21][0]*1000.0<<" ; "<< proton_plus_pars[21][1]*1000.0<<endl;
			//cout<< "rp_025_x ; rp_025_y:"<<proton_plus_pars[25][0]*1000.0<<" ; "<< proton_plus_pars[25][1]*1000.0<<endl;
			fiducial_plus_021025 = ((proton_plus_pars[21][0]*1000.0 > -1.5 && proton_plus_pars[21][0]*1000.0 < 6.5)&&
					(fabs(proton_plus_pars[21][1]*1000.0) > 7.0 && fabs(proton_plus_pars[21][1]*1000.0) < 29.0)&&
					(proton_plus_pars[25][0]*1000.0 > -1.5 && proton_plus_pars[25][0]*1000.0 < 6.5)&&
					(fabs(proton_plus_pars[25][1]*1000.0) > 7.0 && fabs(proton_plus_pars[25][1]*1000.0) < 29.0));}

		//RP Plus
		/*bool proton_plus_rp_accept = (proton_pz_plus > 0.) && 
		  ( ( proton_plus_rp_accept_020 && proton_plus_rp_accept_024 && fiducial_plus_020024 )|| 
		  ( proton_plus_rp_accept_021 && proton_plus_rp_accept_025 && fiducial_plus_021025));*/

		if( proton_minus_rp_accept_120 && proton_minus_rp_accept_124){
			fiducial_minus_120124 = ((proton_minus_pars[120][0]*1000.0 > -1.5 && proton_minus_pars[120][0]*1000.0 < 6.5) &&
					(proton_minus_pars[120][1]*1000.0 > 7.0 && proton_minus_pars[120][1]*1000.0 < 29.0) &&
					(proton_minus_pars[124][0]*1000.0 > -1.5 && proton_minus_pars[124][0]*1000.0 < 6.5) &&
					(proton_minus_pars[124][1]*1000.0 > 7.0 && proton_minus_pars[124][1]*1000.0 < 29.0));}
		else if(( proton_minus_rp_accept_121 && proton_minus_rp_accept_125)){
			fiducial_minus_121125 = ((proton_minus_pars[121][0]*1000.0 > -1.5 && proton_minus_pars[121][0]*1000.0 < 6.5)&&
					(fabs(proton_minus_pars[121][1]*1000.0) > 7.0 && fabs(proton_minus_pars[121][1]*1000.0) < 29.0)&&
					(proton_minus_pars[125][0]*1000.0 > -1.5 && proton_minus_pars[125][0]*1000.0 < 6.5)&&
					(fabs(proton_minus_pars[125][1]*1000.0) > 7.0 && fabs(proton_minus_pars[125][1]*1000.0) < 29.0));}




		// Particle-flow
		double soma1 = 0;
		double soma2 = 0;
		double eta_max=-999.;
		double eta_min=999.;
		/*      double cm = 8000;*/
		double sumEHFMinus = 0.;
		double sumEHFPlus = 0.;

		for(vector<MyPFCand>::iterator it_pfcand = pFlow_coll->begin(); it_pfcand != pFlow_coll->end(); ++it_pfcand){
			int partType = it_pfcand->particleId;
			double eta = it_pfcand->Eta();
			double energy = it_pfcand->Energy();
			double pz = it_pfcand->Pz();

			// HF eta rings 29, 30, 40, 41
			if( ( (fabs(eta) >= 2.866) && (fabs(eta) < 3.152) ) || (fabs(eta) >= 4.730) ) continue;



			////////
			soma1 += (energy + pz);
			soma2 += (energy - pz);

			if (eta > eta_max) {eta_max = eta; PF_eta_max = true;} 
			if (eta < eta_min) {eta_min = eta; PF_eta_min = true;}


			histosTH2F["energyVsEtaAllTypes"]->Fill( eta, energy, event_weight );

			if(partType == MyPFCand::X)
				histosTH2F["energyVsEtaUndefined"]->Fill( eta, energy, event_weight );
			else if(partType == MyPFCand::h)
				histosTH2F["energyVsEtaChargedHadron"]->Fill( eta, energy, event_weight ); 
			else if(partType == MyPFCand::e) 
				histosTH2F["energyVsEtaElectron"]->Fill( eta, energy, event_weight );
			else if(partType == MyPFCand::mu) 
				histosTH2F["energyVsEtaMuon"]->Fill( eta, energy, event_weight );
			else if(partType == MyPFCand::gamma) 
				histosTH2F["energyVsEtaGamma"]->Fill( eta, energy, event_weight );
			else if(partType == MyPFCand::h0) 
				histosTH2F["energyVsEtaNeutralHadron"]->Fill( eta, energy, event_weight );
			else if(partType == MyPFCand::h_HF){ 
				histosTH2F["energyVsEtaHadronHF"]->Fill( eta, energy, event_weight );}
			else if(partType == MyPFCand::egamma_HF) 
				histosTH2F["energyVsEtaEGammaHF"]->Fill( eta, energy, event_weight );
			if((3.0 < eta) && (eta < 5.1) ){
				sumEHFPlus += energy;


			}
			if((-5.1 < eta) && (eta < -3.0) ){
				sumEHFMinus += energy;

			}


		}
		if(!PF_eta_max) continue;
		if(!PF_eta_min) continue;
		weight_total_PF_selected += event_weight;

		++nevents_pf;

		double xi_plus_Reco = soma1/cm;
		double xi_minus_Reco = soma2/cm;
		double delta_eta_maxmin = eta_max - eta_min;  

		double correction = xi_plus_Reco/xi_plus_gen;
		double resolution_before = (xi_plus_gen-xi_plus_Reco)/xi_plus_gen;
		double xi_reconst = xi_plus_Reco/0.8;
		double resolution_after = (xi_plus_gen-xi_reconst)/xi_plus_gen;

		if(((dimuon_mass > 3.05) && (dimuon_mass < 3.15))){
			histosTH1F["jpsi_xi_cms_pfplus"]->Fill(xi_plus_Reco,event_weight);
			histosTH1F["jpsi_xi_cms_pfminus"]->Fill(xi_minus_Reco,event_weight);
			histosTH1F["jpsi_Eta_max"]->Fill( eta_max, event_weight  );
			histosTH1F["jpsi_Eta_min"]->Fill( eta_min, event_weight  );
			histosTH1F["jpsi_Delta_eta_maxmin"]->Fill( delta_eta_maxmin, event_weight  );
			histosTH1F["jpsi_cms_sumEHFminus"]->Fill(sumEHFMinus,event_weight);
			histosTH1F["jpsi_cms_sumEHFplus"]->Fill(sumEHFPlus,event_weight);
		}
		if(verbose)cout<<" Fim do PF"<<endl;
		// smearing /////////////////////////////////////////////////////////////
		//totem proton reconstructed
		float sigma_xi45=0.00714986 - 0.0408903*xi_proton_plus + 0.0965813*xi_proton_plus*xi_proton_plus; // sigma45 vs xi from Hubert
		float sigma_xi56=0.00720615 - 0.0418783*xi_proton_minus + 0.0999515*xi_proton_minus*xi_proton_minus; // sigma56 vs xi from Hubert
		float xi_proton_plus_rec = xi_proton_plus + gRandom->Gaus(0,sigma_xi45);
		float xi_proton_minus_rec = xi_proton_minus + gRandom->Gaus(0,sigma_xi56);

		double sigma_t45=0.233365*t_proton_plus - 0.0975751*t_proton_plus*t_proton_plus; // sigma_t45 vs t from Hubert
		double sigma_t56=0.233365*t_proton_minus - 0.0975751*t_proton_minus*t_proton_minus; // sigma_t56 vs t from Hubert
		double t_proton_plus_rec = t_proton_plus + gRandom->Gaus(0,sigma_t45);
		double t_proton_minus_rec = t_proton_minus + gRandom->Gaus(0,sigma_t56);

		histosTH1F["t_proton_minus_rec"]->Fill( fabs(t_proton_minus_rec) , event_weight );
		histosTH1F["xi_proton_minus_rec"]->Fill( xi_proton_minus_rec , event_weight );
		histosTH1F["t_proton_plus_rec"]->Fill( fabs(t_proton_plus_rec) , event_weight );
		histosTH1F["xi_proton_plus_rec"]->Fill( xi_proton_plus_rec , event_weight );
		if(verbose)cout<<" Fim do smearing"<<endl;
		//rp_accept
		bool proton_minus_rp_accept_mc = ( ( proton_minus_rp_accept_120 && proton_minus_rp_accept_124 && fiducial_minus_120124 )||
				( proton_minus_rp_accept_121 && proton_minus_rp_accept_125 && fiducial_minus_121125 ));

		bool proton_plus_rp_accept_mc =  ( ( proton_plus_rp_accept_020 && proton_plus_rp_accept_024 && fiducial_plus_020024 )|| 
				( proton_plus_rp_accept_021 && proton_plus_rp_accept_025 && fiducial_plus_021025));
		bool proton_minus_kinec_accept_t_mc = fabs(t_proton_minus_rec)>=0.03 && fabs(t_proton_minus_rec)<=1.0;
		bool proton_minus_kinec_accept_xi_mc = xi_proton_minus_rec>-0.04 && xi_proton_minus_rec<0.23;
		bool xi_cms_totem_cut_minus =  xi_minus_Reco-xi_proton_minus_rec<0.009;
		bool xi_cms_totem_cut_minus_halo =  xi_minus_Reco-xi_proton_minus_rec>0.009;
		bool proton_plus_kinec_accept_t_mc = fabs(t_proton_plus_rec)>=0.03 && fabs(t_proton_plus_rec)<=1.0;
		bool proton_plus_kinec_accept_xi_mc = xi_proton_plus_rec>-0.04 && xi_proton_plus_rec<0.23;
		bool xi_cms_totem_cut_plus =  xi_plus_Reco-xi_proton_plus_rec<0.009;
		bool xi_cms_totem_cut_plus_halo =  xi_plus_Reco-xi_proton_plus_rec>0.009;
		if(verbose)cout<<" Fim das combinacoes rp accept"<<endl;

		//Access Zero-bias data
		int i_evt_ZB = 0 + gRandom->Rndm()*(zeroBias.size());
		MyZeroBiasData const & zeroBiasData = zeroBias.at(i_evt_ZB);
		double xi_proton_right_zb = zeroBiasData.proton_rec_right_xi;      
		double t_proton_right_zb = zeroBiasData.proton_rec_right_t;
		bool valid_proton_right_zb = zeroBiasData.proton_rec_right_valid;      
		bool valid_proton_left_zb = zeroBiasData.proton_rec_left_valid;      
		double xi_proton_left_zb = zeroBiasData.proton_rec_left_xi;      
		double t_proton_left_zb = zeroBiasData.proton_rec_left_t;
		bool valid_vtx = zeroBiasData.vtx_valid;
		double ndof_vtx = zeroBiasData.vtx_ndof;
		if (valid_vtx && ndof_vtx>=4) continue; 
		++nZBconditions;
		if(verbose)cout<<"(valid_vtx && ndof_vtx>=4)"<<" "<<"valid_vtx"<<" "<<ndof_vtx<<endl;
		if(verbose)cout<<" Fim do zb condicoes"<<endl;
		//double proton_minus_beta_zb = x_minus/-xi_proton_right_zb;
		//double proton_plus_beta_zb = x_plus/-xi_proton_left_zb;
		//(-0.23 < xi_proton_right)&&( xi_proton_right < 0.04)
		bool proton_right_kinec_xi_zb = xi_proton_right_zb<0.04 && xi_proton_right_zb>-0.23;
		bool proton_right_kinec_t_zb = -t_proton_right_zb>=0.03 && -t_proton_right_zb<=1.0; 
		bool xi_cms_totem_cut_right_zb = xi_minus_Reco+xi_proton_right_zb<0.009;   
		bool xi_cms_totem_cut_right_zb_halo = xi_minus_Reco+xi_proton_right_zb>0.009;    
		bool proton_left_kinec_xi_zb = xi_proton_left_zb<0.04 && xi_proton_left_zb>-0.23;
		bool proton_left_kinec_t_zb = -t_proton_left_zb>=0.03 && -t_proton_left_zb<=1.0; 
		bool xi_cms_totem_cut_left_zb = xi_plus_Reco+xi_proton_left_zb<0.009; 
		bool xi_cms_totem_cut_left_zb_halo = xi_plus_Reco+xi_proton_left_zb>0.009;  
		histosTH1F["xi_right_zb"]->Fill( -xi_proton_right_zb , event_weight );
		//decision where the proton in Totem acceptance came from
		bool signal_right = !valid_proton_right_zb && !valid_proton_left_zb && !proton_plus_rp_accept_mc && proton_minus_rp_accept_mc;
		bool signal_left = !valid_proton_right_zb && !valid_proton_left_zb && !proton_minus_rp_accept_mc && proton_plus_rp_accept_mc;
		bool backg_right = !proton_plus_rp_accept_mc && !proton_minus_rp_accept_mc && !valid_proton_left_zb && valid_proton_right_zb;
		bool backg_left = !proton_plus_rp_accept_mc && !proton_minus_rp_accept_mc && !valid_proton_right_zb && valid_proton_left_zb;

		///proton from MC in the RP acceptance
		if (signal_right){
			++nMCRight;
			histosTH1F["t_proton_signal_right"]->Fill( fabs(t_proton_minus_rec) , event_weight );
			//histosTH1F["beta_proton_minus_signal"]->Fill( proton_minus_beta , event_weight );
			histosTH1F["xi_proton_signal_right"]->Fill( xi_proton_minus_rec , event_weight );
			histosTH1F["xi_cms_totem_right_signal"]->Fill( xi_minus_Reco-xi_proton_minus_rec , event_weight );
			//histosTH1F["proton_right_dimuon_mass"]->Fill( dimuon_mass, event_weight );
			///kinematic region
			if (proton_minus_kinec_accept_t_mc) histosTH1F["xi_proton_signal_right_kint"]->Fill( xi_proton_minus_rec , event_weight );
			if (proton_minus_kinec_accept_t_mc && proton_minus_kinec_accept_xi_mc) histosTH1F["xi_cms_totem_right_signal_kin"]->Fill( xi_minus_Reco-xi_proton_minus_rec , event_weight );
			///xi_cms -xi_totem cut 
			if (proton_minus_kinec_accept_t_mc && xi_cms_totem_cut_minus) histosTH1F["xi_proton_signal_right_kint_cut"]->Fill( xi_proton_minus_rec , event_weight );
			if (proton_minus_kinec_accept_t_mc && proton_minus_kinec_accept_xi_mc && xi_cms_totem_cut_minus){
				++nSignalRegionMCRight;
				histosTH1F["xi_proton_signal_right_kin_cut"]->Fill( xi_proton_minus_rec , event_weight );
				histosTH1F["t_proton_signal_right_kin_cut"]->Fill( fabs(t_proton_minus_rec) , event_weight );
				//histosTH1F["beta_proton_minus_signal_kin_cut"]->Fill( proton_minus_beta , event_weight );
			}
			if (proton_minus_kinec_accept_t_mc && proton_minus_kinec_accept_xi_mc && xi_cms_totem_cut_minus_halo){
				//cout<<"right mc:"<<xi_cms_totem_cut_minus_halo<<endl;
				++nBckgRegionMCRight;
				histosTH1F["xi_proton_signal_right_kin_cut_halo"]->Fill( xi_proton_minus_rec , event_weight );
				histosTH1F["t_proton_signal_right_kin_cut_halo"]->Fill( fabs(t_proton_minus_rec) , event_weight );

			}

		}
		if(verbose)cout<<"proton from MC in the RP acceptance"<<endl;

		if (signal_left){
			++nMCLeft;
			histosTH1F["t_proton_signal_left"]->Fill( fabs(t_proton_plus_rec) , event_weight );
			//histosTH1F["beta_proton_plus_signal"]->Fill( proton_plus_beta , event_weight );
			histosTH1F["xi_proton_signal_left"]->Fill( xi_proton_plus_rec , event_weight );
			histosTH1F["xi_cms_totem_left_signal"]->Fill( xi_plus_Reco-xi_proton_plus_rec , event_weight );
			//histosTH1F["proton_left_dimuon_mass"]->Fill( dimuon_mass, event_weight );cout<< dimuon_mass<<endl;
			///kinematic region
			if(proton_plus_kinec_accept_t_mc) histosTH1F["xi_proton_signal_left_kint"]->Fill( xi_proton_plus_rec , event_weight );
			if(proton_plus_kinec_accept_t_mc && proton_plus_kinec_accept_xi_mc) histosTH1F["xi_cms_totem_left_signal_kin"]->Fill( xi_plus_Reco-xi_proton_plus_rec , event_weight );
			//xi_cms -xi_totem cut 
			if(verbose)cout<<"proton from MC signal left"<<endl;
			if (proton_plus_kinec_accept_t_mc && proton_plus_kinec_accept_xi_mc && xi_cms_totem_cut_plus){ 
				++nSignalRegionMCLeft;
				histosTH1F["xi_proton_signal_left_kin_cut"]->Fill( xi_proton_plus_rec , event_weight );
				histosTH1F["t_proton_signal_left_kin_cut"]->Fill( fabs(t_proton_plus_rec) , event_weight );

				if(verbose)cout<<"proton from MC signal left kin"<<endl;
			}
			if (proton_plus_kinec_accept_t_mc && proton_plus_kinec_accept_xi_mc && xi_cms_totem_cut_plus_halo){ 
				++nBckgRegionMCLeft;
				//cout<<"left mc:"<<xi_cms_totem_cut_plus_halo<<endl;
				histosTH1F["xi_proton_signal_left_kin_cut_halo"]->Fill( xi_proton_plus_rec , event_weight );
				histosTH1F["t_proton_signal_left_kin_cut_halo"]->Fill( fabs(t_proton_plus_rec) , event_weight );

			}
			if(verbose)cout<<"left mc:"<<xi_plus_Reco-xi_proton_plus_rec<<endl;
		}

		///proton from ZB 
		if (backg_right){ 
			++nZBRight;
			histosTH1F["t_proton_backg_right"]->Fill( -t_proton_right_zb , event_weight );
			histosTH1F["xi_proton_backg_right"]->Fill( -xi_proton_right_zb , event_weight );
			//histosTH1F["beta_proton_minus_backg"]->Fill( proton_minus_beta_zb , event_weight );
			histosTH1F["xi_cms_totem_right_zb"]->Fill( xi_minus_Reco+xi_proton_right_zb , event_weight );
			///kinematic region
			if (proton_right_kinec_t_zb) histosTH1F["xi_proton_backg_right_kint"]->Fill( -xi_proton_right_zb , event_weight );
			if (proton_right_kinec_t_zb && proton_right_kinec_xi_zb) histosTH1F["xi_cms_totem_right_zb_kin"]->Fill( xi_minus_Reco+xi_proton_right_zb , event_weight );
			///xi_cms-xi_totem cut
			if (proton_right_kinec_t_zb && xi_cms_totem_cut_right_zb) histosTH1F["xi_proton_backg_right_kint_cut"]->Fill( -xi_proton_right_zb , event_weight );
			if (proton_right_kinec_t_zb && proton_right_kinec_xi_zb && xi_cms_totem_cut_right_zb){
				++nSignalRegionZBRight;
				histosTH1F["xi_proton_backg_right_kin_cut"]->Fill( -xi_proton_right_zb , event_weight );
				histosTH1F["t_proton_backg_right_kin_cut"]->Fill( -t_proton_right_zb , event_weight );
				//histosTH1F["beta_proton_minus_backg_kin_cut"]->Fill( proton_minus_beta_zb , event_weight );
			}
			if (proton_right_kinec_t_zb && proton_right_kinec_xi_zb && xi_cms_totem_cut_right_zb_halo){
				++nBckgRegionZBRight;      
				//cout<<"right zb:"<<xi_minus_Reco+xi_proton_right_zb<<endl;
				histosTH1F["xi_proton_backg_right_kin_cut_halo"]->Fill( -xi_proton_right_zb , event_weight );
				histosTH1F["t_proton_backg_right_kin_cut_halo"]->Fill( -t_proton_right_zb , event_weight );

			}
		}

		if (backg_left){
			++nZBLeft;
			histosTH1F["t_proton_backg_left"]->Fill( -t_proton_left_zb , event_weight );
			histosTH1F["xi_proton_backg_left"]->Fill( -xi_proton_left_zb , event_weight );
			//histosTH1F["beta_proton_plus_backg"]->Fill( proton_plus_beta_zb , event_weight );
			histosTH1F["xi_cms_totem_left_zb"]->Fill( xi_plus_Reco+xi_proton_left_zb , event_weight );
			///kinematic region
			if (proton_left_kinec_t_zb) histosTH1F["xi_proton_backg_left_kint"]->Fill( -xi_proton_left_zb , event_weight );
			if (proton_left_kinec_t_zb && proton_left_kinec_xi_zb) histosTH1F["xi_cms_totem_left_zb_kin"]->Fill( xi_plus_Reco+xi_proton_left_zb , event_weight );
			///xi_cms-xi_totem cut
			if (proton_left_kinec_t_zb && proton_left_kinec_xi_zb && xi_cms_totem_cut_left_zb){
				++nSignalRegionZBLeft;
				histosTH1F["xi_proton_backg_left_kin_cut"]->Fill( -xi_proton_left_zb , event_weight );
				histosTH1F["t_proton_backg_left_kin_cut"]->Fill( -t_proton_left_zb , event_weight );
				//histosTH1F["beta_proton_plus_backg_kin_cut"]->Fill( proton_plus_beta_zb , event_weight );
			}
			if (proton_left_kinec_t_zb && proton_left_kinec_xi_zb && xi_cms_totem_cut_left_zb_halo){
				++nBckgRegionZBLeft;
				//cout<<"left zb:"<<xi_plus_Reco+xi_proton_left_zb<<endl;
				histosTH1F["xi_proton_backg_left_kin_cut_halo"]->Fill( -xi_proton_left_zb , event_weight );
				histosTH1F["t_proton_backg_left_kin_cut_halo"]->Fill( -t_proton_left_zb , event_weight );

			}
		}

		///proton from MC and from ZB -> smallest xi
		if (proton_minus_rp_accept_mc && valid_proton_right_zb) {
			++nBothsRight;
			double t_proton_right_sel =  (xi_proton_minus_rec<-xi_proton_right_zb) ? fabs(t_proton_minus_rec) : -t_proton_right_zb;
			double xi_proton_right_sel =  (xi_proton_minus_rec<-xi_proton_right_zb) ? xi_proton_minus_rec : -xi_proton_right_zb;
			//double beta_proton_right_sel = x_minus/xi_proton_right_sel;
			histosTH1F["t_proton_right_both"]->Fill( t_proton_right_sel , event_weight );
			histosTH1F["xi_proton_right_both"]->Fill( xi_proton_right_sel , event_weight );
			histosTH1F["xi_cms_totem_right_both"]->Fill( xi_minus_Reco-xi_proton_right_sel , event_weight );
			//histosTH1F["beta_proton_both_right"]->Fill( beta_proton_right_sel , event_weight );
			//kinematic region
			if (t_proton_right_sel>=0.03 && t_proton_right_sel<=1.0) histosTH1F["xi_proton_both_right_kint"]->Fill( xi_proton_right_sel , event_weight );
			if (t_proton_right_sel>=0.03 && t_proton_right_sel<=1.0 && xi_proton_right_sel>-0.04 && xi_proton_right_sel<0.23) histosTH1F["xi_cms_totem_right_both_kin"]->Fill( xi_minus_Reco-xi_proton_right_sel , event_weight );
			//if(t_proton_right_sel>=0.03 && t_proton_right_sel<=1.0 ) histosTH1F["xi_cms_totem_right_both_kin"]->Fill( xi_minus_Reco-xi_proton_right_sel , event_weight );
			///xi_cms -xi_totem cut 
			if (t_proton_right_sel>=0.03 && t_proton_right_sel<=1.0 && xi_minus_Reco-xi_proton_right_sel<0.009) histosTH1F["xi_proton_both_right_kint_cut"]->Fill( xi_proton_right_sel , event_weight );
			if (t_proton_right_sel>=0.03 && t_proton_right_sel<=1.0 && xi_proton_right_sel>-0.04 && xi_proton_right_sel<0.23 && xi_minus_Reco-xi_proton_right_sel<0.009){
				++nSignalRegionBothsRight;
				histosTH1F["xi_proton_both_right_kin_cut"]->Fill( xi_proton_right_sel , event_weight );
				histosTH1F["t_proton_both_right_kin_cut"]->Fill( t_proton_right_sel , event_weight );

			}
			if (t_proton_right_sel>=0.03 && t_proton_right_sel<=1.0 && xi_proton_right_sel>-0.04 && xi_proton_right_sel<0.23 && xi_minus_Reco-xi_proton_right_sel>0.009){
				++nBckgRegionBothsRight;
				//cout<<"both R:"<<xi_minus_Reco-xi_proton_right_sel<<endl;
				histosTH1F["xi_proton_both_right_kin_cut_halo"]->Fill( xi_proton_right_sel , event_weight );
				histosTH1F["t_proton_both_right_kin_cut_halo"]->Fill( t_proton_right_sel , event_weight );

			}
			if(verbose)cout<<"proton Right from MC and from ZB -> smallest xi"<<endl;
		}


		if (proton_plus_rp_accept_mc && valid_proton_left_zb) {
			++nBothsLeft;
			double t_proton_left_sel =  (xi_proton_plus_rec<-xi_proton_left_zb) ? fabs(t_proton_plus_rec) : -t_proton_left_zb;
			double xi_proton_left_sel =  (xi_proton_plus_rec<-xi_proton_left_zb) ? xi_proton_plus_rec : -xi_proton_left_zb;
			//double beta_proton_left_sel = x_minus/xi_proton_left_sel;
			histosTH1F["t_proton_left_both"]->Fill( t_proton_left_sel , event_weight );
			histosTH1F["xi_proton_left_both"]->Fill( xi_proton_left_sel , event_weight );
			histosTH1F["xi_cms_totem_left_both"]->Fill( xi_plus_Reco-xi_proton_left_sel , event_weight );
			//histosTH1F["beta_proton_both_right"]->Fill( beta_proton_right_sel , event_weight );
			//kinematic region
			if (t_proton_left_sel>=0.03 && t_proton_left_sel<=1.0) histosTH1F["xi_proton_both_left_kint"]->Fill( xi_proton_left_sel , event_weight );
			//if (t_proton_left_sel>=0.03 && t_proton_left_sel<=1.0) histosTH1F["xi_cms_totem_left_both_kin"]->Fill( xi_plus_Reco-xi_proton_left_sel , event_weight );
			if (t_proton_left_sel>=0.03 && t_proton_left_sel<=1.0 && xi_proton_left_sel>-0.04 && xi_proton_left_sel<0.23) histosTH1F["xi_cms_totem_left_both_kin"]->Fill( xi_plus_Reco-xi_proton_left_sel , event_weight );
			///xi_cms -xi_totem cut 
			if (t_proton_left_sel>=0.03 && t_proton_left_sel<=1.0 && xi_plus_Reco-xi_proton_left_sel<0.009) histosTH1F["xi_proton_both_left_kint_cut"]->Fill( xi_proton_left_sel , event_weight );
			if (t_proton_left_sel>=0.03 && t_proton_left_sel<=1.0 && xi_proton_left_sel>-0.04 && xi_proton_left_sel<0.23 && xi_plus_Reco-xi_proton_left_sel<0.009){
				++nSignalRegionBothsLeft;
				//if(t_proton_left_sel>=0.03 && t_proton_left_sel<=1.0  && xi_plus_Reco-xi_proton_left_sel<=0.009){
				histosTH1F["xi_proton_both_left_kin_cut"]->Fill( xi_proton_left_sel , event_weight );
				histosTH1F["t_proton_both_left_kin_cut"]->Fill( t_proton_left_sel , event_weight );
				//histosTH1F["beta_proton_both_right_kin_cut"]->Fill( beta_proton_right_sel , event_weight );
			}

			if (t_proton_left_sel>=0.03 && t_proton_left_sel<=1.0 && xi_proton_left_sel>-0.04 && xi_proton_left_sel<0.23 && xi_plus_Reco-xi_proton_left_sel>0.009){
				++nBckgRegionBothsLeft;
				if(verbose)cout<<"both:"<<xi_plus_Reco-xi_proton_left_sel<<endl;
				histosTH1F["xi_proton_both_left_kin_cut_halo"]->Fill( xi_proton_left_sel , event_weight );
				histosTH1F["t_proton_both_left_kin_cut_halo"]->Fill( t_proton_left_sel , event_weight );

			}
			if(verbose)cout<<"proton Left from MC and from ZB -> smallest xi"<<endl;
			}
			histosTH2F["xi_plus_reco_gen"]->Fill( xi_plus_gen, xi_plus_Reco, event_weight );
			histosTH2F["xi_minus_reco_gen"]->Fill( xi_minus_gen, xi_minus_Reco, event_weight );
			histosTH2F["logxi_plus_reco_gen"]->Fill( log10(xi_plus_gen), log10(xi_plus_Reco), event_weight );
			histosTH2F["logxi_minus_reco_gen"]->Fill( log10(xi_minus_gen), log10(xi_minus_Reco), event_weight );
			histosTH1F["Eta_max"]->Fill( eta_max, event_weight  );
			histosTH1F["Eta_min"]->Fill( eta_min, event_weight  );
			histosTH1F["Delta_eta_maxmin"]->Fill( delta_eta_maxmin, event_weight  );
			histosTH1F["xi_plus_Reco"]->Fill( xi_plus_Reco, event_weight  );
			histosTH1F["xi_minus_Reco"]->Fill( xi_minus_Reco, event_weight  );
			histosTH1F["logxi_plus"]->Fill( log10(xi_plus_Reco), event_weight  );

			if(verbose)cout<<"variaveis do PF, deltaEta, xi_cms, ..."<<endl;

		}//end loop for events

		cout <<"After PF selection " << nevents_pf << " events "<< endl;
		cout <<"nBothsLeft = "<< nBothsLeft<<endl;
		cout <<"nBckgRegionBothsLeft= "<< nBckgRegionBothsLeft<<endl;
		cout <<"nSignalRegionBothsLeft= "<<nSignalRegionBothsLeft <<endl;

		cout <<"nBothsRight= "<< nBothsRight<<endl;
		cout <<"nBckgRegionBothsRight= "<< nBckgRegionBothsRight<<endl;
		cout <<"nSignalRegionBothsRight= "<<nSignalRegionBothsRight <<endl;

		cout <<"nBckgRegionZBLeft= "<<nBckgRegionZBLeft <<endl;
		cout <<"nSignalRegionZBLeft= "<< nSignalRegionZBLeft<<endl;
		cout <<"nZBLeft= "<<nZBLeft <<endl;

		cout <<"nBckgRegionZBRight= "<<nBckgRegionZBRight <<endl;
		cout <<"nSignalRegionZBRight= "<<nSignalRegionZBRight <<endl;
		cout <<"nZBRight= "<< nZBRight<<endl;

		cout <<"nBckgRegionMCLeft= "<<nBckgRegionMCLeft <<endl;
		cout <<"nSignalRegionMCLeft= "<<nSignalRegionMCLeft <<endl;
		cout <<"nMCLeft= "<< nMCLeft<<endl;

		cout <<"nBckgRegionMCRight= "<< nBckgRegionMCRight<<endl;
		cout <<"nSignalRegionMCRight= "<< nSignalRegionMCRight<<endl;
		cout <<"nMCRight= "<<nMCRight <<endl;

		cout <<"nZBconditions= "<< nZBconditions<<endl;
		cout <<"nJpsi = "<< n_jpsi_mass<<endl;

		file->Close();

		}//end of loop over files

		//output file
		TFile* output = new TFile(outputFileName.c_str(),"RECREATE");
		output->cd();



		histosTH2F["xi_plus_reco_gen"]->SetOption("colz");
		histosTH2F["xi_minus_reco_gen"]->SetOption("colz");
		histosTH2F["logxi_plus_reco_gen"]->SetOption("colz");
		histosTH2F["logxi_minus_reco_gen"]->SetOption("colz");

		for(map<string,TH1F*>::iterator it_histo = histosTH1F.begin();it_histo != histosTH1F.end(); ++it_histo){
			//(*it_histo).second->Scale(scale);
			(*it_histo).second->Write();
		}
		for(map<string,TH2F*>::iterator it_histo = histosTH2F.begin();
				it_histo != histosTH2F.end(); ++it_histo)
			(*it_histo).second->Write();

		output->Close();
	}
