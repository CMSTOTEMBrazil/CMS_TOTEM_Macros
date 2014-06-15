//STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
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
#include "MyFSCHit.h"
#include "MyFSCDigi.h"

// TOTEM data formats
#include "T1Event.h"
#include "T2Event.h"
#include "RPRootDumpReconstructedProton.h"
#include "RPRootDumpReconstructedProtonPair.h"
#include "RPRootDumpTrackInfo.h"
#include "RPRootDumpDigiInfo.h"
#include "RPRootDumpPatternInfo.h"

#include "analysis_tools.h"

//STANDARD C++ INCLUDES
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <fstream>
#define PI 3.141592653589793
#include "deltaPhi.h"
using namespace std;

void data_SDJpsi_CMSTOTEM_MB_Ds(string const& outputFileName = "/afs/cern.ch/work/e/eliza/private/TOTEM_JPSI/CMSTotem/Workspace/Jpsi_SingleArmRecProton/MB_output_semfidWN.root",const Double_t t_proton_down_=0.03, const Double_t t_proton_up_=1.0,const Int_t Bin_mass=200 ,const Int_t nevt_max = -1){
	bool verbose = false;
	string treeName = "cms_totem";
	double ptMax = 9999.0;
	double etaMaxThreshold = 2.0;

	bool selectBunchCrossing = false;
	bool selectVertex = true;
	bool selectTrack = true;
	bool selectMuons = true;
	bool selectEtaMax = false;
	bool selectEtaMin = false;
	bool selectZeroHitsT2Plus = false;
	bool selectZeroHitsT2Minus = false;
	bool selectSingleArmRecProton = true;
	bool selectDoubleArmRecProton = false;
	bool selectElastic = false;
	bool selectNonElastic = false;
	bool signal_left = true;
	bool signal_right = true;

	const Int_t nevt_max_corr = (nevt_max >= 0) ? nevt_max : 99999999;

	vector<string> selectHLTPathNames;
	selectHLTPathNames.push_back("HLT_L1Tech53_MB_2_v1");
	//selectHLTPathNames.push_back("HLT_L1_DoubleMu0");
	//selectHLTPathNames.push_back("HLT_ZeroBias_v7");
	vector<string> hltPathNames;
	hltPathNames.push_back("HLT_L1DoubleEG3_FwdVeto_v1");
	hltPathNames.push_back("HLT_L1DoubleMu0_v1");
	hltPathNames.push_back("HLT_L1DoubleJet20_RomanPotsOR_v1");
	hltPathNames.push_back("HLT_L1DoubleJet20part1_v1");
	hltPathNames.push_back("HLT_L1DoubleJet24_v1");
	hltPathNames.push_back("HLT_L1DoubleJet20part2_v1");
	hltPathNames.push_back("HLT_L1Tech40_BPTXAND_v1");
	hltPathNames.push_back("HLT_L1Tech53_MB_1_v1");
	hltPathNames.push_back("HLT_L1Tech_HF9OR10_v1");
	hltPathNames.push_back("HLT_T1minbias_Tech55_v1");
	hltPathNames.push_back("HLT_L1Tech53_MB_2_v1");
	hltPathNames.push_back("HLT_L1Tech53_MB_3_v1");
	hltPathNames.push_back("HLT_RomanPots_Tech52_v1");
	hltPathNames.push_back("HLT_L1Tech54_ZeroBias_v1");
	hltPathNames.push_back("HLT_ZeroBias_v7");

	// Declaration of histograms
	map<string,TH1F*> histosTH1F;
	map<string,TH2F*> histosTH2F;

	vector<string> selections;
	selections.push_back("All");
	selections.push_back("BunchCrossing");
	selections.push_back("HLT");
	selections.push_back("Vertex");
	selections.push_back("Muons");
	selections.push_back("EtaMax");
	selections.push_back("EtaMin");
	selections.push_back("ZeroHitsT2Plus");
	selections.push_back("ZeroHitsT2Minus");
	selections.push_back("SingleArmRP");
	selections.push_back("DoubleArmRP");
	selections.push_back("Elastic");
	selections.push_back("NonElastic");
	int nBinsEventSelection = selections.size();
	histosTH1F["EventSelection"] = new TH1F("EventSelection","EventSelection",nBinsEventSelection,0,nBinsEventSelection);
	for(size_t k = 0; k < selections.size(); ++k)
		histosTH1F["EventSelection"]->GetXaxis()->SetBinLabel( (k + 1), selections[k].c_str() );

	histosTH1F["bunchCrossingNumber"] = new TH1F("bunchCrossingNumber", "bunchCrossingNumber" , 3900 , 0 , 3900);

	histosTH1F["decisionPhysTrig"] = new TH1F("decisionPhysTrig", "decisionPhysTrig" , 128 , 0 , 128);
	histosTH1F["decisionTechTrig"] = new TH1F("decisionTechTrig", "decisionTechTrig" , 64 , 0 , 64);

	int nBinsHLT = hltPathNames.size(); 
	histosTH1F["hltTrigFired"] = new TH1F("hltTrigFired", "hltTrigFired" , nBinsHLT , 0 , nBinsHLT);
	for(size_t k = 0; k < nBinsHLT; ++k) 
		histosTH1F["hltTrigFired"]->GetXaxis()->SetBinLabel( (k + 1) , hltPathNames[k].c_str() );


	histosTH1F["vtx_zpos"] = new TH1F("vtx_zpos", "z(vtx)" , 150 , -30. , 30.);
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
	histosTH1F["Delta_eta_maxmin"] = new TH1F("Delta_eta_maxmin", "#eta^{max} - #eta^{min}" , 20 , 0,11);
	histosTH1F["xi_plus_Reco"] = new TH1F("xi_plus", "#xi^{+}" , 20 , 0,1);
	histosTH1F["xi_minus_Reco"] = new TH1F("xi_minus", "#xi^{-}" , 20 , 0,1);
	histosTH1F["log_xi_plus_Reco"] = new TH1F("log_xi_plus", "Log #xi^{+}" , 20 , -3,0.5);
	histosTH1F["log_x_plus"] = new TH1F("log_x_plus", "Log x^{+}" , 20 , -4, 0);
	histosTH1F["log_x_minus"] = new TH1F("log_x_minus", "Log x^{-}" , 20 , -4, 0);

	histosTH1F["fscHit_energy"] = new TH1F("fscHit_energy", "FSC hit energy" , 150 , -100. , 200.);
	histosTH1F["fscHit_time"] = new TH1F("fscHit_time", "FSC hit time" , 150 , 0. , 300.);

	histosTH1F["t2_track_chi2Prob_zplus"] = new TH1F("t2_track_chi2Prob_zplus", "#chi^{2}" , 100 , 0. , 1.);
	histosTH1F["t2_track_entryX_zplus"] = new TH1F("t2_track_entryX_zplus", "x_{trk}" , 160 , -160. , 160.);
	histosTH1F["t2_track_entryY_zplus"] = new TH1F("t2_track_entryY_zplus", "y_{trk}" , 160 , -160. , 160.);
	histosTH1F["t2_track_multiplicity_zplus"] = new TH1F("t2_track_multiplicity_zplus", "n tracks" , 100 , 0 , 100);
	histosTH1F["t2_track_chi2Prob_zminus"] = new TH1F("t2_track_chi2Prob_zminus", "#chi^{2}" , 100 , 0. , 1.);
	histosTH1F["t2_track_entryX_zminus"] = new TH1F("t2_track_entryX_zminus", "x_{trk}" , 160 , -160. , 160.);
	histosTH1F["t2_track_entryY_zminus"] = new TH1F("t2_track_entryY_zminus", "y_{trk}" , 160 , -160. , 160.);
	histosTH1F["t2_track_multiplicity_zminus"] = new TH1F("t2_track_multiplicity_zminus", "n tracks" , 100 , 0 , 100);
	////// Dimuon + RP
	histosTH1F["dimuon_proton_right_xi"] = new TH1F("dimuon_proton_right_xi", "#xi" , 200,-1. ,1. );
	histosTH1F["dimuon_proton_right_t"] = new TH1F("dimuon_proton_right_t", "-t" , 100 , 0.,5.);
	histosTH1F["muonDeltaPt_proton_right"] = new TH1F("muonDeltaPt_proton_right", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta_proton_right"] = new TH1F("muonDeltaEta_proton_right", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["muonDeltaPhi_proton_right"] = new TH1F("muonDeltaPhi_proton_right", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);
	histosTH1F["muonDeltaY_proton_right"] = new TH1F("muonDeltaY_proton_right", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);

	histosTH1F["dimuon_mass_proton_right"] = new TH1F("dimuon_mass_proton_right", "mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["dimuon_pt2_proton_right"] = new TH1F("dimuon_pt2_proton_right", "p_{T}^{2}(mu1,mu2)" , 100 , 0. , 1000.);
	histosTH1F["dimuon_pt_proton_right"] = new TH1F("dimuon_pt_proton_right", "p_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["dimuon_eta_proton_right"] = new TH1F("dimuon_eta_proton_right", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["dimuon_rapidity_proton_right"] = new TH1F("dimuon_rapidity_proton_right", "y(mu1,mu2)" , 100 , -15. , 15.);
	histosTH1F["dimuon_multiplicity_proton_right"] = new TH1F("dimuon_multiplicity_proton_right", "n dimuons" , 100 , 0 , 100);
	//_t_cut_proton_right
	//Float_t tbins_matrix[34] = { 0.0144, 0.0256, 0.04, 0.0576, 0.0784, 0.1024, 0.1296, 0.16, 0.1936, 0.2304, 0.2704, 0.3136, 0.36, 0.4096, 0.4624, 0.5184, 0.5778, 0.64, 0.7056, 0.7744, 0.8464, 0.9216, 1., 1.0816, 1.1664, 1.2544, 1.3456, 1.44, 1.5376, 1.6384, 1.7424, 1.8496, 1.96, 2.0736};
	histosTH1F["dimuon_mass_t_cut_proton_right"] = new TH1F("dimuon_mass_t_cut_proton_right", "mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["dimuon_pt2_t_cut_proton_right"] = new TH1F("dimuon_pt2_t_cut_proton_right", "p_{T}^{2}(mu1,mu2)" , 100 , 0. , 1000.);
	histosTH1F["dimuon_pt_t_cut_proton_right"] = new TH1F("dimuon_pt_t_cut_proton_right", "p_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["dimuon_eta_t_cut_proton_right"] = new TH1F("dimuon_eta_t_cut_proton_right", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["dimuon_rapidity_t_cut_proton_right"] = new TH1F("dimuon_rapidity_t_cut_proton_right", "y(mu1,mu2)" , 100 , -15. , 15.);
	histosTH1F["dimuon_multiplicity_t_cut_proton_right"] = new TH1F("dimuon_multiplicity_t_cut_proton_right", "n dimuons" , 100 , 0 , 100);
	histosTH1F["proton_right_xi_t_cut"] = new TH1F("proton_right_xi_t_selected", "#xi" , 200, -1.0 ,1.0);
	histosTH1F["proton_right_t_t_cut"] = new TH1F("proton_right_t_t_selected", "-t" , 100, 0., 5.0);
	//Deltas info
	histosTH1F["muonDeltaPt_t_cut_proton_right"] = new TH1F("muonDeltaPt_t_cut_proton_right", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta_t_cut_proton_right"] = new TH1F("muonDeltaEta_t_cut_proton_right", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["muonDeltaPhi_t_cut_proton_right"] = new TH1F("muonDeltaPhi_t_cut_proton_right", "#Delta#phi(mu1,mu2)" , 100 ,-2.2*M_PI , 2.2*M_PI);
	histosTH1F["muonDeltaY_t_cut_proton_right"] = new TH1F("muonDeltaY_t_cut_proton_right", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);
	//histosTH2F["DeltaPhi_vs_dimuon_pt_rpminus_accept"] = new TH2F("DeltaPhi_vs_dimuon_pt_rpminus_accept", "#Delta#phi_vs_dimuon_p_{T}" , 100, 0.0, 20.0, 100, 0.0,1.2*M_PI );
    histosTH2F["dimuon_y_vs_dimuon_pt_t_cut_proton_right"] = new TH2F("dimuon_y_vs_dimuon_pt_t_cut_proton_right", "y(mu1,mu2) vs pT" , 100 , -15. , 15.,100 , 0. , 100.);
    histosTH2F["dimuon_y_vs_dimuon_pt_t_cut_proton_left"] = new TH2F("dimuon_y_vs_dimuon_pt_t_cut_proton_left", "y(mu1,mu2) vs pT" , 100 , -15. , 15.,100 , 0. , 100.);
    histosTH2F["jpsi_y_vs_pt_t_cut_proton_right"] = new TH2F("jpsi_y_vs_pt_t_cut_proton_right", "y(mu1,mu2) vs pT" , 100 , -15. , 15.,100 , 0. , 100.);
    histosTH2F["jpsi_y_vs_pt_t_cut_proton_left"] = new TH2F("jpsi_y_vs_pt_t_cut_proton_left", "y(mu1,mu2) vs pT" , 100 , -15. , 15.,100 , 0. , 100.);
	//dimuon_proton_left
	histosTH1F["dimuon_proton_left_xi"] = new TH1F("dimuon_proton_left_xi", "#xi" , 200, -1.,1.);
	histosTH1F["dimuon_proton_left_t"] = new TH1F("dimuon_proton_left_t", "-t" , 100, 0., 5.);
	//////////
	//proton_right
	histosTH1F["jpsi_proton_right_xi"] = new TH1F("jpsi_proton_right_xi", "#xi" ,200, -1., 1.);
	histosTH1F["jpsi_proton_right_t"] = new TH1F("jpsi_proton_right_t", "-t" , 100, 0.,5.);
	histosTH1F["jpsi_mass_proton_right"] = new TH1F("jpsi_mass_proton_right", "jpsi_mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["jpsi_pt_proton_right"] = new TH1F("jpsi_pt_proton_right", "jpsi p_{T}(mu1,mu2)" , 100 , 0. , 20.);
	histosTH1F["jpsi_pt2_proton_right"] = new TH1F("jpsi_pt2_proton_right", "jpsi p_{T}^{2}(mu1,mu2)" , 100 , 0. , 1000.);
	histosTH1F["jpsi_eta_proton_right"] = new TH1F("jpsi_eta_proton_right", "jpsi #eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["jpsi_rapidity_proton_right"] = new TH1F("jpsi_rapidity_proton_right", "jpsi y(mu1,mu2)" , 100 , -15. , 15.);
	histosTH1F["jpsi_multiplicity_proton_right"] = new TH1F("jpsi_multiplicity_proton_right", "n jpsi(dimuons)" , 100 , 0 , 100);

	//Deltas info
	histosTH1F["muonDeltaPt_jpsi_proton_right"] = new TH1F("muonDeltaPt_jpsi_proton_right", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta_jpsi_proton_right"] = new TH1F("muonDeltaEta_jpsi_proton_right", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["muonDeltaPhi_jpsi_proton_right"] = new TH1F("muonDeltaPhi_jpsi_proton_right", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);
	histosTH1F["muonDeltaY_jpsi_proton_right"] = new TH1F("muonDeltaY_jpsi_proton_right", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);

	//proton_left
	histosTH1F["jpsi_proton_left_xi"] = new TH1F("jpsi_proton_left_xi", "#xi" , 200, -1.,1. );
	histosTH1F["jpsi_proton_left_t"] = new TH1F("jpsi_proton_left_t", "-t" , 100, 0.,5.);
	histosTH1F["jpsi_mass_proton_left"] = new TH1F("jpsi_mass_proton_left", "jpsi_mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["jpsi_pt_proton_left"] = new TH1F("jpsi_pt_proton_left", "jpsi p_{T}(mu1,mu2)" , 100 , 0. , 20.);
	histosTH1F["jpsi_pt2_proton_left"] = new TH1F("jpsi_pt2_proton_left", "jpsi p_{T}^{2}(mu1,mu2)" , 100 , 0. , 1000.);
	histosTH1F["jpsi_eta_proton_left"] = new TH1F("jpsi_eta_proton_left", "jpsi #eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["jpsi_rapidity_proton_left"] = new TH1F("jpsi_rapidity_proton_left", "jpsi y(mu1,mu2)" , 100 , -15. , 15.);
	histosTH1F["jpsi_multiplicity_proton_left"] = new TH1F("jpsi_multiplicity_proton_left", "n jpsi(dimuons)" , 100 , 0 , 100);

	//Deltas info
	histosTH1F["muonDeltaPt_jpsi_proton_left"] = new TH1F("muonDeltaPt_jpsi_proton_left", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta_jpsi_proton_left"] = new TH1F("muonDeltaEta_jpsi_proton_left", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["muonDeltaPhi_jpsi_proton_left"] = new TH1F("muonDeltaPhi_jpsi_proton_left", "#Delta#phi(mu1,mu2)" , 100 ,-2.2*M_PI , 2.2*M_PI);
	histosTH1F["muonDeltaY_jpsi_proton_left"] = new TH1F("muonDeltaY_jpsi_proton_left", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);        

	//_t_cut_proton_right
	histosTH1F["jpsi_dimuon_mass_t_cut_proton_right"] = new TH1F("jpsi_dimuon_mass_t_cut_proton_right", "jpsi_mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["jpsi_dimuon_pt_t_cut_proton_right"] = new TH1F("jpsi_dimuon_pt_t_cut_proton_right", "p_{T}(mu1,mu2)" , 100 , 0. , 20.);
	histosTH1F["jpsi_dimuon_pt2_t_cut_proton_right"] = new TH1F("jpsi_dimuon_pt2_t_cut_proton_right", "p_{T}^{2}(mu1,mu2)" , 100 , 0. , 1000.);
	histosTH1F["jpsi_dimuon_eta_t_cut_proton_right"] = new TH1F("jpsi_dimuon_eta_t_cut_proton_right", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["jpsi_dimuon_rapidity_t_cut_proton_right"] = new TH1F("jpsi_dimuon_rapidity_t_cut_proton_right", "y(mu1,mu2)" , 100 , -15. , 15.);
	histosTH1F["jpsi_dimuon_multiplicity_t_cut_proton_right"] = new TH1F("jpsi_dimuon_multiplicity_t_cut_proton_right", "n dimuons" , 100 , 0 , 20);

	histosTH1F["muonDeltaPt_jpsi_t_cut_proton_right"] = new TH1F("muonDeltaPt_jpsi_t_cut_proton_right", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta_jpsi_t_cut_proton_right"] = new TH1F("muonDeltaEta_jpsi_t_cut_proton_right", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["muonDeltaPhi_jpsi_t_cut_proton_right"] = new TH1F("muonDeltaPhi_jpsi_t_cut_proton_right", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);
	histosTH1F["muonDeltaY_jpsi_t_cut_proton_right"] = new TH1F("muonDeltaY_jpsi_t_cut_proton_right", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);

	histosTH1F["jpsi_proton_right_xi_t_cut"] = new TH1F("jpsi_proton_right_xi_t_selected", "#xi" , 200, -1.,1.);
	histosTH1F["jpsi_proton_right_t_t_cut"] = new TH1F("jpsi_proton_right_t_t_selected", "-t" , 100, 0.,5.);
	histosTH1F["jpsipfxiMinus_minus_proton_right_xi_t_cut"] = new TH1F("jpsipfxiMinus_minus_proton_right_xi_t_cut", "#xi diff." , 200 , -1. , 1.);
	histosTH1F["jpsiproton_t_cut_right_xi_cut"] = new TH1F("jpsiproton_t_cut_right_xi_cut", "#xi" , 200 , -1. , 1.);
	Float_t tbins[12] = { 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.22, 0.30, 0.40, 0.50, 0.65, 1.};

	histosTH1F["jpsi_Eta_max_t_cut_proton_right"]= new TH1F("Jpsi_Eta_max_t_cut_proton_right", "#eta^{max}" , 82 , etaBinsHCALBoundaries);
	histosTH1F["jpsi_Eta_min_t_cut_proton_right"]= new TH1F("Jpsi_Eta_min_t_cut_proton_right", "#eta^{min}" , 82 , etaBinsHCALBoundaries);


	histosTH1F["jpsiproton_right_t_true"] = new TH1F("jpsiproton_right_t_true", "-t" , 11 , tbins);
	histosTH1F["jpsihalo_right"] = new TH1F("jpsihalo_right", "-t halo" , 11 , tbins);
	histosTH1F["jpsiproton_right_t_signal"] = new TH1F("jpsiproton_right_t_signal", "-t" , 11 , tbins);
	histosTH1F["jpsiproton_right_t_halo"] = new TH1F("jpsiproton_right_t_halo", "-t" , 11 , tbins);
	histosTH1F["jpsi_dimuon_xi_halo_right"] = new TH1F("jpsi_dimuon_xi_halo_right", "#xi" , 200, -1.,1.);
	histosTH1F["jpsi_dimuon_t_halo_right"] = new TH1F("jpsi_dimuon_t_halo_right", "-t" , 100, 0.,5.);
	histosTH1F["jpsi_dimuon_xi_true_right"] = new TH1F("jpsi_dimuon_xi_true_right", "#xi" , 200, -1.,1.);
	histosTH1F["jpsi_dimuon_t_true_right"] = new TH1F("jpsi_dimuon_t_true_right", "-t" , 100, 0.,5.);
	histosTH1F["jpsi_dimuon_mass_halo_right"] = new TH1F("jpsi_dimuon_mass_halo_right", "jpsi_mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["jpsi_dimuon_mass_true_right"] = new TH1F("jpsi_dimuon_mass_true_right", "jpsi_mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["jpsi_dimuon_pt_halo_right"] = new TH1F("jpsi_dimuon_pt_halo_right", "p_{T}(mu1,mu2)" , 100 , 0. , 20.);
	histosTH1F["jpsi_dimuon_pt_true_right"] = new TH1F("jpsi_dimuon_pt_true_right", "p_{T}(mu1,mu2)" , 100 , 0. , 20.);
	histosTH1F["jpsi_dimuon_eta_halo_right"] = new TH1F("jpsi_dimuon_eta_halo_right", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["jpsi_dimuon_eta_true_right"] = new TH1F("jpsi_dimuon_eta_true_right", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["jpsi_dimuon_rapidity_halo_right"] = new TH1F("jpsi_dimuon_rapidity_halo_right", "y(mu1,mu2)" , 100 , -15. , 15.);
	histosTH1F["jpsi_dimuon_rapidity_true_right"] = new TH1F("jpsi_dimuon_rapidity_true_right", "y(mu1,mu2)" , 100 , -15. , 15.);
	//background
	histosTH1F["jpsi_dimuon_mass_t_cut_proton_right_bh"] = new TH1F("jpsi_dimuon_mass_t_cut_proton_right_bh", "jpsi_mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["jpsi_dimuon_pt_t_cut_proton_right_bh"] = new TH1F("jpsi_dimuon_pt_t_cut_proton_right_bh", "p_{T}(mu1,mu2)" , 100 , 0. , 20.);
	histosTH1F["jpsi_dimuon_eta_t_cut_proton_right_bh"] = new TH1F("jpsi_dimuon_eta_t_cut_proton_right_bh", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["jpsi_dimuon_rapidity_t_cut_proton_right_bh"] = new TH1F("jpsi_dimuon_rapidity_t_cut_proton_right_bh", "y(mu1,mu2)" , 100 , -15. , 15.);

	histosTH1F["jpsi_proton_right_xi_t_cut_bh"] = new TH1F("jpsi_proton_right_xi_t_selected_bh", "#xi" , 200, -1.,1.);
	histosTH1F["jpsi_proton_right_t_t_cut_bh"] = new TH1F("jpsi_proton_right_t_t_selected_bh", "-t" , 100, 0.,5.);
	//signal
	histosTH1F["jpsi_dimuon_mass_t_cut_proton_right_sig"] = new TH1F("jpsi_dimuon_mass_t_cut_proton_right_sig", "jpsi_mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["jpsi_dimuon_pt_t_cut_proton_right_sig"] = new TH1F("jpsi_dimuon_pt_t_cut_proton_right_sig", "p_{T}(mu1,mu2)" , 100 , 0. , 20.);
	histosTH1F["jpsi_dimuon_eta_t_cut_proton_right_sig"] = new TH1F("jpsi_dimuon_eta_t_cut_proton_right_sig", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["jpsi_dimuon_rapidity_t_cut_proton_right_sig"] = new TH1F("jpsi_dimuon_rapidity_t_cut_proton_right_sig", "y(mu1,mu2)" , 100 , -15. , 15.);

	histosTH1F["jpsi_proton_right_xi_t_cut_sig"] = new TH1F("jpsi_proton_right_xi_t_selected_sig", "#xi" , 200, -1.,1.);
	histosTH1F["jpsi_proton_right_t_t_cut_sig"] = new TH1F("jpsi_proton_right_t_t_selected_sig", "-t" , 100, 0.,5.);

	//histosTH1F["xi_all_before"] = new TH1F("xi_all_before", "#xi" ,100, -1.,1. ); 

	//_t_cut_proton_left
	histosTH1F["jpsi_dimuon_mass_t_cut_proton_left"] = new TH1F("jpsi_dimuon_mass_t_cut_proton_left", "jpsi_mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["jpsi_dimuon_pt_t_cut_proton_left"] = new TH1F("jpsi_dimuon_pt_t_cut_proton_left", "p_{T}(mu1,mu2)" , 100 , 0. , 20.);
	histosTH1F["jpsi_dimuon_pt2_t_cut_proton_left"] = new TH1F("jpsi_dimuon_pt2_t_cut_proton_left", "p_{T}^{2}(mu1,mu2)" , 100 , 0. , 1000.);
	histosTH1F["jpsi_dimuon_eta_t_cut_proton_left"] = new TH1F("jpsi_dimuon_eta_t_cut_proton_left", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["jpsi_dimuon_rapidity_t_cut_proton_left"] = new TH1F("jpsi_dimuon_rapidity_t_cut_proton_left", "y(mu1,mu2)" , 100 , -15. , 15.);
	histosTH1F["jpsi_dimuon_multiplicity_t_cut_proton_left"] = new TH1F("jpsi_dimuon_multiplicity_t_cut_proton_left", "n jpsi_dimuons" , 100 , 0 , 20);

	histosTH1F["jpsi_proton_left_xi_t_cut"] = new TH1F("jpsi_proton_left_xi_t_selected", "#xi" ,200, -1.,1. );
	histosTH1F["jpsi_proton_left_t_t_cut"] = new TH1F("jpsi_proton_left_t_t_selected", "-t" , 100, 0., 5.);


	histosTH1F["jpsi_Eta_max_t_cut_proton_left"]= new TH1F("Jpsi_Eta_max_t_cut_proton_left", "#eta^{max}" , 82 , etaBinsHCALBoundaries);
	histosTH1F["jpsi_Eta_min_t_cut_proton_left"]= new TH1F("Jpsi_Eta_min_t_cut_proton_left", "#eta^{min}" , 82 , etaBinsHCALBoundaries);


	histosTH1F["jpsipfxiPlus_minus_proton_left_xi_t_cut"] = new TH1F("jpsipfxiPlus_minus_proton_left_xi_t_cut", "#xi diff." , 200 , -1. , 1.);
	histosTH1F["jpsiproton_t_cut_left_xi_cut"] = new TH1F("jpsiproton_t_cut_left_xi_cut", "#xi" , 200 , -1. , 1.);

	histosTH1F["jpsiproton_left_t_true"] = new TH1F("jpsiproton_left_t_true", "-t" , 11 , tbins);
	histosTH1F["jpsihalo_left"] = new TH1F("jpsihalo_left", "-t halo" , 11 , tbins);
	histosTH1F["jpsiproton_left_t_signal"] = new TH1F("jpsiproton_left_t_signal", "-t" , 11 , tbins);
	histosTH1F["jpsiproton_left_t_halo"] = new TH1F("jpsiproton_left_t_halo", "-t" , 11 , tbins);
	histosTH1F["jpsi_dimuon_xi_halo_left"] = new TH1F("jpsi_dimuon_xi_halo_left", "#xi" , 200, -1.,1.);
	histosTH1F["jpsi_dimuon_t_halo_left"] = new TH1F("jpsi_dimuon_t_halo_left", "-t" , 100, 0.,5.);
	histosTH1F["jpsi_dimuon_xi_true_left"] = new TH1F("jpsi_dimuon_xi_true_left", "#xi" , 200, -1.,1.);
	histosTH1F["jpsi_dimuon_t_true_left"] = new TH1F("jpsi_dimuon_t_true_left", "-t" , 100, 0.,5.);
	histosTH1F["jpsi_dimuon_mass_halo_left"] = new TH1F("jpsi_dimuon_mass_halo_left", "jpsi_mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["jpsi_dimuon_mass_true_left"] = new TH1F("jpsi_dimuon_mass_true_left", "jpsi_mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["jpsi_dimuon_pt_halo_left"] = new TH1F("jpsi_dimuon_pt_halo_left", "p_{T}(mu1,mu2)" , 100 , 0. , 20.);
	histosTH1F["jpsi_dimuon_pt_true_left"] = new TH1F("jpsi_dimuon_pt_true_left", "p_{T}(mu1,mu2)" , 100 , 0. , 20.);
	histosTH1F["jpsi_dimuon_eta_halo_left"] = new TH1F("jpsi_dimuon_eta_halo_left", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["jpsi_dimuon_eta_true_left"] = new TH1F("jpsi_dimuon_eta_true_left", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["jpsi_dimuon_rapidity_halo_left"] = new TH1F("jpsi_dimuon_rapidity_halo_left", "y(mu1,mu2)" , 100 , -15. , 15.);
	histosTH1F["jpsi_dimuon_rapidity_true_left"] = new TH1F("jpsi_dimuon_rapidity_true_left", "y(mu1,mu2)" , 100 , -15. , 15.);

	//backgroun
	histosTH1F["jpsi_dimuon_mass_t_cut_proton_left_bh"] = new TH1F("jpsi_dimuon_mass_t_cut_proton_left_bh", "jpsi_mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["jpsi_dimuon_pt_t_cut_proton_left_bh"] = new TH1F("jpsi_dimuon_pt_t_cut_proton_left_bh", "p_{T}(mu1,mu2)" , 100 , 0. , 20.);
	histosTH1F["jpsi_dimuon_eta_t_cut_proton_left_bh"] = new TH1F("jpsi_dimuon_eta_t_cut_proton_left_bh", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["jpsi_dimuon_rapidity_t_cut_proton_left_bh"] = new TH1F("jpsi_dimuon_rapidity_t_cut_proton_left_bh", "y(mu1,mu2)" , 100 , -15. , 15.);
	histosTH1F["jpsi_proton_left_xi_t_cut_bh"] = new TH1F("jpsi_proton_left_xi_t_selected_bh", "#xi" ,200, -1.,1. );
	histosTH1F["jpsi_proton_left_t_t_cut_bh"] = new TH1F("jpsi_proton_left_t_t_selected_bh", "-t" , 100, 0., 5.);
	//signal
	histosTH1F["jpsi_dimuon_mass_t_cut_proton_left_sig"] = new TH1F("jpsi_dimuon_mass_t_cut_proton_left_sig", "jpsi_mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["jpsi_dimuon_pt_t_cut_proton_left_sig"] = new TH1F("jpsi_dimuon_pt_t_cut_proton_left_sig", "p_{T}(mu1,mu2)" , 100 , 0. , 20.);
	histosTH1F["jpsi_dimuon_eta_t_cut_proton_left_sig"] = new TH1F("jpsi_dimuon_eta_t_cut_proton_left_sig", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["jpsi_dimuon_rapidity_t_cut_proton_left_sig"] = new TH1F("jpsi_dimuon_rapidity_t_cut_proton_left_sig", "y(mu1,mu2)" , 100 , -15. , 15.);
	histosTH1F["jpsi_proton_left_xi_t_cut_sig"] = new TH1F("jpsi_proton_left_xi_t_selected_sig", "#xi" ,200, -1.,1. );
	histosTH1F["jpsi_proton_left_t_t_cut_sig"] = new TH1F("jpsi_proton_left_t_t_selected_sig", "-t" , 100, 0., 5.);


	//Deltas info
	histosTH1F["muonDeltaPt_jpsi_t_cut_proton_left"] = new TH1F("muonDeltaPt_jpsi_t_cut_proton_left", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta_jpsi_t_cut_proton_left"] = new TH1F("muonDeltaEta_jpsi_t_cut_proton_left", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["muonDeltaPhi_jpsi_t_cut_proton_left"] = new TH1F("muonDeltaPhi_jpsi_t_cut_proton_left", "#Delta#phi(mu1,mu2)" , 100 ,-2.2*M_PI ,2.2*M_PI);
	histosTH1F["muonDeltaY_jpsi_t_cut_proton_left"] = new TH1F("muonDeltaY_jpsi_t_cut_proton_left", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);
	//t_cut_xi_cut_proton_left
	histosTH1F["jpsi_dimuon_mass_t_cut_xi_cut_proton_left"] = new TH1F("jpsi_dimuon_mass_t_cut_xi_cut_proton_left", "jpsi_mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["jpsi_dimuon_pt_t_cut_xi_cut_proton_left"] = new TH1F("jpsi_dimuon_pt_t_cut_xi_cut_proton_left", "p_{T}(mu1,mu2)" , 100 , 0. , 20.);
	histosTH1F["jpsi_dimuon_pt2_t_cut_xi_cut_proton_left"] = new TH1F("jpsi_dimuon_pt2_t_cut_xi_cut_proton_left", "p_{T}^{2}(mu1,mu2)" , 100 , 0. , 1000.);
	histosTH1F["jpsi_dimuon_eta_t_cut_xi_cut_proton_left"] = new TH1F("jpsi_dimuon_eta_t_cut_xi_cut_proton_left", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["jpsi_dimuon_rapidity_t_cut_xi_cut_proton_left"] = new TH1F("jpsi_dimuon_rapidity_t_cut_xi_cut_proton_left", "y(mu1,mu2)" , 100 , -15. , 15.);
	histosTH1F["jpsi_dimuon_multiplicity_t_cut_xi_cut_proton_left"] = new TH1F("jpsi_dimuon_multiplicity_t_cut_xi_cut_proton_left", "n jpsi_dimuons" , 100 , 0 , 20);
	histosTH1F["jpsi_proton_left_xi_t_cut_xi_cut"] = new TH1F("jpsi_proton_left_xi_selected_t_and_xi ", "#xi" , 200, -1.,1. );
	histosTH1F["jpsi_proton_left_t_t_cut_xi_cut"] = new TH1F("jpsi_proton_left_t_selected_t_and_xi", "-t" ,  100, 0., 5.);

	//Deltas info
	histosTH1F["muonDeltaPt_jpsi_t_cut_xi_cut_proton_left"] = new TH1F("muonDeltaPt_jpsi_t_cut_xi_cut_proton_left", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta_jpsi_t_cut_xi_cut_proton_left"] = new TH1F("muonDeltaEta_jpsi_t_cut_xi_cut_proton_left", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["muonDeltaPhi_jpsi_t_cut_xi_cut_proton_left"] = new TH1F("muonDeltaPhi_jpsi_t_cut_xi_cut_proton_left", "#Delta#phi(mu1,mu2)" , 100 ,-2.2*M_PI ,2.2*M_PI);
	histosTH1F["muonDeltaY_jpsi_t_cut_xi_cut_proton_left"] = new TH1F("muonDeltaY_jpsi_t_cut_xi_cut_proton_left", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);
	//

	histosTH2F["DeltaPhi_vs_dimuon_pt"] = new TH2F("DeltaPhi_vs_dimuon_pt", "#Delta#phi_vs_dimuon_p_{T}" , 100, 0.0, 20.0, 100, 0.0,1.2*M_PI );

	histosTH2F["jpsi_DeltaPhi_vs_dimuon_pt"] = new TH2F("jpsi_DeltaPhi_vs_dimuon_pt", "#Delta#phi_vs_dimuon_p_{T}" , 100, 0.0, 20.0, 100, 0.0,1.2*M_PI );

	histosTH2F["DeltaPhi_vs_dimuon_pt_proton_right"] = new TH2F("DeltaPhi_vs_dimuon_pt_proton_right", "#Delta#phi_vs_dimuon_p_{T}" , 100, 0.0, 20.0, 100, 0.0,1.2*M_PI );

	histosTH2F["DeltaPhi_vs_dimuon_pt_t_cut_proton_right"] = new TH2F("DeltaPhi_vs_dimuon_pt_t_cut_proton_right", "#Delta#phi_vs_dimuon_p_{T}" , 100, 0.0, 20.0, 100, 0.0,1.2*M_PI );

	histosTH2F["DeltaPhi_vs_dimuon_pt_proton_left"] = new TH2F("DeltaPhi_vs_dimuon_pt_proton_left", "#Delta#phi_vs_dimuon_p_{T}" , 100, 0.0, 20.0, 100, 0.0,1.2*M_PI );

	histosTH2F["DeltaPhi_vs_dimuon_pt_t_cut_proton_left"] = new TH2F("DeltaPhi_vs_dimuon_pt_t_cut_proton_left", "#Delta#phi_vs_dimuon_p_{T}" , 100, 0.0, 20.0, 100, 0.0,1.2*M_PI );
	////
	//histosTH2F["jpsi_DeltaPhi_vs_dimuon_pt"] = new TH2F("jpsi_DeltaPhi_vs_dimuon_pt", "#Delta#phi_vs_dimuon_p_{T}" , 100, 0.0, 20.0, 100, 0.0,1.2*M_PI );

	histosTH2F["jpsi_DeltaPhi_vs_dimuon_pt_proton_right"] = new TH2F("jpsi_DeltaPhi_vs_dimuon_pt_proton_right", "#Delta#phi_vs_dimuon_p_{T}" , 100, 0.0, 20.0, 100, 0.0,1.2*M_PI ); 

	histosTH2F["jpsi_DeltaPhi_vs_dimuon_pt_t_cut_proton_right"] = new TH2F("jpsi_DeltaPhi_vs_dimuon_pt_t_cut_proton_right", "#Delta#phi_vs_dimuon_p_{T}" , 100, 0.0, 20.0, 100, 0.0,1.2*M_PI );

	histosTH2F["jpsi_DeltaPhi_vs_dimuon_pt_proton_left"] = new TH2F("jpsi_DeltaPhi_vs_dimuon_pt_proton_left", "#Delta#phi_vs_dimuon_p_{T}" , 100, 0.0, 20.0, 100, 0.0,1.2*M_PI ); 

	histosTH2F["jpsi_DeltaPhi_vs_dimuon_pt_t_cut_proton_left"] = new TH2F("jpsi_DeltaPhi_vs_dimuon_pt_t_cut_proton_left", "#Delta#phi_vs_dimuon_p_{T}" , 100, 0.0, 20.0, 100, 0.0,1.2*M_PI );

	//MC
	histosTH2F["DeltaPhi_vs_dimuon_pt_rpplus_accept"] = new TH2F("DeltaPhi_vs_dimuon_pt_rpplus_accept", "#Delta#phi_vs_dimuon_p_{T}" , 100, 0.0, 20.0, 100, 0.0,1.2*M_PI );
	histosTH2F["DeltaPhi_vs_dimuon_pt_rpplus_accept_tsel"] = new TH2F("DeltaPhi_vs_dimuon_pt_rpplus_accept_tsel", "#Delta#phi_vs_dimuon_p_{T}" , 100, 0.0, 20.0, 100, 0.0,1.2*M_PI );

	histosTH2F["DeltaPhi_vs_dimuon_pt_rpminus_accept"] = new TH2F("DeltaPhi_vs_dimuon_pt_rpminus_accept", "#Delta#phi_vs_dimuon_p_{T}" , 100, 0.0, 20.0, 100, 0.0,1.2*M_PI );
	histosTH2F["DeltaPhi_vs_dimuon_pt_rpminus_accept_tsel"] = new TH2F("DeltaPhi_vs_dimuon_pt_rpminus_accept_tsel", "#Delta#phi_vs_dimuon_p_{T}" , 100, 0.0, 20.0, 100, 0.0,1.2*M_PI );

	histosTH2F["jpsi_DeltaPhi_vs_dimuon_pt_rpplus_accept"] = new TH2F("jpsi_DeltaPhi_vs_dimuon_pt_rpplus_accept", "#Delta#phi_vs_dimuon_p_{T}" , 100, 0.0, 20.0, 100, 0.0,1.2*M_PI );
	histosTH2F["jpsi_DeltaPhi_vs_dimuon_pt_rpplus_accept_tsel"] = new TH2F("jpsi_DeltaPhi_vs_dimuon_pt_rpplus_accept_tsel", "#Delta#phi_vs_dimuon_p_{T}" , 100, 0.0, 20.0, 100, 0.0,1.2*M_PI );

	histosTH2F["jpsi_DeltaPhi_vs_dimuon_pt_rpminus_accept"] = new TH2F("jpsi_DeltaPhi_vs_dimuon_pt_rpminus_accept", "#Delta#phi_vs_dimuon_p_{T}" , 100, 0.0, 20.0, 100, 0.0,1.2*M_PI );
	histosTH2F["jpsi_DeltaPhi_vs_dimuon_pt_rpminus_accept_tsel"] = new TH2F("jpsi_DeltaPhi_vs_dimuon_pt_rpminus_accept_tsel", "#Delta#phi_vs_dimuon_p_{T}" , 100, 0.0, 20.0, 100, 0.0,1.2*M_PI );

	histosTH1F["pf_EPlusPz"] = new TH1F("pf_EPlusPz","sum(E + pz)",24,binningEPlusPz);
	histosTH1F["pf_EMinusPz"] = new TH1F("pf_EMinusPz","sum(E - pz)",24,binningEPlusPz);
    ////////////////////////////////////
    //dimuon
    

    histosTH1F["dimuon_mass_t_cut_proton_right_halo"]= new TH1F("dimuon_mass_t_cut_proton_right_halo", "mass(mu1,mu2)" , Bin_mass , 0. , 10.);
    histosTH1F["dimuon_pt_t_cut_proton_right_halo"] = new TH1F("dimuon_pt_t_cut_proton_right_halo", "p_{T}(mu1,mu2)" , 100 , 0. , 100.);  
    histosTH1F["dimuon_pt2_t_cut_proton_right_halo"] = new TH1F("dimuon_pt2_t_cut_proton_right_halo", "p_{T}^{2}(mu1,mu2)" , 100 , 0. ,   1000.); 
    histosTH1F["dimuon_eta_t_cut_proton_right_halo"] = new TH1F("dimuon_eta_t_cut_proton_right_halo", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
    histosTH1F["dimuon_rapidity_t_cut_proton_right_halo"] = new TH1F("dimuon_rapidity_t_cut_proton_right_halo", "y(mu1,mu2)" , 100 , -  15. , 15.);

    histosTH1F["dimuon_mass_t_cut_proton_right_signal"]= new TH1F("dimuon_mass_t_cut_proton_right_signal", "mass(mu1,mu2)" ,  Bin_mass , 0. , 10.);
    histosTH1F["dimuon_pt_t_cut_proton_right_signal"] = new TH1F("dimuon_pt_t_cut_proton_right_signal", "p_{T}(mu1,mu2)" , 100 ,    0. , 100.);  
    histosTH1F["dimuon_pt2_t_cut_proton_right_signal"] = new TH1F("dimuon_pt2_t_cut_proton_right_signal", "p_{T}^{2}(mu1,mu2)" ,    100 , 0. ,   1000.); 
    histosTH1F["dimuon_eta_t_cut_proton_right_signal"] = new TH1F("dimuon_eta_t_cut_proton_right_signal", "#eta(mu1,mu2)" , 100 , - 5.2 , 5.2);                                                                                                        
    histosTH1F["dimuon_rapidity_t_cut_proton_right_signal"] = new TH1F("dimuon_rapidity_t_cut_proton_right_signal", "y(mu1,mu2)" ,  100 , -  15. , 15.);
    //left


    histosTH1F["dimuon_mass_t_cut_proton_left_halo"] = new TH1F("dimuon_mass_t_cut_proton_left_halo", "mass(mu1,mu2)" , Bin_mass , 0. , 10.);
    histosTH1F["dimuon_pt_t_cut_proton_left_halo"] = new TH1F("dimuon_pt_t_cut_proton_left_halo", "p_{T}(mu1,mu2)" , 100 , 0. , 100. );
    histosTH1F["dimuon_pt2_t_cut_proton_left_halo"] = new TH1F("dimuon_pt2_t_cut_proton_left_halo", "p_{T}^{2}(mu1,mu2)" , 100 , 0. , 1000.);
    histosTH1F["dimuon_eta_t_cut_proton_left_halo"] = new TH1F("dimuon_eta_t_cut_proton_left_halo", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
    histosTH1F["dimuon_rapidity_t_cut_proton_left_halo"] = new TH1F("dimuon_rapidity_t_cut_proton_left_halo", "y(mu1,mu2)" , 100 , -15. , 15.);

    histosTH1F["dimuon_mass_t_cut_proton_left_signal"] = new TH1F("dimuon_mass_t_cut_proton_left_signal", "mass(mu1,mu2)" ,    Bin_mass , 0. ,10.);
    histosTH1F["dimuon_pt_t_cut_proton_left_signal"] = new TH1F("dimuon_pt_t_cut_proton_left_signal", "p_{T}(mu1,mu2)" , 100 , 0. , 100. );
    histosTH1F["dimuon_pt2_t_cut_proton_left_signal"] = new TH1F("dimuon_pt2_t_cut_proton_left_signal", "p_{T}^{2}(mu1,mu2)" , 100 , 0. , 1000.);
    histosTH1F["dimuon_eta_t_cut_proton_left_signal"] = new TH1F("dimuon_eta_t_cut_proton_left_signal", "#eta(mu1,mu2)" , 100 ,-5.2 , 5.2);
    histosTH1F["dimuon_rapidity_t_cut_proton_left_signal"] = new TH1F("dimuon_rapidity_t_cut_proton_left_signal", "y(mu1,mu2)" , 100 , -15. , 15.);
    //jpsi
    histosTH1F["jpsi_dimuon_eta_t_cut_proton_right_sel"] = new TH1F("jpsi_dimuon_eta_t_cut_proton_right_sel", "#eta(mu1,mu2)" , 100 ,-5.2 , 5.2);
    histosTH1F["jpsi_dimuon_eta_t_cut_proton_right_halosel"] = new TH1F("jpsi_dimuon_eta_t_cut_proton_right_halosel", "#eta(mu1,mu2)" , 100 ,-5.2 , 5.2);
    histosTH1F["jpsi_dimuon_eta_t_cut_proton_left_halosel"] = new TH1F("jpsi_dimuon_eta_t_cut_proton_left_halosel", "#eta(mu1,mu2)" , 100 ,-5.2 , 5.2);
    histosTH1F["jpsi_dimuon_eta_t_cut_proton_left_sel"] = new TH1F("jpsi_dimuon_eta_t_cut_proton_left_sel", "#eta(mu1,mu2)" , 100 ,-5.2 , 5.2);
   
    //TH2
    histosTH2F["etaMuonsVsetaJpsi_pright_Signal"]= new TH2F("etaMuonsVsetaJpsi_pright_Signal","dimuon_#eta_vs_jpsi_#eta" ,100 , -5.2 , 5.2, 100, -5.2,5.2 );
    histosTH2F["etaMuonsVsetaJpsi_pright_Halo"]= new TH2F("etaMuonsVsetaJpsi_pright_Halo","dimuon_#eta_vs_jpsi_#eta" , 100 , -5.2 ,5.2, 100, -5.2,5.2 ); 
    histosTH2F["etaMuonsVsetaJpsi_pleft_Signal"]= new TH2F("etaMuonsVsetaJpsi_pleft_Signal","dimuon_#eta_vs_jpsi_#eta" ,100 , -5.2 , 5.2, 100, -5.2,5.2 );
    histosTH2F["etaMuonsVsetaJpsi_pleft_Halo"]= new TH2F("etaMuonsVsetaJpsi_pleft_Halo","dimuon_#eta_vs_jpsi_#eta" ,100 , -5.2 , 5.2, 100, -5.2,5.2 );

    ///////////////////////////////////

	histosTH1F["dimuon_mass_proton_left"] = new TH1F("dimuon_mass_proton_left", "mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["dimuon_pt_proton_left"] = new TH1F("dimuon_pt_proton_left", "p_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["dimuon_pt2_proton_left"] = new TH1F("dimuon_pt2_proton_left", "p_{T}^{2}(mu1,mu2)" , 100 , 0. , 1000.);
	histosTH1F["dimuon_eta_proton_left"] = new TH1F("dimuon_eta_proton_left", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["dimuon_rapidity_proton_left"] = new TH1F("dimuon_rapidity_proton_left", "y(mu1,mu2)" , 100 , -15. , 15.);
	histosTH1F["dimuon_multiplicity_proton_left"] = new TH1F("dimuon_multiplicity_proton_left", "n dimuons" , 100 , 0 , 100);
	histosTH1F["muonDeltaPt_proton_left"] = new TH1F("muonDeltaPt_proton_left", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta_proton_left"] = new TH1F("muonDeltaEta_proton_left", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["muonDeltaPhi_proton_left"] = new TH1F("muonDeltaPhi_proton_left", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);
	histosTH1F["muonDeltaY_proton_left"] = new TH1F("muonDeltaY_proton_left", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);

	//_t_cut_proton_left
	histosTH1F["dimuon_mass_t_cut_proton_left"] = new TH1F("dimuon_mass_t_cut_proton_left", "mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["dimuon_pt_t_cut_proton_left"] = new TH1F("dimuon_pt_t_cut_proton_left", "p_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["dimuon_pt2_t_cut_proton_left"] = new TH1F("dimuon_pt2_t_cut_proton_left", "p_{T}^{2}(mu1,mu2)" , 100 , 0. , 1000.);
	histosTH1F["dimuon_eta_t_cut_proton_left"] = new TH1F("dimuon_eta_t_cut_proton_left", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["dimuon_rapidity_t_cut_proton_left"] = new TH1F("dimuon_rapidity_t_cut_proton_left", "y(mu1,mu2)" , 100 , -15. , 15.);
	histosTH1F["dimuon_multiplicity_t_cut_proton_left"] = new TH1F("dimuon_multiplicity_t_cut_proton_left", "n dimuons" , 100 , 0 , 100);
	histosTH1F["proton_left_xi_t_cut"] = new TH1F("proton_left_xi_t_selected", "#xi" ,  200, -1.,1.);
	histosTH1F["proton_left_t_t_cut"] = new TH1F("proton_left_t_t_selected", "-t" , 100 , 0., 5.0);

	//Deltas info
	histosTH1F["muonDeltaPt_t_cut_proton_left"] = new TH1F("muonDeltaPt_t_cut_proton_left", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta_t_cut_proton_left"] = new TH1F("muonDeltaEta_t_cut_proton_left", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["muonDeltaPhi_t_cut_proton_left"] = new TH1F("muonDeltaPhi_t_cut_proton_left", "#Delta#phi(mu1,mu2)" , 100 ,-2.2*M_PI , 2.2*M_PI);
	histosTH1F["muonDeltaY_t_cut_proton_left"] = new TH1F("muonDeltaY_t_cut_proton_left", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["jpsi_cms_sumEHFminus_rpminus_accept_tsel"] = new TH1F("jpsi_cms_sumEHFminus_rpminus_accept_tsel","sumEHF-",10,0.,200.);
	histosTH1F["jpsi_cms_sumEHFplus_rpminus_accept_tsel"] = new TH1F("jpsi_cms_sumEHFplus_rpminus_accept_tsel","sumEHF+",10, 0.,200.);
	histosTH1F["jpsi_cms_sumEHFminus_rpplus_accept_tsel"] = new TH1F("jpsi_cms_sumEHFminus_rpplus_accept_tsel","sumEHF-",10, 0.,200.);
	histosTH1F["jpsi_cms_sumEHFplus_rpplus_accept_tsel"] = new TH1F("jpsi_cms_sumEHFplus_rpplus_accept_tsel","sumEHF+",10, 0.,200.);
	histosTH1F["jpsi_cms_sumEHFminus_proton_left_t_t_cut_sig"] = new TH1F("jpsi_cms_sumEHFminus_proton_left_t_t_cut_sig","sumEHF-",10, 0.,200.);
	histosTH1F["jpsi_cms_sumEHFplus_proton_left_t_t_cut_sig"] = new TH1F("jpsi_cms_sumEHFplus_proton_left_t_t_cut_sig","sumEHF+",10, 0.,200.);
	histosTH1F["jpsi_cms_sumEHFminus_proton_left_t_t_cut_bh"] = new TH1F("jpsi_cms_sumEHFminus_proton_left_t_t_cut_bh","sumEHF-",10, 0.,200.);
	histosTH1F["jpsi_cms_sumEHFplus_proton_left_t_t_cut_bh"] = new TH1F("jpsi_cms_sumEHFplus_proton_left_t_t_cut_bh","sumEHF+",10, 0.,200.);
	histosTH1F["jpsi_cms_sumEHFminus_proton_right_t_t_cut_sig"] = new TH1F("jpsi_cms_sumEHFminus_proton_right_t_t_cut_sig","sumEHF-",10, 0.,200.);
	histosTH1F["jpsi_cms_sumEHFplus_proton_right_t_t_cut_sig"] = new TH1F("jpsi_cms_sumEHFplus_proton_right_t_t_cut_sig","sumEHF+",10, 0.,200.);

	histosTH1F["jpsi_cms_sumEHFminus_proton_right_t_t_cut_bh"] = new TH1F("jpsi_cms_sumEHFminus_proton_right_t_t_cut_bh","sumEHF-",10, 0.,200.);
	histosTH1F["jpsi_cms_sumEHFplus_proton_right_t_t_cut_bh"] = new TH1F("jpsi_cms_sumEHFplus_proton_right_t_t_cut_bh","sumEHF+",10, 0.,200.);
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
	histosTH1F["rp_track_posx_020"] = new TH1F("rp_track_posx_020", "x(RP track)" , 200, -10., 10.);
	histosTH1F["rp_track_posy_020"] = new TH1F("rp_track_posy_020", "y(RP track)" , 500, -50., 50.);
	histosTH1F["rp_track_posx_021"] = new TH1F("rp_track_posx_021", "x(RP track)" , 200, -10., 10.);
	histosTH1F["rp_track_posy_021"] = new TH1F("rp_track_posy_021", "y(RP track)" , 500, -50., 50.);
	histosTH1F["rp_track_posx_022"] = new TH1F("rp_track_posx_022", "x(RP track)" , 200, -10., 10.);
	histosTH1F["rp_track_posy_022"] = new TH1F("rp_track_posy_022", "y(RP track)" , 500, -50., 50.);
	histosTH1F["rp_track_posx_023"] = new TH1F("rp_track_posx_023", "x(RP track)" , 200, -10., 10.);
	histosTH1F["rp_track_posy_023"] = new TH1F("rp_track_posy_023", "y(RP track)" , 500, -50., 50.);
	histosTH1F["rp_track_posx_024"] = new TH1F("rp_track_posx_024", "x(RP track)" , 200, -10., 10.);
	histosTH1F["rp_track_posy_024"] = new TH1F("rp_track_posy_024", "y(RP track)" , 500, -50., 50.);
	histosTH1F["rp_track_posx_025"] = new TH1F("rp_track_posx_025", "x(RP track)" , 200, -10., 10.);
	histosTH1F["rp_track_posy_025"] = new TH1F("rp_track_posy_025", "y(RP track)" , 500, -50., 50.);

	histosTH1F["rp_track_posx_120"] = new TH1F("rp_track_posx_120", "x(RP track)" , 200, -10., 10.);
	histosTH1F["rp_track_posy_120"] = new TH1F("rp_track_posy_120", "y(RP track)" , 500, -50., 50.);
	histosTH1F["rp_track_posx_121"] = new TH1F("rp_track_posx_121", "x(RP track)" , 200, -10., 10.);
	histosTH1F["rp_track_posy_121"] = new TH1F("rp_track_posy_121", "y(RP track)" , 500, -50., 50.);
	histosTH1F["rp_track_posx_122"] = new TH1F("rp_track_posx_122", "x(RP track)" , 200, -10., 10.);
	histosTH1F["rp_track_posy_122"] = new TH1F("rp_track_posy_122", "y(RP track)" , 500, -50., 50.);
	histosTH1F["rp_track_posx_123"] = new TH1F("rp_track_posx_123", "x(RP track)" , 200, -10., 10.);
	histosTH1F["rp_track_posy_123"] = new TH1F("rp_track_posy_123", "y(RP track)" , 500, -50., 50.);
	histosTH1F["rp_track_posx_124"] = new TH1F("rp_track_posx_124", "x(RP track)" , 200, -10., 10.);
	histosTH1F["rp_track_posy_124"] = new TH1F("rp_track_posy_124", "y(RP track)" , 500, -50., 50.);
	histosTH1F["rp_track_posx_125"] = new TH1F("rp_track_posx_125", "x(RP track)" , 200, -10., 10.);
	histosTH1F["rp_track_posy_125"] = new TH1F("rp_track_posy_125", "y(RP track)" , 500, -50., 50.);

	histosTH2F["Deltaxleft_x020x024_vs_posx_020"] =  new TH2F("Deltaxleft_vs_rp_track_posx_020","#Deltax_vs_x020",  200, -20., 20.,100,-10,10);
	histosTH2F["Deltaxright_x120x124_vs_posx_120"] =  new TH2F("Deltaxright_vs_rp_track_posx_120","#Deltax_vs_x120",  200, -20., 20.,100,-10,10);  
	histosTH2F["Deltayleft_y020y024_vs_posy_020"] =  new TH2F("Deltayleft_vs_rp_track_posy_020","#Deltay_vs_y020",  500, -50., 50.,500,-50,50);
	histosTH2F["Deltayright_y120y124_vs_posy_120"] =  new TH2F("Deltayright_vs_rp_track_posy_120","#Deltay_vs_y120",  500, -50., 50.,500,-50,50);


	histosTH2F["proton_left_xi_vs_rp_track_posx_020"] =  new TH2F("proton_left_xi_vs_rp_track_posx_020","#x_vs_xi left",  200, -10., 10.,33,-0.3,0.04); 
	histosTH2F["proton_left_xi_vs_rp_track_posx_021"] =  new TH2F("proton_left_xi_vs_rp_track_posx_021","#x_vs_xi left",  200, -10., 10.,33,-0.3,0.04); 
	histosTH2F["proton_left_xi_vs_rp_track_posx_022"] =  new TH2F("proton_left_xi_vs_rp_track_posx_022","#x_vs_xi left",  200, -10., 10.,33,-0.3,0.04);         
	histosTH2F["proton_left_xi_vs_rp_track_posx_023"] =  new TH2F("proton_left_xi_vs_rp_track_posx_023","#x_vs_xi left",  200, -10., 10.,33,-0.3,0.04);   
	histosTH2F["proton_left_xi_vs_rp_track_posx_024"] =  new TH2F("proton_left_xi_vs_rp_track_posx_024","#x_vs_xi left",  200, -10., 10.,33,-0.3,0.04);
	histosTH2F["proton_left_xi_vs_rp_track_posx_025"] =  new TH2F("proton_left_xi_vs_rp_track_posx_025","#x_vs_xi left",  200, -10., 10.,33,-0.3,0.04); 

	histosTH2F["proton_right_xi_vs_rp_track_posx_120"] =  new TH2F("proton_right_xi_vs_rp_track_posx_120","#x_vs_xi rigth",  200, -10., 10.,33,-0.3,0.04); 
	histosTH2F["proton_right_xi_vs_rp_track_posx_121"] =  new TH2F("proton_right_xi_vs_rp_track_posx_121","#x_vs_xi rigth",  200, -10., 10.,33,-0.3,0.04);   
	histosTH2F["proton_right_xi_vs_rp_track_posx_122"] =  new TH2F("proton_right_xi_vs_rp_track_posx_122","#x_vs_xi rigth",  200, -10., 10.,33,-0.3,0.04);   
	histosTH2F["proton_right_xi_vs_rp_track_posx_123"] =  new TH2F("proton_right_xi_vs_rp_track_posx_123","#x_vs_xi rigth",  200, -10., 10.,33,-0.3,0.04);   
	histosTH2F["proton_right_xi_vs_rp_track_posx_124"] =  new TH2F("proton_right_xi_vs_rp_track_posx_124","#x_vs_xi rigth",  200, -10., 10.,33,-0.3,0.04);
	histosTH2F["proton_right_xi_vs_rp_track_posx_125"] =  new TH2F("proton_right_xi_vs_rp_track_posx_125","#x_vs_xi rigth",  200, -10., 10.,33,-0.3,0.04);

	histosTH1F["proton_right_xi"] = new TH1F("proton_right_xi", "#xi" , 200 , -1. , 1.);
	histosTH1F["proton_right_xi_cut"] = new TH1F("proton_right_xi_cut", "#xi" , 200 , -1. , 1.);
	histosTH1F["proton_left_xi_totem"] = new TH1F("proton_left_xi_totem", "#xi" , 200 , -1. , 1.);
	histosTH1F["proton_right_xi_totem"] = new TH1F("proton_right_xi_totem", "#xi" , 200 , -1. , 1.);
	histosTH1F["proton_right_logXi"] = new TH1F("proton_right_logXi","log(#xi)",200,-5.,0.);
	histosTH1F["proton_right_t"] = new TH1F("proton_right_t", "-t" , 100 , 0. , 5.);
	histosTH1F["proton_right_chi2"] = new TH1F("proton_right_chi2", "#chi^{2}" , 100 , 0. , 100.);
	histosTH1F["proton_left_xi"] = new TH1F("proton_left_xi", "#xi" , 200 , -1. , 1.);
	histosTH1F["proton_left_xi_cut"] = new TH1F("proton_left_xi_cut", "#xi" , 200 , -1. , 1.);
	histosTH1F["proton_left_logXi"] = new TH1F("proton_left_logXi","log(#xi)",200,-5.,0.);
	histosTH1F["proton_left_t"] = new TH1F("proton_left_t", "-t" , 100 , 0. , 5.);
	histosTH1F["proton_left_chi2"] = new TH1F("proton_left_chi2", "#chi^{2}" , 100 , 0. , 100.);

	//Float_t tbins[13] = { 0.0,0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.22, 0.30, 0.40, 0.50, 0.65, 1.};
	//===================================================
	histosTH1F["proton_right_t_sel"] = new TH1F("proton_right_t_sel", "-t" , 11 , tbins);
	histosTH1F["proton_right_t_kint"] = new TH1F("proton_right_t_kint", "-t" , 11 , tbins);
	histosTH1F["proton_right_xi_kint"] = new TH1F("proton_right_xi_kint", "#xi Right RPs" , 20 , -0.1, 0.3);
	histosTH1F["proton_right_xi_sel"] = new TH1F("proton_right_xi_sel", "#xi Right RPs" , 20 , -0.1, 0.3);
	histosTH1F["proton_right_t_halosel"] = new TH1F("proton_right_t_halosel", "-t" , 11 , tbins);
	histosTH1F["proton_right_xi_halosel"] = new TH1F("proton_right_xi_halosel", "#xi Right RPs" , 20 , -0.1, 0.3);
        histosTH2F["proton_right_xi_vs_t_halosel"] = new TH2F("proton_right_xi_vs_t_halosel", "#xi_vs_t Right RPs" , 20 , -0.1, 0.3,11 , 0.03, 1.0);
        histosTH2F["proton_right_xi_vs_t_sel"] = new TH2F("proton_right_xi_vs_t_sel", "#xi_vs_t Right RPs" , 20 , -0.1, 0.3,11 , 0.03, 1.0);
        histosTH2F["proton_right_xi_vs_t_tsel"] = new TH2F("proton_right_xi_vs_t_tsel", "#xi_vs_t Right RPs" , 200 , -1.0, 1.0, 100 , 0. , 5.);
        histosTH2F["proton_right_xi_vs_t_kin"] = new TH2F("proton_right_xi_vs_t_kin", "#xi_vs_t Right RPs" , 20 , -0.1, 0.3,11 , 0.03, 1.0);
        histosTH2F["proton_right_xi_vs_t"] = new TH2F("proton_right_xi_vs_t", "#xi_vs_t Right RPs" , 200 , -1. , 1.,100 ,0., 5.0);

	histosTH1F["proton_left_t_sel"] = new TH1F("proton_Left_t_sel", "-t" , 11 , tbins);
	histosTH1F["proton_left_t_kint"] = new TH1F("proton_left_t_kint", "-t" , 11 , tbins);
	histosTH1F["proton_left_t_halosel"] = new TH1F("proton_Left_t_halosel", "-t" , 11 , tbins);
	histosTH1F["proton_left_xi_halosel"] = new TH1F("proton_Left_xi_halosel", "#xi Left RPs" , 20 , -0.1, 0.3);
	histosTH1F["proton_left_xi_kint"] = new TH1F("proton_Left_xi_kint", "#xi Left RPs" , 20 , -0.1, 0.3);
	histosTH1F["proton_left_xi_sel"] = new TH1F("proton_Left_xi_sel", "#xi Left RPs" , 20 , -0.1, 0.3);
        histosTH2F["proton_left_xi_vs_t_halosel"] = new TH2F("proton_left_xi_vs_t_halosel", "#xi_vs_t Left RPs" , 20 , -0.1, 0.3,11 , 0.03, 1.0);
        histosTH2F["proton_left_xi_vs_t_sel"] = new TH2F("proton_left_xi_vs_t_sel", "#xi_vs_t Left RPs" , 20 , -0.1, 0.3,11 , 0.03, 1.0);
        histosTH2F["proton_left_xi_vs_t_tsel"] = new TH2F("proton_left_xi_vs_t_tsel", "#xi_vs_t Left RPs" , 200 , -1.0, 1.0, 100 , 0. , 5.);
        histosTH2F["proton_left_xi_vs_t_kin"] = new TH2F("proton_left_xi_vs_t_kin", "#xi_vs_t Left RPs" , 20 , -0.1, 0.3,11 , 0.03, 1.0);
        histosTH2F["proton_left_xi_vs_t"] = new TH2F("proton_left_xi_vs_t", "#xi_vs_t Left RPs" , 200 , -1. , 1.,100 ,0., 5.0);


	//histosTH1F["xitotem_xicms_rightRPs"] = new TH1F("xitotem_xicms_rightRPs", "Right RPs" , 20 , -0.4 , 0.4);
	histosTH1F["xitotem_xicms_rightRPs_kin"] = new TH1F("xitotem_xicms_rightRPs_kin", "Right RPs" , 20 , -0.4 , 0.4);
	histosTH1F["xitotem_xicms_leftRPs_kin"] = new TH1F("xitotem_xicms_leftRPs_kin", "Left RPs" , 20 , -0.4 , 0.4);
	//===================================================
	histosTH1F["proton_right_tbins"] = new TH1F("proton_right_tbins", "-t" , 11 , tbins);
	histosTH1F["proton_right_t_cut"] = new TH1F("proton_right_t_cut", "-t" , 11 , tbins);
	histosTH1F["proton_right_t_signal"] = new TH1F("proton_right_t_signal", "-t" , 11 , tbins);
	histosTH1F["proton_right_t_true"] = new TH1F("proton_right_t_true", "-t" , 11 , tbins);
	histosTH1F["proton_right_t_true_constbin"] = new TH1F("proton_right_t_true_constbin", "-t" , 20 , 0, 1);
	histosTH1F["proton_right_t_signal_constbin"] = new TH1F("proton_right_t_sigal_constbin", "-t" , 20 , 0, 1);
	histosTH1F["proton_right_t_signal_eff"] = new TH1F("proton_right_t_signal_eff", "-t" , 11 , tbins);
	histosTH1F["proton_right_t_signal_effweight"] = new TH1F("proton_right_t_signal_effweight", "-t" , 11 , tbins);
	histosTH1F["proton_right_t_signal_averagept_eff"] = new TH1F("proton_right_t_signal_averagept_eff", "-t" , 11 , tbins);
	histosTH1F["proton_right_t_signal_averagept_effweight"] = new TH1F("proton_right_t_signal_averagept_effweight", "-t" , 11 , tbins);
	histosTH1F["proton_right_t_halo"] = new TH1F("proton_right_t_halo", "-t" , 11 , tbins);
	histosTH1F["proton_right_t_halo_constbin"] = new TH1F("proton_right_t_halo_constbin", "-t" , 20 , 0, 1);
	histosTH1F["halo_right"] = new TH1F("halo_right", "-t halo" , 11 , tbins);
	histosTH1F["halo_right_constbin"] = new TH1F("halo_right_constbin", "-t halo" , 20 , 0, 1);
	histosTH1F["proton_right_xi_halo"] = new TH1F("proton_right_xi_halo", "#xi Right RPs" , 200 , -1., 1.0);
	histosTH1F["proton_right_xi_signal"] = new TH1F("proton_right_xi_signal", "#xi Right RPs" , 200 , -1., 1.0);
	histosTH1F["proton_right_xi_bin"] = new TH1F("proton_right_xi_bin", "#xi Right RPs" , 20 , -0.1, 0.3);
	histosTH1F["proton_right_beta"] = new TH1F("proton_right_beta", "#beta Right RPs" , 20 , 0, 1);

	histosTH1F["proton_left_tbins"] = new TH1F("proton_left_tbins", "-t" , 11 , tbins);
	histosTH1F["proton_left_t_cut"] = new TH1F("proton_left_t_cut", "-t" , 11 , tbins);
	histosTH1F["proton_left_t_signal"] = new TH1F("proton_left_t_signal", "-t" , 11 , tbins);
	histosTH1F["proton_left_t_true"] = new TH1F("proton_left_t_true", "-t" , 11 , tbins);
	histosTH1F["proton_left_t_true_constbin"] = new TH1F("proton_left_t_true_constbin", "-t" , 20 , 0, 1);
	histosTH1F["proton_left_t_signal_constbin"] = new TH1F("proton_left_t_sigal_constbin", "-t" , 20 , 0, 1);
	histosTH1F["proton_left_t_signal_eff"] = new TH1F("proton_left_t_signal_eff", "-t" , 11 , tbins);
	histosTH1F["proton_left_t_signal_effweight"] = new TH1F("proton_left_t_signal_effweight", "-t" , 11 , tbins);
	histosTH1F["proton_left_t_signal_averagept_eff"] = new TH1F("proton_left_t_signal_averagept_eff", "-t" , 11 , tbins);
	histosTH1F["proton_left_t_signal_averagept_effweight"] = new TH1F("proton_left_t_signal_averagept_effweight", "-t" , 11 , tbins);
	histosTH1F["proton_left_t_halo"] = new TH1F("proton_left_t_halo", "-t" , 11 , tbins);
	histosTH1F["proton_left_t_halo_constbin"] = new TH1F("proton_left_t_halo_constbin", "-t" , 20 , 0, 1);
	histosTH1F["halo_left"] = new TH1F("halo_left", "-t halo" , 11 , tbins);
	histosTH1F["halo_left_constbin"] = new TH1F("halo_left_constbin", "-t halo" , 20 , 0, 1);

	histosTH1F["proton_left_xi_signal"] = new TH1F("proton_left_xi_signal", "#xi Left RPs" , 200 , -1.0, 1.0);
	histosTH1F["proton_left_xi_halo"] = new TH1F("proton_left_xi_halo", "#xi Left RPs" , 200 , -1.0, 1.0);
	histosTH1F["proton_left_xi_bin"] = new TH1F("proton_left_xi_bin", "#xi Left RPs" , 200 , -1.0, 1.0);
	histosTH1F["proton_left_beta"] = new TH1F("proton_left_beta", "#beta Left RPs" , 200 , -1.0, 1.0);
	//histosTH1F["proton_left_xi_cut"] = new TH1F("proton_left_xi_tcut", "#xi Left RPs" , 250 , -1.0, 1.0);

	//histosTH1F["xitotem_xicms_leftRPs"] = new TH1F("xitotem_xicms_leftRPs", "Left RPs" , 250 , -1.0 , 1.0);
	//histosTH1F["xitotem_xicms_rightRPs"] = new TH1F("xitotem_xicms_rightRPs", "Right RPs" , 250 , -1.0 , 1.0);
	//histosTH1F["xitotem_xicms_rightRPs"] = new TH1F("xitotem_xicms_rightRPs", "Right RPs" , 20 , -0.4 , 0.4);
	//histosTH1F["xitotem_xicms_leftRPs"] = new TH1F("xitotem_xicms_leftRPs", "Left RPs" , 20 , -0.4 , 0.4);
	//histosTH1F["xitotem_xicms_rightRPs_tcut"] = new TH1F("xitotem_xicms_rightRPs_tcut", "Right RPs" , 20 , -0.4 , 0.4);
	//histosTH1F["xitotem_xicms_rightRPs_cut"] = new TH1F("xitotem_xicms_rightRPs_cut", "Right RPs" , 20 , -0.4 , 0.4);
	histosTH1F["xi_cms_totem_background_simulated"] = new TH1F("xitotem_xicms_rightRPs_simulated", "Right RPs" , 200 , -1.0 , 1.0);
	histosTH1F["xi_cms_totem_background_simulatedleft"] = new TH1F("xitotem_xicms_leftRPs_simulated", "Left RPs" , 200 , -1.0 , 1.0);

	histosTH1F["proton_pair_right_xi"] = new TH1F("proton_pair_right_xi", "#xi" , 200 , -1. , 1.);
	histosTH1F["proton_pair_right_logXi"] = new TH1F("proton_pair_right_logXi","log(#xi)",200,-5.,0.);
	histosTH1F["proton_pair_right_t"] = new TH1F("proton_pair_right_t", "-t" , 100 , 0. , 5.);
	histosTH1F["proton_pair_left_xi"] = new TH1F("proton_pair_left_xi", "#xi" , 200 , -1. , 1.);
	histosTH1F["proton_pair_left_logXi"] = new TH1F("proton_pair_left_logXi","log(#xi)",200,-5.,0.);
	histosTH1F["proton_pair_left_t"] = new TH1F("proton_pair_left_t", "-t" , 100 , 0. , 5.);
	histosTH1F["proton_pair_chi2"] = new TH1F("proton_pair_chi2", "#chi^{2}" , 100 , 0. , 100.);

	histosTH1F["pf_xiPlus_minus_proton_left_xi"] = new TH1F("pf_xiPlus_minus_proton_left_xi", "#xi diff." , 200 , -1. , 1.);
	histosTH1F["pf_xiMinus_minus_proton_right_xi"] = new TH1F("pf_xiMinus_minus_proton_right_xi", "#xi diff." , 200 , -1. , 1.);
	histosTH1F["pfxiMinus_minus_proton_right_xi"] = new TH1F("pfxiMinus_minus_proton_right_xi", "#xi diff." , 200 , -1. , 1.);
	histosTH1F["pfxiPlus_minus_proton_left_xi"] = new TH1F("pfxiPlus_minus_proton_left_xi", "#xi diff." , 200 , -1. , 1.);
	histosTH1F["jpsipfxiMinus_minus_proton_right_xi"] = new TH1F("jpsipfxiMinus_minus_proton_right_xi", "#xi diff." , 200 , -1. , 1.);
	histosTH1F["jpsipfxiPlus_minus_proton_left_xi"] = new TH1F("jpsipfxiPlus_minus_proton_left_xi", "#xi diff." , 200 , -1. , 1.);
	histosTH1F["random_pfxiMinus_minus_proton_right_xi"] = new TH1F("random_pfxiMinus_minus_proton_right_xi", "#xi diff." , 200 , -1. , 1.);
	histosTH1F["random_pfxiPlus_minus_proton_left_xi"] = new TH1F("random_pfxiPlus_minus_proton_left_xi", "#xi diff." , 200 , -1. , 1.);

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

	//////
	double energyMin = -10.;
	double energyMax = 190.;
	int nBinsEnergy = 1000;
	histosTH2F["energyVsEtaAllTypes"] = new TH2F("energyVsEtaAllTypes","energy Vs Eta AllTypes",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
	histosTH2F["energyVsEtaUndefined"] = new TH2F("energyVsEtaUndefined","energy Vs Eta Undefined",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
	histosTH2F["energy Vs Eta Charged Hadron"] = new TH2F("energyVsEtaChargedHadron","energy Vs Eta Charged Hadron",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
	histosTH2F["energyVsEtaElectron"] = new TH2F("energyVsEtaElectron","energy Vs Eta Electron",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
	histosTH2F["energyVsEtaMuon"] = new TH2F("energyVsEtaMuon","energy Vs Eta Muon",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
	histosTH2F["energyVsEtaGamma"] = new TH2F("energyVsEtaGamma","energy Vs Eta Gamma",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
	histosTH2F["energyVsEtaNeutralHadron"] = new TH2F("energyVsEtaNeutralHadron","energy Vs Eta Neutral Hadron",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
	histosTH2F["energyVsEtaHadronHF"] = new TH2F("energyVsEtaHadronHF","energy Vs Eta HadronHF",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
	histosTH2F["energyVsEtaHadronHFEcalEnergy"] = new TH2F("energyVsEtaHadronHFEcalEnergy","energyVsEtaHadronHFEcalEnergy",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
	histosTH2F["energyVsEtaHadronHFNoEcalEnergy"] = new TH2F("energyVsEtaHadronHFNoEcalEnergy","energyVsEtaHadronHFNoEcalEnergy",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
	histosTH2F["energyVsEtaEGammaHF"] = new TH2F("energyVsEtaEGammaHF","energy Vs Eta GammaHF",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
	histosTH2F["xi+Vseta_max"] = new TH2F("xi+Vseta_max","#xi^{+] Vs #eta_{max}",200,0,5.2,nBinsEnergy,0,0.5);

	for(map<string,TH1F*>::const_iterator it = histosTH1F.begin(); it != histosTH1F.end(); ++it)
		it->second->Sumw2();
	for(map<string,TH2F*>::const_iterator it = histosTH2F.begin(); it != histosTH2F.end(); ++it)
		it->second->Sumw2();
	//===================
	//int i_tot = 0 , nevt_tot = 0;
	//const char *ext=".root";

	//vector<TString>* vfiles = new vector<TString>; 
	//for(size_t idx_file = 0; idx_file < fileNames.size(); ++idx_file) vfiles->push_back( fileNames[idx_file] );
 
  int i_tot = 0 , nevt_tot = 0;
	const char *ext=".root";

	vector<TString>* vdirs = new vector<TString>; 
	vdirs->push_back("/storage1/eliza/TOTEM/MinimumBias/MinBias1-198902/");
	vdirs->push_back("/storage1/eliza/TOTEM/MinimumBias/MinBias1-198903/");
	//vdirs->push_back("/storage1/eliza/TOTEM/MinimumBias/");
	vector<TString>* vfiles = new vector<TString>; 
	//for(size_t idx_file = 0; idx_file < fileNames.size(); ++idx_file) vfiles->push_back( fileNames[idx_file] );
	for(vector<TString>::iterator itdirs = vdirs->begin(); itdirs != vdirs->end(); ++itdirs){
		TString& dirname = *itdirs;
		//vector<TString>* vfiles = new vector<TString>; 
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
	}

	//Declaration of tree and its branches variables
	TTree* tree = NULL;
	MyEvtId*           evtId        = NULL;
	MyL1TrigOld*       l1Trig       = NULL;  
	MyHLTrig*          hltTrig      = NULL;
	//vector<MyGenPart>* genPart      = NULL;
	vector<MyTracks>*  track_coll   = NULL;
	vector<MyVertex>*  vertex_coll  = NULL;
	vector<MyMuon>*    muon_coll   = NULL;
	vector<MyPFCand>*  pFlow_coll   = NULL;
	vector<MyFSCHit>*  fscHits_coll = NULL;
	vector<MyFSCDigi>* fscDigis_coll = NULL;
	//===================
	T2Event* t2_event = NULL;
	RPRootDumpReconstructedProton* rec_proton_left  = NULL;
	RPRootDumpReconstructedProton* rec_proton_right = NULL;
	RPRootDumpReconstructedProtonPair* rec_proton_pair  = NULL;
	map<unsigned int, RPRootDumpTrackInfo*> rp_track_info;
	map<unsigned int, RPRootDumpDigiInfo*> rp_digi_info;
	map<unsigned int, RPRootDumpPatternInfo*> rp_par_patterns_info;
	map<unsigned int, RPRootDumpPatternInfo*> rp_nonpar_patterns_info;
	map<unsigned int, std::vector<RPRootDumpTrackInfo>*> rp_multi_track_info;

	TString outtxt_right = outputFileName;
	TString outtxt_left = outputFileName;

	outtxt_right.ReplaceAll("root","txt_right");
	outtxt_left.ReplaceAll("root","txt_left");  

	ofstream outstring_left(outtxt_left); 
	ofstream outstring_right(outtxt_right);
	//===================  
	int n_vertices_selected = 0;
	int n_select_Vertex_After_vtx_cut =0;
	int n_evt_selectTrack = 0;
	int n_muon_selected = 0;
	int n_dimuons_selected = 0;
	int n_dimuons_t_selected_proton_left = 0;
	int n_dimuons_t_selected_proton_right = 0;
	int n_evt_PF = 0;
	int n_evt_passed_threshold = 0; 
	int n_jpsi_selected = 0;
    int n_jpsi_mass =0;
    int n_jpsi_mass_p =0;
	int n_jpsi_t_selected_proton_left = 0;
	int n_jpsi_t_selected_proton_right = 0;
	int n_jpsi_proton_left = 0;
	int n_jpsi_proton_right=0;
	int n_dimuons_proton_left = 0;
	int n_dimuons_proton_right=0;
	int nevtxisignalleft = 0;
	int nevtxisignalright = 0;
	int nevtxibkcgleft = 0;
	int nevtxibkcgright = 0;
	int nevtxibhleftjpsi = 0;
	int nevtxibhrightjpsi = 0;
	int nevtxisignalleftjpsi = 0;
	int nevtxisignalrightjpsi = 0; 

	//starting Loop over files, stops at end of list of files or when reached nevt_max
	for(vector<TString>::iterator itfiles = vfiles->begin(); itfiles != vfiles->end() && i_tot < nevt_max_corr; ++itfiles){

		TFile* file = TFile::Open(*itfiles,"READ");

		//getting the tree form the current file
		tree = (TTree*) file->Get( treeName.c_str() );

		//ofstream ofs;
		//ofs.Open ("eff_ptJet2.txt");

		//Getting number of events
		int nev = int(tree->GetEntriesFast());
		nevt_tot += nev;
		cout <<"The current file has " << nev << " entries : "<< endl << *itfiles << endl;
		//adding branches to the tree ----------------------------------------------------------------------
		tree->SetBranchAddress("cmsEvtUA",&evtId);
		tree->SetBranchAddress("cmsTrigUA",&l1Trig);
		tree->SetBranchAddress("cmsHLTTrigUA",&hltTrig);
		tree->SetBranchAddress("cmsTracksUA",&track_coll);
		tree->SetBranchAddress("cmsVerticesUA",&vertex_coll);
		tree->SetBranchAddress("cmsMuonsUA",&muon_coll);
		tree->SetBranchAddress("cmsParticleFlowUA",&pFlow_coll);
		//tree->SetBranchAddress("cmsFSCHitsUA",&fscHits_coll);
		//tree->SetBranchAddress("cmsFSCDigisUA",&fscDigis_coll);
		tree->SetBranchAddress("branchT2EV.",&t2_event);
		tree->SetBranchAddress("rec_prot_left.",&rec_proton_left);
		tree->SetBranchAddress("rec_prot_right.",&rec_proton_right);
		tree->SetBranchAddress("rec_prot_pair.",&rec_proton_pair);

		std::vector<unsigned int> rp_list;
		rp_list.push_back(20); rp_list.push_back(21); rp_list.push_back(22); rp_list.push_back(23); rp_list.push_back(24); rp_list.push_back(25);
		rp_list.push_back(120); rp_list.push_back(121); rp_list.push_back(122); rp_list.push_back(123); rp_list.push_back(124); rp_list.push_back(125);
		char br_name[200];
		for (unsigned int a = 0; a < 2; ++a) {
			int s = 2;
			for (unsigned int r = 0; r < 6; r++) {
				unsigned int id = 100 * a + 10 * s + r;
				if( std::find(rp_list.begin(), rp_list.end(), id) == rp_list.end() ) continue;

				sprintf(br_name, "track_rp_%u.", id);
				std::cout << br_name << std::endl;
				tree->SetBranchAddress(br_name, &rp_track_info[id]);
			}
		} 

		//starting loop over events, stops when reached end of file or nevt_max
		for(int i_evt = 0; i_evt < nev && i_tot < nevt_max_corr; ++i_evt , ++i_tot){

			//printing the % of events done every 10k evts
			if( ((i_tot+1) % 10000) == 0) cout <<int(double(i_tot+1)/1000)<<"k done"<<endl;

			//Filling the variables defined setting branches
			tree->GetEntry(i_evt);
			double event_weight = 1.;
			bool passed_HLTMuon= false;
			string HLT_muon = "HLT_L1Tech53_MB_2_v1";
			//string HLT_muon = "HLT_L1DoubleMu0_v1"; 
			//string HLT_muon = "HLT_ZeroBias_v7";
			bool PF_eta_max = false;
			bool PF_eta_min = false;
			bool mu1_selected = false;
			bool mu2_selected = false;
			bool select_Track = false;
			bool jpsi_mass=false;
			//bool select_OneVertex = false;
			histosTH1F["EventSelection"]->Fill( "All", event_weight );

			//AT THIS POINT ON, CAN START USING THE VARIABLES LIKE TRACKS, VERTEX ETC !
			// Event selection w/o RP
			//-------------------------------------------------------------------------------------------------
			// HLT

			for (int itrig = 0 ; itrig < 128 ; ++itrig){
				if( l1Trig->PhysTrigWord[itrig] == 1) 
					histosTH1F["decisionPhysTrig"]->Fill( itrig, event_weight );
			}

			for (int itrig = 0 ; itrig < 64 ; ++itrig){
				if( l1Trig->TechTrigWord[itrig] == 1 )
					histosTH1F["decisionTechTrig"]->Fill( itrig, event_weight );
			}

			map<string,bool>::iterator it_hlt = (*hltTrig).HLTmap.begin();
			map<string,bool>::iterator it_hlt_end = (*hltTrig).HLTmap.end();
		   for(; it_hlt != it_hlt_end; ++it_hlt){
           string const& hltName = it_hlt->first;
           vector<string>::const_iterator it_pos = find(hltPathNames.begin(),hltPathNames.end(),hltName);
           if(it_pos != hltPathNames.end()){
            size_t idx = it_pos - hltPathNames.begin();//cout <<hltName<<endl;
            if( hltName == "HLT_L1Tech53_MB_1_v1" || hltName == "HLT_L1Tech53_MB_2_v1" || hltName == "HLT_L1Tech53_MB_3_v1"){
	      if( it_hlt->second == true ){
	         passed_HLTMuon = true; 
		 histosTH1F["hltTrigFired"]->Fill( idx, event_weight );
	      }
	    }
	 }   
	         /*for(int ibin = 1; ibin <= histosTH1F["hltTrigFired"]->GetNbinsX(); ++ibin){
            if( hltName.c_str() != histosTH1F["hltTrigFired"]->GetXaxis()->GetBinLabel(ibin) ) continue;
            
            if( it_hlt->second ) 
               histosTH1F["hltTrigFired"]->Fill( histosTH1F["hltTrigFired"]->GetBinCenter( ibin ) );
         }*/ 
      }
			/*for(; it_hlt != it_hlt_end; ++it_hlt){
				string const& hltName = it_hlt->first;
				vector<string>::const_iterator it_pos = find(hltPathNames.begin(),hltPathNames.end(),hltName);
				if(it_pos != hltPathNames.end()){
					// if( it_hlt->second ) histosTH1F["hltTrigFired"]->Fill( hltName.c_str(), event_weight );

					if( hltName == HLT_muon){ 
						passed_HLTMuon = true;

						if( it_hlt->second ) histosTH1F["hltTrigFired"]->Fill( hltName.c_str(), event_weight );
					}
				}


			}*/

			if(!passed_HLTMuon) continue;

			histosTH1F["EventSelection"]->Fill( "HLT", event_weight );
			//-------------------------------------------------------------------------------------------------
			// Vertices

			double xpos; double ypos; double zpos;
			if(verbose)cout<<"MyVertex"<<endl;
			// Find number of good vertices
			int nGoodVertices = 0;
			for(vector<MyVertex>::iterator it_vtx = vertex_coll->begin() ; it_vtx != vertex_coll->end() ; ++it_vtx){

				if( it_vtx->fake ) continue;
				if( !it_vtx->validity ) continue;
				++nGoodVertices;
				if(verbose)cout<<" nGoodVertices : "<<nGoodVertices<<endl;

				zpos = it_vtx->z;
				ypos = it_vtx->y;
				xpos = it_vtx->x;
				histosTH1F["vtx_zpos"]->Fill( it_vtx->z, event_weight );
				histosTH1F["vtx_xpos"]->Fill( it_vtx->x, event_weight );
				histosTH1F["vtx_ypos"]->Fill( it_vtx->y, event_weight );
				histosTH1F["vtx_ndof"]->Fill( it_vtx->ndof, event_weight );
				histosTH1F["vtx_chi2"]->Fill( it_vtx->chi2, event_weight );
			}

			//histosTH1F["vertex_multiplicity"]->Fill(nGoodVertices, event_weight );
			//if(verbose)cout<<"Before nGoodVertices: "<<nGoodVertices<<endl;
			//bool select_OneVertex =(nGoodVertices > 0 && nGoodVertices <= 1);
			//if (selectVertex && !select_OneVertex) continue;
			//if(verbose)cout<<"After nGoodVertices= 1 : "<<nGoodVertices<<endl;
			++n_vertices_selected;

			//MyVertex& primaryVertex = vertex_coll->at(0);
			//histosTH1F["prim_vtx_zpos"]->Fill( primaryVertex.z, event_weight );
			//histosTH1F["prim_vtx_xpos"]->Fill( primaryVertex.x, event_weight );
			//histosTH1F["prim_vtx_ypos"]->Fill( primaryVertex.y, event_weight );
			//double prim_vtx_r = sqrt( primaryVertex.x*primaryVertex.x + primaryVertex.y*primaryVertex.y );
			//bool select_Vertex = ( !primaryVertex.fake && primaryVertex.validity &&
			//		primaryVertex.ndof > 4 && fabs( primaryVertex.z ) < 15.0 && prim_vtx_r < 2.0);

			//if(selectVertex && !select_Vertex) continue;

			++n_select_Vertex_After_vtx_cut;


			/*histosTH1F["vertex_multiplicity_after_vtx_sel"]->Fill( nGoodVertices, event_weight );  
			histosTH1F["prim_vtx_zpos_after_vtx_sel"]->Fill( primaryVertex.z, event_weight );
			histosTH1F["prim_vtx_xpos_after_vtx_sel"]->Fill( primaryVertex.x, event_weight );
			histosTH1F["prim_vtx_ypos_after_vtx_sel"]->Fill( primaryVertex.y, event_weight );

			histosTH1F["prim_vtx_ndof_after_vtx_sel"]->Fill( primaryVertex.ndof, event_weight );
			histosTH1F["prim_vtx_chi2_after_vtx_sel"]->Fill( primaryVertex.chi2, event_weight );
			histosTH1F["prim_vtx_chi2n_after_vtx_sel"]->Fill( primaryVertex.chi2n(), event_weight );
			histosTH1F["prim_vtx_ntracks_after_vtx_sel"]->Fill( primaryVertex.ntracks, event_weight );
			histosTH1F["prim_vtx_sumpt_after_vtx_sel"]->Fill( primaryVertex.SumPtTracks, event_weight );

			histosTH1F["EventSelection"]->Fill( "Vertex", event_weight );

			int prim_vtx_id = primaryVertex.id;*/
			//-------------------------------------------------------------------------------------------------     
			// Tracks
			/*int n_tracks_selected = 0;
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
			if(selectTrack && !select_Track) continue;*/
			++n_evt_selectTrack;
			//-------------------------------------------------------------------------------------------------
			// Muons variables 
			/*int n_muons = 0; int n_dimuons = 0.;
			double chmu1 = 0.; double chmu2 = 0.;
			double phimu1 = 0.; double ptmu1 = 0.;double etamu1 = 0.;double ymu1 = 0.;
			double phimu2 = 0.; double ptmu2 = 0.;double etamu2 = 0.;double ymu2 = 0.;
			double deltaphi = 0.; double deltaeta = 0.; double deltapt = 0.; double deltay = 0.; double Dphi = 0.;
			double dimuon_mass = 0.; double dimuon_pt=0.; double dimuon_pt2=0.; double dimuon_eta =0.;
			double dimuon_rapidity = 0.;  double jpsi_pt=0.;
			double jpsi_pt2=0.; double jpsi_eta =0.; double jpsi_rapidity = 0.; double dphijpsi = 0.;
            double muons_eta_t_cut_proton_right = 0.;
            double muons_rapidity_t_cut_proton_right = 0.;

            double muons_eta_t_cut_proton_right_halo = 0.;
            double muons_rapidity_t_cut_proton_right_halo = 0.;

            double muons_eta_t_cut_proton_right_signal = 0.;
            double muons_rapidity_t_cut_proton_right_signal = 0.;

            double muons_eta_t_cut_proton_left = 0.;
            double muons_rapidity_t_cut_proton_left = 0.;

            double muons_eta_t_cut_proton_left_halo = 0.;
            double muons_rapidity_t_cut_proton_left_halo = 0.;

            double muons_eta_t_cut_proton_left_signal = 0.;
            double muons_rapidity_t_cut_proton_left_signal = 0.;

            double jpsi_eta_t_cut_proton_right_halo = 0.;
            double jpsi_rapidity_t_cut_proton_right_halo = 0.;

            double jpsi_eta_t_cut_proton_right_signal = 0.;
            double jpsi_rapidity_t_cut_proton_right_signal = 0.;

            double jpsi_eta_t_cut_proton_left_halo = 0.;
            double jpsi_rapidity_t_cut_proton_left_halo = 0.;

            double jpsi_eta_t_cut_proton_left_signal = 0.;
            double jpsi_rapidity_t_cut_proton_left_signal = 0.; 
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

				histosTH1F["muon_multiplicity"]->Fill( n_muons, event_weight );
			}
			bool select_Muons = ( muons_selected.size() >= 2 );
			if(selectMuons && !select_Muons) continue;
			histosTH1F["EventSelection"]->Fill( "Muons", event_weight );
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
					histosTH2F["DeltaPhi_vs_dimuon_pt"]->Fill(dimuon_pt,deltaphi, event_weight ); 
					//histosTH1F["muonDphi"]->Fill(Dphi, event_weight );
					histosTH1F["lorentzdphi"]->Fill(dphijpsi, event_weight );



				}
			}
			if(!mu2_selected)continue;
			//Dimuons
			if(((dimuon_mass > 3.05) && (dimuon_mass < 3.15))){
				++n_jpsi_selected;
				jpsi_mass=true;
			    jpsi_eta = dimuon_eta;
                //double Dphi_jpsi = std::fabs(deltaPhi(phimu1 ,phimu2));
				histosTH1F["jpsi_mass"]->Fill( dimuon_mass, event_weight );
				histosTH1F["jpsi_pt"]->Fill( dimuon_pt, event_weight );
				histosTH1F["jpsi_pt2"]->Fill( dimuon_pt2, event_weight );
				histosTH1F["jpsi_eta"]->Fill( dimuon_eta, event_weight );
				histosTH1F["jpsi_rapidity"]->Fill( dimuon_rapidity, event_weight );
				histosTH1F["muonDeltaPt_jpsi"]->Fill(deltapt, event_weight );
				histosTH1F["muonDeltaEta_jpsi"]->Fill(deltaeta, event_weight );	
				histosTH1F["muonDeltaPhi_jpsi"]->Fill(deltaphi, event_weight );
				//histosTH1F["muonDphi_jpsi"]->Fill(Dphi, event_weight );
				//histosTH1F["lorentzdphi_jpsi"]->Fill(dphijpsi, event_weight );
				histosTH1F["muonDeltaY_jpsi"]->Fill(deltay, event_weight ); 
				histosTH2F["jpsi_DeltaPhi_vs_dimuon_pt"]->Fill(dimuon_pt,deltaphi, event_weight );

				if(verbose)cout<<"finalizando jpsi mass cut"<<endl; 
			}   */
			//if(!jpsi_mass)continue; cout<<dimuon_mass<<endl;
                        // ++n_jpsi_mass;
			//-----------------------------------------------------------------------------------------
			// Particle-flow
			double pfEPlusPz = 0.;
			double pfEMinusPz = 0.;
			double eta_max=-999.;
			double eta_min=999;
			//sum_E_HF
			double sumEHFMinus = 0.;
			double sumEHFPlus = 0.;
			int NMuonsPF = 0;
			//double cm = 8000;
			for(vector<MyPFCand>::iterator it_pfcand = pFlow_coll->begin(); it_pfcand != pFlow_coll->end(); ++it_pfcand){
				int partType = it_pfcand->particleId;
				double eta = it_pfcand->Eta();
				double energy = it_pfcand->Energy();
				double pz = it_pfcand->Pz();

				// Apply thresholds
				// HF eta rings 29, 30, 40, 41
				if( ( (fabs(eta) >= 2.866) && (fabs(eta) < 3.152) ) || (fabs(eta) >= 4.730) ) continue;
				/*bool Barrel_region = false;
				  bool Endcap_region = false;
				  bool Transition_region = false;
				  bool Forward_region = false;
				  bool passed_threshold = true;
				  if( eta>=-1.4 && eta<=1.4 ) Barrel_region = true;
				  else if( (eta>=-2.6 && eta<-1.4) || (eta>1.4 && eta<=2.6) ) Endcap_region = true;
				  else if( (eta>=-3.2 && eta<-2.6) || (eta>2.6 && eta<=3.2) ) Transition_region = true;
				  else if( eta<-3.2 || eta>3.2 ) Forward_region = true;
				  else cout << "ERROR!!!!!!!!!" << endl;
				//Applying threshold
				if (Barrel_region == true){
				//if(partType == MyPFCand::h0 && energy>=1.4) Treshold_Barrel1==true;
				if(partType == MyPFCand::h0 && energy < 1.4) passed_threshold = false; 
				if(partType == MyPFCand::gamma && energy<0.9) passed_threshold = false; 
				}  
				if (Endcap_region == true){
				if(partType == MyPFCand::h0 && energy<2.7) passed_threshold = false; 
				if(partType == MyPFCand::gamma && energy<2.5) passed_threshold = false; 
				}  
				if (Transition_region == true){
				if(partType == MyPFCand::h0 && energy<3.8) passed_threshold = false; 
				if(partType == MyPFCand::gamma && energy<2.5) passed_threshold = false; 
				if(partType == MyPFCand::h_HF && energy<4) passed_threshold = false; 
				if(partType == MyPFCand::egamma_HF && energy<3.5) passed_threshold = false; 
				}  
				if (Forward_region == true){
				if(partType == MyPFCand::h_HF && energy<4) passed_threshold = false; 
				if(partType == MyPFCand::egamma_HF && energy<3.5) passed_threshold = false; 
				}  

				if( !passed_threshold ) continue;*/
				// soma1 += (energy + pz);
				// soma2 += (energy - pz);
				
				//if(partType == MyPFCand::mu) ++NMuonsPF;
				
				if((3.0 < eta) && (eta < 5.1) ){
					sumEHFPlus += energy;


				}
				if((-5.1 < eta) && (eta < -3.0) ){
					sumEHFMinus += energy;

				}

				pfEPlusPz  += (energy + pz);
				pfEMinusPz += (energy - pz);

				if (eta > eta_max) {eta_max = eta; PF_eta_max = true;} 
				if (eta < eta_min) {eta_min = eta; PF_eta_min = true;}

				/*histosTH2F["energyVsEtaAllTypes"]->Fill( eta, energy, event_weight );

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
					histosTH2F["energyVsEtaEGammaHF"]->Fill( eta, energy, event_weight );*/

			}

 //MyPFCand const& pf1m=NULL;
  //MyPFCand const& pf2m=NULL;

  //int pfsize = PFCandidates->size();
  
			if(!PF_eta_max) continue;
			if(!PF_eta_min) continue;
			++n_evt_PF;

			// Xi (CMS)
			double delta_eta_maxmin = eta_max - eta_min;
			//double pfXiPlusReco = xiCorrFactor*pfEPlusPz/8000.;
			//double pfXiMinusReco = xiCorrFactor*pfEMinusPz/8000.;
			double pfXiPlusReco = pfEPlusPz/8000.;
			double pfXiMinusReco = pfEMinusPz/8000.;
			histosTH1F["xi_cms_pfplus"]->Fill(pfXiPlusReco,event_weight);
			histosTH1F["xi_cms_pfminus"]->Fill(pfXiMinusReco,event_weight);
			histosTH1F["cms_sumEHFminus"]->Fill(sumEHFMinus,event_weight);
			histosTH1F["cms_sumEHFplus"]->Fill(sumEHFPlus,event_weight);
			histosTH1F["xi_cms_pfminus"]->Fill(pfXiMinusReco,event_weight);
			histosTH1F["Eta_max"]->Fill( eta_max, event_weight  );
			histosTH1F["Eta_min"]->Fill( eta_min, event_weight  );
			histosTH1F["Delta_eta_maxmin"]->Fill( delta_eta_maxmin, event_weight  );

			/*if(((dimuon_mass > 3.05) && (dimuon_mass < 3.15))){
			  ++n_jpsi_mass;
			  histosTH1F["jpsi_xi_cms_pfplus"]->Fill(pfXiPlusReco,event_weight);
			  histosTH1F["jpsi_xi_cms_pfminus"]->Fill(pfXiMinusReco,event_weight);
			  histosTH1F["jpsi_Eta_max"]->Fill( eta_max, event_weight  );
			  histosTH1F["jpsi_Eta_min"]->Fill( eta_min, event_weight  );
			  histosTH1F["jpsi_Delta_eta_maxmin"]->Fill( delta_eta_maxmin, event_weight  );
			  histosTH1F["jpsi_cms_sumEHFminus"]->Fill(sumEHFMinus,event_weight);
			  histosTH1F["jpsi_cms_sumEHFplus"]->Fill(sumEHFPlus,event_weight);
			}*/
			//Find eta_min & eta_max
			histosTH1F["EventSelection"]->Fill( "EtaMax", event_weight );
			histosTH1F["EventSelection"]->Fill( "EtaMin", event_weight );
            //cout<< "fim do PF"<<endl;
			//-----------------------------------------------------------------------------------------
			// TOTEM T2
			int n_t2_tracks_selected = 0;
			int n_t2_tracks_selected_zplus = 0;
			int n_t2_tracks_selected_zminus = 0;

			vector<double> const& t2_trk_entryX = t2_event->TrkEntryX;
			vector<double> const& t2_trk_entryY = t2_event->TrkEntryY;
			vector<double> const& t2_trk_entryZ =  t2_event->TrkEntryZ;
			vector<double> const& t2_trk_chiProb =  t2_event->TrkChiProb;

			size_t n_t2_tracks = t2_trk_chiProb.size();
			for(size_t i_t2_trk = 0; i_t2_trk < n_t2_tracks; ++i_t2_trk){
				double trk_entryZ = t2_trk_entryZ[i_t2_trk];
				int zside = ( trk_entryZ >= 0. ) ? 1 : -1;
				if( zside > 0 )
					histosTH1F["t2_track_chi2Prob_zplus"]->Fill( t2_trk_chiProb[i_t2_trk], event_weight );
				else
					histosTH1F["t2_track_chi2Prob_zminus"]->Fill( t2_trk_chiProb[i_t2_trk], event_weight );

				// Select tracks
				if( t2_trk_chiProb[i_t2_trk] < 0.2 ) continue;

				++n_t2_tracks_selected;
				if( zside > 0 ) ++n_t2_tracks_selected_zplus;
				else            ++n_t2_tracks_selected_zminus;

				if( zside > 0 ){
					histosTH1F["t2_track_entryX_zplus"]->Fill( t2_trk_entryX[i_t2_trk], event_weight );
					histosTH1F["t2_track_entryY_zplus"]->Fill( t2_trk_entryY[i_t2_trk], event_weight );
					histosTH2F["t2_track_entryY_vs_entryX_zplus"]->Fill( t2_trk_entryX[i_t2_trk], t2_trk_entryY[i_t2_trk], event_weight );
				} else{
					histosTH1F["t2_track_entryX_zminus"]->Fill( t2_trk_entryX[i_t2_trk], event_weight );
					histosTH1F["t2_track_entryY_zminus"]->Fill( t2_trk_entryY[i_t2_trk], event_weight );
					histosTH2F["t2_track_entryY_vs_entryX_zminus"]->Fill( t2_trk_entryX[i_t2_trk], t2_trk_entryY[i_t2_trk], event_weight );
				}
			}

			if( selectZeroHitsT2Plus && (n_t2_tracks_selected_zplus > 0) ) continue;
			histosTH1F["EventSelection"]->Fill( "ZeroHitsT2Plus", event_weight );

			if( selectZeroHitsT2Minus && (n_t2_tracks_selected_zminus > 0) ) continue;
			histosTH1F["EventSelection"]->Fill( "ZeroHitsT2Minus", event_weight );
			//=============================================================================	
			//bool proton_right_valid = false;
			//bool proton_left_valid = false;

			if(verbose)cout<<"before - select proton"<<endl;
			bool proton_right_valid = rec_proton_right->valid;
			bool proton_left_valid = rec_proton_left->valid;
			bool double_arm_rec_proton = (proton_right_valid && proton_left_valid);
			bool single_arm_rec_proton = (proton_right_valid || proton_left_valid) && !double_arm_rec_proton;

			if( selectSingleArmRecProton && !single_arm_rec_proton ) continue;
			histosTH1F["EventSelection"]->Fill( "SingleArmRP", event_weight );

			if( selectDoubleArmRecProton && !double_arm_rec_proton ) continue;
			histosTH1F["EventSelection"]->Fill( "DoubleArmRP", event_weight );

			bool tag_elastic_top45_bot56 = elastic_top45_bot56(rp_track_info);      
			bool tag_elastic_bot45_top56 = elastic_bot45_top56(rp_track_info);      
			if( selectElastic && !(tag_elastic_top45_bot56 || tag_elastic_bot45_top56) ) continue;
			histosTH1F["EventSelection"]->Fill( "Elastic", event_weight );

			if( selectNonElastic && (tag_elastic_top45_bot56 || tag_elastic_bot45_top56) ) continue;
			histosTH1F["EventSelection"]->Fill( "NonElastic", event_weight );

			// RP Track protons
			if(verbose)cout<<"RP_track_info"<<endl;
			RPRootDumpTrackInfo* rp_track_info_020 = rp_track_info[20];
			RPRootDumpTrackInfo* rp_track_info_021 = rp_track_info[21];
			RPRootDumpTrackInfo* rp_track_info_022 = rp_track_info[22];
			RPRootDumpTrackInfo* rp_track_info_023 = rp_track_info[23];
			RPRootDumpTrackInfo* rp_track_info_024 = rp_track_info[24];
			RPRootDumpTrackInfo* rp_track_info_025 = rp_track_info[25];
			if(verbose)cout<<"RP_track_info left"<<endl;
			RPRootDumpTrackInfo* rp_track_info_120 = rp_track_info[120];
			RPRootDumpTrackInfo* rp_track_info_121 = rp_track_info[121];
			RPRootDumpTrackInfo* rp_track_info_122 = rp_track_info[122];
			RPRootDumpTrackInfo* rp_track_info_123 = rp_track_info[123];
			RPRootDumpTrackInfo* rp_track_info_124 = rp_track_info[124];
			RPRootDumpTrackInfo* rp_track_info_125 = rp_track_info[125];
			if(verbose) cout<<"RP_track_info right"<<endl;
			double xi_left = rec_proton_left->xi;                       
			double xi_right = rec_proton_right->xi;
			double deltax020024 = 0.;
			double deltax120124 = 0.;
			double deltay020024 = 0.;
			double deltay120124 = 0.;

			bool rp_track_valid_020 = rp_track_info_020->valid && proton_left_valid;
			if( rp_track_valid_020 ){
				double rp_track_posx_020 = rp_track_info_020->x;
				double rp_track_posy_020 = rp_track_info_020->y;
				histosTH1F["rp_track_posx_020"]->Fill( rp_track_posx_020, event_weight );
				histosTH1F["rp_track_posy_020"]->Fill( rp_track_posy_020, event_weight );
				histosTH2F["rp_track_pos_y_vs_x_020"]->Fill( rp_track_posx_020, rp_track_posy_020, event_weight );
				histosTH2F["proton_left_xi_vs_rp_track_posx_020"]->Fill( rp_track_posx_020, xi_left, event_weight );
				if(verbose) cout<<"RP_track_info 020"<<endl;

			}
			bool rp_track_valid_021 = rp_track_info_021->valid && proton_left_valid;
			if( rp_track_valid_021 ){
				double rp_track_posx_021 = rp_track_info_021->x;
				double rp_track_posy_021 = rp_track_info_021->y;
				histosTH1F["rp_track_posx_021"]->Fill( rp_track_posx_021, event_weight );
				histosTH1F["rp_track_posy_021"]->Fill( rp_track_posy_021, event_weight );
				histosTH2F["rp_track_pos_y_vs_x_021"]->Fill( rp_track_posx_021, rp_track_posy_021, event_weight );
				histosTH2F["proton_left_xi_vs_rp_track_posx_021"]->Fill( rp_track_posx_021, xi_left, event_weight );
				if(verbose) cout<<"RP_track_info 021"<<endl;
			}
			bool rp_track_valid_022 = rp_track_info_022->valid && proton_left_valid;
			if( rp_track_valid_022 ){
				double rp_track_posx_022 = rp_track_info_022->x;
				double rp_track_posy_022 = rp_track_info_022->y;
				histosTH1F["rp_track_posx_022"]->Fill( rp_track_posx_022, event_weight );
				histosTH1F["rp_track_posy_022"]->Fill( rp_track_posy_022, event_weight );
				histosTH2F["rp_track_pos_y_vs_x_022"]->Fill( rp_track_posx_022, rp_track_posy_022, event_weight );
				histosTH2F["proton_left_xi_vs_rp_track_posx_022"]->Fill( rp_track_posx_022, xi_left, event_weight );
				if(verbose) cout<<"RP_track_info 022"<<endl;
			}
			bool rp_track_valid_023 = rp_track_info_023->valid && proton_left_valid;
			if( rp_track_valid_023 ){
				double rp_track_posx_023 = rp_track_info_023->x;
				double rp_track_posy_023 = rp_track_info_023->y;
				histosTH1F["rp_track_posx_023"]->Fill( rp_track_posx_023, event_weight );
				histosTH1F["rp_track_posy_023"]->Fill( rp_track_posy_023, event_weight );
				histosTH2F["rp_track_pos_y_vs_x_023"]->Fill( rp_track_posx_023, rp_track_posy_023, event_weight );
				histosTH2F["proton_left_xi_vs_rp_track_posx_023"]->Fill( rp_track_posx_023, xi_left, event_weight );
				if(verbose) cout<<"RP_track_info 023"<<endl;
			}
			bool rp_track_valid_024 = rp_track_info_024->valid && proton_left_valid;
			//double xi_left = rec_proton_left->xi;
			if( rp_track_valid_024 ){
				double rp_track_posx_024 = rp_track_info_024->x;
				double rp_track_posy_024 = rp_track_info_024->y;
				histosTH1F["rp_track_posx_024"]->Fill( rp_track_posx_024, event_weight );
				histosTH1F["rp_track_posy_024"]->Fill( rp_track_posy_024, event_weight );
				histosTH2F["rp_track_pos_y_vs_x_024"]->Fill( rp_track_posx_024, rp_track_posy_024, event_weight );
				//histosTH2F["proton_left_xi_vs_rp_track_posx_024"]->Fill( rp_track_posx_024, xi_left, event_weight );
				if( rp_track_valid_020 ){
					deltay020024 = ( rp_track_info_020->y -  rp_track_info_024->y);
					deltax020024 = ( rp_track_info_020->x -  rp_track_info_024->x);
					histosTH2F["Deltaxleft_x020x024_vs_posx_020"]->Fill( deltax020024,  rp_track_info_020->x, event_weight );
					histosTH2F["Deltayleft_y020y024_vs_posy_020"]->Fill( deltay020024,  rp_track_info_020->y, event_weight );
					histosTH2F["proton_left_xi_vs_rp_track_posx_024"]->Fill( rp_track_posx_024, xi_left, event_weight );
					if(verbose) cout<<"RP_track_info deltax020024"<<endl;
				}
				if(verbose) cout<<"RP_track_info 024"<<endl;
			}
			bool rp_track_valid_025 = rp_track_info_025->valid && proton_left_valid;
			if( rp_track_valid_025 ){
				double rp_track_posx_025 = rp_track_info_025->x;
				double rp_track_posy_025 = rp_track_info_025->y;
				histosTH1F["rp_track_posx_025"]->Fill( rp_track_posx_025, event_weight );
				histosTH1F["rp_track_posy_025"]->Fill( rp_track_posy_025, event_weight );
				histosTH2F["rp_track_pos_y_vs_x_025"]->Fill( rp_track_posx_025, rp_track_posy_025, event_weight );
				histosTH2F["proton_left_xi_vs_rp_track_posx_025"]->Fill( rp_track_posx_025, xi_left, event_weight );
				if(verbose) cout<<"RP_track_info 025"<<endl;
			}

			bool rp_track_valid_120 = rp_track_info_120->valid && proton_right_valid;
			if( rp_track_valid_120 ){
				double rp_track_posx_120 = rp_track_info_120->x;
				double rp_track_posy_120 = rp_track_info_120->y;
				histosTH1F["rp_track_posx_120"]->Fill( rp_track_posx_120, event_weight );
				histosTH1F["rp_track_posy_120"]->Fill( rp_track_posy_120, event_weight );
				histosTH2F["rp_track_pos_y_vs_x_120"]->Fill( rp_track_posx_120, rp_track_posy_120, event_weight );
				histosTH2F["proton_right_xi_vs_rp_track_posx_120"]->Fill( rp_track_posx_120, xi_right, event_weight );
				if(verbose) cout<<"RP_track_info 120"<<endl;
			}
			bool rp_track_valid_121 = rp_track_info_121->valid && proton_right_valid;
			if( rp_track_valid_121 ){
				double rp_track_posx_121 = rp_track_info_121->x;
				double rp_track_posy_121 = rp_track_info_121->y;
				histosTH1F["rp_track_posx_121"]->Fill( rp_track_posx_121, event_weight );
				histosTH1F["rp_track_posy_121"]->Fill( rp_track_posy_121, event_weight );
				histosTH2F["rp_track_pos_y_vs_x_121"]->Fill( rp_track_posx_121, rp_track_posy_121, event_weight );
				histosTH2F["proton_right_xi_vs_rp_track_posx_121"]->Fill( rp_track_posx_121, xi_right, event_weight );
				if(verbose) cout<<"RP_track_info 121"<<endl;
			}
			bool rp_track_valid_122 = rp_track_info_122->valid && proton_right_valid;
			if( rp_track_valid_122 ){
				double rp_track_posx_122 = rp_track_info_122->x;
				double rp_track_posy_122 = rp_track_info_122->y;
				histosTH1F["rp_track_posx_122"]->Fill( rp_track_posx_122, event_weight );
				histosTH1F["rp_track_posy_122"]->Fill( rp_track_posy_122, event_weight );
				histosTH2F["rp_track_pos_y_vs_x_122"]->Fill( rp_track_posx_122, rp_track_posy_122, event_weight );
				histosTH2F["proton_right_xi_vs_rp_track_posx_122"]->Fill( rp_track_posx_122, xi_right, event_weight );
				if(verbose) cout<<"RP_track_info 122"<<endl;
			}
			bool rp_track_valid_123 = rp_track_info_123->valid && proton_right_valid;
			if( rp_track_valid_123 ){
				double rp_track_posx_123 = rp_track_info_123->x;
				double rp_track_posy_123 = rp_track_info_123->y;
				histosTH1F["rp_track_posx_123"]->Fill( rp_track_posx_123, event_weight );
				histosTH1F["rp_track_posy_123"]->Fill( rp_track_posy_123, event_weight );
				histosTH2F["rp_track_pos_y_vs_x_123"]->Fill( rp_track_posx_123, rp_track_posy_123, event_weight );
				histosTH2F["proton_right_xi_vs_rp_track_posx_123"]->Fill( rp_track_posx_123, xi_right, event_weight );
				if(verbose) cout<<"RP_track_info 123"<<endl;
			}
			bool rp_track_valid_124 = rp_track_info_124->valid && proton_right_valid;
			//double xi_right = rec_proton_right->xi;
			if( rp_track_valid_124 ){
				double rp_track_posx_124 = rp_track_info_124->x;
				double rp_track_posy_124 = rp_track_info_124->y;
				histosTH1F["rp_track_posx_124"]->Fill( rp_track_posx_124, event_weight );
				histosTH1F["rp_track_posy_124"]->Fill( rp_track_posy_124, event_weight );
				histosTH2F["rp_track_pos_y_vs_x_124"]->Fill( rp_track_posx_124, rp_track_posy_124, event_weight );
				//histosTH2F["proton_right_xi_vs_rp_track_posx_124"]->Fill( rp_track_posx_124, xi_right, event_weight );
				if( rp_track_valid_120 ){
					deltax120124 = ( rp_track_info_120->x - rp_track_posx_124);
					deltay120124 = ( rp_track_info_120->y - rp_track_posy_124);
					histosTH2F["Deltaxright_x120x124_vs_posx_120"]->Fill( deltax120124,rp_track_info_120->x, event_weight );
					histosTH2F["Deltayright_y120y124_vs_posy_120"]->Fill( deltay120124,rp_track_info_120->y, event_weight );
					histosTH2F["proton_right_xi_vs_rp_track_posx_124"]->Fill( rp_track_posx_124, xi_right, event_weight );
				}
				if(verbose) cout<<"RP_track_info 124"<<endl;
			}
			bool rp_track_valid_125 = rp_track_info_125->valid && proton_right_valid;
			if( rp_track_valid_125 ){
				double rp_track_posx_125 = rp_track_info_125->x;
				double rp_track_posy_125 = rp_track_info_125->y;
				histosTH1F["rp_track_posx_125"]->Fill( rp_track_posx_125, event_weight );
				histosTH1F["rp_track_posy_125"]->Fill( rp_track_posy_125, event_weight );
				histosTH2F["rp_track_pos_y_vs_x_125"]->Fill( rp_track_posx_125, rp_track_posy_125, event_weight );
				histosTH2F["proton_right_xi_vs_rp_track_posx_125"]->Fill( rp_track_posx_125, xi_right, event_weight );
				if(verbose) cout<<"RP_track_info 125"<<endl;
			}

			//RP Track protons
			bool proton_right_rp_accept = ( ( rp_track_valid_120 && rp_track_valid_124 )|| 
					( rp_track_valid_121 && rp_track_valid_125) );
			
			bool fiducial_right_xcut120124 = ((rp_track_info_120->x > -1.5 && rp_track_info_120->x < 6.5) && (rp_track_info_124->x > -1.5 && rp_track_info_124->x < 6.5));
			bool fiducial_right_xcut121125 = ((rp_track_info_121->x > -1.5 && rp_track_info_121->x < 6.5) && (rp_track_info_125->x > -1.5 && rp_track_info_125->x < 6.5));
            
			bool fiducial_right_ycut120124 = ((rp_track_info_120->y > 7.0 && rp_track_info_120->y < 29.0) && (rp_track_info_124->y > 7.0 && rp_track_info_124->x < 29.0));
			//bool fiducial_right_ycut121125 = ((fabs(rp_track_info_121->y) > 7.0 && fabs(rp_track_info_121->y) < 29.0) && (fabs(rp_track_info_125->y) > 7.0 && fabs(rp_track_info_125->y) < 29.0));
			bool fiducial_right_ycut121125 = (((rp_track_info_121->y) < -7.0 && (rp_track_info_121->y) > -29.0) && ((rp_track_info_125->y) < -7.0 && (rp_track_info_125->y) > -29.0));
            
            
			bool proton_left_rp_accept = ( ( rp_track_valid_020 && rp_track_valid_024 )|| 
					( rp_track_valid_021 && rp_track_valid_025) );

			bool fiducial_left_xcut020024 = ((rp_track_info_020->x > -1.5 && rp_track_info_020->x < 6.5) && (rp_track_info_024->x > -1.5  && rp_track_info_024->x < 6.5) );                                
			bool fiducial_left_xcut021025 = ((rp_track_info_021->x > -1.5 && rp_track_info_021->x < 6.5) && (rp_track_info_025->x > -1.5 && rp_track_info_025->x < 6.5));

			bool fiducial_left_ycut020024 = ((rp_track_info_020->y > 7.0 && rp_track_info_020->y < 29.0) && (rp_track_info_024->y > 7.0  && rp_track_info_024->x < 29.0)); 
			bool //fiducial_left_ycut021025 = ((fabs(rp_track_info_021->y) > 7.0 && fabs(rp_track_info_021->y) < 29.0) && (fabs(rp_track_info_025->y) > 7.0 && fabs(rp_track_info_025->y) < 29.0));
			fiducial_left_ycut021025 =( ((rp_track_info_021->y) < -7.0 && (rp_track_info_021->y) > -29.0) && ((rp_track_info_025->y) < -7.0 && (rp_track_info_025->y) > -29.0));
			/////////////////////////////////////////////
			double chi2_proton_right = rec_proton_right->chi2;
			double xi_proton_right = rec_proton_right->xi;
			double t_proton_right = rec_proton_right->t;
			//bool good_proton_right =( proton_right_valid && ((-0.23 < xi_proton_right)&&( xi_proton_right < 0.04)) && proton_right_rp_accept && ((fiducial_right_xcut120124 && fiducial_right_ycut120124) ||(fiducial_right_xcut121125 && fiducial_right_ycut121125)));
			bool good_proton_right = ( proton_right_valid && proton_right_rp_accept && ((fiducial_right_xcut120124 &&  fiducial_right_ycut120124) ||(fiducial_right_xcut121125 && fiducial_right_ycut121125)));
			double chi2_proton_left = rec_proton_left->chi2;
			double xi_proton_left = rec_proton_left->xi;
			double t_proton_left = rec_proton_left->t;
			//bool good_proton_left = (proton_left_valid && ((-0.23 < xi_proton_left)&&( xi_proton_left < 0.04))&& proton_left_rp_accept && ((fiducial_left_xcut020024 && fiducial_left_ycut020024) || (fiducial_left_xcut021025 && fiducial_left_ycut021025)));
            bool good_proton_left = ( proton_left_valid && proton_left_rp_accept && ((fiducial_left_xcut020024 &&         fiducial_left_ycut020024) || (fiducial_left_xcut021025 && fiducial_left_ycut021025)) );
                      
			bool t_region_right = -t_proton_right>=0.03 && -t_proton_right<=1;
			bool t_region_left = -t_proton_left>=0.03 && -t_proton_left<=1;
			bool runs_cms_902 = (evtId->Run==198902);
			bool runs_cms_903 = (evtId->Run==198903);
            bool evet_1 =  (evtId->Evt=317193);
            bool evet_2 =  (evtId->Evt=3197060);
            bool evet_3 =  (evtId->Evt=3427673);
            bool evet_4 =  (evtId->Evt=6577058);
            bool evet_5 =  (evtId->Evt=9906633);
            bool evet_6 =  (evtId->Evt=10146018); 
            bool evet_7 =  (evtId->Evt=12335728);
            bool evet_8 =  (evtId->Evt=16583902);
            bool evet_9 =  (evtId->Evt=18985082);
            bool evet_10 = (evtId->Evt=3450535);
            bool evet_11 = (evtId->Evt=3482302);
			bool list_evts = (evtId->Evt=317193)||(evtId->Evt=3197060)||(evtId->Evt=3427673)||(evtId->Evt=6577058)||(evtId->Evt=9906633)||(evtId->Evt=10146018)||(evtId->Evt==12335728)||(evtId->Evt=16583902)||(evtId->Evt=18985082)||(evtId->Evt=3450535)||(evtId->Evt=3482302);
			//outstring_left << evtId->Run << ":"<< evtId->LumiSect << ":"<< evtId->Evt << endl;
//=============================================

                    cout<< evtId->Run << ":"<< evtId->LumiSect << ":"<< evtId->Evt << endl; 
       
//=============================================       
//if(evet_1 || evet_2 || evet_3 || evet_4 || evet_5 || evet_6 || evet_7 || evet_8 || evet_9 || evet_10 || evet_11){ 
            if(list_evts){
            //if(jpsi_mass){
              // cout<<dimuon_mass<<endl;
               //++n_jpsi_mass_p;
               cout<<"after: "<< evtId->Run << ":"<< evtId->LumiSect << ":"<< evtId->Evt << endl; 
			   if( good_proton_right ){
				 cout<<"-good proton rigth"<<endl;
                 cout<<"after right: "<< evtId->Run << ":"<< evtId->LumiSect << ":"<< evtId->Evt << endl; 

				histosTH1F["proton_right_chi2"]->Fill( chi2_proton_right, event_weight ); 
				histosTH1F["proton_right_xi"]->Fill(- xi_proton_right, event_weight );
				histosTH1F["proton_right_xi_bin"]->Fill(- xi_proton_right, event_weight );
				histosTH1F["proton_right_t"]->Fill( -t_proton_right, event_weight );
				histosTH1F["proton_right_tbins"]->Fill( -t_proton_right, event_weight );
				histosTH2F["proton_right_xi_vs_pf_xiPlus"]->Fill( pfXiPlusReco, -xi_proton_right, event_weight );
				histosTH2F["proton_right_xi_vs_pf_xiMinus"]->Fill( pfXiMinusReco, -xi_proton_right, event_weight );
                histosTH2F["proton_right_xi_vs_t"]->Fill( -xi_proton_right,-t_proton_right, event_weight );
				//histosTH1F["proton_right_dimuon_mass"]->Fill( dimuon_mass, event_weight );
				if (pfXiMinusReco>0.05){ 
					histosTH1F["proton_right_xi_cut"]->Fill(-xi_proton_right, event_weight );
					if(verbose) cout<<"proton_right_xi_cut"<<-xi_proton_right<<endl;
				} 


				histosTH1F["pfxiMinus_minus_proton_right_xi"]->Fill( (pfXiMinusReco + xi_proton_right), event_weight );
				if (t_region_right) {

					n_jpsi_t_selected_proton_right++;
					/*histosTH1F["jpsi_dimuon_mass_t_cut_proton_right"]->Fill( dimuon_mass, event_weight );
					histosTH1F["jpsi_dimuon_pt_t_cut_proton_right"]->Fill( dimuon_pt, event_weight );
					histosTH1F["jpsi_dimuon_pt2_t_cut_proton_right"]->Fill( dimuon_pt2, event_weight );
					histosTH1F["jpsi_dimuon_eta_t_cut_proton_right"]->Fill( dimuon_eta, event_weight );
					histosTH1F["jpsi_dimuon_rapidity_t_cut_proton_right"]->Fill( dimuon_rapidity, event_weight );
					histosTH1F["muonDeltaPt_jpsi_t_cut_proton_right"]->Fill(deltapt, event_weight );
					histosTH1F["muonDeltaEta_jpsi_t_cut_proton_right"]->Fill(deltaeta, event_weight );
					histosTH1F["muonDeltaPhi_jpsi_t_cut_proton_right"]->Fill(deltaphi, event_weight );
					histosTH1F["muonDeltaY_jpsi_t_cut_proton_right"]->Fill(deltay, event_weight );
					histosTH2F["jpsi_DeltaPhi_vs_dimuon_pt_t_cut_proton_right"]->Fill(dimuon_pt,deltaphi, event_weight );
					histosTH2F["jpsi_y_vs_pt_t_cut_proton_right"]->Fill(dimuon_rapidity,dimuon_pt, event_weight );
*/
					histosTH1F["proton_right_xi_kint"]->Fill(-xi_proton_right, event_weight );
					histosTH1F["proton_right_t_kint"]->Fill(-t_proton_right, event_weight );
					histosTH1F["xitotem_xicms_rightRPs_kin"]->Fill( pfXiMinusReco+xi_proton_right, event_weight ); 
                    histosTH2F["proton_right_xi_vs_t_kin"]->Fill( -xi_proton_right,-t_proton_right, event_weight );       
                    histosTH2F["proton_right_xi_vs_t_tsel"]->Fill( -xi_proton_right,-t_proton_right, event_weight );
					if(verbose) cout<<"pfxiMinus_minus_proton_right_xi" << (pfXiMinusReco + xi_proton_right)<<endl;
					
                        if(pfXiMinusReco+xi_proton_right>0.009){
					   ++nevtxibkcgright;
					  //cout<<"(pfXiMinusReco+xi_proton_right>0.009) = "<<(pfXiMinusReco+xi_proton_right)<<endl;
                        histosTH2F["proton_right_xi_vs_t_halosel"]->Fill( -xi_proton_right,-t_proton_right, event_weight );
						histosTH1F["proton_right_t_halo"]->Fill(-t_proton_right, event_weight );
						histosTH1F["proton_right_t_halo_constbin"]->Fill(-t_proton_right, event_weight );
						histosTH1F["proton_right_xi_halo"]->Fill(-xi_proton_right, event_weight );
	     		        histosTH1F["proton_right_xi_halosel"]->Fill( -xi_proton_right, event_weight );
						histosTH1F["proton_right_t_halosel"]->Fill( -t_proton_right, event_weight );
                        //jpsi_eta_t_cut_proton_right_halo = dimuon_eta;
                        //jpsi_rapidity_t_cut_proton_right_halo = dimuon_rapidity; 
                        //histosTH1F["jpsi_dimuon_eta_t_cut_proton_right_halosel"]->Fill( dimuon_eta, event_weight );
                        if(verbose)cout<<"good_proton_right, Jpsi halo"<<endl;
					}

					if( (pfXiMinusReco + xi_proton_right) < 0.009){
						//cout<<"(pfXiMinusReco+xi_proton_right<0.009) = "<<(pfXiMinusReco+xi_proton_right)<<endl;
						++nevtxisignalright;
						histosTH1F["proton_right_t_signal"]->Fill( -t_proton_right, event_weight );
						histosTH1F["proton_right_t_signal_constbin"]->Fill( -t_proton_right, event_weight );
						histosTH1F["proton_right_xi_signal"]->Fill(-xi_proton_right, event_weight );
						histosTH1F["proton_right_xi_sel"]->Fill( -xi_proton_right, event_weight );
						histosTH1F["proton_right_t_sel"]->Fill( -t_proton_right, event_weight );
                        histosTH2F["proton_right_xi_vs_t_sel"]->Fill( -xi_proton_right,-t_proton_right, event_weight );                                         
                        //jpsi_eta_t_cut_proton_right_signal = dimuon_eta;
                        //jpsi_rapidity_t_cut_proton_right_signal = dimuon_rapidity; 
                        //histosTH1F["jpsi_dimuon_eta_t_cut_proton_right_sel"]->Fill( dimuon_eta, event_weight );
                        if(verbose)cout<<"good_proton_right, Jpsi signal"<<endl;
					}
				}

			}//t region
			if( good_proton_left ){
				cout<<"-good proton left"<<endl;
                 cout<<"after left: "<< evtId->Run << ":"<< evtId->LumiSect << ":"<< evtId->Evt << endl; 
				histosTH1F["proton_left_chi2"]->Fill( chi2_proton_left, event_weight );
				histosTH1F["proton_left_xi"]->Fill( -xi_proton_left, event_weight );
				histosTH1F["proton_left_t"]->Fill( -t_proton_left, event_weight );
				histosTH1F["proton_left_xi_bin"]->Fill( -xi_proton_left, event_weight );
				histosTH1F["proton_left_tbins"]->Fill( -t_proton_left, event_weight );
				histosTH2F["proton_left_xi_vs_pf_xiPlus"]->Fill( pfXiPlusReco, -xi_proton_left, event_weight );
				histosTH2F["proton_left_xi_vs_pf_xiMinus"]->Fill( pfXiPlusReco, -xi_proton_left, event_weight );
				//histosTH1F["proton_left_dimuon_mass"]->Fill( dimuon_mass, event_weight );
                histosTH2F["proton_left_xi_vs_t"]->Fill( -xi_proton_left,-t_proton_left, event_weight );
                
				if (pfXiPlusReco>0.05){ 
					//cout<<"xi_totem_left_cut:"<<-xi_proton_left<<endl;
					histosTH1F["proton_left_xi_cut"]->Fill(-xi_proton_left, event_weight );
				} 
				histosTH1F["pfxiPlus_minus_proton_left_xi"]->Fill((pfXiPlusReco+xi_proton_left), event_weight );
				
				if (t_region_left) {

                   	n_jpsi_t_selected_proton_left++;
					/*histosTH1F["jpsi_dimuon_mass_t_cut_proton_left"]->Fill( dimuon_mass, event_weight );
					histosTH1F["jpsi_dimuon_pt_t_cut_proton_left"]->Fill( dimuon_pt, event_weight );
					histosTH1F["jpsi_dimuon_pt2_t_cut_proton_left"]->Fill( dimuon_pt2, event_weight );
					histosTH1F["jpsi_dimuon_eta_t_cut_proton_left"]->Fill( dimuon_eta, event_weight );
					histosTH1F["jpsi_dimuon_rapidity_t_cut_proton_left"]->Fill( dimuon_rapidity, event_weight );
					histosTH1F["muonDeltaPt_jpsi_t_cut_proton_left"]->Fill(deltapt, event_weight );
					histosTH1F["muonDeltaEta_jpsi_t_cut_proton_left"]->Fill(deltaeta, event_weight );
					histosTH1F["muonDeltaPhi_jpsi_t_cut_proton_left"]->Fill(deltaphi, event_weight );
					histosTH1F["muonDeltaY_jpsi_t_cut_proton_left"]->Fill(deltay, event_weight );
					histosTH2F["jpsi_DeltaPhi_vs_dimuon_pt_t_cut_proton_left"]->Fill(dimuon_pt,deltaphi, event_weight );
					histosTH2F["jpsi_y_vs_pt_t_cut_proton_left"]->Fill(dimuon_rapidity,dimuon_pt, event_weight );
*/
					histosTH1F["xitotem_xicms_leftRPs_kin"]->Fill( pfXiPlusReco+xi_proton_left, event_weight );
					histosTH1F["proton_left_t_kint"]->Fill(-t_proton_left, event_weight );
					histosTH1F["proton_left_xi_kint"]->Fill(-xi_proton_left, event_weight );
                    histosTH2F["proton_left_xi_vs_t_tsel"]->Fill( -xi_proton_left,-t_proton_left, event_weight );
                    histosTH2F["proton_left_xi_vs_t_kin"]->Fill( -xi_proton_left,-t_proton_left, event_weight );

					if(pfXiPlusReco+xi_proton_left>0.009){
					    ++nevtxibkcgleft;
						histosTH1F["proton_left_t_halo"]->Fill(-t_proton_left, event_weight );
						histosTH1F["proton_left_t_halo_constbin"]->Fill(-t_proton_left, event_weight );
						histosTH1F["proton_left_xi_halo"]->Fill(-xi_proton_left, event_weight );
		      	        histosTH1F["proton_left_xi_halosel"]->Fill( -xi_proton_left, event_weight );
						histosTH1F["proton_left_t_halosel"]->Fill( -t_proton_left, event_weight );
                        histosTH2F["proton_left_xi_vs_t_halosel"]->Fill( -xi_proton_left,-t_proton_left, event_weight );
                        //jpsi_eta_t_cut_proton_left_halo = dimuon_eta;
                        //jpsi_rapidity_t_cut_proton_left_halo = dimuon_rapidity; 
						//histosTH1F["jpsi_dimuon_eta_t_cut_proton_left_halosel"]->Fill( dimuon_eta, event_weight );
					}

					if( (pfXiPlusReco+xi_proton_left) < 0.009 ){
						//select_proton_plus = true;
						++nevtxisignalleft;
						histosTH1F["proton_left_t_signal"]->Fill( -t_proton_left, event_weight );
						histosTH1F["proton_left_t_signal_constbin"]->Fill( -t_proton_left, event_weight );
						histosTH1F["proton_left_xi_signal"]->Fill(-xi_proton_left, event_weight );
						histosTH1F["proton_left_xi_sel"]->Fill( -xi_proton_left, event_weight );
						histosTH1F["proton_left_t_sel"]->Fill( -t_proton_left, event_weight );
                        histosTH2F["proton_left_xi_vs_t_sel"]->Fill( -xi_proton_left,-t_proton_left, event_weight );
			            //jpsi_eta_t_cut_proton_left_signal = dimuon_eta;
                        //jpsi_rapidity_t_cut_proton_left_signal = dimuon_rapidity;
                        //histosTH1F["jpsi_dimuon_eta_t_cut_proton_left_sel"]->Fill( dimuon_eta, event_weight );
					}
				}
			} //t region
			bool proton_pair_valid = rec_proton_pair->valid;
			double chi2_proton_pair = rec_proton_pair->chi2;
			//double chindf_proton_pair = rec_proton_pair->chindf;
			double xi_proton_pair_right = rec_proton_pair->xir;
			double xi_proton_pair_left = rec_proton_pair->xil;
			//bool good_proton_pair = proton_pair_valid && (chi2_proton_pair/chindf_proton_pair > 2);
			bool good_proton_pair = proton_pair_valid; //&& (-0.23< xi_proton_right < 0.04) && (-0.23< xi_proton_left < 0.04);
			if( good_proton_pair ){
				histosTH1F["proton_pair_chi2"]->Fill( chi2_proton_pair, event_weight );
				//double xi_proton_pair_right = rec_proton_pair->xir;
				double t_proton_pair_right = rec_proton_pair->tr;
				histosTH1F["proton_pair_right_xi"]->Fill( -xi_proton_pair_right, event_weight );
				histosTH1F["proton_pair_right_t"]->Fill( -t_proton_pair_right, event_weight );
				if(-xi_proton_pair_right > 0.)
					histosTH1F["proton_pair_right_logXi"]->Fill( log10(-xi_proton_pair_right), event_weight );

				double t_proton_pair_left = rec_proton_pair->tl;
				histosTH1F["proton_pair_left_xi"]->Fill( -xi_proton_pair_left, event_weight );
				histosTH1F["proton_pair_left_t"]->Fill( -t_proton_pair_left, event_weight );
				if(-xi_proton_pair_left > 0.){
					histosTH1F["proton_pair_left_logXi"]->Fill( log10(-xi_proton_pair_left), event_weight );  }

			}
//}//jpsimass
}//lista de eventos

     
				//-------------------
				// Detector-level distributions
				//-------------------


                //histosTH2F["etaMuonsVsetaJpsi_pright_Signal"]->Fill(muons_eta_t_cut_proton_right_signal,jpsi_eta_t_cut_proton_right_signal , event_weight );
                //histosTH2F["etaMuonsVsetaJpsi_pright_Halo"]->Fill(muons_eta_t_cut_proton_right_halo, jpsi_eta_t_cut_proton_right_halo , event_weight );
                //histosTH2F["etaMuonsVsetaJpsi_pleft_Signal"]->Fill( muons_eta_t_cut_proton_left_signal,jpsi_eta_t_cut_proton_left_signal, event_weight );
                //histosTH2F["etaMuonsVsetaJpsi_pleft_Halo"]->Fill(muons_eta_t_cut_proton_left_halo,jpsi_eta_t_cut_proton_left_halo , event_weight );
				histosTH1F["pf_EPlusPz"]->Fill( pfEPlusPz, event_weight );
				histosTH1F["pf_EMinusPz"]->Fill( pfEMinusPz, event_weight );
				histosTH1F["pf_xiPlus"]->Fill( pfXiPlusReco, event_weight );
				histosTH1F["pf_xiMinus"]->Fill( pfXiMinusReco, event_weight );
				histosTH1F["pf_logXiPlus"]->Fill( log10(pfXiPlusReco), event_weight );
				histosTH1F["pf_logXiMinus"]->Fill( log10(pfXiMinusReco), event_weight );

				histosTH1F["t2_track_multiplicity_zplus"]->Fill( n_t2_tracks_selected_zplus, event_weight );
				histosTH1F["t2_track_multiplicity_zminus"]->Fill( n_t2_tracks_selected_zminus, event_weight );
				//histosTH2F["t2_track_multiplicity_vs_track_multiplicity"]->Fill( n_tracks_selected, n_t2_tracks_selected, event_weight );





				///////////////////////////////////////////////////////////////////////////////


			} // End of loop over events in a file

			cout<<"Total of evts="<< nev << endl << *itfiles << endl;

			//cout<<"n_select_Vertex_After_vtx_cut="<< n_select_Vertex_After_vtx_cut<<endl;
			//cout<< "n_evt_selectTrack"<<n_evt_selectTrack<<endl;
			//cout<<"n_muon_selected="<< n_muon_selected<<endl;
			//cout<<"n_dimuons_selected= "<< n_dimuons_selected<<endl;
			//cout<< "n_evt_passed_threshold "<<n_evt_passed_threshold<<endl;
            //cout<<"n_jpsi_selected= "<< n_jpsi_selected<<endl;
            //cout<<"n_jpsi_mass= "<< n_jpsi_mass<<endl;
            //cout<<"n_jpsi_mass_p= "<< n_jpsi_mass_p<<endl;
            //cout<<"nevtxibkcgright= "<<nevtxibkcgright<<endl;
            //cout<<"nevtxisignalright= "<<nevtxisignalright<<endl;
            //cout<<"nevtxibkcgleft= "<<nevtxibkcgleft<<endl;
            //cout<<"nevtxisignalleft= "<<nevtxisignalleft<<endl;
            //cout<<"n_jpsi_t_selected_proton_left="<< n_jpsi_t_selected_proton_left<<endl;
			//cout<<"n_jpsi_t_selected_proton_right="<< n_jpsi_t_selected_proton_right<<endl;
            //cout<<"n_dimuons_t_selected_proton_left="<< n_dimuons_t_selected_proton_left<<endl;
			//cout<<"n_dimuons_t_selected_proton_right="<< n_dimuons_t_selected_proton_right<<endl;
			/*cout<<"n_dimuons_proton_left="<<n_dimuons_proton_left<<endl;
			cout<<"n_dimuons_proton_right="<<n_dimuons_proton_right<<endl;
			cout<<"n_dimuons_t_selected_proton_left="<< n_dimuons_t_selected_proton_left<<endl;
			cout<<"n_dimuons_t_selected_proton_right="<< n_dimuons_t_selected_proton_right<<endl;
			cout<<"n_jpsi_selected="<< n_jpsi_selected<<endl;
			cout<<"n_jpsi_proton_left="<<n_jpsi_proton_left<<endl;
			cout<<"n_jpsi_proton_right="<<n_jpsi_proton_right<<endl;
			cout<<"n_jpsi_t_selected_proton_left="<< n_jpsi_t_selected_proton_left<<endl;
			cout<<"n_jpsi_t_selected_proton_right="<< n_jpsi_t_selected_proton_right<<endl;
			cout<<"nr pf="<<n_evt_PF<<endl;
			cout<<"nevtxisignalleft="<<nevtxisignalleft<<endl;
			cout<<"nevtxisignalright="<<nevtxisignalright<<endl;
			cout<<"nevtxisignalleftjpsi="<<nevtxisignalleftjpsi<<endl;
			cout<<"nevtxisignalrightjpsi="<<nevtxisignalrightjpsi<<endl;
			cout<<"nevtxibhleftjpsi="<<nevtxibhleftjpsi<<endl;
			cout<<"nevtxibhrightjpsi="<<nevtxibhrightjpsi<<endl;*/



			// Close current file
			file->Close();

		} // End of loop over files
		//==========

		//==========
		//output file
		TFile* output = new TFile(outputFileName.c_str(),"RECREATE");
		output->cd();

		for(map<string,TH1F*>::iterator it_histo = histosTH1F.begin();
				it_histo != histosTH1F.end(); ++it_histo)
			(*it_histo).second->Write();
		for(map<string,TH2F*>::iterator it_histo = histosTH2F.begin();
				it_histo != histosTH2F.end(); ++it_histo)
			(*it_histo).second->Write();


		output->Close();

		}




