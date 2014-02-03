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
#include <TRandom.h>
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

#include "EventMetaData.h"
#include "NtupleInfo.h"
#include "RPEvent.h"
#include "RPRootDumpDigiInfo.h"
#include "RPRootDumpPattern.h"
#include "RPRootDumpPatternInfo.h"
#include "RPRootDumpReconstructedProton.h"
#include "RPRootDumpReconstructedProtonPair.h"
#include "RPRootDumpTrackInfo.h"
#include "T1Event.h"
#include "T2Event.h"
#include "TriggerData.h"

#include "analysis_tools.h"
#include "rp_aperture_config.h"
#include "deltaPhi.h"
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

using namespace std;

void DimuonTriggerEff(string const& outputFileName = "data_DimuonTriggerEff_28Jan2014.root",const Double_t t_proton_down_=0.0, const Double_t t_proton_up_=1.0,const Int_t Bin_mass=200 , const Int_t nevt_max = -1){
	//void ana_psel_DataMC_Jpsi_halo(vector<string> const& fileNames, string const& outputFileName = "output.root",const Double_t t_proton_down_=0.0, const Double_t t_proton_up_=1.0,const Int_t Bin_mass=200 ,const Int_t nevt_max = -1){

	bool isMC  = false;
	bool verbose = false;
	bool Vertex = false;
	//string treeName ="evt";
	//string treeName ="cms_totem";
	string treeName = (!isMC) ? "cms_totem" : "evt";
	//double etaMaxThreshold = 2.0;
	double ptMax = 9999.0;
	bool selectBunchCrossing = false;
	bool selectVertex = true;
	bool selectMuons = true;
	bool selectTrack = true;
	//bool selectEtaMax = false;
	//bool selectEtaMin = false;
	bool selectZeroHitsT2Plus = false;
	bool selectZeroHitsT2Minus = false;
	bool selectSingleArmRecProton = true;
	//bool selectSingleArmRecProton = false;
	bool selectDoubleArmRecProton = false;
	//bool selectDoubleArmRecProton = true;
	bool selectElastic = false;
	bool selectNonElastic = false;
	//bool selectNonElastic = true;
	//bool selectRPProton = true;
	//MC
	bool selectRPPlusAccept = false;
	bool selectRPMinusAccept = true;
	bool sdplus = false;
	bool sdminus = true;

	// Declaration of histograms
	map<string,TH1F*> histosTH1F;
	map<string,TH2F*> histosTH2F;

	vector<int> bunchCrossingList;
	bunchCrossingList.push_back(27);
	bunchCrossingList.push_back(649);
	bunchCrossingList.push_back(2991);

	//vector<string> selectHLTPathNames;
	//selectHLTPathNames.push_back("HLT_L1_DoubleMu0");
	//selectHLTPathNames.push_back("HLT_ZeroBias_v7");

	ThresholdsPerRegion thresholdsPFlow;
	thresholdsPFlow[Barrel] = ThresholdsPerType(); 
	thresholdsPFlow[Endcap] = ThresholdsPerType(); 
	thresholdsPFlow[Transition] = ThresholdsPerType(); 
	thresholdsPFlow[Endcap] = ThresholdsPerType(); 
	resetPFThresholds(thresholdsPFlow[Barrel]);
	resetPFThresholds(thresholdsPFlow[Endcap]);
	resetPFThresholds(thresholdsPFlow[Transition]);
	resetPFThresholds(thresholdsPFlow[Forward]);

	thresholdsPFlow[Barrel][MyPFCand::h0]            = make_pair(-1.,1.4);
	thresholdsPFlow[Barrel][MyPFCand::gamma]         = make_pair(-1.,0.9);
	thresholdsPFlow[Endcap][MyPFCand::h0]            = make_pair(-1.,2.7);
	thresholdsPFlow[Endcap][MyPFCand::gamma]         = make_pair(-1.,2.5);
	thresholdsPFlow[Transition][MyPFCand::h0]        = make_pair(-1.,3.8);
	thresholdsPFlow[Transition][MyPFCand::gamma]     = make_pair(-1.,2.5);
	thresholdsPFlow[Transition][MyPFCand::h_HF]      = make_pair(-1.,4.0);
	thresholdsPFlow[Transition][MyPFCand::egamma_HF] = make_pair(-1.,3.5);
	thresholdsPFlow[Forward][MyPFCand::h_HF]         = make_pair(-1.,4.0);
	thresholdsPFlow[Forward][MyPFCand::egamma_HF]    = make_pair(-1.,3.5);

	ThresholdsPerType::const_iterator pfThreshold = thresholdsPFlow[Barrel].begin();
	ThresholdsPerType::const_iterator pfThresholds_end = thresholdsPFlow[Barrel].end(); 
	ostringstream oss;
	oss << "Using the following PF thresholds:\n";
	for(; pfThreshold != pfThresholds_end; ++pfThreshold){
		int key = pfThreshold->first;    
		oss << "  " << key << ": "
			<< "(" << thresholdsPFlow[Barrel][key].first
			<< "," << thresholdsPFlow[Barrel][key].second << ")  "
			<< "(" << thresholdsPFlow[Endcap][key].first
			<< "," << thresholdsPFlow[Endcap][key].second << ")  "
			<< "(" << thresholdsPFlow[Transition][key].first
			<< "," << thresholdsPFlow[Transition][key].second << ")  "
			<< "(" << thresholdsPFlow[Forward][key].first
			<< "," << thresholdsPFlow[Forward][key].second << ")\n";   
	}
	cout << oss.str();
	//==============================
	const Int_t nevt_max_corr = (nevt_max >= 0) ? nevt_max : 99999999;
	//if(verbose)cout<<"pass MC 1"<<endl;
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
	selections.push_back("RPProton");
	selections.push_back("MCRPPlusAccept");
	selections.push_back("MCRPMinusAccept");
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
		histosTH1F["hltTrigFired"]->GetXaxis()->SetBinLabel( (k + 1), hltPathNames[k].c_str() );

	histosTH1F["vtx_zpos"] = new TH1F("vtx_zpos", "z(vtx)" , 150 , -30. , 30.);
	histosTH1F["vtx_xpos"] = new TH1F("vtx_xpos", "x(vtx)" , 150 , -1.5 , 1.5);
	histosTH1F["vtx_xpos_random"] = new TH1F("vtx_xpos_random", "x(vtx)" , 150 , -1.5 , 1.5);
	histosTH1F["vtx_ypos"] = new TH1F("vtx_ypos", "y(vtx)" , 150 , -1.5 , 1.5);
	histosTH1F["vtx_ypos_random"] = new TH1F("vtx_ypos_random", "y(vtx)" , 150 , -1.5 , 1.5);
	histosTH1F["vtx_ndof"] = new TH1F("vtx_ndof", "ndof(vtx)" , 100 , 0. , 15.);
	histosTH1F["vtx_chi2"] = new TH1F("vtx_chi2", "chi2(vtx)" , 100 , 0. , 10.);
	histosTH1F["vtx_xpos_after"] = new TH1F("vtx_xpos_after", "x(vtx)" , 150 , -1.5 , 1.5);
	histosTH1F["vtx_ypos_after"] = new TH1F("vtx_ypos_after", "y(vtx)" , 150 , -1.5 , 1.5);

	histosTH1F["vertex_multiplicity"] = new TH1F("vertex_multiplicity", "n vertices" , 30 , 0 , 30);
	histosTH1F["vertex_multiplicity_after_vtx_sel"] = new TH1F("vertex_multiplicity_after_vtx_sel", "n vertices after vtx sel" , 30 , 0 , 30);
	histosTH1F["prim_vtx_zpos"] = new TH1F("prim_vtx_zpos", "z(vtx)" , 150 , -30. , 30.);
	histosTH1F["prim_vtx_xpos"] = new TH1F("prim_vtx_xpos", "x(vtx)" , 150 , -1.5 , 1.5);
	histosTH1F["prim_vtx_ypos"] = new TH1F("prim_vtx_ypos", "y(vtx)" , 150 , -1.5 , 1.5);
	histosTH1F["prim_vtx_xpos_random"] = new TH1F("prim_vtx_xpos_random", "x(vtx)" , 150 , -1.5 , 1.5);
	histosTH1F["prim_vtx_ypos_random"] = new TH1F("prim_vtx_ypos_random", "y(vtx)" , 150 , -1.5 , 1.5);
	histosTH1F["prim_vtx_zpos_after"] = new TH1F("prim_vtx_zpos_after", "z(vtx)" , 150 , -30. , 30.);
	histosTH1F["prim_vtx_xpos_after"] = new TH1F("prim_vtx_xpos_after", "x(vtx)" , 150 , -1.5 , 1.5);

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
	histosTH1F["track_pt"] = new TH1F("track_pt", "p_{T}(trk)" , 100 , 0. , 15.);
	histosTH1F["track_eta"] = new TH1F("track_eta", "#eta(trk)" , 100 , -5.2 , 5.2);
	histosTH1F["track_rapidity"] = new TH1F("track_rapidity", "#rapidity(trk)" , 100 , -5.2 , 5.2);
	histosTH1F["track_phi"] = new TH1F("track_phi", "#phi(trk)" , 100 , -M_PI , M_PI);
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
	//Muon1 info

	histosTH1F["muon2_pt"] = new TH1F("muon2_pt", "p_{T}(mu2)" , 100 , 0. , 100.);
	histosTH1F["muon2_eta"] = new TH1F("muon2_eta", "#eta(mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["muon2_phi"] = new TH1F("muon2_phi", "#phi(muon2)" , 100 , -1.2*M_PI ,1.2*M_PI);
	histosTH1F["muon2_rapidity"] = new TH1F("muon2_rapidity", "y(mu2)" , 100 , -15. , 15.);


	histosTH1F["leadingMuon_selected_pt15"] = new TH1F("leadingMuon_selected_pt15", "p_{T}(mu1)" , 100 , 0. , 100.);
	histosTH1F["leadingMuon_selected_pt20"] = new TH1F("leadingMuon_selected_pt20", "p_{T}(mu1)" , 100 , 0. , 100.);
	histosTH1F["secondMuon_leadingMuon15"] = new TH1F("secondMuon_leadingMuon15", "p_{T}(mu2)" , 100 , 0. , 100.);
	histosTH1F["secondMuon_leadingMuon20"] = new TH1F("secondMuon_leadingMuon20", "p_{T}(mu2)" , 100 , 0. , 100.);
	histosTH1F["secondMuon_selected_Muon1pt15"] = new TH1F("leadingMuon_selected_Muon1pt15", "p_{T}(mu2)" , 100 , 0. , 100.);
	histosTH1F["secondMuon_selected_Muon1pt20"] = new TH1F("leadingMuon_selected_Muon1pt20", "p_{T}(mu2)" , 100 , 0. , 100.);
	//Deltas info
	histosTH1F["muonDeltaPt"] = new TH1F("muonDeltaPt", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta"] = new TH1F("muonDeltaEta", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);	
	histosTH1F["muonDeltaPhi"] = new TH1F("muonDeltaPhi", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);
	histosTH1F["muonDeltaY"] = new TH1F("muonDeltaY", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["lorentzdphi"] = new TH1F("lorentzDPhi", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);    
	//Dimuon
	histosTH1F["dimuon_mass"] = new TH1F("dimuon_mass", "mass(mu1,mu2)" , Bin_mass , 0. , 10.);
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
	//Dimuon
	histosTH1F["lorentzdphi_jpsi"] = new TH1F("lorentzDPhi_jpsi", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);

	//
	histosTH2F["DeltaPhi_vs_dimuon_pt"] = new TH2F("DeltaPhi_vs_dimuon_pt", "#Delta#phi_vs_dimuon_p_{T}" , 100, 0.0, 20.0, 100, 0.0,1.2*M_PI );
	histosTH2F["jpsi_DeltaPhi_vs_dimuon_pt"] = new TH2F("jpsi_DeltaPhi_vs_dimuon_pt", "#Delta#phi_vs_dimuon_p_{T}" , 100, 0.0, 20.0, 100, 0.0,1.2*M_PI );
	//
	histosTH1F["pf_etaMax"] = new TH1F("pf_etaMax","#eta^{max}",82,etaBinsHCALBoundaries);
	histosTH1F["pf_etaMin"] = new TH1F("pf_etaMin","#eta^{min}",82,etaBinsHCALBoundaries);
	histosTH1F["pf_deltaEta"] = new TH1F("pf_deltaEta","#Delta#eta",100,0.,10.);
	histosTH1F["pf_EPlusPz"] = new TH1F("pf_EPlusPz","sum(E + pz)",24,binningEPlusPz);
	histosTH1F["pf_EMinusPz"] = new TH1F("pf_EMinusPz","sum(E - pz)",24,binningEPlusPz);
	//histosTH1F["pf_EPlusPz"] = new TH1F("pf_EPlusPz","sum(E + pz)",200 , -1.,1.);
	//histosTH1F["pf_EMinusPz"] = new TH1F("pf_EMinusPz","sum(E - pz)",200 , -1.,1.);
	//histosTH1F["random_xi_cms_plus"] = new TH1F("random_xi_cms_plus", "#xi^{+} CMS" , 250 , -1.,1.);
	//histosTH1F["random_xi_cms_minus"] = new TH1F("random_xi_cms_minus", "#xi^{-} CMS" , 250 , -1.,1.);
	//histosTH1F["random_xi_totem_right"] = new TH1F("random_xi_totem_right","#xi TOTEM",250,-1.,1.);
	//histosTH1F["random_xi_totem_left"] = new TH1F("random_xi_totem_left","#xi TOTEM",250,-1.,1.);
	histosTH1F["pf_xiPlus"] = new TH1F("pf_xiPlus","#xi^{+}",200,-1.,1.);
	histosTH1F["pf_xiMinus"] = new TH1F("pf_xiMinus","#xi^{-}",200,-1.,1.);
	histosTH1F["xi_cms_pfplus"] =  new TH1F("xi_cms_pfplus","#xi^{+}",200,-1.,1.);
	histosTH1F["xi_cms_pfminus"] =  new TH1F("xi_cms_pfminus","#xi^{-}",200,-1.,1.);
	histosTH1F["jpsi_xi_cms_pfplus"] =  new TH1F("jpsi_xi_cms_pfplus","#xi^{+}",200,-1.,1.);
	histosTH1F["jpsi_xi_cms_pfminus"] =  new TH1F("jpsi_xi_cms_pfminus","#xi^{-}",200,-1.,1.);
	//histosTH1F["xi_cms_pfplus"] =  new TH1F("xi_cms_pfplus","#xi^{+}",20,0.,1.);
	//histosTH1F["xi_cms_pfminus"] =  new TH1F("xi_cms_pfminus","#xi^{-}",20,0.,1.);
	histosTH1F["pf_logXiPlus"] = new TH1F("pf_logXiPlus","log(#xi^{+})",20,-4.,0.);
	histosTH1F["pf_logXiMinus"] = new TH1F("pf_logXiMinus","log(#xi^{-})",20,-4.,0.);


	histosTH1F["fscHit_energy"] = new TH1F("fscHit_energy", "FSC hit energy" , 150 , -100. , 200.);
	histosTH1F["fscHit_time"] = new TH1F("fscHit_time", "FSC hit time" , 150 , 0. , 300.);


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
	//===========MC

	//============

	for(map<string,TH1F*>::const_iterator it = histosTH1F.begin(); it != histosTH1F.end(); ++it)
		it->second->Sumw2();
	for(map<string,TH2F*>::const_iterator it = histosTH2F.begin(); it != histosTH2F.end(); ++it)
		it->second->Sumw2();
	//===================
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
	//TTree* tree = new TTree(treeName.c_str(),"");
	TTree* tree = NULL;
	MyEvtId*           evtId        = NULL;
	MyL1TrigOld*       l1Trig       = NULL;  
	MyHLTrig*          hltTrig      = NULL;
	vector<MyGenPart>* genPart      = NULL;
	vector<MyTracks>*  track_coll   = NULL;
	vector<MyVertex>*  vertex_coll  = NULL;
	//vector<MyPFJet>*   pfJet_coll   = NULL;
	vector<MyMuon>*    muon_coll   = NULL;
	vector<MyPFCand>*  pFlow_coll   = NULL;
	//vector<MyFSCHit>*  fscHits_coll = NULL;
	//vector<MyFSCDigi>* fscDigis_coll = NULL;
	MyGenKin*          genKin        = NULL;
	//===================
	T1Event* t1_event = NULL;
	T2Event* t2_event = NULL; 
	RPRootDumpReconstructedProton* rec_proton_left  = NULL;
	RPRootDumpReconstructedProton* rec_proton_right = NULL;
	RPRootDumpReconstructedProtonPair* rec_proton_pair  = NULL;
	//RPRootDumpReconstructedProton* sim_proton_right = NULL;
	//RPRootDumpReconstructedProton* sim_proton_left = NULL;
	map<unsigned int, RPRootDumpTrackInfo*> rp_track_info;
	map<unsigned int, RPRootDumpDigiInfo*> rp_digi_info;
	map<unsigned int, RPRootDumpPatternInfo*> rp_par_patterns_info;
	map<unsigned int, RPRootDumpPatternInfo*> rp_nonpar_patterns_info;
	map<unsigned int, std::vector<RPRootDumpTrackInfo>*> rp_multi_track_info;
	//===================
	//if(!isMC){
	/*TString outtxt_right = outputFileName;
	  TString outtxt_left = outputFileName;

	  outtxt_right.ReplaceAll("root","txt_right");
	  outtxt_left.ReplaceAll("root","txt_left");  

	//const string left_ = "left_"+ outtxt_left ;
	//std::string right_ = "right_"+ outtxt_right ;
	ofstream outstring_left(outtxt_left); 
	ofstream outstring_right(outtxt_right);*/
	//}
	//===================
	//===================
	if(isMC) rp_aperture_config();

	//if(verbose)cout<<"rp_aperture_config"<<endl;
	int n_vertices_selected = 0;
	int n_select_Vertex_After_vtx_cut =0;
	int n_tracks_selected = 0;
	int n_muon_selected = 0;
	int n_dimuons_selected = 0;

	int n_jpsi_selected = 0;


	//int i_tot = 0 , nevt_tot = 0;
	if(verbose)cout<<"i_tot"<<endl;
	//starting Loop over files, stops at end of list of files or when reached nevt_max
	for(vector<TString>::iterator itfiles = vfiles->begin(); itfiles != vfiles->end() && i_tot < nevt_max_corr; ++itfiles){

		TFile* file = TFile::Open(*itfiles,"READ");

		//getting the tree form the current file
		tree = (TTree*) file->Get( treeName.c_str() );

		//Getting number of events
		int nev = int(tree->GetEntriesFast());

		nevt_tot += nev;
		cout <<"The current file has " << nev << " entries : "<< endl << *itfiles << endl;

		// Add branches to TTree ----------------------------------------------------------------------
		if(verbose)cout<<"branches"<<endl;

		if(isMC){
			tree->SetBranchAddress("evtId",&evtId);
			tree->SetBranchAddress("generalTracks",&track_coll); 
			tree->SetBranchAddress("offlinePrimaryVertices",&vertex_coll);
			tree->SetBranchAddress("muons",&muon_coll);
			tree->SetBranchAddress("particleFlow",&pFlow_coll);
			tree->SetBranchAddress("genKin",&genKin);
			tree->SetBranchAddress("genPart",&genPart);
			if(verbose)cout<<"add branches to ttree"<<endl;

		}else{ 
			tree->SetBranchAddress("cmsEvtUA",&evtId);
			tree->SetBranchAddress("cmsTrigUA",&l1Trig);
			tree->SetBranchAddress("cmsHLTTrigUA",&hltTrig);
			tree->SetBranchAddress("cmsTracksUA",&track_coll); 
			tree->SetBranchAddress("cmsVerticesUA",&vertex_coll);
			//tree->SetBranchAddress("cmsPFJetsUA",&pfJet_coll);
			//tree->SetBranchAddress("cmsak5PFJetsUA",&pfJet_coll);
			tree->SetBranchAddress("cmsMuonsUA",&muon_coll);
			tree->SetBranchAddress("cmsParticleFlowUA",&pFlow_coll);
			tree->SetBranchAddress("branchT1EV.",&t1_event);
			tree->SetBranchAddress("branchT2EV.",&t2_event);
			tree->SetBranchAddress("rec_prot_left.",&rec_proton_left);
			tree->SetBranchAddress("rec_prot_right.",&rec_proton_right);
			tree->SetBranchAddress("rec_prot_pair.",&rec_proton_pair);
			//tree->SetBranchAddress("sim_prot_left.",&sim_proton_left);
			//tree->SetBranchAddress("sim_prot_right.",&sim_proton_right);//}
			std::vector<unsigned int> rp_list;
			rp_list.push_back(20); rp_list.push_back(21); rp_list.push_back(22); rp_list.push_back(23); rp_list.push_back(24); rp_list.push_back(25);
			rp_list.push_back(120); rp_list.push_back(121); rp_list.push_back(122); rp_list.push_back(123); rp_list.push_back(124); rp_list.push_back(125);
			char br_name[200];
			//char name[200];
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
	}


	if(verbose)cout<<"pass MC 2"<<endl;




	//starting loop over events, stops when reached end of file or nevt_max
	for(int i_evt = 0; i_evt < nev && i_tot < nevt_max_corr; ++i_evt , ++i_tot){

		//printing the % of events done every 10k evts
		if( ((i_tot+1) % 5000) == 0) cout <<int(double(i_tot+1)/1000)<<"k done"<<endl;
		//if( ((i_tot+1) % 100) == 0 ) cout << (i_tot+1) << " done" << endl;

		//Filling the variables defined setting branches
		tree->GetEntry(i_evt);

		//continue;

		//AT THIS POINT ON, CAN START USING THE VARIABLES LIKE TRACKS, VERTEX ETC !
		double event_weight = 1.;
		bool passed_HLTMuon= false;
		//string HLT_muon = "HLT_L1DoubleMu0_v1"; 
		bool passed_leadingMuon_pt15 = false;
		bool passed_leadingMuon_pt20 = false;
		//string HLT_muon = "HLT_ZeroBias_v7";
		bool passedHLT_Reference = false;
		//bool passed_vrt= false;      
		//bool passed_muon = false;
		bool mu1_selected = false;
		bool mu2_selected = false;
		//bool chmuons = false;
		bool select_Track = false;
		histosTH1F["EventSelection"]->Fill( "All", event_weight );
		//========================================================

		bool proton_plus_rp_accept = false;
		bool proton_minus_rp_accept = false;
		bool proton_plus_xi_range = false;
		bool proton_plus_t_range = false;
		bool proton_minus_xi_range = false;
		bool proton_minus_t_range = false;
		double xi_proton_plus = -1.;
		double xi_proton_minus = -1.;
		double t_proton_plus = 0.;
		double t_proton_minus = 0.;
		float xi_proton_plus_rec = -1.;
		float xi_proton_minus_rec = -1.;
		double t_proton_plus_rec = 0.;
		double t_proton_minus_rec = 0.;

		// Event selection
		//-------------------------------------------------------------------------------------------------
		// Bunch crossing & HLT
		if(!isMC){
			//==========================================================
			//int orbitNumber = evtId->Orbit;
			//int orbitNumber = evtId->Orbit;
			if(verbose)cout<<"1347- bunchCrossingNumber"<<endl;
			int bunchCrossingNumber = evtId->Bunch;
			histosTH1F["bunchCrossingNumber"]->Fill( bunchCrossingNumber, event_weight );

			if( selectBunchCrossing &&
					find(bunchCrossingList.begin(),bunchCrossingList.end(),bunchCrossingNumber) == bunchCrossingList.end() )
				continue;

			histosTH1F["EventSelection"]->Fill( "BunchCrossing", event_weight );

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

			/*for(; it_hlt != it_hlt_end; ++it_hlt){
			  string const& hltName = it_hlt->first;
			  vector<string>::const_iterator it_pos = find(hltPathNames.begin(),hltPathNames.end(),hltName);
			  if(it_pos != hltPathNames.end()){
			  size_t idx = it_pos - hltPathNames.begin();//cout <<hltName<<endl;
			// if( it_hlt->second ) histosTH1F["hltTrigFired"]->Fill( hltName.c_str(), event_weight );
			if( hltName == "HLT_L1Tech53_MB_1_v1" || hltName == "HLT_L1Tech53_MB_2_v1" || hltName == "HLT_L1Tech53_MB_3_v1"){
			if( it_hlt->second == true ){
			passedHLT_Reference = true;
			histosTH1F["hltTrigFired"]->Fill( idx, event_weight );
			}
			}
			if( hltName == HLT_muon){ 
			passed_HLTMuon = true;

			if( it_hlt->second ) histosTH1F["hltTrigFired"]->Fill( hltName.c_str(), event_weight );
			}
			}


			}

			//if(!passed_HLTMuon) continue;
			if(!passedHLT_Reference) continue;
			histosTH1F["EventSelection"]->Fill( "HLT", event_weight );
			}*/
			for(; it_hlt != it_hlt_end; ++it_hlt){
				string const& hltName = it_hlt->first;
				vector<string>::const_iterator it_pos = find(hltPathNames.begin(),hltPathNames.end(),hltName);
				if(it_pos != hltPathNames.end()){
					size_t idx = it_pos - hltPathNames.begin(); //cout <<hltName<<endl;
					// if( hltName == "HLT_ZeroBias_v7" || hltName == "HLT_L1Tech54_ZeroBias_v1");{
					if( hltName == "HLT_L1Tech53_MB_1_v1" || hltName == "HLT_L1Tech53_MB_2_v1" || hltName == "HLT_L1Tech53_MB_3_v1"){
						if( it_hlt->second == true ){
							passedHLT_Reference = true;
							histosTH1F["hltTrigFired"]->Fill( idx, event_weight );
						}
					}

					if( hltName == "HLT_L1DoubleMu0_v1"){
						if( it_hlt->second == true ){
							passed_HLTMuon = true;
							histosTH1F["hltTrigFired"]->Fill( idx, event_weight );
						}
					}
				}
				}
				if(!passedHLT_Reference) continue;
			}
			//-------------------------------------------------------------------------------------------------
			//-------------------------------------------------------------------------------------------------
			/*int idx_vtx_max_sumpt = -1;
			  double sumpt_max = 0.;*/
			double xpos; double ypos; 
			double zpos;
			if(verbose)cout<<"MyVertex"<<endl;
			for(vector<MyVertex>::iterator it_vtx = vertex_coll->begin() ; it_vtx != vertex_coll->end() ; ++it_vtx){
				//int idx_vtx = it_vtx - vertex_coll->begin();
				//if( it_vtx->SumPtTracks > sumpt_max ){ idx_vtx_max_sumpt = idx_vtx; sumpt_max = it_vtx->SumPtTracks; }

				if( it_vtx->fake ) continue;
				if( !it_vtx->validity ) continue;
				++n_vertices_selected;
				zpos = it_vtx->z;
				ypos = it_vtx->y;
				xpos = it_vtx->x;
				histosTH1F["vtx_zpos"]->Fill( it_vtx->z, event_weight );
				histosTH1F["vtx_xpos"]->Fill( it_vtx->x, event_weight );
				histosTH1F["vtx_ypos"]->Fill( it_vtx->y, event_weight );
				histosTH1F["vtx_ndof"]->Fill( it_vtx->ndof, event_weight );
				histosTH1F["vtx_chi2"]->Fill( it_vtx->chi2, event_weight );
			}
			//histosTH1F["vtx_sumpt_max"]->Fill( idx_vtx_max_sumpt, event_weight );
			histosTH1F["vertex_multiplicity"]->Fill( n_vertices_selected, event_weight );

			//vertex-simulated-xy from data fit
			Double_t mean_pxpos = 0.07280;
			Double_t sigma_pxpos = 0.0115;
			Double_t mean_pypos = 0.0687554;
			Double_t sigma_pypos = 0.0107789;
			Double_t mean_xpos = 0.0769463;
			Double_t sigma_xpos = 0.0187961;
			Double_t mean_ypos = 0.0680937;
			Double_t sigma_ypos = 0.0157877;
			Double_t xpos_random = gRandom->Gaus(mean_xpos, sigma_xpos)/100;
			Double_t ypos_random = gRandom->Gaus(mean_ypos, sigma_ypos)/100;
			Double_t pxpos_random = gRandom->Gaus(mean_pxpos, sigma_pxpos)/100;
			Double_t pypos_random = gRandom->Gaus(mean_pypos, sigma_pypos)/100;
			if(isMC){histosTH1F["vtx_xpos_random"]->Fill(xpos_random*100, event_weight );
				histosTH1F["vtx_ypos_random"]->Fill( ypos_random*100, event_weight );
			}
			//MyVertex const& primaryVertex = vertex_coll->at(0);
			MyVertex& primaryVertex = vertex_coll->at(0);
			histosTH1F["prim_vtx_zpos"]->Fill( primaryVertex.z, event_weight );
			histosTH1F["prim_vtx_xpos"]->Fill( primaryVertex.x, event_weight );
			histosTH1F["prim_vtx_ypos"]->Fill( primaryVertex.y, event_weight );

			if(isMC){
				histosTH1F["prim_vtx_xpos_random"]->Fill( pxpos_random*100, event_weight );
				histosTH1F["prim_vtx_ypos_random"]->Fill( pypos_random*100, event_weight );}

				histosTH1F["prim_vtx_ndof"]->Fill( primaryVertex.ndof, event_weight );
				histosTH1F["prim_vtx_chi2"]->Fill( primaryVertex.chi2, event_weight );
				histosTH1F["prim_vtx_chi2n"]->Fill( primaryVertex.chi2n(), event_weight );
				histosTH1F["prim_vtx_ntracks"]->Fill( primaryVertex.ntracks, event_weight );
				histosTH1F["prim_vtx_sumpt"]->Fill( primaryVertex.SumPtTracks, event_weight );

				double prim_vtx_r = sqrt( primaryVertex.x*primaryVertex.x + primaryVertex.y*primaryVertex.y );
				bool select_Vertex = ( !primaryVertex.fake && primaryVertex.validity &&
						primaryVertex.ndof > 4 && fabs( primaryVertex.z ) < 15.0 && prim_vtx_r < 2.0);
				// Include Muon selection REf CMS AN AN 12-067 
				//bool select_Vertex = ( !primaryVertex.fake && primaryVertex.validity &&
				//                        primaryVertex.ndof < 10 && fabs( primaryVertex.z ) < 15.0 && prim_vtx_r < 2.0); 

				if(selectVertex && !select_Vertex) continue;

				++n_select_Vertex_After_vtx_cut;


				histosTH1F["vertex_multiplicity_after_vtx_sel"]->Fill( n_select_Vertex_After_vtx_cut, event_weight );  
				histosTH1F["prim_vtx_zpos_after_vtx_sel"]->Fill( primaryVertex.z, event_weight );
				histosTH1F["prim_vtx_xpos_after_vtx_sel"]->Fill( primaryVertex.x, event_weight );
				histosTH1F["prim_vtx_ypos_after_vtx_sel"]->Fill( primaryVertex.y, event_weight );

				histosTH1F["prim_vtx_ndof_after_vtx_sel"]->Fill( primaryVertex.ndof, event_weight );
				histosTH1F["prim_vtx_chi2_after_vtx_sel"]->Fill( primaryVertex.chi2, event_weight );
				histosTH1F["prim_vtx_chi2n_after_vtx_sel"]->Fill( primaryVertex.chi2n(), event_weight );
				histosTH1F["prim_vtx_ntracks_after_vtx_sel"]->Fill( primaryVertex.ntracks, event_weight );
				histosTH1F["prim_vtx_sumpt_after_vtx_sel"]->Fill( primaryVertex.SumPtTracks, event_weight );

				histosTH1F["EventSelection"]->Fill( "Vertex", event_weight );

				int prim_vtx_id = primaryVertex.id;
				if(verbose)cout<<"pass MyVertex"<<endl;
				// Tracks
				//int n_tracks_selected = 0;
				for(vector<MyTracks>::iterator it_trk = track_coll->begin() ; it_trk != track_coll->end() ; ++it_trk){
					histosTH1F["track_pt"]->Fill( it_trk->Pt(), event_weight );
					histosTH1F["track_eta"]->Fill( it_trk->Eta(), event_weight );
					histosTH1F["track_rapidity"]->Fill( it_trk->Rapidity(), event_weight );
					histosTH1F["track_phi"]->Fill( it_trk->Phi(), event_weight );

					if( it_trk->Pt() < 0.5 ) continue;
					if( fabs( it_trk->Eta() ) > 2.5 ) continue;
					if( ( it_trk->dz / it_trk->edz ) > 5. ) continue;
					if( ( it_trk->d0 / it_trk->ed0 ) > 5. ) continue;

					/*
					   outtrack.quality[0] = intrack.quality(TrackBase::qualityByName("loose"));
					   outtrack.quality[1] = intrack.quality(TrackBase::qualityByName("tight"));
					   outtrack.quality[2] = intrack.quality(TrackBase::qualityByName("highPurity"));
					   outtrack.quality[3] = intrack.quality(TrackBase::qualityByName("confirmed"));
					   outtrack.quality[4] = intrack.quality(TrackBase::qualityByName("goodIterative"));
					 */ 
					if( !it_trk->quality[2] ) continue;

					select_Track = true;
					++n_tracks_selected;

				}
				histosTH1F["track_multiplicity"]->Fill( n_tracks_selected, event_weight );
				if(selectTrack && !select_Track) continue;
				if(verbose)cout<<"pass MyTrack"<<endl;



				if(verbose)cout<<"starting Particle-flow"<<endl;
				//-----------------------------------------------------------------------------------------

				//-------------------------------------------------------------------------------------------------

				// Muons variables
				//int n_muon_selected = 0;
				double chmu1 = 0.; double chmu2 = 0.;
				double phimu1 = 0.; double ptmu1 = 0.;double etamu1 = 0.;double ymu1 = 0.;
				double phimu2 = 0.; double ptmu2 = 0.;double etamu2 = 0.;double ymu2 = 0.;
				double deltaphi = 0.; double deltaeta = 0.; double deltapt = 0.; double deltay = 0.; double Dphi = 0.;
				double dimuon_mass = 0.; double dimuon_pt=0.; double dimuon_pt2=0.; double dimuon_eta =0.;
				double dimuon_rapidity = 0.; double jpsi_mass = 0.; double jpsi_pt=0.;
				double jpsi_pt2=0.; double jpsi_eta =0.; double jpsi_rapidity = 0.; double dphijpsi = 0.;
				// Muons
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

					++n_muon_selected;

					histosTH1F["muon_pt"]->Fill( it_muon->Pt(), event_weight );
					histosTH1F["muon_eta"]->Fill( it_muon->Eta(), event_weight );
					histosTH1F["muon_phi"]->Fill( it_muon->Phi(), event_weight );

					muons_selected.push_back( *it_muon );

					histosTH1F["muon_multiplicity"]->Fill( n_muon_selected, event_weight );
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
					if(etamu1< 2.45){
						if(ptmu1 > 1.5) {
							if(verbose)cout<<"ptmu1 > 1.5 :"<<ptmu1<<endl;
							passed_leadingMuon_pt15 = true;
							histosTH1F["leadingMuon_selected_pt15"]->Fill(it_mu1->Pt(), event_weight);
						}
						if(ptmu1 > 2.0) {
							if(verbose)cout<<"ptmu1 > 2.0 :"<<ptmu1<<endl;
							passed_leadingMuon_pt20 = true;
							histosTH1F["leadingMuon_selected_pt20"]->Fill(it_mu1->Pt(), event_weight);
						}

					}
					for(vector<MyMuon>::iterator it_mu2 = muons_selected.begin() ;
							it_mu2 != muons_selected.end() ; ++it_mu2){
						bool os_muons = ( it_mu1->charge*it_mu2->charge < 0. );
						if( !os_muons ) continue;
						++n_dimuons_selected;
						mu2_selected = true;
						if(verbose)cout<<"mu2_selected :"<<mu2_selected<<endl;
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
						//dphijpsi = dimuon_lorentz.Phi();
						//cout<<"Delta Definitions :"<<endl;		
						//deltapt = fabs(muon1_lorentz.Pt() - muon2_lorentz.Pt());
						//deltaeta = fabs(muon1_lorentz.Eta() - muon2_lorentz.Eta());
						//deltaphi = fabs(phimu1 - phimu2);
						//if(deltaphi > M_PI)deltaphi = (2*M_PI - deltaphi);
						//deltay = fabs(muon1_lorentz.Rapidity() - muon2_lorentz.Rapidity());
						//Dphi = std::fabs(deltaPhi(phimu1 ,phimu2));  
						//cout<<"Delta Phi :"<<deltaphi<<endl;
						//histosTH1F["muonDeltaPt"]->Fill(deltapt, event_weight );
						//histosTH1F["muonDeltaEta"]->Fill(deltaeta, event_weight );	
						//histosTH1F["muonDeltaPhi"]->Fill(deltaphi, event_weight );
						//histosTH1F["muonDeltaY"]->Fill(deltay, event_weight ); 
						//histosTH2F["DeltaPhi_vs_dimuon_pt"]->Fill(dimuon_pt,deltaphi, event_weight ); 
						//histosTH1F["muonDphi"]->Fill(Dphi, event_weight );
						//histosTH1F["lorentzdphi"]->Fill(dphijpsi, event_weight );
						//cout<<"mass_dimuon = "<< dimuon_lorentz.M() << endl;
						////////////////////////////////////////////////////////////////////////////////////////////
						if( etamu2<2.45){
							if(passed_leadingMuon_pt15==true) {
								//passed_leadingMuon_pt15 = true;
								if(verbose)cout<<"ptmu1 > 1.5 "<<endl;
								histosTH1F["secondMuon_leadingMuon15"]->Fill(it_mu2->Pt(), event_weight);
								if(passed_HLTMuon == true){
									if(verbose)cout<<"passed_HLTMuon; ptmu1 > 1.5 ::"<<passed_HLTMuon <<endl;
									histosTH1F["secondMuon_selected_Muon1pt15"]->Fill(it_mu2->Pt(), event_weight);
								}
							}
							if(passed_leadingMuon_pt20==true) {
								histosTH1F["secondMuon_leadingMuon20"]->Fill(it_mu2->Pt(), event_weight);
								//passed_leadingMuon_pt20 = true;
								if(verbose)cout<<"ptmu1 > 2.0 "<<endl;
								if(passed_HLTMuon == true){
									if(verbose)cout<<"passed_HLTMuon; ptmu1 > 2.0 ::"<<passed_HLTMuon <<endl;
									histosTH1F["secondMuon_selected_Muon1pt20"]->Fill(it_mu2->Pt(), event_weight);
								}
							}
						}
						//double Jpsi_mass = dimuon_lorentz.M();
						//cout<<"jpsi_mass_dimuon = "<< Jpsi_mass << endl;                              
						/*if(((dimuon_mass > 3.0) && (dimuon_mass < 3.2))){
						//cout<<"jpsi_mass = "<< dimuon_mass << endl;
						++n_jpsi_selected;
						//double Dphi_jpsi = std::fabs(deltaPhi(phimu1 ,phimu2));
						//histosTH1F["jpsi_multiplicity"]->Fill(n_dimuons_selected, event_weight );
						histosTH1F["jpsi_mass"]->Fill( dimuon_lorentz.M(), event_weight );
						histosTH1F["jpsi_pt"]->Fill( dimuon_lorentz.Pt(), event_weight );
						histosTH1F["jpsi_pt2"]->Fill( (dimuon_lorentz.Pt()*dimuon_lorentz.Pt()), event_weight );
						histosTH1F["jpsi_eta"]->Fill( dimuon_lorentz.Eta(), event_weight );
						histosTH1F["jpsi_rapidity"]->Fill( dimuon_lorentz.Rapidity(), event_weight );
						histosTH1F["muonDeltaPt_jpsi"]->Fill(deltapt, event_weight );
						histosTH1F["muonDeltaEta_jpsi"]->Fill(deltaeta, event_weight );	
						histosTH1F["muonDeltaPhi_jpsi"]->Fill(deltaphi, event_weight );
						//histosTH1F["muonDphi_jpsi"]->Fill(Dphi, event_weight );
						histosTH1F["lorentzdphi_jpsi"]->Fill(dphijpsi, event_weight );
						histosTH1F["muonDeltaY_jpsi"]->Fill(deltay, event_weight ); 
						histosTH2F["jpsi_DeltaPhi_vs_dimuon_pt"]->Fill(dimuon_pt,deltaphi, event_weight );
						jpsi_mass = dimuon_lorentz.M();  
						jpsi_pt=dimuon_lorentz.Pt(); 
						jpsi_pt2=dimuon_lorentz.Pt()*dimuon_lorentz.Pt(); 
						jpsi_eta =dimuon_lorentz.Eta();
						jpsi_rapidity = dimuon_lorentz.Rapidity();
						if(verbose)cout<<"finalizando jpsi mass cut"<<endl; 
						} */  
						//cout<<"final do loop mass range"<<endl;
					}
				}
				//if(!mu2_selected)continue;

				if(verbose)cout<<"::FIM ::"<<endl;

				//===========================================================================
				//-------------------
				// After selection 
				//-------------------

				//-------------------
				// Generator-level proton distributions
				//-------------------

				//	}

				//-------------------
				// Detector-level distributions
				//-------------------
				/*histosTH1F["track_multiplicity"]->Fill( n_tracks_selected, event_weight );

				  histosTH1F["pf_etaMax"]->Fill( pfEtaMax, event_weight ); 
				  histosTH1F["pf_etaMin"]->Fill( pfEtaMin, event_weight );

				  double pfDeltaEta = pfEtaMax - pfEtaMin;
				  histosTH1F["pf_deltaEta"]->Fill( pfDeltaEta, event_weight ); 
				  histosTH1F["vtx_xpos_after"]->Fill( xpos, event_weight );
				  histosTH1F["vtx_ypos_after"]->Fill( ypos, event_weight );

				  histosTH1F["pf_EPlusPz"]->Fill( pfEPlusPz, event_weight );
				  histosTH1F["pf_EMinusPz"]->Fill( pfEMinusPz, event_weight );
				  histosTH1F["pf_xiPlus"]->Fill( pfXiPlusReco, event_weight );
				  histosTH1F["pf_xiMinus"]->Fill( pfXiMinusReco, event_weight );
				  histosTH1F["pf_logXiPlus"]->Fill( log10(pfXiPlusReco), event_weight );
				  histosTH1F["pf_logXiMinus"]->Fill( log10(pfXiMinusReco), event_weight );
				 */


				/*			for(vector<MyPFCand>::iterator it_pfcand = pFlow_coll->begin();
							it_pfcand != pFlow_coll->end(); ++it_pfcand){
							int partType = it_pfcand->particleId;
							double eta = it_pfcand->Eta();
							double energy = it_pfcand->Energy();
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
							histosTH2F["energyVsEtaHadronHF"]->Fill( eta, energy, event_weight );
				//if( part->ecalEnergy() > 0. ) histosTH2F["energyVsEtaHadronHFEcalEnergy"]->Fill( eta, energy, event_weight );
				//else                          histosTH2F["energyVsEtaHadronHFNoEcalEnergy"]->Fill( eta, energy, event_weight );
				} else if(partType == MyPFCand::egamma_HF) 
				histosTH2F["energyVsEtaEGammaHF"]->Fill( eta, energy, event_weight );
				}*/

				///////////////////////////////////////////////////////////////////////////////

} // End of loop over events in a file

cout<<"Total of evts="<< nev << endl << *itfiles << endl;
cout<<"n_vertices_selected ="<< n_vertices_selected << endl;
cout<<"n_select_Vertex_After_vtx_cut="<< n_select_Vertex_After_vtx_cut<<endl;
cout<<"n_tracks_selected="<< n_tracks_selected<<endl;
cout<<"n_muon_selected="<< n_muon_selected<<endl;
cout<<"n_dimuons_selected="<< n_dimuons_selected<<endl;

//cout<<"n_jpsi_selected="<< n_jpsi_selected<<endl;

//  }	
// Close current file
file->Close();

} // End of loop over files
//==========

//}
//==========
// Output file
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
//}
}
/*TFile* output = new TFile(outputFileName.c_str(),"RECREATE");
  output->cd();



  for(map<string,TH1F*>::iterator it_histo = histosTH1F.begin();
  it_histo != histosTH1F.end(); ++it_histo)
  (*it_histo).second->Write();
  for(map<string,TH2F*>::iterator it_histo = histosTH2F.begin();
  it_histo != histosTH2F.end(); ++it_histo)
  (*it_histo).second->Write();

  output->Close();
  }*/
