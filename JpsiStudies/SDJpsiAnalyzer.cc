////////// Header section /////////////////////////////////////////////
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "TH1F.h"

class SDJpsiAnalyzer: public edm::EDAnalyzer {
	public:
		/// Constructor
		SDJpsiAnalyzer(const edm::ParameterSet& pset);

		/// Destructor
		virtual ~SDJpsiAnalyzer();

		// Operations

		void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);

		//virtual void beginJob(const edm::EventSetup& eventSetup) ;
		virtual void beginJob(); 
		virtual void endJob() ;

	private:
		// Input from cfg file
		edm::InputTag genParticlesTag_;
		int particle1Id_;
		int particle2Id_;

		// Histograms
		TH1F* hPart1Pt;
		TH1F* hPart1Eta;
		TH1F* hPart1Phi;
		TH1F* hPart2Pt;
		TH1F* hPart2Eta;
		TH1F* hPart2Phi;
		TH1F* hMesonPt;
		TH1F* hMesonEta;
		TH1F* hMesonPhi;
		TH1F* hMesonM;
		TH1F* hMesonMaftercuts;
		TH1F* hMesonEtaaftercuts;
		TH1F* hMesonPt2;

		TH1F* hEnergyvsEta;
		TH1F* hXiGen;
		TH1F* hProtonPt2;	

		TH1F* hXiGenPlus;
		TH1F* hXiGenMinus;
		TH1F* hTGenPlus;
		TH1F* hTGenMinus;
		TH1F* hTDirGenPlus;
		TH1F* hTDirGenMinus;
		TH1F* hTPlus;
		TH1F* hTMinus;


		TH1F* hXiApGenPlus;
		TH1F* hXiApGenMinus;

		TH1F* hDeltaPt;
		TH1F* hDeltaEta;
		TH1F* hDeltaPhi;

		int nevents;
		bool debug;
		double Ebeam;
		int neventsaftercuts;
};

////////// Source code ////////////////////////////////////////////////
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

/// Constructor
SDJpsiAnalyzer::SDJpsiAnalyzer(const edm::ParameterSet& pset)
{
	genParticlesTag_ = pset.getParameter<edm::InputTag>("GenParticleTag");
	particle1Id_ = pset.getParameter<int>("Particle1Id");
	particle2Id_ = pset.getParameter<int>("Particle2Id");

	debug = pset.getUntrackedParameter<bool>("debug",false);
	if(debug){
		std::cout << ">>> First particle Id: " << particle1Id_ << std::endl;
		std::cout << ">>> Second particle Id: " << particle2Id_ << std::endl;
	}
}

/// Destructor
SDJpsiAnalyzer::~SDJpsiAnalyzer(){
}

void SDJpsiAnalyzer::beginJob(){
	edm::Service<TFileService> fs;
	//TH1::SetDefaultSumw2(true);

	hPart1Pt = fs->make<TH1F>("hPart1Pt","hPart1Pt",100,0.,100.);
	hPart1Eta = fs->make<TH1F>("hPart1Eta","hPart1Eta",100,-5.,5.);
	hPart1Phi = fs->make<TH1F>("hPart1Phi","hPart1Phi",100,-1.2*M_PI , 1.2*M_PI);
	hPart2Pt = fs->make<TH1F>("hPart2Pt","hPart2Pt",100,0.,100.);
	hPart2Eta = fs->make<TH1F>("hPart2Eta","hPart2Eta",100,-5.,5.);
	hPart2Phi = fs->make<TH1F>("hPart2Phi","hPart2Phi",100,-1.2*M_PI , 1.2*M_PI);
	hMesonPt = fs->make<TH1F>("hMesonPt","hMesonPt",100,0.,50.);
	hMesonEta = fs->make<TH1F>("hMesonEta","hMesonEta",100,-5.,5.);
	hMesonPhi = fs->make<TH1F>("hMesonPhi","hMesonPhi",100,-1.2*M_PI , 1.2*M_PI);
	hMesonM = fs->make<TH1F>("hMesonM","hMesonM",100,0.,5.);
	hMesonMaftercuts = fs->make<TH1F>("hMesonMaftercuts","hMesonM",100,0.,5.);
	hMesonEtaaftercuts = fs->make<TH1F>("hMesonEtaaftercuts","hMesonEta",100,-5.,5.);

	hMesonPt2 = fs->make<TH1F>("hMesonPt2","hMesonPt2",100,0.,3.0);

	hXiGenPlus = fs->make<TH1F>("hXiGenPlus","hXiGenPlus",100,0.,1.);
	hXiGenMinus = fs->make<TH1F>("hXiGenMinus","hXiGenMinus",100,0.,1.);
	hTGenPlus = fs->make<TH1F>("hTGenPlus","hTGenPlus",100,0.,5.);
	hTGenMinus = fs->make<TH1F>("hTGenMinus","hTGenMinus",100,0.,5.);
	hTMinus = fs->make<TH1F>("hTMinus","hTMinus",100,0.,5.);
	hTPlus = fs->make<TH1F>("hTPlus","hTPlus",100,0.,5.);
	hTDirGenPlus = fs->make<TH1F>("hTDirGenPlus","hTDirGenPlus",100,0.,5.);
	hTDirGenMinus = fs->make<TH1F>("hTDirGenMinus","hTDirGenMinus",100,0.,5.);

	hXiApGenPlus = fs->make<TH1F>("hXiApGenPlus","hXiApGenPlus",100,0.,1.);
	hXiApGenMinus = fs->make<TH1F>("hXiApGenMinus","hXiApGenMinus",100,0.,1.);

	hDeltaPt = fs->make<TH1F>("hDeltaPt", "#Delta p_{T}(#mu)" , 150 , 0. , 150.);
	hDeltaEta = fs->make<TH1F>("hDeltaEta", "#Delta#eta(#mu)" , 200 , 0 , 5.2);
	hDeltaPhi = fs->make<TH1F>("hDeltaPhi", "#Delta#phi(#mu)" , 200 , 0. , 1.2*M_PI);

	hEnergyvsEta = fs->make<TH1F>("hEnergyvsEta","hEnergyvsEta",100,-15.0,15.0); 		
	hXiGen = fs->make<TH1F>("hXiGen","hXiGen",100,0.,0.21);
	hProtonPt2 = fs->make<TH1F>("hProtonPt2","hProtonPt2",100,0.,3.0);

	nevents = 0;
	Ebeam = 4000.;//Fix get the Ebeam from the event
}

void SDJpsiAnalyzer::endJob(){
	std::cout<<" N events:"<< nevents << std::endl;
	hEnergyvsEta->Scale(1/(float)nevents);
	std::cout<<" N events after cuts:"<< neventsaftercuts << std::endl;	
}

void SDJpsiAnalyzer::analyze(const edm::Event & ev, const edm::EventSetup&){
	nevents++; 
	//std::cout<<" N events:"<< nevents << std::endl;
	// Generator Information
	edm::Handle<reco::GenParticleCollection> genParticles;
	ev.getByLabel(genParticlesTag_, genParticles);
	double pz1max = 0.;
	double pz2min = 0.;
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
	double xi_plus_gen;
	double xi_minus_gen;
	double proton_pf_plus = 0.;
	double proton_pf_minus= 0.;
	double proton_pi_plus = 0.;
	reco::GenParticleCollection::const_iterator proton1 = genParticles->end();
	reco::GenParticleCollection::const_iterator proton2 = genParticles->end();
	reco::GenParticleCollection::const_iterator particle1 = genParticles->end();
	reco::GenParticleCollection::const_iterator particle2 = genParticles->end();

	for(reco::GenParticleCollection::const_iterator genpart = genParticles->begin(); genpart != genParticles->end(); ++genpart){
		//std::cout << ">>>>>>> pid,status,px,py,px,e= "  << genpart->pdgId() << " , " << genpart->status() << " , " << genpart->px() << " , " << genpart->py() << " , " << genpart->pz() << " , " << genpart->energy() << std::endl;	
		if(genpart->status() != 1) continue;

		hEnergyvsEta->Fill(genpart->eta(),genpart->energy());

		energy_gen = genpart->energy();
		px_gen = genpart->px();
		py_gen = genpart->py();
		pz_gen = genpart->pz();
		double pz = genpart->pz();
		double pz_cut = 0.7*proton_pi;
		//proton_pi_plus = sqrt(6500*4000+4000*4000);
		//proton_pf_plus = sqrt(px_gen*px_gen+py_gen*py_gen+pz_gen*pz_gen+energy_gen*energy_gen);
		//if((genpart->pdgId() == 2212)&&(pz > 0.7*Ebeam)){
		if((genpart->pdgId() == 2212)){
			//if(pz > pz1max){proton1 = genpart;pz1max=pz;
			if(pz > 0.){proton1 = genpart;pz1max=pz;
				proton_pz_plus = pz_gen; proton_energy_plus = energy_gen;
				proton_px_plus = px_gen; proton_py_plus = py_gen;    
				proton_pf_plus = sqrt(px_gen*px_gen+py_gen*py_gen+pz_gen*pz_gen);
				//TLorentzVector p3(proton_px_plus,proton_py_plus,proton_pz_plus,proton_energy_plus); // 4 momentum of p3   
			}
		}

		if((genpart->pdgId() == 2212)){ 
                //if((genpart->pdgId() == 2212)&&(pz < -0.7*Ebeam)){
			//if(pz < pz2min){proton2 = genpart;pz2min=pz;
			if(pz < 0.){proton2 = genpart;pz2min=pz;
				proton_pz_minus = pz_gen; proton_energy_minus = energy_gen;
				proton_px_minus = px_gen; proton_py_minus = py_gen;
				proton_pf_minus = sqrt(px_gen*px_gen+py_gen*py_gen+pz_gen*pz_gen);
			}

		}

		if (fabs(pz_gen) < pz_cut){
			genEPlusPz += (energy_gen + pz_gen);
			genEMinusPz += (energy_gen - pz_gen);   	
		}
		//Fix add constraint on mother/daughter relation
		if((particle1 == genParticles->end())&&(abs(genpart->pdgId()) == abs(particle1Id_))) {particle1 = genpart;continue;}
		if((particle2 == genParticles->end())&&(abs(genpart->pdgId()) == abs(particle2Id_))) {particle2 = genpart;continue;}
		}


		//double xi_plus_gen = genEPlusPz/cm; std::cout<<xi_plus_gen<<std::endl;
		//double xi_minus_gen = genEMinusPz/cm;
		double xi_proton_plus = -1.;
		double xi_proton_minus = -1.;
		double t_proton_plus = 0.;
		double t_proton_minus = 0.;
		double t_dir_plus = 0.;
		double t_dir_minus = 0.;
		double t_plus = 0.;
		double t_minus = 0.;
		double t = 0.; 

		//if(pz_gen > 0.){
			if(proton1 != genParticles->end()){
				if(debug) std::cout << "Proton 1: " << proton1->pt() << "  " << proton1->eta() << "  " << proton1->phi() << std::endl;
				double xigen1 = 1 - proton1->pz()/Ebeam;
				hXiGen->Fill(fabs(xigen1));
				hProtonPt2->Fill(proton1->pt()*proton1->pt());
				if(proton_pz_plus > 0.){
					xi_proton_plus =  ( 1 - (proton_pz_plus/proton_pi) );
					//math::XYZTLorentzVector p2(0,0,6500,6500);
					math::XYZTLorentzVector p2(0,0,4000,4000);
					math::XYZTLorentzVector p3(proton_px_plus,proton_py_plus,proton_pz_plus,proton_energy_plus); // 4 momentum of p3
					math::XYZTLorentzVector vec_t = (p3 - p2);
					t_plus=vec_t.M2();
					if(debug) std::cout << ">>> t_plus: " << t_plus << ">>> xi_plus: "<<xi_proton_plus <<" "<< xigen1 <<std::endl;
				}

			}
		//}	
		//if(pz_gen < 0.){
			if(proton2 != genParticles->end()){
				if(debug) std::cout << "Proton 2: " << proton2->pt() << "  " << proton2->eta() << "  " << proton2->phi() << std::endl;	
				double xigen2 = 1 + proton2->pz()/Ebeam;
				hXiGen->Fill(xigen2);
				hProtonPt2->Fill(proton2->pt()*proton2->pt());
				if(proton_pz_minus < 0.){
					xi_proton_minus = (proton_pz_minus < 0.) ? ( 1 + (proton_pz_minus/proton_pi) ) : -1.;
					//math::XYZTLorentzVector p1(0,0,-6500,6500);
					math::XYZTLorentzVector p1(0,0,-4000,4000);
					math::XYZTLorentzVector pm(proton_px_minus,proton_py_minus,proton_pz_minus,proton_energy_minus); // 4 momentum of p3
					math::XYZTLorentzVector vec_t = (pm - p1);
					t_minus=vec_t.M2();
				}

			}
		//}

		if((particle1 != genParticles->end())&&(particle2 != genParticles->end())){
			if(debug) std::cout << ">>> particle 1 pt,eta: " << particle1->pt() << " , " << particle1->eta() << std::endl;
			hPart1Pt->Fill(particle1->pt());
			hPart1Eta->Fill(particle1->eta());
			hPart1Phi->Fill(particle1->phi());

			if(debug) std::cout << ">>> particle 2 pt,eta: " << particle2->pt() << " , " << particle2->eta() << std::endl;
			hPart2Pt->Fill(particle2->pt());
			hPart2Eta->Fill(particle2->eta());
			hPart2Phi->Fill(particle2->phi());

			double deltapt = fabs(particle1->pt() - particle2->pt());
			double deltaeta = fabs(particle1->eta() - particle2->eta());
			double deltaphi = fabs(particle1->phi() - particle2->phi());
			if(deltaphi > M_PI){deltaphi = 2.0*M_PI - deltaphi;}

			math::XYZTLorentzVector mymeson(particle1->px() + particle2->px(),
					particle1->py() + particle2->py(),
					particle1->pz() + particle2->pz(),
					particle1->energy() + particle2->energy());

			xi_plus_gen = genEPlusPz/cm; std::cout<<xi_plus_gen<<std::endl;
			xi_minus_gen = genEMinusPz/cm;
			hMesonPt->Fill(mymeson.pt());
			hMesonEta->Fill(mymeson.eta());
			hMesonPhi->Fill(mymeson.phi());
			hMesonM->Fill(mymeson.M());
			hMesonPt2->Fill(mymeson.pt()*mymeson.pt());
			hDeltaPt->Fill(deltapt);
			hDeltaEta->Fill(deltaeta);
			hDeltaPhi->Fill(deltaphi);
			hXiGenPlus->Fill(fabs(xi_proton_plus));
			hXiGenMinus->Fill(xi_proton_minus);
			hXiApGenPlus->Fill(xi_plus_gen);
			hXiApGenMinus->Fill(xi_minus_gen);
			/*hTGenPlus->Fill(fabs(t_proton_plus));
			  hTGenMinus->Fill(fabs(t_proton_minus));
			  hTDirGenPlus->Fill(fabs(t_dir_plus));
			  hTDirGenMinus->Fill(fabs(t_dir_minus));*/
			hTPlus->Fill(fabs(t_plus));
			hTMinus->Fill(fabs(t_minus));

			if(particle1->pt() > 0.0 && particle2->pt() > 0.0){
				std::cout << "pTmu1:" <<particle1->pt()<<"pTmu2:" <<particle2->pt()<<std::endl;
				if(fabs(particle1->eta()) < 2.45 && fabs(particle2->eta() < 2.45)){
				if(mymeson.M()> 2.9 && mymeson.M()< 3.2){
				// if(fabs(t_plus)> 0.03 && fabs(t_plus)< 1.0 && fabs(xi_proton_plus) < 0.23){
				if(fabs(xi_proton_plus) < 0.04){
                                //if(fabs(xi_proton_minus) < 0.04){
					++neventsaftercuts;
					hMesonMaftercuts->Fill(mymeson.M());
					hMesonEtaaftercuts->Fill(mymeson.eta());
				}//xi
				}//mass
				}//eta
			}

			}	
		}

		DEFINE_FWK_MODULE(SDJpsiAnalyzer);

