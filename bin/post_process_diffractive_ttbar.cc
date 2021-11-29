#include <iostream>
#include <math.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include "TApplication.h"
#include "TMatrixD.h"
#include "TF1.h"
#include "TTree.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TMath.h"

using namespace std;


int main(int argc, char* argv[]){
//int  main(TString input="", bool keepProtons = true){

	// Check arguments
	/*if (argc < 2 || argc > 4) {
		cout << "ERROR: missing arguments" << endl
			<< "Usage: " << argv[0] << " inputMCFileName keepProton=true" << endl;
		cout << "Exampe: " << argv[0] << " /eos/cms/store/group/phys_top/TTbarCentralExclProd/ntuples/mc/excl_ttbar_semilep_QED_xa120_era2017_preTS2.root" << endl;
		return 1;
	}*/


	// Check input files
	if(gSystem->AccessPathName(argv[1])){
		cout << "ERROR: missing input file " << argv[1] << endl;
		return 0;
	}
	
	
	bool keepProtons = true;
	if(argc>2) keepProtons=keepProtons && TString(argv[2])==TString("true");

	cout<<"Keep proton tree ? "<<(keepProtons==true)<<endl;

	TString infile(argv[1]);
	bool isData = infile.Contains("Data") || infile.Contains("SingleElectron") || infile.Contains("SingleMuon");

	//Preparing output file
	TString outfile = TString(infile.Tokenize("/")->Last()->GetName()).Prepend("slimmed_");

		

	TFile file(argv[1]);
	TTree* tree = (TTree*) file.Get("tree");
	TFile output(outfile,"RECREATE","");
	TTree* newTree = tree->CloneTree(0);
	
	if(isData && keepProtons){
	TTree* protons = (TTree*) file.Get("protons");
	TTree* newProtonTree = protons->CloneTree();
	}	
	
	TH1F* evt_count = (TH1F*)file.Get("evt_count");
	TH1F* ch_tag = (TH1F*)file.Get("ch_tag");
	TH1F* nbjets = (TH1F*)file.Get("nbjets");
	TH1F* njets = (TH1F*)file.Get("njets");
	TH1F* nvtx = (TH1F*)file.Get("nvtx");
	evt_count->Write();
	ch_tag->Write();
	nbjets->Write();
	njets->Write();
	nvtx->Write();

	float acoplanarity, xi_tt, Xi, weight, xi_diff,lightJetsTotalEnergy, JetsMass, SumEnergy, xi_PF;
	TBranch *b_acoplanarity = newTree->Branch("acoplanarity",&acoplanarity,"acoplanarity/F");
	TBranch *b_xitt = newTree->Branch("xitt",&xi_tt,"xitt/F");
	TBranch *b_xiPF = newTree->Branch("xiPF",&xi_PF,"xiPF/F");
	TBranch *b_xiDiff = newTree->Branch("xiDiff",&xi_diff,"xiDiff/F");
	TBranch *b_LightJetsTotalEnergy = newTree->Branch("lightJetsTotalNergy",&lightJetsTotalEnergy,"lightJetsTotalEnergy/F");
	TBranch *b_JetsMass = newTree->Branch("JetsMass",&JetsMass,"JetsMass/F");
	TBranch *b_SumEnergy = newTree->Branch("SumEnergy",&SumEnergy,"SumEnergy/F");



	int nevents = tree->GetEntries();

	for(int h=0; h < nevents; h++){ 

		tree->GetEvent(h);

		//	weight = 1.0;
		//	if(!isData)weight=tree->GetLeaf("weight")->GetValue(0);

		float cat = tree->GetLeaf("cat")->GetValue(0);

		float Xip1 = tree->GetLeaf("p1_xi")->GetValue(0);
		float Xip2 = tree->GetLeaf("p2_xi")->GetValue(0);

		TLorentzVector top1 = TLorentzVector();
		TLorentzVector top2 = TLorentzVector();

		float Etap3 = tree->GetLeaf("t_eta")->GetValue(0); //CAREFULLLLLL THIS IS A "RAPIDITY"
		float Ptp3 = tree->GetLeaf("t_pt")->GetValue(0);
		float Phip3 = tree->GetLeaf("t_phi")->GetValue(0);
		float Mp3 = tree->GetLeaf("t_m")->GetValue(0);


		float Etap4 = tree->GetLeaf("tbar_eta")->GetValue(0);
		float Ptp4 = tree->GetLeaf("tbar_pt")->GetValue(0);
		float Phip4 = tree->GetLeaf("tbar_phi")->GetValue(0);
		float Mp4 = tree->GetLeaf("tbar_m")->GetValue(0);



		float top1X=Ptp3*TMath::Cos(Phip3);
		float top1Y=Ptp3*TMath::Sin(Phip3);
		float top1Z=TMath::SinH(Etap3)*TMath::Sqrt(Ptp3*Ptp3+Mp3*Mp3);

		float top2X=Ptp4*TMath::Cos(Phip4);
		float top2Y=Ptp4*TMath::Sin(Phip4);
		float top2Z=TMath::SinH(Etap4)*TMath::Sqrt(Ptp4*Ptp4+Mp4*Mp4);

		top1.SetXYZM(top1X,top1Y,top1Z,Mp3);
		top2.SetXYZM(top2X,top2Y,top2Z,Mp4);

		TLorentzVector pair = TLorentzVector();

		float EtaPair = tree->GetLeaf("ttbar_eta")->GetValue(0);
		float PtPair = tree->GetLeaf("ttbar_pt")->GetValue(0);
		float PhiPair = tree->GetLeaf("ttbar_phi")->GetValue(0);
		float MPair = tree->GetLeaf("ttbar_m")->GetValue(0);
		float EPair = tree->GetLeaf("ttbar_E")->GetValue(0);

		float PairX=PtPair*TMath::Cos(PhiPair);
		float PairY=PtPair*TMath::Sin(PhiPair);
		float PairZ=TMath::SinH(EtaPair)*TMath::Sqrt(PtPair*PtPair+MPair*MPair);


		pair.SetXYZM(PairX,PairY,PairZ,MPair);

		acoplanarity=1.-fabs(fmod(((top1.Phi()-top2.Phi())/TMath::Pi()),2.)) ;


		double m_b1=tree->GetLeaf("bJet0_m")->GetValue(0);
		double m_b2=tree->GetLeaf("bJet1_m")->GetValue(0);
		double m_b3=tree->GetLeaf("bJet2_m")->GetValue(0);
		double m_b4=tree->GetLeaf("bJet3_m")->GetValue(0);

		double E_b1=tree->GetLeaf("bJet0_E")->GetValue(0);
		double E_b2=tree->GetLeaf("bJet1_E")->GetValue(0);
		double E_b3=tree->GetLeaf("bJet2_E")->GetValue(0);
		double E_b4=tree->GetLeaf("bJet3_E")->GetValue(0);

		double m_q1= tree->GetLeaf("lightJet0_m")->GetValue(0);
		double m_q2= tree->GetLeaf("lightJet1_m")->GetValue(0);
		double m_q3=tree->GetLeaf("lightJet2_m")->GetValue(0);
		double m_q4=tree->GetLeaf("lightJet3_m")->GetValue(0);

		double E_q1= tree->GetLeaf("lightJet0_E")->GetValue(0);
		double E_q2= tree->GetLeaf("lightJet1_E")->GetValue(0);
		double E_q3=tree->GetLeaf("lightJet2_E")->GetValue(0);
		double E_q4=tree->GetLeaf("lightJet3_E")->GetValue(0);

		double E_l=tree->GetLeaf("l_E")->GetValue(0);
		double met_pt=tree->GetLeaf("met_pt")->GetValue(0);

		double sumPVChPz = tree->GetLeaf("sumPVChPz")->GetValue(0);
		double sumPVChEn = tree->GetLeaf("sumPVChEn")->GetValue(0);

		lightJetsTotalEnergy = E_q1+E_q2+E_q3+E_q4;
		JetsMass = m_b1+m_b2+m_b3+m_b4+m_q1+m_q2+m_q3+m_q4;
		SumEnergy = fabs(E_b1+E_b2+E_b3+E_b4+E_q1+E_q2+E_q3+E_q4+E_l+met_pt-pair.E());

		if ((bool(Xip1>0.0) ^ bool(Xip2>0.0)) && (cat==5. || cat==4.)){

			float xi_tt_plus = (1./(13000.))*( (pair).E()+(pair).Pz());
			float xi_tt_minus = (1./(13000.))*( (pair).E()-(pair).Pz());


			float xi_PF_plus = (1./(13000.))*( sumPVChEn+sumPVChPz);
			float xi_PF_minus = (1./(13000.))*( sumPVChEn-sumPVChPz);			



			if(Xip1 > 0.0){
				xi_diff=xi_tt_plus-Xip1;
				xi_tt = xi_tt_plus;
				xi_PF = xi_PF_plus;
			}
			if(Xip2 > 0.0){
				xi_diff=xi_tt_minus-Xip2;
				xi_tt = xi_tt_minus;
				xi_PF = xi_PF_minus;
			}

			newTree->Fill();
		}

	}

	newTree->Write();
	output.Close();		

	return 0;
}
