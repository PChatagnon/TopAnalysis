#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TF1.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH3F.h"
#include "THStack.h"
#include <iostream>
#include <fstream>
using namespace std;



int Plot()
{
	gROOT->SetBatch(kTRUE);

	gStyle->SetPaintTextFormat("4.1f");
	gStyle->SetOptStat(1);
	gStyle->SetPalette(55);
	gStyle->SetLabelSize(.04, "xyz");
	gStyle->SetTitleSize(.04, "xyz");
	gStyle->SetTitleSize(.07, "t");
	gStyle->SetFrameLineWidth(2);
	gStyle->SetLineWidth(2);
	gStyle->SetHistLineWidth(2);
	gStyle->SetMarkerStyle(13);
	gStyle->SetTitleW(0.8); //per cent of the pad width
	gStyle->SetTitleH(0.1); //per cent of the pad height

	//TFile *BG_file = new TFile("merged_slimmed_BG_October.root");
	//TFile *Signal_file = new TFile("merged_slimmed_Signal_October.root");
	TFile *BG_file = new TFile("merged_BG_November2021_enriched.root");//"merged_slimmed_BG_XiPF_October21.root");
	TFile *Signal_file = new TFile("merged_slimmed_Signal_XiPF_October21.root");
	TFile *Data_file = new TFile("slimmed_merge_ntuple_CERT_Sept2021_V2.root");
	
	TTree *BG_tree = (TTree*)BG_file->Get("tree");
	TTree *Signal_tree = (TTree*)Signal_file->Get("tree");
	TTree *Data_tree = (TTree*)Data_file->Get("tree");
	
	std::vector<TString> labels {"xiDiff","maxEtaJets","minEtaJets","ttbar_m","ttbar_eta","ttbar_phi",
	"p1_xi","p2_xi","nJets","nBjets","gen_ttbar_phi","acoplanarity",
	"lightJet0_pt","lightJet0_eta","nvtx"};
	//TString labels[nblabels] = {"maxEtaJets"};
	
	
	//Xsec BG
	float ttbar_xsec = 832.0;
	
	//Xsec Signal
	float strenght_signal = 200.;
	float ttdiff_xsec = 3.27*strenght_signal;
	
	//Lumi Data
	float lumi_eraH = 223.502;
	
	float Trigger_factor = 49.1/3.53; //Ratio of events with trigger tag wrt total nb of events (trigger+ no trigger tag)
	
	int nbEvents_BG = BG_tree->GetEntries()*Trigger_factor;
	int nbEvents_Signal = Signal_tree->GetEntries()*Trigger_factor/0.04;
	
	cout<<"nbEvents_BG "<<nbEvents_BG<<endl;
	
	for(int j=0;j<2;j++){
	TCut cut = "1";
	TCut cut_1 = "p1_xi>0 && cat==5";
	TCut cut_2 = "p2_xi>0 && cat==5";
	if(j==0){cut=cut_1;}
	else {cut=cut_2;}
	TString cut_string=cut.GetTitle();
	
	for(int i=0;i<labels.size();i++){
	TString label = labels[i];//"maxEtaJets";
	float maximum = (Data_tree->GetMaximum(label));
	float minimum = (Data_tree->GetMinimum(label));
	int nBins = 25;
	TH1F *Data_hist = new TH1F("Data_hist","",nBins,minimum-0.2*(maximum-minimum),maximum+0.2*(maximum-minimum));
	Data_tree->Draw(label+">>Data_hist",cut);
	
	float essai = Data_tree->GetMaximum(label);
	cout<<"ici "<<essai<<endl;
	
	int nbBins = Data_hist->GetNbinsX();
	float min_histo = Data_hist->GetXaxis()->GetXmin();
 	float max_histo = Data_hist->GetXaxis()->GetXmax();
	cout<<nbBins<<" "<<min_histo<<" "<<max_histo<<endl;
	
	TH1F *BG_hist = new TH1F("BG_hist","",nbBins,min_histo,max_histo);
	TH1F *Signal_hist = new TH1F("Signal_hist","",nbBins,min_histo,max_histo);
	TCut weight = Form("%s*%f*%f/%i","weight",ttbar_xsec,lumi_eraH,nbEvents_BG);
	//TCut weight = Form("%f*%f",ttbar_xsec,lumi_eraH);
	TCut weight_Signal = Form("%s*%f*%f/%i","weight",ttdiff_xsec,lumi_eraH,nbEvents_Signal);
	BG_tree->Draw(label+">>BG_hist",cut*weight);
	Signal_tree->Draw(label+">>Signal_hist",cut*weight_Signal);


	Signal_hist->SetLineWidth(2);
	Signal_hist->SetLineColor(kRed);
	BG_hist->SetLineWidth(1);
	BG_hist->SetLineColor(kBlack);
	BG_hist->SetMarkerSize(0);
	//BG_hist->SetFillStyle(3001);
	BG_hist->SetFillColor(865);

	
	Data_hist->SetLineWidth(2);
	Data_hist->SetLineColor(kBlack);
	Data_hist->SetMarkerColor(kBlack);
	Data_hist->SetMarkerSize(2);
	Data_hist->SetMarkerStyle(20);
	Data_hist->SetTitle(";"+label+";Events");

	Signal_hist->SetStats(kFALSE);
	BG_hist->SetStats(kFALSE);
	Data_hist->SetStats(kFALSE);

	

	THStack *hs = new THStack("hs", "");
	hs->SetTitle(";"+label+";Events");
	hs->Add(BG_hist);
	hs->Add(Signal_hist);
	
	
	
	
	
	//Ratio plots and uncertainty
	TH1F *ratio_hist = (TH1F*)Data_hist->Clone("Data_hist");
	ratio_hist->Divide(BG_hist);
	
	TH1F *ratio_hist_uncertainty = (TH1F*)Data_hist->Clone("Data_hist");
	for(int ii=1;ii<nBins+1;ii++){
          ratio_hist_uncertainty->SetBinContent(ii,1);
          if(BG_hist->GetBinContent(ii)>0) {
            ratio_hist_uncertainty->SetBinError(ii,BG_hist->GetBinError(ii)/BG_hist->GetBinContent(ii));
          }
          else {ratio_hist_uncertainty->SetBinError(ii,1);}
        }


	
	TCanvas *cancG0 = new TCanvas("", "can0",1500,1000);
	cancG0->cd();
	
	
	 // Upper plot will be in pad1
   	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
   	//pad1->SetGridx();         // Vertical grid
  	//pad1->SetGridy();
  	pad1->SetBottomMargin(0.);
  	pad1->SetTicks(1,1);
  	
  	//pad1->SetLogy();
  	hs->SetMaximum((Data_hist->GetMaximum())*1.5);//*20.);
  	
  	
	pad1->Draw();             // Draw the upper pad: pad1
	pad1->cd();               // pad1 becomes the current pad
	
	
	
	hs->Draw("nostack   hist");
	
	TH1F* CLoned_BG_hist = (TH1F*)BG_hist->Clone("clone_gb");
	CLoned_BG_hist->SetFillColor(kBlack);
	CLoned_BG_hist->SetFillStyle(3144);
	CLoned_BG_hist->Draw(" same  e2");
	//BG_hist->Draw("same hist");
	
	Data_hist->Draw("e same");
	
	
	TText *t = new TText();
	t->SetTextAlign(22);
	t->SetTextColorAlpha(kGray, 0.50);
	t->SetTextFont(40);
	t->SetTextSize(0.25);
	t->SetTextAngle(25);
	

	auto legend = new TLegend(0.54, 0.87, 0.87, 0.67);
	legend->AddEntry(Data_hist, Form("Data (%3.1f)",Data_hist->Integral()), "lp");	
	legend->AddEntry(BG_hist, Form("t#bar{t}X (%3.1f)",BG_hist->Integral()), "f1");
	legend->AddEntry(Signal_hist, Form("Signal (%3.1f), #sigma=%3.1f pb",Signal_hist->Integral(),ttdiff_xsec), "l");
	legend->SetFillStyle(0);
	legend->SetLineWidth(0);
	legend->Draw("same ");
	
	//Labels
	double x_top_label = 0.80;
	
	TPaveText * CMS_Internal = new TPaveText(0.20,x_top_label,0.388191,x_top_label+0.1,"NDC");
	/*lumi->SetFillColor(0); */
	CMS_Internal->SetFillStyle(4050);
	CMS_Internal->SetLineColor(0);
	CMS_Internal->SetTextFont(42); CMS_Internal->SetTextSize(0.0599401); CMS_Internal->SetBorderSize(0);
	CMS_Internal->AddText(Form("CMS Internal"));
	CMS_Internal->Draw();
	
	TPaveText * lumi = new TPaveText(0.20,x_top_label-0.2,0.388191,x_top_label,"NDC");
	/*lumi->SetFillColor(0); */
	lumi->SetFillStyle(4050);
	lumi->SetLineColor(0);
	lumi->SetTextFont(42); lumi->SetTextSize(0.0599401); lumi->SetBorderSize(0);
	lumi->AddText(Form("\\sqrt{s}= 13 TeV, \\mathscr{L} = %1.1f pb^{-1}",lumi_eraH));
	//lumi->AddText("#mu+jets channel, ".append(CutString));
	lumi->AddText("#mu+jets channel"+cut);
	lumi->Draw();
	
	TPaveText * Warning = new TPaveText(0.20,x_top_label-0.25,0.388191,x_top_label-0.15,"NDC");
	Warning->SetFillStyle(4050);
	Warning->SetLineColor(0);
	Warning->SetTextFont(42); Warning->SetTextSize(0.052); Warning->SetBorderSize(0);
	Warning->AddText("Extra /10 factor");
	Warning->Draw();
	
	t->DrawTextNDC(.5, .53, "Preliminary");

	// lower plot will be in pad
  	cancG0->cd();          // Go back to the main canvas before defining pad2
   	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.25);
   	pad2->SetTopMargin(0);
   	pad2->SetTicks(1,1);
   	pad2->SetGridy();

   	
   	pad2->SetBottomMargin(0.3);
   	//pad2->SetGridx(); // vertical grid
   	pad2->Draw();
   	pad2->cd();       // pad2 becomes the current pad

   	
   	ratio_hist_uncertainty->SetFillColor(kYellow-2);
   	//ratio_hist_uncertainty->SetFillStyle(3001);
        ratio_hist_uncertainty->SetLineColor(1);
        ratio_hist_uncertainty->SetLineWidth(2);
        ratio_hist_uncertainty->SetMarkerSize(0);
        
        
        ratio_hist_uncertainty->SetMaximum(1.5);
        ratio_hist_uncertainty->SetMinimum(0.5);
	ratio_hist_uncertainty->Draw("e2");
	ratio_hist->Draw("ep same");
	pad2->Update();
	/*TH2C *hgrid = new TH2C("hgrid","",nBins,min_histo,max_histo,14,1.7,0.3);   
   	hgrid->GetYaxis()->SetNdivisions(14);
  	hgrid->Draw("same text");*/

   	// Ratio plot (h3) settings
   	ratio_hist_uncertainty->SetTitle(""); // Remove the ratio title

   	// Y axis ratio plot settings
   	ratio_hist_uncertainty->GetYaxis()->SetTitle("Data/MC ratio");
   	ratio_hist_uncertainty->GetYaxis()->SetNdivisions(505);
   	ratio_hist_uncertainty->GetYaxis()->SetTitleSize(30);
   	ratio_hist_uncertainty->GetYaxis()->SetTitleFont(43);
   	ratio_hist_uncertainty->GetYaxis()->SetTitleOffset(1.55);
   	ratio_hist_uncertainty->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   	ratio_hist_uncertainty->GetYaxis()->SetLabelSize(30);

   	// X axis ratio plot settings
   	ratio_hist_uncertainty->GetXaxis()->SetTitleSize(30);
   	ratio_hist_uncertainty->GetXaxis()->SetTitleFont(43);
   	ratio_hist_uncertainty->GetXaxis()->SetTitleOffset(4.);
   	ratio_hist_uncertainty->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   	ratio_hist_uncertainty->GetXaxis()->SetLabelSize(30);

	auto legendR = new TLegend(0.81, 0.95, 0.9, 0.75);
	legendR->AddEntry(ratio_hist_uncertainty, "MC Uncert.", "f1");
	legendR->Draw("same ");
	legendR->SetFillStyle(0);
	legendR->SetLineWidth(0);

	//cancG0->SaveAs(Form("pdfs/%s.pdf",label));
	cancG0->SaveAs("pdfs/"+label+cut_string+".pdf");
	cancG0->SaveAs("pngs/"+label+cut_string+".png");
	
	//Data_hist->SaveAs("essai.root");
	
	}
	}
	

	//gApplication->Terminate();

	return 0;
}
