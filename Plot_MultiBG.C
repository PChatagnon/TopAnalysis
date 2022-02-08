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



int Plot_MultiBG()
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


	std::vector<vector<TString>> samples
	{
		{"slimmed_merge_Wt_BG_2_enriched.root","1","34","Wt","35.24"},
		{"slimmed_merge_Wtbar__BG_enriched.root","1","94","Wtbar","35.24"},
		{"slimmed_merged_BG_January_2022_with_HLT_enriched_MVA.root","832","865","ttX",""},	
		{"slimmed_merge_January_2021_Signal_2_enriched_MVA.root","3","99","Signal",""}
	};


	TFile *Data_file = new TFile("slimmed_merge_January_2022_5_Data_MVA.root");
	TTree *Data_tree = (TTree*)Data_file->Get("tree");

	std::vector<vector<TString>> labels
	{
		/*{"weight_BDT","#w_{BDT}","no_arrow","no_log","unblind","default_range", "0","0","default_binning","0","no_cut","tree"},
			{"xiDiff","#xi Diff.","no_arrow","no_log","unblind","default_range","0","0","default_binning","0","cut","tree"},
			{"maxEtaJets","Max. #eta Jets","arrow","no_log","unblind","default_range","0","0","default_binning","0","cut","tree"},
			{"minEtaJets","Min #eta Jets","arrow","no_log","unblind","default_range","0","0","default_binning","0","cut","tree"},
			{"ttbar_m","M_{t#bar{t}} (GeV)","no_arrow","no_log","unblind","default_range","0","0","default_binning","0","no_cut","tree"},
			{"ttbar_eta","#eta_{t#bar{t}}","arrow","no_log","unblind","default_range","0","0","default_binning","0","cut","tree"},
			{"ttbar_phi","#phi_{t#bar{t}}","no_arrow","no_log","unblind","default_range","0","0","custom_binning","10","no_cut","tree"},
			{"p1_xi","#xi (P1/+z)","no_arrow","no_log","unblind","custom_range","0.01","0.15","custom_binning","10","cut","tree"},
			{"p2_xi","#xi (P2/-z)","no_arrow","no_log","unblind","custom_range","0.01","0.15","custom_binning","10","cut","tree"},
			{"nJets","Nb. Jets","no_arrow","no_log","unblind","custom_range","4","10","custom_binning","6","no_cut","tree"},
			{"nBjets","Nb. B jets","no_arrow","no_log","unblind","custom_range","1","5","custom_binning","4","no_cut","tree"},
			{"gen_ttbar_phi","gen_ttbar_phi","no_arrow","no_log","unblind","default_range","0","0","default_binning","0","cut","tree"},*/
			{"acoplanarity","Acoplanarity","no_arrow","log","unblind","custom_range","-0.6","0.6","custom_binning","40","no_cut","tree"},
			{"lightJet0_pt","Leading light jet Pt (GeV)","no_arrow","no_log","unblind","default_range","0","0","default_binning","0","no_cut","tree"},
			{"lightJet0_eta","Leading light jet #eta","arrow","no_log","unblind","default_range","0","70","default_binning","0","cut","tree"},/*
			{"nvtx","Nb vertex","no_arrow","no_log","unblind","custom_range","0","12","custom_binning","12","no_cut","tree"},
			{"l_pt","l_pt","no_arrow","no_log","unblind","default_range","0","0","default_binning","0","no_cut","tree"},
			{"nchPV","nchPV","no_arrow","no_log","unblind","default_range","0","0","default_binning","00","no_cut","tree"},
			{"xiPF-xitt","xiPF-xitt","no_arrow","no_log","blind","custom_range","-0.16","0.02","default_binning","0","cut","tree"},
			{"xiPF","xi","no_arrow","no_log","unblind","default_range","0","0","default_binning","0","cut","tree"},
			{"xitt","xi","no_arrow","no_log","unblind","default_range","0","0","default_binning","0","cut","tree"},
			{"xiPF+xiDiff-xitt","#xi_{PF} Diff.","no_arrow","no_log","blind","custom_range","-0.140","0.08","default_binning","0","cut","tree"},
			{"sumPVChPz","sumPVChPz","no_arrow","no_log","unblind","default_range","0","0","default_binning","0","no_cut","tree"},
			{"sumPVChEn","sumPVChEn","no_arrow","no_log","unblind","default_range","0","0","default_binning","0","no_cut","tree"},*/

			{"evt_count","Selection Stage","no_arrow","log","unblind","default_range","0","0","default_binning","0","no_cut","hist"}/*,
			{"nbjets","Nb Jets","no_arrow","log","unblind","default_range","0","0","default_binning","0","no_cut","hist"},
			{"nvtx","Nb vtx","no_arrow","log","unblind","default_range","0","0","default_binning","0","no_cut","hist"},
			{"njets","Nb Jets","no_arrow","log","unblind","default_range","0","0","default_binning","0","no_cut","hist"}*/

	};

	
	//Xsec Signal
	float strenght_signal = 40.;
	
	//Lumi Data
	float lumi_eraH = 196.5;//223.502; //As in the exclusive approach

	float Norm_factor = 1.0;//0.2;//1.0;//


	for(int i=0;i<labels.size();i++){


		for(int k=0;k<2;k++){
			TCut cut = "1";
			TCut cut_1 = "p1_xi>0";
			TCut cut_2 = "p2_xi>0";
			if(k==0){cut=cut_1;}
			else{cut=cut_2;}


			//////////////////////////////////////////////////
			//Set options for each label
			TString label = labels[i][0];
			TString xAxis_label = labels[i][1];
			TString arrow_option = labels[i][2];
			TString log_option = labels[i][3];
			TString blind_option = labels[i][4];
			TString range_option = labels[i][5];
			TString binning_option = labels[i][8];
			TString no_cut = labels[i][10];
			TString hist_mode = labels[i][11];
			//////////////////////////////////////////////////

			if(no_cut=="no_cut"){
				cut="1";
			}

			TString cut_string=cut.GetTitle();

			cout<<"cut option "<<no_cut<<" "<<(no_cut=="no_cut")<<"  "<<cut_string<<endl;

			float maximum = (Data_tree->GetMaximum(label));
			float minimum = (Data_tree->GetMinimum(label));
			float min_histo_ini = minimum-0.2*(maximum-minimum);
			float max_histo_ini = maximum+0.2*(maximum-minimum);
			int nBins = 20;

			if(range_option=="custom_range"){
				cout<<"here "<<endl;
				min_histo_ini = stof((string)labels[i][6].Data());
				max_histo_ini = stof((string)labels[i][7].Data());
			}
			if(binning_option=="custom_binning"){
				cout<<"here "<<labels[i][9].Data()<<endl;
				nBins=stoi((string)labels[i][9].Data());
			}

			TH1F *Data_hist = new TH1F("Data_hist","",nBins,min_histo_ini,max_histo_ini);
			Data_tree->Draw(label+">>Data_hist",cut);

			if(hist_mode=="hist"){Data_hist = (TH1F*)Data_file->Get(label);}

			int nbBins = Data_hist->GetNbinsX();
			float min_histo = Data_hist->GetXaxis()->GetXmin();
			float max_histo = Data_hist->GetXaxis()->GetXmax();
			cout<<"Bining and range "<<nbBins<<" "<<min_histo<<" "<<max_histo<<endl;

			double r = 0.0;
			r=Data_hist->GetBinContent(1)+Data_hist->GetBinContent(0); Data_hist->SetBinContent(1,r);
			r=Data_hist->GetBinContent(nBins)+Data_hist->GetBinContent(nBins+1); Data_hist->SetBinContent(nBins,r);
			
			Data_hist->SetLineWidth(2);
			Data_hist->SetLineColor(kBlack);
			Data_hist->SetMarkerColor(kBlack);
			Data_hist->SetMarkerSize(2);
			Data_hist->SetMarkerStyle(20);
			Data_hist->SetTitle(";"+label+";Events");
			Data_hist->SetStats(kFALSE);
			
			auto legend = new TLegend(0.54, 0.87, 0.87, 0.67);
			legend->AddEntry(Data_hist, Form("Data (%3.1f)",Data_hist->Integral()), "lp");	

			THStack *hs = new THStack("hs", "");
			TH1F* Signal_hist = new TH1F("Signal_hist","",nBins,min_histo_ini,max_histo_ini);

			TH1F* BG_hist = new TH1F("BG_hist","",nBins,min_histo_ini,max_histo_ini);

			for(int j=0;j<samples.size();j++){

				TFile *sample_file = new TFile(samples[j][0]);
				
				int nbEvents_sample = ((TH1F*)sample_file->Get("evt_count"))->GetBinContent(1)/Norm_factor;//BG_tree->GetEntries()*Trigger_factor;

				TTree *sample_tree = (TTree*)sample_file->Get("tree");
				
				//Xsec BG
				double xsec = stof((string)samples[j][1].Data());
				
				//add xsec scale
				if(samples[j][3]=="Signal"){xsec = xsec*strenght_signal;}
				
				TString hist_name = Form("sample_hist_%s",samples[j][3].Data());

				TH1F *sample_hist = new TH1F(hist_name,"",nbBins,min_histo,max_histo);
				
				
				TCut weight = Form("%s*%f*%f/%i","weight",xsec,lumi_eraH,nbEvents_sample);
				


				sample_tree->Draw(label+">>"+hist_name,cut*weight);

				if(hist_mode=="hist"){
					sample_hist = (TH1F*)sample_file->Get(label);
				}

				cout<<cut*weight<<endl;

				//UnderFlow-OverFlow
				
				r=sample_hist->GetBinContent(1)+sample_hist->GetBinContent(0); sample_hist->SetBinContent(1,r);
				r=sample_hist->GetBinContent(nBins)+sample_hist->GetBinContent(nBins+1); sample_hist->SetBinContent(nBins,r);
				
				
				
				if(samples[j][3]!="Signal"){
					sample_hist->SetLineWidth(0);
					sample_hist->SetLineColor(kBlack);
					sample_hist->SetMarkerSize(0);
					cout<<std::stoi(samples[j][2].Data())<<endl;
					sample_hist->SetFillColor(std::stoi(samples[j][2].Data()));
				}
				
				if(samples[j][3]=="Signal"){
					sample_hist->SetLineWidth(2);
					sample_hist->SetLineColor(kRed);
				}

				sample_hist->SetStats(kFALSE);
				
				
				hs->SetTitle(";"+label+";Events");
				
				
				if(samples[j][3]!="Signal"){
				hs->Add(sample_hist);
				BG_hist->Add(sample_hist);
				}
				else Signal_hist = sample_hist;
				
				
				if(samples[j][3]!="Signal" && samples[j][3]!="Wt" && samples[j][3]!="Wtbar"){
					legend->AddEntry(sample_hist, Form("%s (%3.1f), #sigma=%3.1f pb", samples[j][3].Data(),sample_hist->Integral(),stof(samples[j][1].Data())), "f1" );
				}
				else if(samples[j][3]!="Signal"){
					legend->AddEntry(sample_hist, Form("%s (%3.1f), #sigma=%3.1f pb", samples[j][3].Data(),sample_hist->Integral(),stof(samples[j][4].Data())), "f1" );
				}
				else legend->AddEntry(Signal_hist, Form("Signal (%3.1f), #sigma=%3.1f (SM) x %3.1f pb",Signal_hist->Integral(),stof(samples[j][1].Data()),strenght_signal), "l");

			}

			

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

			if(log_option=="log")pad1->SetLogy();
			float max_display = (Data_hist->GetMaximum())*1.5;
			if(blind_option!="blind") hs->SetMaximum(max_display);//*20.);


			pad1->Draw();             // Draw the upper pad: pad1
			pad1->cd();               // pad1 becomes the current pad


			
			hs->Draw("hist");//nostack
			Signal_hist->Draw("hist same");

			/*TH1F* CLoned_BG_hist = (TH1F*)BG_hist->Clone("clone_gb");*/
			BG_hist->SetFillColor(kBlack);
			BG_hist->SetFillStyle(3144);
			BG_hist->Draw(" same  e2");
			//BG_hist->Draw("same hist");

			if(blind_option!="blind")Data_hist->Draw("e same");


			TText *t = new TText();
			t->SetTextAlign(22);
			t->SetTextColorAlpha(kGray, 0.50);
			t->SetTextFont(40);
			t->SetTextSize(0.25);
			t->SetTextAngle(25);


			
			
			legend->SetFillStyle(0);
			legend->SetLineWidth(0);
			legend->Draw("same ");

			//Labels
			double x_top_label = 0.90;

			TPaveText * CMS_Internal = new TPaveText(0.10,x_top_label,0.288191,x_top_label+0.1,"NDC");
			/*lumi->SetFillColor(0); */
			CMS_Internal->SetFillStyle(4050);
			CMS_Internal->SetLineColor(0);
			CMS_Internal->SetTextFont(42); CMS_Internal->SetTextSize(0.0599401); CMS_Internal->SetBorderSize(0);
			CMS_Internal->AddText(Form("CMS Internal"));
			CMS_Internal->Draw();

			TPaveText * lumi = new TPaveText(0.20,x_top_label-0.22,0.388191,x_top_label-0.05,"NDC");
			/*lumi->SetFillColor(0); */
			lumi->SetFillStyle(4050);
			lumi->SetLineColor(0);
			lumi->SetTextFont(42); lumi->SetTextSize(0.0599401); lumi->SetBorderSize(0);
			lumi->AddText(Form("#sqrt{s}= 13 TeV, %1.1f pb^{-1}",lumi_eraH));
			//lumi->AddText("#mu+jets channel, ".append(CutString));
			lumi->AddText("#mu+jets channel"+cut);
			lumi->Draw();

			TPaveText * Warning = new TPaveText(0.20,x_top_label-0.26,0.388191,x_top_label-0.16,"NDC");
			Warning->SetFillStyle(4050);
			Warning->SetLineColor(0);
			Warning->SetTextFont(42); Warning->SetTextSize(0.052); Warning->SetBorderSize(0);
			Warning->AddText(Form("Extra %3.1f factor",Norm_factor));
			//Warning->Draw();

			t->DrawTextNDC(.5, .53, "Preliminary");

			//Add an arrow to indicate the intact proton in the case the option is selected
			///////////////////////////////////////////////////////////////////////////////
			if(arrow_option=="arrow"){
				float x_min_plot = hs->GetXaxis()->GetXmin();
				float x_max_plot = hs->GetXaxis()->GetXmax();

				float center_arrow_x = x_min_plot+ (x_max_plot-x_min_plot)*0.81;
				float center_arrow_y = max_display*0.65;
				float size_x_arrow = (x_max_plot-x_min_plot)*0.1;


				float x_min_arrow = 0.0;
				float x_max_arrow = 0.0;

				if(k==0){ //j gives the direction of the intact proton
					x_min_arrow = (center_arrow_x-size_x_arrow/2.);
					x_max_arrow = (center_arrow_x+size_x_arrow/2.);
				}else{
					x_min_arrow = (center_arrow_x+size_x_arrow/2.);
					x_max_arrow = (center_arrow_x-size_x_arrow/2.);
				}

				TArrow *ar2 = new TArrow(x_min_arrow,center_arrow_y,x_max_arrow,center_arrow_y,0.02,">");
				ar2->SetAngle(40);
				ar2->SetLineWidth(7);
				ar2->SetLineColor(12);
				ar2->SetFillColor(12);
				ar2->Draw();

				TPaveText * Proton_text = new TPaveText(0.60,0.48,0.9,0.5,"NDC");
				Proton_text->SetFillStyle(4050);
				Proton_text->SetLineColor(0);
				Proton_text->SetTextFont(42); Proton_text->SetTextSize(0.052); Proton_text->SetBorderSize(0);
				Proton_text->AddText("Intact Proton dir.");
				Proton_text->Draw();
			}
			///////////////////////////////////////////////////////////////////////////////



			// lower plot will be in pad
			///////////////////////////////////////////////////////////////////////////////
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


			ratio_hist_uncertainty->SetMaximum(2.0);
			ratio_hist_uncertainty->SetMinimum(0.0);
			ratio_hist_uncertainty->Draw("e2");
			if(blind_option!="blind")ratio_hist->Draw("ep same");
			pad2->Update();


			// Ratio plot (h3) settings
			ratio_hist_uncertainty->SetTitle(""); // Remove the ratio title

			// Y axis ratio plot settings
			ratio_hist_uncertainty->GetYaxis()->SetTitle("Data/MC ratio");
			ratio_hist_uncertainty->GetXaxis()->SetTitle(xAxis_label);
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
			///////////////////////////////////////////////////////////////////////////////


			cancG0->SaveAs("pdfs/"+label+cut_string+".pdf");
			cancG0->SaveAs("pngs/"+label+cut_string+".png");
			cancG0->SaveAs("root/essai.root");



		}
	}


	//gApplication->Terminate();

	return 0;
}
