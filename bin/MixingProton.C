#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <TDirectory.h>
#include <TBranch.h>
#include <TRandom3.h>
#include <TMath.h>
#include "TSystem.h"

#include "TopLJets2015/TopAnalysis/interface/protonTrackRatios.h"
#include "TopLJets2015/TopAnalysis/interface/PPSEff.h"

#include <time.h>
using namespace std;

int main(int argc, char *argv[])
{
    // CMSSW settings
    const char *CMSSW_BASE = getenv("CMSSW_BASE");
    TString data_path = Form("%s/src/TopLJets2015/TopAnalysis/data/era2017/", CMSSW_BASE);

    cout << data_path << endl;

    // Parameters
    int rndSeed = 1234567890;
    int nEventsToMix = -1;
    int nMCEventsToSkip = 0;

    TRandom3 *rand_gen = new TRandom3(rndSeed);

    string inMCFileName = argv[1];
    TString name_mu = argv[2];
    string outFileName = inMCFileName.substr(inMCFileName.find_last_of('/') + 1, inMCFileName.find_last_of('.') - inMCFileName.find_last_of('/') - 1) + "_enriched.root";
    bool isSignal = TString(inMCFileName.c_str()).Contains("signal") || TString(inMCFileName.c_str()).Contains("Signal");

    if (argc >= 4)
    {
        TString arg3 = TString(argv[3]);
        if (arg3.IsDec())
            nEventsToMix = arg3.Atoi();
        else
            cout << "Argument 3 must be an integer number. All input events will be mixed." << endl;
    }

    cout << "Check config: Input MC File " << inMCFileName << " Data File: " << name_mu << " output: " << outFileName << endl;

    // create proton pools Prefilter the proton tree to remove no protons events
    TChain *_ch = new TChain("protons");
    _ch->Add(name_mu);
    TTree *PUpr = (_ch->CopyTree("p1_xi>0 || p2_xi>0.0"))->CloneTree(); //_ch->CloneTree();

    // variables to be used to read from PU pools
    unsigned int poll_run;
    float poll_p1_xi, poll_p2_xi;
    PUpr->SetBranchAddress("run", &poll_run);
    PUpr->SetBranchAddress("p1_xi", &poll_p1_xi);
    PUpr->SetBranchAddress("p2_xi", &poll_p2_xi);
    cout << "done " << endl;

    // read MC file
    TFile *oldfile = new TFile(inMCFileName.c_str());
    TTree *chMCEvents = (TTree *)oldfile->Get("tree");
    cout << "Store MC nvtx distribution " << endl;
    TH1F *mc_pu = new TH1F("mc_pu", ";nvtx;w", 100, 0, 100);
    chMCEvents->Draw("nvtx>>mc_pu", "", "norm");
    TH1F *evt_count = (TH1F *)oldfile->Get("evt_count"); // event counter with SumWeights

    //list of branches to update:
    unsigned int run;
    int nvtx;
    float beamXangle, p1_xi, p2_xi, weight, ppsSF_wgt, ppsSF_wgt_err, pu_wgt, ptag_wgt, ptag_wgt_err;
    float p1_x = 0, p1_y = 0, p2_x = 0, p2_y = 0;
    chMCEvents->SetBranchAddress("run", &run);
    chMCEvents->SetBranchAddress("beamXangle", &beamXangle);
    chMCEvents->SetBranchAddress("p1_xi", &p1_xi);
    chMCEvents->SetBranchAddress("p2_xi", &p2_xi);
    chMCEvents->SetBranchAddress("p1_x", &p1_x);
    chMCEvents->SetBranchAddress("p1_y", &p1_y);
    chMCEvents->SetBranchAddress("p2_x", &p2_x);
    chMCEvents->SetBranchAddress("p2_y", &p2_y);
    chMCEvents->SetBranchAddress("nvtx", &nvtx);
    chMCEvents->SetBranchAddress("weight", &weight);
    chMCEvents->SetBranchAddress("pu_wgt", &pu_wgt);
    chMCEvents->SetBranchAddress("ppsSF_wgt", &ppsSF_wgt);
    chMCEvents->SetBranchAddress("ppsSF_wgt_err", &ppsSF_wgt_err);
    chMCEvents->SetBranchAddress("ptag_wgt", &ptag_wgt);
    chMCEvents->SetBranchAddress("ptag_wgt_err", &ptag_wgt_err);
    chMCEvents->SetBranchAddress("nvtx", &nvtx);
    chMCEvents->SetBranchStatus("*", 1); // activate all branches to copy

    int nMCEntries = chMCEvents->GetEntries();
    cout << nMCEntries << " events read from file(s) " << inMCFileName << endl;
    if (nEventsToMix == -1)
        nEventsToMix = nMCEntries;

    // Proton efficiency class
    PPSEff *MultiRP_eff = new PPSEff(Form("%s/pixelEfficiencies_multiRP.root", data_path.Data()));
    PPSEff *Strip_eff = new PPSEff(Form("%s/PreliminaryEfficiencies_March302021_1D2DMultiTrack.root", data_path.Data()));

    // Get pileup proton ratios
    protonTrackRatios ptr;
    int nLines = ptr.readFromFile(Form("%s/protonRatios_2017.dat", data_path.Data()));
    if (nLines <= 0)
    {
        cout << "No true-zero-track ratio read from file! No mixing performed." << endl;
        return 4;
    }

    //Probability of proton config in PPS
    //no proton, proton in arm 1, proton in arm 2, proton in both
    //To calculate from Zero bias events
    double proton_config_probability[4] = {0.25, 0.25, 0.25, 0.25};
    double proton_config_probability_error[4] = {0.25, 0.25, 0.25, 0.25};

    // Read 2 proton hits and total number of events
    int n_p2_m = int(((TH1F *)TFile::Open(name_mu)->Get("pn_count"))->GetBinContent(4)); //ERROR HERE in the original code I think
    int n_m = int(((TH1F *)TFile::Open(name_mu)->Get("evt_count"))->GetBinContent(3));
    int n_p2 = (n_p2_m), n = (n_m);

    // Get systematics for n_proton==2 event fraction from sub-selection of nBjet>0 events:
    int n_m_sys = int(((TH1F *)TFile::Open(name_mu)->Get("pn_count"))->GetBinContent(5));
    TChain *_ch_protons = new TChain("protons");
    _ch_protons->Add(name_mu);
    int n_sys = (n_m_sys), n_p2_sys = _ch_protons->GetEntries("nBjets>0");

    // Get 1 proton hits (depend on the arm)
    int n_p1_mRP0 = int(((TH1F *)TFile::Open(name_mu)->Get("pn_count"))->GetBinContent(2));
    int n_p1_mRP1 = int(((TH1F *)TFile::Open(name_mu)->Get("pn_count"))->GetBinContent(3));
    int n_p1_RP0 = (n_p1_mRP0), n_p1_RP1 = (n_p1_mRP1);

    // events with exactly 0 pu tracks
    int n_p0 = n - n_p1_RP0 - n_p1_RP1 - n_p2;

    int n_p0_bis = int(((TH1F *)TFile::Open(name_mu)->Get("pn_count"))->GetBinContent(1));
    cout<<"Check probabilities"<<endl;
    cout<<n_p0<<" "<<n_p0_bis<<endl;

    //Contruct PPS configurations probabilities

    // probabilities for 0 track in signal events
    proton_config_probability[0] = (n_p0) / float(n); // probability of 0 tracks in both arms
    proton_config_probability_error[0] = 0.95;        // 5% flat

    // probabilities for 1 track in signal events
    proton_config_probability[1] = n_p1_RP0 / float(n); // probability of 0 tracks in RP0
    proton_config_probability_error[1] = 0.95;          // 5% flat

    proton_config_probability[2] = n_p1_RP1 / float(n); // probability of 0 tracks in RP1
    proton_config_probability_error[2] = 0.95;          // 5% flat

    //Two protons events
    proton_config_probability[3] = n_p2 / float(n); // probability of 2 tracks
    proton_config_probability_error[3] = 0.95;      //(n_p2_sys / float(n_sys)) / norm_weight[i_era * n_xa + i_xa];

    cout << "Check probabilities" << endl;
    cout << "0p: " << proton_config_probability[0] << " 1pRP1: " << proton_config_probability[1] << " 1pRP2: " << proton_config_probability[2] << " 2p: " << proton_config_probability[3] << endl;
    cout << "Sum proba. " << (proton_config_probability[0] + proton_config_probability[1] + proton_config_probability[2] + proton_config_probability[3]) << endl;
    cout << "Is Signal :" << isSignal << endl;
    // Open output mixed file
    TFile *fOutMixed = new TFile(outFileName.c_str(), "RECREATE");
    if (!fOutMixed->IsOpen())
    {
        cout << "Error opening output file mixedPUProtons. Abort." << endl;
        return 5;
    }

    // Create new output with updated PPS values:
    TTree *tMCMixed = chMCEvents->CloneTree(0);

    int signal_protons = 0; // additional variables used for signal to indecate if all/part/none of the protons are from pileup
    tMCMixed->Branch("signal_protons", &signal_protons);

    int times, timed;
    times = time(NULL);
    cout << "Loop over MC entries and add pileup protons" << endl;
    // Loop over MC entries and mix pileup protons
    int iMCEntry;
    for (iMCEntry = 0; (iMCEntry < nEventsToMix) && (iMCEntry + nMCEventsToSkip < nMCEntries); iMCEntry++)
    {
        if (iMCEntry % 1000 == 0)
            printf("\r [%3.0f%%] done", 100. * (float)iMCEntry / (float)nEventsToMix);

        //Some events in the proton pool have
        int i_event = rand_gen->Rndm() * PUpr->GetEntries();
        ;
        float p1_xi_in_pool = 0.0;
        float p2_xi_in_pool = 0.0;
        while (p1_xi_in_pool == 0.0 || p2_xi_in_pool == 0.0)
        {
            //cout << i_event << endl;
            PUpr->GetEntry(i_event);
            i_event++;
            p1_xi_in_pool = (p1_xi_in_pool > 0.0) ? p1_xi_in_pool : poll_p1_xi;
            p2_xi_in_pool = (p2_xi_in_pool > 0.0) ? p2_xi_in_pool : poll_p2_xi;
        }

        //cout << "Sampled Xi " << p1_xi_in_pool << " " << p2_xi_in_pool << endl;

        // asign the protons to the MC event
        chMCEvents->GetEntry(iMCEntry + nMCEventsToSkip);

        // proton efficiency implementation (xi of protons that fail reco. will be set to zero)
        double ppsSF_wgt = 1.;
        double ppsSF_wgt_err = 0;

        if (isSignal)
        {
            if (p1_xi > 0)
            {
                float SF = Strip_eff->getEff(p1_x, p1_y, 0, run) * MultiRP_eff->getEff(p1_x, p1_y, 0, run);
                if (rand_gen->Rndm() > SF) p1_xi = 0;
                else ppsSF_wgt *= SF;
                ppsSF_wgt_err += MultiRP_eff->getRelEffErrSq(p1_x, p1_y, 0, run);
            }
            if (p2_xi > 0)
            {
                float SF = Strip_eff->getEff(p2_x, p2_y, 1, run) * MultiRP_eff->getEff(p2_x, p2_y, 1, run);
                if (rand_gen->Rndm() > SF) p2_xi = 0;
                else ppsSF_wgt *= SF;
                ppsSF_wgt_err += MultiRP_eff->getRelEffErrSq(p2_x, p2_y, 1, run);
            }
            ppsSF_wgt_err = sqrt(ppsSF_wgt_err);
        }

        //To finish the signal part
        if (isSignal && p1_xi == 0 && p2_xi > 0)
        { // one signal proton in arm1 + 1 PU proton in each arm
            p1_xi = p1_xi_in_pool;
            p2_xi = 0.0;
            ptag_wgt = proton_config_probability[3];
            ptag_wgt *= ptr.trueZeroTracksRatio(run, beamXangle, 1);
            ptag_wgt_err = proton_config_probability_error[3];
            weight *= ptag_wgt;
            signal_protons = 01;
        }
        else if (isSignal && p1_xi > 0 && p2_xi == 0)
        { // one signal proton in arm0 + 1 PU proton in each arm
            p1_xi = 0.0;
            p2_xi = p2_xi_in_pool;
            ptag_wgt = proton_config_probability[3];
            ptag_wgt *= ptr.trueZeroTracksRatio(run, beamXangle, 0);
            ptag_wgt_err = proton_config_probability_error[3];
            weight *= ptag_wgt;
            signal_protons = 10;
        }

        //Process Background (ttX) events
        else if (!isSignal)
        {
            // BG event with no signal protons
            // First select one arm of PPS
            int arm_PPS = rand_gen->Rndm() > 0.5 ? 1 : 2;
            if (arm_PPS == 1){
                p1_xi = p1_xi_in_pool;
                p2_xi = 0.0;
            }
            if (arm_PPS == 2){
                p1_xi = 0.0;
                p2_xi = p2_xi_in_pool;
            }
            ptag_wgt = proton_config_probability[arm_PPS];
            //ptag_wgt *= ptr.trueZeroTracksRatio(run, beamXangle, 0);
            ptag_wgt_err = proton_config_probability_error[arm_PPS];

            weight *= ptag_wgt;

            //cout << arm_PPS << " " << ptag_wgt << " " << p1_xi << " " << p2_xi << endl;

            signal_protons = arm_PPS == 1 ? 10 : 01;
        }

        tMCMixed->Fill();
    }

    cout << "\nDone. " << iMCEntry << " events processed." << endl;
    timed = time(NULL);
    times = timed - times;

    cout << "time from start to end = " << times << endl;

    // Write mixed tree to file and close everything
    fOutMixed->cd();
    cout << "Writes " << fOutMixed->GetName() << endl;
    tMCMixed->Write();
    /*if (nMCEventsToSkip == 0)
        evt_count->Write();*/
    fOutMixed->Close();

    return iMCEntry;
}
