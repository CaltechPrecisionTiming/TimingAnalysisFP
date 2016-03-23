#include <iostream>
#include <fstream>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TF1.h"
#include <stdexcept>

/*
 * For all specified files, calculate the average charge deposit and compare
 * with beam energy
 */

/* argv = {file, energy in GeV}*
 * argc = odd
 */

#define TREENAME "tree"
#define NCHANS 4
#define MAINCHAN 0 // main mcp (photonis) channel
#define REFCHAN 1 // reference mcp (photek) channel
#define BINNING 50
#define SCALEFACTOR 4. // conversion from 200 mV*ps -> fC, assuming 50 ohms
#define ATTENUATION 3.16227766 // 10 dB
#define ATTENERGY 8 // Energy above which signal was attenuated
#define NRGYERR 0

int main (int argc, char **argv) {
    // usage
    if (argc == 1)
        std::cout << "usage: chargeVSenergy {file, energy in GeV}*" << std::endl;
    // argc must be odd
    if (argc % 2 != 1) {
        std::cout << "Wrong number of inputs" << std::endl;
        return 0;
    }

    int nDataPoints = argc / 2;

    // interpret beam energy values
    float *beamE = new float[nDataPoints];
    try {
        for (int ii = 0; ii < nDataPoints; ii++)
            beamE[ii] = std::stof(argv[2+2*ii]);
    }
    catch (const std::invalid_argument& ia) {
        std::cout << "Cannot understand energy values" << std::endl;
        return 0;
    }

    // Check existence of all files
    for (int ii = 0; ii < nDataPoints; ii++) {
        ifstream f(argv[1+2*ii]);
        if (f.good())
            f.close();
        else {
            f.close();
            std::cout << "Bad input file: " << argv[1+2*ii] << std::endl;
            return 0;
        }
    }

    // make a plot for signal int by beam energy
    TGraphErrors *QvsE = new TGraphErrors(nDataPoints);

    // make histograms to hold charge info for each file 
    TH1F **intHist = new TH1F*[nDataPoints];
    for (int ii = 0; ii < nDataPoints; ii++)
        intHist[ii] = new TH1F(("Int_Q_" + std::to_string(int(beamE[ii])) + "GeV").c_str(),
                               ("Signal Int " + std::to_string(ii+1)).c_str(),
                    /* nbins */BINNING,
                    /* Q min */0,
                    /* Q max */25);
    float intValue[NCHANS];
    int quality[NCHANS];

    // Open files, open trees and calculate average charge
    for (int iFile = 0; iFile < nDataPoints; iFile++) {
        TFile *file = new TFile(argv[1+2*iFile]);
        TTree *tree = (TTree *) file->Get(TREENAME);
        tree->SetBranchAddress("Int", &intValue[MAINCHAN]);
        tree->SetBranchAddress("QualityBit", &quality);

        // Read all entries and extract values to fill the histograms
        Long64_t nentries = tree->GetEntries();
        for (Long64_t iEntry = 0; iEntry < nentries; iEntry++) {
            tree->GetEntry(iEntry);
            float factor = SCALEFACTOR * (beamE[iFile] > ATTENERGY ? ATTENUATION : 1);
            if (!quality[MAINCHAN] && !(quality[REFCHAN] & 0x18)) { // ref, only care if exactly 1 pulse
                intHist[iFile]->Fill(intValue[MAINCHAN] * factor);
            }
        }
        // check number of data points and fit limits
        std::cout << "Histogram " << iFile << " # Data Points: " << intHist[iFile]->GetEntries() << std::endl;
        std::cout << "Mean, RMS: " << intHist[iFile]->GetMean() << ", " << intHist[iFile]->GetRMS() << std::endl;

        // calculate average signal integration
        TF1 *IntGausPeak = new TF1("intpeak", "gaus",
                intHist[iFile]->GetMean() -     intHist[iFile]->GetRMS(),
                intHist[iFile]->GetMean() + 1.5*intHist[iFile]->GetRMS());//neg skewed
        intHist[iFile]->Fit(IntGausPeak, "LMQR");

        // Place into graph
        QvsE->SetPoint(iFile, beamE[iFile], IntGausPeak->GetParameter(1));
        QvsE->SetPointError(iFile, NRGYERR, IntGausPeak->GetParameter(2));

        // Print data points
        std::cout << beamE[iFile] << ", " << IntGausPeak->GetParameter(1) << std::endl;

        file->Close();
        delete file, IntGausPeak;
    }
    // Sort graph by x-axis values
    QvsE->Sort();

    // Save graph
    TFile *file = new TFile("charge_vs_energy.root", "RECREATE");
    QvsE->Write();

    for (int ii = 0; ii < argc/2; ii++) {
        intHist[ii]->Write(); // save histograms+fits for debugging
        delete intHist[ii];
    }
    file->Close();
    delete beamE, intHist, QvsE, file;
}

