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
 * For all specified files, calculate the time resolution and compare
 * with beam energy
 */

/* argv = {file, energy in GeV}*
 * argc = odd
 */

#define TREENAME "tree"
#define NCHANS 4
#define MAINCHAN 0 // main mcp (photonis) channel
#define REFCHAN 1 // reference mcp (photek) channel
#define SCALEFACTOR 1 // conversion from 200 mV*ps -> fC, assuming 50 ohms
#define NRGYERR 0
#define BINNING 50

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

    // make a plot for time value by beam energy
    TGraphErrors *TvsE = new TGraphErrors(nDataPoints);

    // make histograms to hold time info for each file 
    TH1F **dtHist = new TH1F*[nDataPoints];
    for (int ii = 0; ii < nDataPoints; ii++)
        dtHist[ii] = new TH1F(("DT_" + std::to_string(int(beamE[ii])) + "GeV").c_str(),
                               ("Time resolution " + std::to_string(ii+1)).c_str(),
                   /* nbins  */BINNING,
                   /* DT min */-1.2,
                   /* DT max */ -.4);
    float timeval[NCHANS];
    float amplitudes[NCHANS];
    int quality[NCHANS];

    // Open files, open trees and calculate average charge
    for (int iFile = 0; iFile < nDataPoints; iFile++) {
        TFile *file = new TFile(argv[1+2*iFile]);
        TTree *tree = (TTree *) file->Get(TREENAME);
        tree->SetBranchAddress("tgausroot", &timeval);
        tree->SetBranchAddress("Amplitude", &amplitudes);
        tree->SetBranchAddress("QualityBit", &quality);

        // Read all entries and extract values to fill the histograms
        Long64_t nentries = tree->GetEntries();
        for (Long64_t iEntry = 0; iEntry < nentries; iEntry++) {
            tree->GetEntry(iEntry);
            // We perform a cut few cuts in addition to the quality bit to remove observed background
            if (!(quality[MAINCHAN] || quality[REFCHAN]) &&
                    (amplitudes[REFCHAN] > 0.12 && amplitudes[REFCHAN] < 0.49 && amplitudes[MAINCHAN] > 0.02)) {
                
                dtHist[iFile]->Fill(timeval[MAINCHAN] - timeval[REFCHAN]);
            }
        }
        // check number of data points and fit limits
        std::cout << "Histogram " << iFile << " # Data Points: " << dtHist[iFile]->GetEntries() << std::endl;
        std::cout << "Mean, RMS: " << dtHist[iFile]->GetMean() << ", " << dtHist[iFile]->GetRMS() << std::endl;

        // fit histogram with gaussian
        TF1 *DTGausPeak = new TF1("intpeak", "gaus",
                dtHist[iFile]->GetMean() - 1.7*dtHist[iFile]->GetRMS(),
                dtHist[iFile]->GetMean() +     dtHist[iFile]->GetRMS()); // tails are skewed
        dtHist[iFile]->Fit(DTGausPeak, "LMQR");

        // Place into graph
        TvsE->SetPoint(iFile, beamE[iFile], DTGausPeak->GetParameter(2) * SCALEFACTOR);
        TvsE->SetPointError(iFile, NRGYERR, DTGausPeak->GetParError(2) * SCALEFACTOR);

        // Print data points
        std::cout << beamE[iFile] << ", " << DTGausPeak->GetParameter(2) * SCALEFACTOR << std::endl;

        file->Close();
        delete file, DTGausPeak;
    }
    // Sort graph by x-axis values
    TvsE->Sort();

    // Save graph
    TFile *file = new TFile("timeres_vs_energy.root", "RECREATE");
    TvsE->Write();

    for (int ii = 0; ii < argc/2; ii++) {
        dtHist[ii]->Write(); // save histogram+fits for debugging
        delete dtHist[ii];
    }
    file->Close();
    delete beamE, dtHist, TvsE, file;
}

