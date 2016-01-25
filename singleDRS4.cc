#include <iostream>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "Math/Interpolator.h"

#define NSAMPLES 1024
#define NCHANS 4

#define SPLINE 1
#define MINPULSE 0.02
#define INTRANGE 4

enum PulseQuality {
    kNegativePolarity         = 0x000001, // bit 0
    kSuddenJump               = 0x000002, // bit 1
    kFlatTop                  = 0x000004, // bit 2
    kSecondPulse              = 0x000008, // bit 3
    kNoPulse                  = 0x000010, // bit 4
    kLargeNegativeAmplitude   = 0x000020, // bit 5
    kSaturated                = 0x000040  // bit 6
};

TH1F *InterpolateWaveform(int nsamples, float *outputwaveform, float *inputwaveform, int splineBinFactor, std::string name);
int FindMin(int n, int splineBinFactor, float *a);
int FindMax(int n, int splineBinFactor, float *a);
float LinearFit_Baseline(TH1F *pulse, const int index_min, const int range);
float ChannelIntegral(float *a, int peak);
unsigned int CheckPulseQuality(int binMin, int binMax, float *a, float offset);
float GausFit_MeanTime(TH1F *pulse, const int index_first, const int index_last);
/*
int FindRisingEdge(int n, int binMax, float *a);
int FindFirstPulsePeak(int n, float *a);
float LED(TH1F *pulse, double threshold, int nsamples, int splineBinFactor);
void FitRisingEdge(TH1F* pulse, int nbinsL, int nbinsH, float &THM, float &risetime, float base);
*/

int main (int argc, char **argv) {
    TFile *f;

    if (argc >= 3) {
            f = new TFile(argv[1]);
            std::cout << ">> Opening file " << argv[1] << " ......" << std::endl;
            // terminate if the file can't be opened
            if (!f->IsOpen()) {
                std::cerr << "!! File open error:" << argv[1] << std::endl;
                return 1;
            }
    }
    // terminate if there is no input file or more than 1 input file
    else {
        std::cerr << "!! No input file" << std::endl;
        return 1;
    }
    bool includePulseshapeInOutput = true;
    if (argc >= 4) includePulseshapeInOutput = bool(atoi(argv[3]));

    // Create the output file with a TTree
    TFile *fout;
    if (strncmp(argv[2], "same", 5) == 0) {
        std::string fn(argv[1]);
        int pf = fn.find("_pulse.root");
        int pi = fn.rfind("/") + 1;
        fn = fn.substr(pi, pf-pi) + "_anal.root";
        std::cout << "fname: " << fn << std::endl;
        fout = new TFile(fn.c_str(), "recreate");
    }
    else {
        fout = new TFile(argv[2], "recreate");
    }

    // time array 0-1024
    int t_[NSAMPLES];
    for (int i = 0; i < 1024; i++) t_[i] = i;
    // Create Waveform arrays
    float VoltagesRaw_[NCHANS][NSAMPLES];
    float Voltages_[NCHANS][NSAMPLES*SPLINE];

    TTree *t1 = (TTree*)f->Get("T");
    for (int a = 0; a < NCHANS; a++)
            t1->SetBranchAddress(("c" + std::to_string(a+1)).c_str(), VoltagesRaw_[a]);


    // Create Array of histogram pointers
    TH1F *CHPulsesRaw[NCHANS];
    TH1F *CHPulses[NCHANS];

    /*
    TH1F *CH1Amp = new TH1F("CH1Amp", "CH1Amp", 40, -0.3, -0.6);
    TH1F *CH2Amp = new TH1F("CH2Amp", "CH2Amp", 40, -0.3, -0.6);
    TH1F *GausPeak_CH12_dt  = new TH1F("GausPeak_CH12_dt", "GausPeak_CH12_dt; t1-t2 [ns]; Events", 4000, -4, 4);
    TH1F *GausPeak_CH34_dt  = new TH1F("GausPeak_CH34_dt", "GausPeak_CH34_dt; t3-t4 [ns]; Events", 4000, -4, 4);
    TH1F *GausPeak_TOF_CH13 = new TH1F("GausPeak_TOF_CH13", "GausPeak_TOF_CH13; t1-t3 [ns]; Events", 4000, -4, 4);
    TH1F *GausPeak_TOF_CH14 = new TH1F("GausPeak_TOF_CH14", "GausPeak_TOF_CH14; t1-t4 [ns]; Events", 4000, -4, 4);
    */

    TTree *treeOut = new TTree("tree", "tree");
    
    unsigned int eventNumber = 0;
    float time_gausfitroot[NCHANS];
    float amplitude[NCHANS];
    float thm[NCHANS];
    float risetime[NCHANS];
    float tff[NCHANS];
    float tff_v2[NCHANS];
    float bl[NCHANS];
    float aff[NCHANS];
    float integral[NCHANS];
    unsigned int QualityBit[NCHANS];
    
    float chisq[NCHANS];
    for (int a = 0; a < NCHANS; a++)
            chisq[a] = -1;
    
    treeOut->Branch("event", &eventNumber, "event/i");
    treeOut->Branch("tgausroot", &time_gausfitroot, TString::Format("tgausroot[%d]/F", NCHANS));
    treeOut->Branch("Amplitude", &amplitude, TString::Format("Amplitude[%d]/F", NCHANS));
    treeOut->Branch("THM", &thm, TString::Format("THM[%d]/F", NCHANS));
    treeOut->Branch("Risetime", &risetime, TString::Format("Risetime[%d]/F", NCHANS));
    treeOut->Branch("BL", &bl, TString::Format("BL[%d]/F", NCHANS));
    treeOut->Branch("TFF", &tff, TString::Format("TFF[%d]/F", NCHANS));
    treeOut->Branch("TFF_v2", &tff_v2, TString::Format("TFF_v2[%d]/F", NCHANS));
    treeOut->Branch("AFF", &aff, TString::Format("AFF[%d]/F", NCHANS));
    treeOut->Branch("QualityBit", &QualityBit, TString::Format("QualityBit[%d]/i", NCHANS));
    treeOut->Branch("Int", &integral, TString::Format("Int[%d]/F", NCHANS));
    treeOut->Branch("chisq", &chisq, TString::Format("chisq[%d]/F", NCHANS));
    
    // Initialize raw pulse histograms to empty
    for (int a = 0; a < NCHANS; a++)
            CHPulsesRaw[a] = new TH1F(("PulseRawCh" + std::to_string(a)).c_str(), \
                ("PulseRawCh" + std::to_string(a)).c_str(), NSAMPLES, 0, NSAMPLES);

    if (includePulseshapeInOutput) {
        treeOut->Branch("chnls", VoltagesRaw_, TString::Format("chnls[%d][%d]/F", NCHANS, NSAMPLES));
        treeOut->Branch("t", t_, TString::Format("t[%d]/I", NSAMPLES));
    }

    // Read all entries and fill the histograms
    Long64_t nentries = t1->GetEntries();
    for (Long64_t iEntry = 0; iEntry < nentries; iEntry++) {
        if (iEntry % 100 == 0)
            std::cout << "Processing Event: " << iEntry << " out of: " << nentries << std::endl;
        
        t1->GetEntry(iEntry);
        eventNumber = iEntry + 1;

        // Convert to Volts and flip signal
        for (int a = 0; a < NCHANS; a++)
                for (int c = 0; c < NSAMPLES; c++)
                    VoltagesRaw_[a][c] *= -0.001;

        // Initialize histograms
        for (int a = 0; a < NCHANS; a++)
                for (int c = 0; c < NSAMPLES; c++) {
                    CHPulsesRaw[a]->SetBinContent(c+1, VoltagesRaw_[a][c]); // change sign over here
                    CHPulsesRaw[a]->SetBinError(c+1, 0.001);
                }

        // Do spline to add more points to the waveform
        for (int a = 0; a < NCHANS; a++)
                CHPulses[a] = InterpolateWaveform(NSAMPLES, Voltages_[a], VoltagesRaw_[a], SPLINE,\
                    "Pulse_CH" + std::to_string(a));

        // Find Min and Max of the Channel data (Voltage)
        int index_min[NCHANS];
        int index_max[NCHANS];
        for (int a = 0; a < NCHANS; a++) {
                index_min[a] = FindMin(NSAMPLES, SPLINE, Voltages_[a]); // NB: pulse is already flipped
                index_max[a] = FindMax(NSAMPLES, SPLINE, Voltages_[a]);
            }

        /*
        // Find the rising edge on CH1
        int fbin1[NCHANS];
        for (int a = 0; a < NCHANS; a++)
                fbin1[a] = FindRisingEdge(NSAMPLES, index_max[a], Voltages_[a]);
        */

        /*
        // For the first version of Photonis data the pulse has many peaks: find the first one
        // this is useful ONLY for Photonis MCP
        int index_firstPulse[NCHANS];
        for (int a = 0; a < NCHANS; a++)
                index_firstPulse[a] = FindFirstPulsePeak(NSAMPLES, Voltages_[a]);
        */

        ////////////////////////
        // Done with setup, start the fits
        ////////////////////////
        
        // Fit the baseline
        float base[NCHANS];
        // If not initialized, defaults to 0.0
        for (int a = 0; a < NCHANS; a++)
                base[a] = LinearFit_Baseline(CHPulses[a], index_min[a], 10);
        
        // Find the amplitudes
        for (int a = 0; a < NCHANS; a++)
                amplitude[a] = Voltages_[a][(index_max[a])] - base[a];

        for (int a = 0; a < NCHANS; a++)
                integral[a] = ChannelIntegral(Voltages_[a], index_max[a]) - (2*INTRANGE+1) * base[a];

        // Find the quality of the pulse
        for (int a = 0; a < NCHANS; a++)
                QualityBit[a] = CheckPulseQuality(index_min[a], index_max[a], VoltagesRaw_[a], base[a]);
    
        ///////////////////
        // Gaussian fit of peak
        ///////////////////

        float timepeak[NCHANS];
        for (int a = 0; a < NCHANS; a++) {
                timepeak[a] = GausFit_MeanTime(CHPulses[a],\
                    index_max[a] - 4*SPLINE, index_max[a] + 5*SPLINE);
                // for CFD
                /* timepeak[a] = LED(CHPulses[a], 0.50 * amplitudes[a], NSAMPLES, SPLINE); */
                time_gausfitroot[a] = timepeak[a] * 0.2 / SPLINE; // 0.2 conversion factor
            }

        // Fit Rising Edge
        /*
        for (int a = 0; a < NCHANS; a++)
                FitRisingEdge(CHPulses, -1, 0, thm[a], risetime[a], base[a]);
        */
        
        //Fill the tree
        treeOut->Fill();

        if (iEntry < nentries-1)
            for (int a = 0; a < NCHANS; a++)
                    delete CHPulses[a];

    }
    /*
    // Save last pulses as examples
    for (int a = 0; a < NCHANS; a++) {
            CHPulses[a]->Write();
            CHPulsesRaw[a]->Write();
            delete CHPulses[a], CHPulsesRaw[a];
        }
    */
    /* 
    GausPeak_CH12_dt->Write();
    GausPeak_CH34_dt->Write();
    GausPeak_TOF_CH13->Write();
    GausPeak_TOF_CH14->Write();
    CH1Amp->Write();
    CH2Amp->Write();
    */

    treeOut->Write();

    fout->Close();
}

/*
 * Performa a spline to interpolate
 */
TH1F *InterpolateWaveform(int nsamples, float *outputwaveform, float *inputwaveform, int splineBinFactor, std::string name) {

    ROOT::Math::Interpolator cspline(nsamples, ROOT::Math::Interpolation::kCSPLINE);    
    TH1F *pulse = new TH1F(name.c_str(),name.c_str(),nsamples * splineBinFactor,0,nsamples*splineBinFactor);
    
    if (splineBinFactor != 1) {
        double bins[nsamples];
        double w[nsamples];
        for (int i = 0; i < nsamples; i++) { 
            bins[i] = i + 1;
            w[i] = inputwaveform[i];
        }

        cspline.SetData(nsamples, bins, w);

        // Do spline
        for (int i = 0; i < nsamples * splineBinFactor; i++) {
            if (i > splineBinFactor)
                outputwaveform[i] = cspline.Eval( pulse->GetXaxis()->GetBinCenter(i+1) / splineBinFactor );
            else
                outputwaveform[i] = 0;
        }
        for (int i = 1; i <= pulse->GetXaxis()->GetNbins(); i++) {
            pulse->SetBinContent(i, outputwaveform[i-1]);
            pulse->SetBinError(i, 0.001);        
        }
    }

    else {
        for (int ii = 0; ii < nsamples; ii++) {
            outputwaveform[ii] = inputwaveform[ii];
        }
        for (int ii = 0; ii < nsamples; ii++) { 
            pulse->SetBinContent(ii+1, outputwaveform[ii]);
            pulse->SetBinError(ii+1, 0.001);
        }
    }

    return pulse;
}

/*
 * Find minimum of the pulse
 */
int FindMin(int n, int splineBinFactor, float *a) {
    if (n <= 0 || !a)
        return -1;
    float xmin = a[5];
    int loc = 0;
    for (int i = 5*splineBinFactor; i < (n-5)*splineBinFactor; i++)
        if (xmin > a[i] && a[i+1] < 0.5*a[i]) {
            xmin = a[i];
            loc = i;
        }
    return loc;
}

/*
 * Find maximum of the pulse 
 */
int FindMax(int n, int splineBinFactor, float *a) {
    if (n <= 0 || !a)
        return -1;
    float xmax = a[0];
    int loc = 0;
    for (int i = 0; i < n*splineBinFactor; i++)
        if (xmax < a[i]) {
            xmax = a[i];
            loc = i;
        }
    return loc;
}

/*
 * Find the constant offset of a signal
 */
float LinearFit_Baseline(TH1F *pulse, const int index_min, const int range) {
    TF1 *fBaseline = new TF1("fBaseline", "pol0", 10, index_min-range);
    pulse->Fit("fBaseline", "Q", "", 10, index_min-range);
    float base = fBaseline->GetParameter(0);
    delete fBaseline;
    return base;
}

/*
 * Assign pulse quality score based on different criteria
 */
unsigned int CheckPulseQuality(int binMin, int binMax, float *a, float offset) {
    unsigned int answer = 0;

    // Check if there is real pulse in the event
    if (a[binMax] - offset < MINPULSE)
        answer |= kNoPulse;

    // Check pulses with large opposite amplitude
    if (a[binMin] - offset < -0.1)
        answer |= kLargeNegativeAmplitude;

    // Check pulses with saturation
    if (a[binMax] == a[binMax+1] || a[binMax] == a[binMax-1])
        answer |= kSaturated;

    // Check that no points near the peak has negative polarity
    for (int i = binMax-3; i < binMax+4; i++)
        if (a[i] < 0)
            answer |= kNegativePolarity;

    // Check that no points near the peak has sudden jump
    bool hasSuddenJump = false;
    if (!(a[binMax] > a[binMax+1] && a[binMax+1] > a[binMax+2] && a[binMax+2] > a[binMax+3] && a[binMax+3] > a[binMax+4]))
        hasSuddenJump = true;
    if (!(a[binMax] > a[binMax-1] && a[binMax-1] > a[binMax-2] && a[binMax-2] > a[binMax-3] && a[binMax-3] > a[binMax-4]))
        hasSuddenJump = true;

    if (hasSuddenJump)
        answer |= kSuddenJump;
    
    // Check that second point away from peak is at least 10% lower
    // prevents strange pulses where there is flattish (but not completely flat) top
    bool hasFlatTop = false;
    if (!(fabs(a[binMax] - a[binMax+2])/fabs(a[binMax]) > 0.05))
        hasFlatTop = true;
    if (!(fabs(a[binMax] - a[binMax-2])/fabs(a[binMax]) > 0.05))
        hasFlatTop = true;

    if (hasFlatTop)
        answer |= kFlatTop;

    // Check for presence of 2nd pulse
    bool hasSecondPulse = false;
    int secondPulseIndex = -99;
    float secondPulseHeight = 0;
    for (int i = 20; i < 1000; i++) {
        if (secondPulseHeight > a[i]
            && fabs(a[i] - a[i-1]) / fabs(a[i]) < 0.5
            && fabs(a[i] - a[i+1]) / fabs(a[i]) < 0.5
            && abs(binMax - i) > 20
            && fabs(fabs(a[i]) - fabs(a[binMax])) / fabs(a[binMax]) < 0.15 ) {
                secondPulseHeight = a[i];
                secondPulseIndex = i;
                hasSecondPulse = true;
        }
    }

    if (hasSecondPulse)
        answer |= kSecondPulse;
    
    return answer; 
}

/*
 * Perform rectangular integration near peak of pulse
 */
float ChannelIntegral(float *a, int peak) {
    float integral = 0.;

    // For sharp gaussian peak type pulses
    for (int i = -INTRANGE; i <= INTRANGE; i++)
        integral += a[peak+i];

    //for scintillation type pulses
    /*for (int i = std::max(peak - 100, 2); i < std::min(peak + 800, 1023); i++)
        integral += a[i];
    */
    return integral;
}

/*
 * Fit Gaussian function to peak of signal and obtain time coordinate
 */
float GausFit_MeanTime(TH1F* pulse, const int index_first, const int index_last) {
    TF1 *fpeak = new TF1("fpeak", "gaus", index_first, index_last);
    pulse->Fit("fpeak", "Q", "", index_first, index_last);
    float timepeak = fpeak->GetParameter(1);
    delete fpeak;
    
    return timepeak;
}


/*
// Find rising edge of the pulse
int FindRisingEdge( int n, int binMax, float *a) {
    if (n <= 0 || !a)
        return -1;
    float xmin = a[0];
    int loc = -99;
    for (int i = binMax-10; i <binMax; i++) { // sometimes there is noise, check that it is rising in three bins
        if ( a[i] < -0.01 && a[i+1] < a[i] && a[i+2] < a[i+1] ) {
        // if ( Channel1Voltages_[i+2] < 0.3*Channel1Voltages_[index_min1])
            loc = i; 
            break;
        }
    }    
    return loc;
}

// For photonis: find the first pulse 
int FindFirstPulsePeak( int n, float *a) {

    if (n <= 0 || !a) return -1;

    int loc = 0;
    for    (int i = 20; i < 1000; i++) {
        if ( a[i] < -0.015
                 && a[i] <= a[i-1] && a[i] <= a[i+1]
                 && a[i-1] <= a[i-2]
                 && a[i-2] <= a[i-3]
            ) {

            loc = i;
            break;
        }
    }
    return loc;

}

#define FITPARAM "RQ"
#define LOW
#deinfe HIGH
void FitRisingEdge(TH1F *pulse, int nbinsL, int nbinsH, float &THM, float &risetime) {
    int bM = pulse->FindFirstBinAbove(HIGH * pulse->GetMaximum());
    int bL = pulse->FindFirstBinAbove(LOW * pulse->GetMaximum());
    TF1* f = new TF1("f", "[0]+x*[1]", pulse->GetBinCenter(bL+nbinsL), pulse->GetBinCenter(bM+nbinsH));
    pulse->Fit(f, FITPARAM);
    float m = f->GetParameter(1);
    float b = f->GetParameter(0);
    delete f;
    // 0.2 -> convert to picoseconds
    THM = 0.2 * 0.5 * (pulse->GetMaximum() - b) / m;
    risetime = 0.2 * (HIGH * pulse->GetMaximum() - b) / m - THM;
}

float LED(TH1F *pulse, double threshold, int nsamples, int splineBinFactor) {
    double crosstime;
    int bin = 0;
    
    // first sample above thresh (in absolute value)
    while(pulse->GetBinContent(bin) < threshold && bin< nsamples*splineBinFactor)
        bin++;
        // linear interpolation
        if (bin < nsamples * splineBinFactor) {
            double inf = pulse->GetBinContent(bin-1);
            double sup = pulse->GetBinContent(bin);     
            crosstime = pulse->GetXaxis()->GetBinCenter(bin-1) +\
                ((threshold-inf) / (sup-inf)) * (pulse->GetXaxis()->GetBinCenter(bin)-pulse->GetXaxis()->GetBinCenter(bin-1));
        }
        else {
            crosstime = 1000000.;
            std::cerr << "signal does not reach threshold\n";
        }
    return crosstime;
}
*/