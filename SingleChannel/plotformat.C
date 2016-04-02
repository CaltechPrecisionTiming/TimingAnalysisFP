#include <iostream>
#include <cmath>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TArrow.h>
#include <TF1.h>
#include <TGraphErrors.h>

//Axis
#define axisTitleSize 0.06
#define axisTitleOffset .8

//Margins
#define leftMargin    0.15
#define rightMargin   0.05
#define topMargin     0.07
#define bottomMargin  0.12

#define bordersize    2
#define markerstyle   20

//TLatex
#define font          42
#define textsize      0.045

void TimeResolution(TString _rootFileName = "", TString outName = "time_resolution") {
    
    std::cout << "loading file" << std::endl;

    // Obtain file
    TFile *f = new TFile(_rootFileName);
    if(f->IsOpen()) {
        std::cout << "[INFO]: opening file: " << _rootFileName << std::endl;
    }
    else {
        std::cerr << "[ERROR]: could not open file: " << _rootFileName << std::endl;
        std::cerr << "[ERROR]: exiting!" << std::endl;
        return;
    }
    TCanvas *c = new TCanvas( "c", "c", 2119, 33, 800, 700 );

    std::cout << "preparing canvas" << std::endl;
    
    // set canvas parameters
    c->SetHighLightColor(2);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(bordersize);
    c->SetLeftMargin(leftMargin);
    c->SetRightMargin(rightMargin);
    c->SetTopMargin(topMargin);
    c->SetBottomMargin(bottomMargin);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);

    std::cout << "obtaining graph" << std::endl;

    // obtain graph
    TGraphErrors *graph = (TGraphErrors *) f->Get("Graph");
    
    std::cout << "formatting graph" << std::endl;

    // Set graph properties
    graph->SetTitle("");
    graph->GetYaxis()->SetTitle("Time Resolution [ns]");
    graph->GetXaxis()->SetTitle("Beam Energy [GeV]");
    graph->GetXaxis()->SetTitleSize(axisTitleSize);
    graph->GetXaxis()->SetTitleOffset(axisTitleOffset);
    graph->GetYaxis()->SetTitleSize(axisTitleSize);
    graph->GetYaxis()->SetTitleOffset(axisTitleOffset);
    graph->SetMarkerStyle(markerstyle);
    graph->Draw("AP");

    std::cout << "saving" << std::endl;

    // save plot
    c->SaveAs( outName + ".pdf" );
    // c->SaveAs( outName + ".C" );
    return;
}



void ChargeResolution(TString _rootFileName = "", TString outName = "charge_dependence") {
    
    // Obtain file
    TFile *f = new TFile(_rootFileName);
    if(f->IsOpen()) {
        std::cout << "[INFO]: opening file: " << _rootFileName << std::endl;
    }
    else {
        std::cerr << "[ERROR]: could not open file: " << _rootFileName << std::endl;
        std::cerr << "[ERROR]: exiting!" << std::endl;
        return;
    }
    TCanvas *c = new TCanvas( "c", "c", 2119, 33, 800, 700 );    


    // set canvas parameters
    c->SetHighLightColor(2);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(bordersize);
    c->SetLeftMargin(leftMargin);
    c->SetRightMargin(rightMargin);
    c->SetTopMargin(topMargin);
    c->SetBottomMargin(bottomMargin);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);


    // obtain graph
    TGraphErrors *graph = (TGraphErrors *) f->Get("Graph");
    
    // make a fitting function
    TF1 *line = new TF1("fit", "pol1");
    graph->Fit(line, "Q");
    line->Draw();

    // Set graph properties
    graph->SetTitle("");
    graph->GetYaxis()->SetTitle("Charge [fC]");
    graph->GetXaxis()->SetTitle("Beam Energy [GeV]");
    graph->GetXaxis()->SetTitleSize(axisTitleSize);
    graph->GetXaxis()->SetTitleOffset(axisTitleOffset);
    graph->GetYaxis()->SetTitleSize(axisTitleSize);
    graph->GetYaxis()->SetTitleOffset(axisTitleOffset);
    graph->SetMarkerStyle(markerstyle);
    graph->Draw("AP");

    // Draw fit constants
    TLatex *tex1 = new TLatex();
    tex1->SetTextColor(kBlack);
    tex1->SetTextFont(font);
    tex1->SetTextSize(textsize);

    tex1->DrawLatex(4, 17,
        Form("#splitline{p0:  %2.2f mm}{p1:  %2.2f mm}", line->GetParameter(0), line->GetParameter(1)));


    // save plot
    c->SaveAs( outName + ".pdf" );
    // c->SaveAs( outName + ".C" );
    return;
}
