#ifndef EXOBackgroundAnalysis_hh_
#define EXOBackgroundAnalysis_hh_
#include "EXOAnalysisManager/EXOAnalysisModule.hh"

#ifndef NOMYSQL
#include "EXOCalibUtilities/EXODriftVelocityCalib.hh"
#endif

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h>
#include <TSpectrum.h>
#include <TMath.h>

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <string.h>
#include <math.h>

#include "EXOUtilities/EXOChargeCluster.hh"
#include "EXOUtilities/EXOScintillationCluster.hh"
#include "EXOUtilities/EXOWaveformData.hh"
#include "EXOUtilities/EXOWaveform.hh"

class EXOBackgroundAnalysis : public EXOAnalysisModule 
{
  public:
    EXOBackgroundAnalysis() {}
    ~EXOBackgroundAnalysis() { ; }
    int Initialize();
    EventStatus BeginOfRun(EXOEventData *ED);
    EventStatus ProcessEvent(EXOEventData *ED);
    int TalkTo(EXOTalkToManager *tm);
    int ShutDown();
    void SetOutputFile(std::string aval) {oName = aval;}
    void SetDriftVelocity(double  aval) {driftVelocity = aval;}
    void SetELifetime(double aval) {elife = aval;}
    void SetWaveformFile(std::string aval) {WFName = aval;}
    void SetFullDataFile(std::string aval) {FullDataFile = aval;}

  private:
    std::string oName;
    TFile *oFile;
    TTree *oTree;

    int nED;
    int nEntriesOld;
    double runTime;

    int runID;

    std::string WFName;
    TFile *WFFile;
    TTree *WFTree;
    EXOEventData *WFED;
    bool WFFileLoaded;

    std::string FullDataFile;
    TFile *FDFile;
    TTree *FDTree;
    EXOEventData *FDED;
    bool FDFileLoaded;

    TH1F *hAPDSumWaveform;
    TH1F *hAPDSumWaveformBLSub;
    TH1F *hWireSumWaveform;
    TH1F *hWireSumWaveformBLSub;
    TH1F *hWireSumWaveformInv;
    TH1F *hWireSumWaveformBLSubInv;

    bool hasWaveforms;
    int n_sample;
    int nWaveFormPeaks;
    int WaveFormPeakX[100];

    int nNoiseEvents;

    double driftVelocity;
    double elife;
    double elifetime;

    int fnr;
    int fne;
    double fELifeTime;
    double fdR;
    double fdt;
    double fercl1;
    double fercl2;
    double feccl1;
    double feccl2;
    double fepcl1;
    double fepcl2;
    double fcsc1;
    double fcsc2;
    double fCL1;
    double fCL2;
    double fpCL1;
    double fpCL2;
    double fPeakRatio;
    double fx1;
    double fy1;
    double fz1;
    double fx2;
    double fy2;
    double fz2;
    int fAlg;
    EXOChargeCluster *fBetaCC;
    EXOChargeCluster *fAlphaCC;
    EXOScintillationCluster *fBetaSC;
    EXOScintillationCluster *fAlphaSC;

    bool first;
    int RunStart;
    int RunEnd;

    bool IsNoiseEvent(EXOEventData *ED);
    void SubtractBaseline(int *DataSamples, double *bl, double *stdev);
    bool FirstAttempt(EXOEventData *ED);
    bool SecondAttempt(EXOEventData *ED);
    bool ThirdAttempt(EXOEventData *ED);
    bool FourthAttempt(EXOEventData *ED);
    bool FirstSubAttempt(EXOEventData *ED);
    bool SecondSubAttempt(EXOEventData *ED);

  DEFINE_EXO_ANALYSIS_MODULE( EXOBackgroundAnalysis )
}; 

#include "EXOUtilities/EXOPluginUtilities.hh"
EXO_DEFINE_PLUGIN_MODULE(EXOBackgroundAnalysis, "bga")
#endif
