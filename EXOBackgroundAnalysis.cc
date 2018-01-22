// ------------------------------------------------------------------------------------------------
// EXOBackgroundAnalysis
//
// This module searches for beta - alpha coincidences. The search is performed in a single trigger
// window which is appropriate for the Bi214 - Po214 decay coincidence. The ideal candidate event
// contains a weak scintillation signal with an associated strong charge cluster (beta decay)
// followed by a strong scintillation signal with an associated weak charge cluster (alpha decay).
// Deviations from this ideal case are also cinsidered and ordered in 4 search algorithms:
//
//  Alg1: two scintillation cluster and two charge cluster
//  Alg2: first scintillation signal is missing
//  Alg3: second charge cluster is missing (can be promoted to case1 event by reassigning the charge
//        cluster
//  Alg4: first scintillation signal and second charge cluster are missing
//
// The analysis relys on APD data and is therfore able to filter out noise events by summing up the
// APD and wire waveforms. We do not expect to see any peaks in the some of the wire waveforms. If
// there are peaks, the event is considered as noise event.

#include "EXOBackgroundAnalysis.hh"
#include "EXOUtilities/EXOErrorLogger.hh"
#include "EXOUtilities/EXOTalkToManager.hh"
#include "EXOUtilities/EXOEventData.hh"
#include "EXOUtilities/EXOChannelMap.hh"
#include "EXOAnalysisManager/EXOAnalysisManager.hh"
#include "EXOCalibUtilities/EXOCalibManager.hh"
#include "EXOCalibUtilities/EXODriftVelocityCalib.hh"
#include <iostream>

int EXOBackgroundAnalysis::Initialize()
{
  // Initialize histograms and output file
  
  std::cout << "Initializing module bga..." << std::endl;

  hAPDSumWaveform = new TH1F("hAPDSumWaveform", "APD Sum Waveform",2048,0,2048);
  hAPDSumWaveformBLSub = new TH1F("hAPDSumWaveformBLSub", "APD Sum Waveform - baseline > threshdold",2048,0,2048);
  hWireSumWaveform = new TH1F("hWireSumWaveform", "Wire Sum Waveform",2048,0,2048);
  hWireSumWaveformBLSub = new TH1F("hWireSumWaveformBLSub", "Wire Sum Waveform - baseline > threshdold",2048,0,2048);
  hWireSumWaveformInv = new TH1F("hWireSumWaveformInv", "Inverted Wire Sum Waveform",2048,0,2048);
  hWireSumWaveformBLSubInv = new TH1F("hWireSumWaveformBLSubInv", "Inverted Wire Sum Waveform - baseline > threshdold",2048,0,2048);

  oTree = new TTree("t","coincidences");

  oTree->Branch("fnr",&fnr,"fnr/I");
  oTree->Branch("fne",&fne,"fne/I");
  oTree->Branch("fELifeTime",&fELifeTime,"fELifeTime/D");
  oTree->Branch("fdR",&fdR,"fdR/D");
  oTree->Branch("fdt",&fdt,"fdt/D");
  oTree->Branch("fercl1",&fercl1,"fercl1/D");
  oTree->Branch("fercl2",&fercl2,"fercl2/D");
  oTree->Branch("feccl1",&feccl1,"feccl1/D");
  oTree->Branch("feccl2",&feccl2,"feccl2/D");
  oTree->Branch("fepcl1",&fepcl1,"fepcl1/D");
  oTree->Branch("fepcl2",&fepcl2,"fepcl2/D");
  oTree->Branch("fcsc1",&fcsc1,"fcsc1/D");
  oTree->Branch("fcsc2",&fcsc2,"fcsc2/D");
  oTree->Branch("fCL1",&fCL1,"fCL1/D");
  oTree->Branch("fCL2",&fCL2,"fCL2/D");
  oTree->Branch("fpCL1",&fpCL1,"fpCL1/D");
  oTree->Branch("fpCL2",&fpCL2,"fpCL2/D");
  oTree->Branch("fPeakRatio",&fPeakRatio,"fPeakRatio/D");
  oTree->Branch("fx1",&fx1,"fx1/D");
  oTree->Branch("fy1",&fy1,"fy1/D");
  oTree->Branch("fz1",&fz1,"fz1/D");
  oTree->Branch("fx2",&fx2,"fx2/D");
  oTree->Branch("fy2",&fy2,"fy2/D");
  oTree->Branch("fz2",&fz2,"fz2/D");
  oTree->Branch("fAlg",&fAlg,"fAlg/I");
  //oTree->Branch("fBetaChargeCluster",&fBetaCC);
  //oTree->Branch("fAlphaChargeCluster",&fAlphaCC);
  //oTree->Branch("fBetaScintCluster",&fBetaSC);
  //oTree->Branch("fAlphaScintCluster",&fAlphaSC);

  //driftVelocity = 0.18; // [mm/us]
  //elifetime = 290000; // [ns]

  nNoiseEvents = 0;

  first = true;

  // Register the relevant handlers for fetching drift velocity from the database
  //analysisManager->GetCalibManager()->registerHandler(dynamic_cast <EXOCalibHandlerBase*>(new EXODriftVelocityHandler(analysisManager->GetCalibManager())),EXOCalib::CTYPEdrift);

  WFED = 0;
  FDED = 0;
  nED = 0;
  nEntriesOld = 0;
  runTime = 0.0;

  return 0;
}

EXOAnalysisModule::EventStatus EXOBackgroundAnalysis::BeginOfRun(EXOEventData *ED)
{
  // retrieve drift velocity from database
  if (driftVelocity == -1) {
     #ifndef NOMYSQL
     std::cout << "EXOBackgroundAnalysis-> Trying to get drift velocity from database" << std::endl;

     EXODriftVelocityCalib *current_driftvelocity = GetCalibrationFor(EXODriftVelocityCalib, EXODriftVelocityHandler, "vanilla", ED->fEventHeader);

     if(current_driftvelocity) {
        double driftVelocity_TPC1 = current_driftvelocity->get_drift_velocity_TPC1();
        double driftVelocity_TPC2 = current_driftvelocity->get_drift_velocity_TPC2();

        // average the drift velocities in the two TPC halves
        driftVelocity = (driftVelocity_TPC1 + driftVelocity_TPC2) / 2.0;
        std::cout << "EXOBackgroundAnalysis-> Successfully retrieved drift velocity: " << driftVelocity << " mm/ns" << std::endl;
     }
     else {
        // set default drift velocity if database read failed
        driftVelocity = 0.0018;
        std::cout << "EXOBackgroundAnalysis-> Failed to read drift velocity from database. Setting default value: " << driftVelocity << " mm/ns" << std::endl;
     }

     #else
     driftVelocity = 0.0018;
     std::cout << "EXOBackgroundAnalysis-> Setting default drift velocity: " << driftVelocity << " mm/ns" << std::endl;
     #endif
  }
  else {std::cout << "EXOBackgroundAnalysis-> Set drift velocity is: " << driftVelocity << std::endl;}

  // Try to get root file that contains the waveforms
  if (WFName != "") {
     WFFile = new TFile(WFName.c_str(),"READ");
     if (WFFile->IsZombie()) {
        std::cout << "EXOBackgroundAnalysis-> Failed to open waveform file!" << std::endl;
        WFFileLoaded = false;
     }

     else {
        std::cout << "EXOBackgroundAnalysis-> Successfully opened waveform file with name " << WFName << std::endl;

        WFTree = (TTree*)WFFile->Get("tree");
        WFTree->SetBranchAddress("EventBranch",&WFED);
        WFTree->BuildIndex("EventBranch.fRunNumber","EventBranch.fEventNumber");

        WFFileLoaded = true;
     }
  }
  else {WFFileLoaded = false;}

  // Try to get root file that contains the waveforms
  if (FullDataFile != "") {
     FDFile = new TFile(FullDataFile.c_str(),"READ");
     if (FDFile->IsZombie()) {
        std::cout << "EXOBackgroundAnalysis-> Failed to open full data file! Assuming non masked data." << std::endl;
        FDFileLoaded = false;
     }

     else {
        std::cout << "EXOBackgroundAnalysis-> Successfully opened full data file with name " << FullDataFile << std::endl;

        FDTree = (TTree*)FDFile->Get("tree");
        FDTree->SetBranchAddress("EventBranch",&FDED);
        FDTree->BuildIndex("EventBranch.fRunNumber","EventBranch.fEventNumber");

        FDFileLoaded = true;
     }
  }
  else {std::cout << "EXOBackgroundAnalysis-> Failed to open full data file! Assuming non masked data." << std::endl; FDFileLoaded = false;}

  runID = ED->fRunNumber;
  //nED = 0;
  
  if (FDFileLoaded) {
    int CurrentRun = ED->fRunNumber;
    char cmd[50];
    sprintf(cmd,"fRunNumber == %i",CurrentRun);
    int CurrentEntries = FDTree->GetEntries(cmd);
    FDTree->GetEntry(nEntriesOld);
    double tStart = FDED->fEventHeader.fTriggerSeconds;
    FDTree->GetEntry(CurrentEntries + nEntriesOld - 1);
    double tEnd = FDED->fEventHeader.fTriggerSeconds;
    runTime += tEnd - tStart;
    nEntriesOld += CurrentEntries;
  }

  return kOk;
}

EXOAnalysisModule::EventStatus EXOBackgroundAnalysis::ProcessEvent(EXOEventData *ED)
{
  // First check if this event is a noise event. If no noise event the four search algorithm are
  // called one after the other until a canidiate event was found. If no event was found it just
  // continues with the next event

  if (first) {RunStart = ED->fEventHeader.fTriggerSeconds; first = false;}

  if (ED->fEventHeader.fTriggerSeconds != 0) {RunEnd = ED->fEventHeader.fTriggerSeconds;}

  runID = ED->fRunNumber;
  nED++;

  // get current lifetime
  if (elife == -1) {
     //double purFitP0 = 286;
     //double purFitP1 = -1.427;
     double purFitP0;
     double purFitP1;
     double purFitP2;
     double purFitP3;
     double purFitP4;

     double purTime = double(ED->fEventHeader.fTriggerSeconds - 1304146800.0) / 3600.0 / 24.0;

     if (purTime < 58) {
           purFitP0 = -284.596;
           purFitP1 = 53.6978;
           purFitP2 = -1.88664;
           purFitP3 = 0.0269101;
           purFitP4 = -0.000133772;
        }
        if (purTime >= 58 && purTime < 81.6) {
           purFitP0 = 14068.5;
           purFitP1 = -908.011;
           purFitP2 = 21.8864;
           purFitP3 = -0.230994;
           purFitP4 = 0.00090631;
        }
        if (purTime >= 81.6 && purTime < 94.0) {
           purFitP0 = -9011.55;
           purFitP1 = 115.417;
           purFitP2 = 0.0;
           purFitP3 = 0.0;
           purFitP4 = 0.0;
        }
        if (purTime >= 94.0 && purTime < 102.5) {
           purFitP0 = 2000.0;
           purFitP1 = 0.0;
           purFitP2 = 0.0;
           purFitP3 = 0.0;
           purFitP4 = 0.0;
        }
        if (purTime >= 102.5 && purTime < 113.0) {
           purFitP0 = -1208000.0;
           purFitP1 = 34380.0;
           purFitP2 = -325.9;
           purFitP3 = 1.03;
           purFitP4 = 0.0;
        }
        if (purTime >= 113.0 && purTime < 129.6) {
           purFitP0 = -48740.0;
           purFitP1 = 805.0;
           purFitP2 = -3.259;
           purFitP3 = 0.0;
           purFitP4 = 0.0;
        }
        if (purTime >= 129.6 && purTime < 142.0) {
           purFitP0 = -29510.0;
           purFitP1 = 230.1;
           purFitP2 = 0.0;
           purFitP3 = 0.0;
           purFitP4 = 0.0;
        }
        if (purTime >= 142.0) {
           purFitP0 = 5000.0;
           purFitP1 = 0.0;
           purFitP2 = 0.0;
           purFitP3 = 0.0;
           purFitP4 = 0.0;
        }

     elifetime = (purFitP4*purTime*purTime*purTime*purTime + purFitP3*purTime*purTime*purTime + purFitP2*purTime*purTime + purFitP1*purTime + purFitP0) * 1000.0;
  }
  else {elifetime = elife;}

  fELifeTime = elifetime;

  if (IsNoiseEvent(ED)) {nNoiseEvents++; return kOk;}

  bool result;
  if (result = FirstAttempt(ED)) {oTree->Fill(); return kOk;}
  if (result = SecondAttempt(ED)) {oTree->Fill(); return kOk;}
  if (result = ThirdAttempt(ED)) {oTree->Fill(); return kOk;}
  if (result = FourthAttempt(ED)) {oTree->Fill(); return kOk;}
  if (result = FirstSubAttempt(ED)) {oTree->Fill(); return kOk;}
  if (result = SecondSubAttempt(ED)) {oTree->Fill(); return kOk;}

  return kOk;
}

int EXOBackgroundAnalysis::TalkTo(EXOTalkToManager *tm)
{
  // Set input commands for this module. The output file, drift velocity and electron lifetime can
  // be set
  
  tm->CreateCommand("/bga/file","name of output file",this,"output.root",&EXOBackgroundAnalysis::SetOutputFile);
  if ( oName == NULL ) {
    LogEXOMsg("Unable to create output file command",EECritical);
  }
  
  tm->CreateCommand("/bga/DriftVelocity","drift velocity in mm/ns",this,0.0018,&EXOBackgroundAnalysis::SetDriftVelocity);
  if ( driftVelocity == 0.0 ) {
    LogEXOMsg("Unable to set drift velocity",EECritical);
  }
  
  tm->CreateCommand("/bga/ElectronLifetime","electron lifetime in ns",this,290000.0,&EXOBackgroundAnalysis::SetELifetime);
  if ( elife == 0.0 ) {
    LogEXOMsg("Unable to set electron lifetime",EECritical);
  }

  tm->CreateCommand("/bga/WaveformData","name of waveform data file",this,"",&EXOBackgroundAnalysis::SetWaveformFile);
  if ( WFName.empty() ) {
    LogEXOMsg("No waveform data file specified",EECritical);
  }

  tm->CreateCommand("/bga/FullDataFile"," name of full data file",this,"",&EXOBackgroundAnalysis::SetFullDataFile);
  if ( FullDataFile.empty() ) {
    LogEXOMsg("No full data file specified",EECritical);
  }

  return 0;
}

int EXOBackgroundAnalysis::ShutDown()
{
  // At shutdown the output tree is written to the file

  std::cout << "EXOBackgroundAnalysis-> Opening the output file to write background analysis results..." << std::endl;
  oFile = new TFile(oName.c_str(),"RECREATE");

  oTree->Write();

  std::cout << "EXOBackgroundAnalysis-> Closing file..." << std::endl;
  oFile->Close();

  std::cout << "EXOBackgroundAnalysis-> Calculating 214Bi rate and 222Rn content in Xenon..." << std::endl;

  double eff = 0.5652; //0.8884; // overall efficiency
  double halflife = 330350.4; // half life of Rn222
  double ln2 = TMath::Log(2); // ln(2)

  double dt = 1.0;
  
  if (FDFileLoaded) {
    int nFD = FDTree->GetEntries();
    FDTree->GetEntry(0);
    double tStart = FDED->fEventHeader.fTriggerSeconds;
    FDTree->GetEntry(nFD - 1);
    double tEnd = FDED->fEventHeader.fTriggerSeconds;
    
    double dtTMP = runTime / 3600.0;
    double subSec = (dtTMP - TMath::Floor(dtTMP)) * 60.0;
    if (subSec < 10.0) {dt = TMath::Floor(dtTMP) * 1200.0 + subSec * 60.0;}
    if (subSec >= 10.0 && subSec < 30.0) {dt = TMath::Floor(dtTMP) * 1200.0 + 600.0;}
    if (subSec >= 30.0 && subSec < 40.0) {dt = TMath::Floor(dtTMP) * 1200.0 + 600.0 + (subSec - 30.0) * 60.0;}
    if (subSec >= 40.0) {dt = TMath::Floor(dtTMP) * 1200.0 + 1200.0;}
    
    std::cout << "dt = " << tEnd - tStart << "    dt (h) = " << dtTMP << "    dt (s) = " << dt << std::endl;
  }
  else {dt = double(RunEnd - RunStart);}

  //double scale = double(nFD)/double(nED);
  //eff /= scale;

  //std::cout << "nFD = " << nFD << "    nED = " << nED << "    scale = " << scale << "    dt = " << dt << std::endl;

  double N_Bi = double(oTree->GetEntries("fAlg != 4")); // total number of Bi-Po events
  double nrTotal = double(oTree->GetEntries("fz2 != -999 && fAlg != 4")); // total number of Bi-Po with reconstructed position
  double nrCathode = double(oTree->GetEntries("fz2 > -5 && fz2 < 5 && fz2 != -999 && fAlg != 4")); // number of events on the cathode
  double nrAnode = double(oTree->GetEntries("TMath::Abs(fz2) > 190 && fz2 != -999")); // number of events at the anode
  double frac = 0.0;
  double frac2 = 0.0;
  double frac_errTMP = 0.0;
  double frac_err = 0.0;
  double frac2_err = 0.0;
  if (nrTotal > 0.0) {
     frac = nrCathode / nrTotal;
     frac_errTMP = TMath::Sqrt(1.0/(nrTotal*nrTotal)*nrCathode + (nrCathode/(nrTotal*nrTotal))*(nrCathode/(nrTotal*nrTotal))*nrTotal);
  }
  frac_err = TMath::Sqrt(frac_errTMP*frac_errTMP + 0.04*0.04); // add 4% systematic error due to the z-cut

  if (nrTotal > 0.0) {
     frac2 = nrAnode / nrTotal;
     frac2_err = TMath::Sqrt(1.0/(nrTotal*nrTotal)*nrAnode + (nrAnode/(nrTotal*nrTotal))*(nrAnode/(nrTotal*nrTotal))*nrTotal);
  }

  double e1 = 2.1317; // 1/efficiency for detecting the Bi-Po event on the cathode
  double e2 = 2.0; // 1/efficiency for detecting the Bi-Po event at the anode

  double N_bulk = N_Bi - N_Bi*frac - N_Bi*frac2; // number of bulk events
  double N_cathode = N_Bi*frac; // number of cathode events
  double N_anode = N_Bi*frac2; // number of anode events
  double N = (N_bulk + e1*N_cathode + e2*N_anode) / eff; // total number of events
  double N_err2 = ((1-frac*(1-e1)-frac2*(1-e2))*(1-frac*(1-e1)-frac2*(1-e2))*N_Bi + (e1 - 1)*(e1 - 1)*N_Bi*N_Bi*frac_err*frac_err + (e2 - 1)*(e2 - 1)*N_Bi*N_Bi*frac2_err*frac2_err) / (eff*eff);

  double A = N / dt * 1000; // activity in mBq
  double A_err = TMath::Sqrt(N_err2) / dt * 1000;

  double m = 122.595; // mass of the active xenon
  double m_err = 12.2595;
  double A_norm = A / m * 1000.0; // activity in uBq/kg
  double A_norm_err = TMath::Sqrt(1/(m*m)*A_err*A_err + (A/(m*m))*(A/(m*m))*m_err*m_err) * 1000;

  double N_Rn = A * halflife / ln2 / 1000;
  double N_Rn_err = A_err * halflife / ln2 / 1000;

  std::cout << "EXOBackgroundAnalysis-> Printing summary..." << std::endl;
  std::cout << "EXOBackgroundAnalysis-> Number noise events: " << nNoiseEvents << std::endl;
  std::cout << "Run start: " << RunStart - 1304146800 << std::endl;
  std::cout << "**************************************************************************************************" << std::endl;
  std::cout << "Nr. Events\tN_Bi\tfrac\t\tExposure\tA_Bi\t\t\tN_Rn" << std::endl;
  std::cout << "--------------------------------------------------------------------------------------------------" << std::endl;
  std::cout << oTree->GetEntries() << "\t\t" << N_Bi << "\t" << frac << "\t" << dt << "\t\t" << A << " +- " << A_err << "\t" << N_Rn << " +- " << N_Rn_err << std::endl;
  std::cout << "**************************************************************************************************" << std::endl;
  std::cout << "**** Formated string **************************************************************" << std::endl;
  //std::cout << runID << "\t" << printf("%.4f",RunStart+dt/2.0) << "\t" << printf("%.4f",dt/2.0) << "\t" << printf("%.5f",A) << "\t" << printf("%.5f",A_err) << "\t" << printf("%.3f",N_Rn) << "\t" << printf("%.3f",N_Rn_err) << std::endl;
  char FormatedOutput[10000];
  sprintf(FormatedOutput,"%i\t%.4f\t%.4f\t%.5f\t%.5f\t%.3f\t%.3f\t%.4f\t%.4f\n",runID,(RunStart - 1304146800 + dt/2.0) / 3600.0,dt/2.0/3600.0,A,A_err,N_Rn,N_Rn_err,frac,frac_err);

  //std::cout << runID << "\t" << (RunStart - 1304146800 + dt/2.0) / 3600.0 << "\t" << dt/2.0/3600.0 << "\t" << A << "\t" << A_err << "\t" << N_Rn << "\t" << N_Rn_err << "\t" << frac << std::endl;
  std::cout << FormatedOutput << std::endl;
  std::cout << "Open file to write string..." << std::endl;
  FILE *ResultFile;
  ResultFile = fopen("RadonData.dat","a");
  fputs(FormatedOutput,ResultFile);
  fclose(ResultFile);

  return 0;
}

bool EXOBackgroundAnalysis::IsNoiseEvent(EXOEventData *ED)
{
  // This function checks if the current event is a noise event. It sums the APD and wire waveforms
  // and searches for peaks. We do not expect large peaks in the summed wire waveform. If there are
  // peaks the event is considered noise event. It also stores the peak positions of the APD signal.
  // This information is later used to match the positions from waveforms and the positions found by
  // reconstruction in case reconstruction finds more peaks than expected.
  
  // check if there are waveforms in the event data
  int nWaveforms = ED->GetWaveformData()->GetNumWaveforms();
  if (nWaveforms > 0) {hasWaveforms = true;}
  else {hasWaveforms = false;}

  if (!hasWaveforms && !WFFileLoaded) {
     LogEXOMsg("No waveforms in event data",EECritical);
     hasWaveforms = false;
     nWaveFormPeaks = 0;
     
     // if no waveforms were found we consider the event noisy if the number of charge clusters is
     // larger than 20
     int ncl = ED->GetNumChargeClusters();
     if (ncl > 8) {return true;}
     
     return false;
  }

  int EvtID = ED->fEventNumber;
  EXOWaveformData* wf_data = 0;
  if (WFFileLoaded) {
     int id = WFTree->GetEntryNumberWithIndex(runID,EvtID);
     WFTree->GetEntry(id);

     wf_data = WFED->GetWaveformData();
     wf_data->Decompress();
  }

  if (hasWaveforms) {n_sample = ED->GetWaveformData()->fNumSamples;}
  else {n_sample = WFED->GetWaveformData()->fNumSamples;}

  hAPDSumWaveform->Reset();
  hAPDSumWaveformBLSub->Reset();
  hWireSumWaveform->Reset();
  hWireSumWaveformBLSub->Reset();
  hWireSumWaveformInv->Reset();
  hWireSumWaveformBLSubInv->Reset();

  double ChargeSigmaThreshold = 4.0;
  double ChargeThreshold = 500;
  double ChargeThresholdInv = 500;
  double LightSigmaThreshold = 4.0;
  double LightThreshold = 300;

  // loop over all channels
  int *sumAPD = new int[n_sample];
  int *sumWire = new int[n_sample];
  for (int k = 0; k < n_sample; k++) {sumAPD[k] = 0; sumWire[k] = 0;}

  int nrSig = 0;
  if (hasWaveforms) {nrSig = ED->GetWaveformData()->GetNumWaveforms();}
  else {nrSig = WFED->GetWaveformData()->GetNumWaveforms();}

  for (int i = 0; i < nrSig; i++) {
     EXOWaveform *wf = 0;
     if (hasWaveforms) {wf = (EXOWaveform*)ED->GetWaveformData()->GetWaveform(i);}
     else {wf = (EXOWaveform*)WFED->GetWaveformData()->GetWaveform(i);}

     int chID = wf->fChannel;

     if (chID >= 152) {
        for (int k = 0; k < n_sample; k++) {sumAPD[k] += wf->At(k);}
     }
     else {
        for (int k = 0; k < n_sample; k++) {sumWire[k] += wf->At(k);}
     }
  }

  for (int k = 0; k < n_sample; k++) {
     hAPDSumWaveform->SetBinContent(k+1,sumAPD[k]);
     hWireSumWaveform->SetBinContent(k+1,sumWire[k]);
     hWireSumWaveformInv->SetBinContent(k+1,-1 * sumWire[k]);
  }

  double bl_APD;
  double bl_Wire;
  double stdev_APD;
  double stdev_Wire;

  SubtractBaseline(sumAPD,&bl_APD,&stdev_APD);
  SubtractBaseline(sumWire,&bl_Wire,&stdev_Wire);

  double noiseThreshold_APD = LightSigmaThreshold * stdev_APD;
  double noiseThreshold_Wire = ChargeSigmaThreshold * stdev_Wire;

  int *sumAPDBLSub = new int[n_sample];
  int *sumWireBLSub = new int[n_sample];
  for (int k = 0; k < n_sample; k++) {
     if (sumAPD[k] > LightThreshold) {hAPDSumWaveformBLSub->SetBinContent(k+1,sumAPD[k]); sumAPDBLSub[k] = sumAPD[k];}
     else {hAPDSumWaveformBLSub->SetBinContent(k+1,0); sumAPDBLSub[k] = 0;}

     if (sumWire[k] > ChargeThreshold) {hWireSumWaveformBLSub->SetBinContent(k+1,sumWire[k]); sumWireBLSub[k] = sumWire[k];}
     else {hWireSumWaveformBLSub->SetBinContent(k+1,0); sumWireBLSub[k] = 0;}

     if (-1 * sumWire[k] > ChargeThresholdInv) {hWireSumWaveformBLSubInv->SetBinContent(k+1,-1 * sumWire[k]);}
     else {hWireSumWaveformBLSubInv->SetBinContent(k+1,0);}
  }

  TSpectrum *specAPD = new TSpectrum();
  specAPD->Search(hAPDSumWaveformBLSub, 2, "goff");

  int nPeaksAPD = specAPD->GetNPeaks();
  float *peakPositionAPDX = specAPD->GetPositionX();
  float *peakPositionAPDY = specAPD->GetPositionY();

  TSpectrum *specWire = new TSpectrum();
  specWire->Search(hWireSumWaveformBLSub, 2, "goff");

  int nPeaksWire = specWire->GetNPeaks();
  float *peakPositionWireX = specWire->GetPositionX();
  float *peakPositionWireY = specWire->GetPositionY();

  TSpectrum *specWireInv = new TSpectrum();
  specWireInv->Search(hWireSumWaveformBLSubInv, 2, "goff");

  int nPeaksWireInv = specWireInv->GetNPeaks();
  float *peakPositionWireInvX = specWireInv->GetPositionX();
  float *peakPositionWireInvY = specWireInv->GetPositionY();

  // get number of waveform peaks and their positions
  if (nPeaksAPD < 100) {
     nWaveFormPeaks = nPeaksAPD;
     for (int i = 0; i < nPeaksAPD; i++) {
        int PeakX = int(peakPositionAPDX[i]) - 1;
        int PeakY;
        while (1) {
           if (PeakX >= 0) {PeakY = sumAPD[PeakX];}
           else {break;}
           if (PeakY < stdev_APD) {break;}
           PeakX--;
        }

        WaveFormPeakX[i] = PeakX;
     }
  }
  else {nWaveFormPeaks = 100;}

  // find matches of peak positions in APD and wire spectrum
  int nMatchesP = 0;
  for (int i = 0; i < nPeaksWire; i++) {
     for (int k = 0; k < nPeaksAPD; k++) {if (TMath::Abs(peakPositionAPDX[k]-peakPositionWireX[i]) <= 5) {nMatchesP++; break;}}
  }

  int nMatchesN = 0;
  for (int i = 0; i < nPeaksWireInv; i++) {
     for (int k = 0; k < nPeaksAPD; k++) {if (TMath::Abs(peakPositionAPDX[k]-peakPositionWireInvX[i]) <= 5) {nMatchesN++; break;}}
  }

  // clean up memory
  delete sumAPD;
  delete sumWire;
  delete specAPD;
  delete specWire;
  delete specWireInv;

  // apply cuts
  //if (nPeaksWire > 0 || nPeaksWireInv > 0) {return true;}
  //if (nMatchesP > 0 || nMatchesN > 0) {return true;}
  if (nMatchesN > 0) {return true;}

  if (WFFileLoaded) {hasWaveforms = true;}

  return false;
}

void EXOBackgroundAnalysis::SubtractBaseline(int *DataSamples, double *bl, double *stdev)
{
  double blTMP = 0.0;
  int nSumSample = n_sample / 6;
  double *DataMean = new double[nSumSample];

  for (int i = 0; i < nSumSample; i++) {
     if (nSumSample > 0) {blTMP += double(DataSamples[i]) / nSumSample;}
     DataMean[i] = double(DataSamples[i]);
  }

  // standard deviation of first nsample/6 samples
  double stdevTMP = 0.0;
  for (int i = 0; i < nSumSample; i++) {
     if (nSumSample > 0) {stdevTMP += (DataMean[i] - blTMP)*(DataMean[i] - blTMP) / nSumSample;}
  }

  stdevTMP = TMath::Sqrt(stdevTMP);

  for (int i = 0; i < n_sample; i++) {DataSamples[i] = DataSamples[i] - int(blTMP);}

  *bl = blTMP;
  *stdev = stdevTMP;

  // clean up memory
  delete DataMean;

  return;
}

// Require two charge cluster at the same position (within R) separated in time and two scintillation cluster.
// Calculate charge/light and peak ratio.
bool EXOBackgroundAnalysis::FirstAttempt(EXOEventData *ED)
{
  double R = 20; // radius within the position must match [mm]
  double tclDiffTresh = 5; // time difference that the two charge cluster must have [us]

  int ncl = int(ED->GetNumChargeClusters());
  int nsc = int(ED->GetNumScintillationClusters());

  if (ncl < 2) {return false;} // require at least two charge cluster
  if (nsc < 2) {return false;} // require at least two scintillation cluster

  EXOScintillationCluster *FirstSC = 0; // this is the scintillation cluster of the beta decay
  EXOScintillationCluster *SecondSC = 0; // this is the scintillation cluster of the alhpa decay

  if (nWaveFormPeaks < 2 && hasWaveforms) {return false;}
  if (nsc == 2) {
     double TimeSC1 = ED->GetScintillationCluster(0)->fTime;
     double TimeSC2 = ED->GetScintillationCluster(1)->fTime;

     if (TimeSC2 > TimeSC1) {
        FirstSC = ED->GetScintillationCluster(0);
        SecondSC = ED->GetScintillationCluster(1);
     }
     else {
        FirstSC = ED->GetScintillationCluster(1);
        SecondSC = ED->GetScintillationCluster(0);
     }
  }
  else {
     if (!hasWaveforms) {return false;}
     
     // if there are more than two scintillation cluster, compare with sum APD peaks
     if (nsc > 2 && nWaveFormPeaks > 2) {return false;} // there are definitely more than two scintillation cluster
     if (nsc > 2 && nWaveFormPeaks == 2) {

        // find first waveform peak in reconstructed cluster
        int m = 0;
        double TimeSC1 = 0.0;
        bool foundFirst = false;
        while (m < nsc) {
           TimeSC1 = ED->GetScintillationCluster(m)->fTime;
           if (WaveFormPeakX[0] > TMath::Nint(TimeSC1/1000) - 2 && WaveFormPeakX[0] < TMath::Nint(TimeSC1/1000) + 2) {foundFirst = true; break;}
           m++;
        }

        // find second waveform peak in reconstructed cluster
        int n = 0;
        double TimeSC2 = 0.0;
        bool foundSecond = false;
        while (n < nsc) {
           TimeSC2 = ED->GetScintillationCluster(n)->fTime;
           if (WaveFormPeakX[1] > TMath::Nint(TimeSC2/1000) - 2 && WaveFormPeakX[1] < TMath::Nint(TimeSC2/1000) + 2) {foundSecond = true; break;}
           n++;
        }

        if (foundFirst && foundSecond) {
           if (TimeSC2 > TimeSC1) {
              FirstSC = ED->GetScintillationCluster(m);
              SecondSC = ED->GetScintillationCluster(n);
           }
           else {
              FirstSC = ED->GetScintillationCluster(n);
              SecondSC = ED->GetScintillationCluster(m);
           }
        }
        else {return false;}
     }
  }

  int k = 0;
  int i = 0;
  double dist = 0.0;
  double tclDiff = 0.0;
  bool found = false;

  double CL1_X = 0.0;
  double CL1_Y = 0.0;
  double CL1_Z = 0.0;
  double CL1_Time = 0.0;

  double CL2_X = 0.0;
  double CL2_Y = 0.0;
  double CL2_Z = 0.0;
  double CL2_Time = 0.0;

  while (k < ncl - 1) {
     CL1_X = ED->GetChargeCluster(k)->fX;
     CL1_Y = ED->GetChargeCluster(k)->fY;
     CL1_Z = ED->GetChargeCluster(k)->fZ;
     CL1_Time = ED->GetChargeCluster(k)->fCollectionTime;

     for (i = k+1; i < ncl; i++) {
        CL2_X = ED->GetChargeCluster(i)->fX;
        CL2_Y = ED->GetChargeCluster(i)->fY;
        CL2_Z = ED->GetChargeCluster(i)->fZ;
        CL2_Time = ED->GetChargeCluster(i)->fCollectionTime;

        dist = TMath::Sqrt((CL1_X - CL2_X)*(CL1_X - CL2_X) + (CL1_Y - CL2_Y)*(CL1_Y - CL2_Y) + (CL1_Z - CL2_Z)*(CL1_Z - CL2_Z));
        tclDiff = TMath::Abs(CL1_Time - CL2_Time) / 1000;
        if (dist < R && tclDiff > tclDiffTresh) {found = true; break;}
     }

     if (found) {break;}
     k++;
  }

  if (!found) {return false;}

  EXOChargeCluster *FirstCC; // this is the charge cluster of the beta decay
  EXOChargeCluster *SecondCC; // this is the charge cluster of the alpha decay

  if (CL2_Time > CL1_Time) {
     FirstCC = ED->GetChargeCluster(k);
     SecondCC = ED->GetChargeCluster(i);
  }
  else {
     FirstCC = ED->GetChargeCluster(i);
     SecondCC = ED->GetChargeCluster(k);
  }

  fnr = ED->fRunNumber;
  fne = ED->fEventNumber;
  fdR = dist;
  fdt = SecondSC->fTime/1000 - FirstSC->fTime/1000;
  if (TMath::Abs(fdt - tclDiff) > 10) {return false;}
  if (fdt < 110 || fdt > 980) {return false;} // apply time window cut (110us > dt < 980)

  fercl1 = FirstCC->fRawEnergy;
  fercl2 = SecondCC->fRawEnergy;
  feccl1 = FirstCC->fCorrectedEnergy;
  feccl2 = SecondCC->fCorrectedEnergy;
  fepcl1 = FirstCC->fCorrectedEnergy * TMath::Exp(FirstCC->fDriftTime / elifetime);
  fepcl2 = SecondCC->fCorrectedEnergy * TMath::Exp(SecondCC->fDriftTime / elifetime);
  fcsc1 = (FirstSC->fCountsSumOnAPDPlaneOne + FirstSC->fCountsSumOnAPDPlaneTwo);
  fcsc2 = (SecondSC->fCountsSumOnAPDPlaneOne + SecondSC->fCountsSumOnAPDPlaneTwo);

  // update purity corrected energy if values are set to zero
  if (FirstCC->fPurityCorrectedEnergy == 0) {FirstCC->fPurityCorrectedEnergy = fepcl1;}
  if (SecondCC->fPurityCorrectedEnergy == 0) {SecondCC->fPurityCorrectedEnergy = fepcl2;}

  double sumEnergy = 0.0;
  double sumPEnergy = 0.0;
  for (int i = 0; i < int(FirstSC->GetNumChargeClusters()); i++) {
     sumEnergy += FirstSC->GetChargeClusterAt(i)->fCorrectedEnergy;
     sumPEnergy += FirstSC->GetChargeClusterAt(i)->fCorrectedEnergy * TMath::Exp(FirstSC->GetChargeClusterAt(i)->fDriftTime / elifetime);
  }

  fCL1 = sumEnergy / fcsc1;
  fCL2 = fercl2 / fcsc2;
  fpCL1 = sumPEnergy / fcsc1;
  fpCL2 = fepcl2 / fcsc2;
  fPeakRatio = fcsc2 / fcsc1;
  fx1 = FirstCC->fX;
  fy1 = FirstCC->fY;
  fz1 = FirstCC->fZ;
  fx2 = SecondCC->fX;
  fy2 = SecondCC->fY;
  fz2 = SecondCC->fZ;
  fAlg = 1;

  // return false if the charge clusters were not correctly assigned to the scintillation signals
  if (fz1 == -100000 || fz2 == -100000) {return false;}

  fBetaCC = FirstCC;
  fAlphaCC = SecondCC;
  fBetaSC = FirstSC;
  fAlphaSC = SecondSC;

  return true;
}

// Require two charge cluster at same x-y position (within R) separated in time. z position can not be compared
// because the first scintillaion cluster is missing. There must be only one scintillation cluster not at t0 and
// one charge cluster at t0.
bool EXOBackgroundAnalysis::SecondAttempt(EXOEventData *ED)
{
  double R = 15; // radis within the position must match [mm]
  double tclDiffTresh = 5; // time difference that the two charge cluster must have [us]
  double cscThresh = 10000; // the scintillation cluster must at least have cscThresh photon counts

  int ncl = int(ED->GetNumChargeClusters());
  int nsc = int(ED->GetNumScintillationClusters());

  if (ncl < 2) {return false;} // require at least two charge cluster
  if (nsc == 0) {return false;} // require one scintillation cluster

  EXOScintillationCluster *FirstSC = 0; // this is the scintillation cluster of the beta decay
  EXOScintillationCluster *SecondSC = 0; // this is the scintillation cluster of the alpha decay

  if (nWaveFormPeaks == 0 && hasWaveforms) {return false;}
  if (nsc == 1) {SecondSC = ED->GetScintillationCluster(0);}
  else {
     if (!hasWaveforms) {return false;}
     
     // if there are more than one scintillation cluster, compare with sum APD peaks
     if (nsc > 1 && nWaveFormPeaks > 1) {return false;} // there are definitely more than one scintillation cluster
     if (nsc > 1 && nWaveFormPeaks == 1) { // find the one waveform peak in reconstructed cluster
        int k = 0;
        bool foundFirst = false;
        while (k < nsc) {
           double TimeSC = ED->GetScintillationCluster(k)->fTime;
           if (WaveFormPeakX[0] > TMath::Nint(TimeSC/1000) - 2 && WaveFormPeakX[0] < TMath::Nint(TimeSC/1000) + 2) {foundFirst = true; break;}
           k++;
        }

        if (!foundFirst) {return false;}
        SecondSC = ED->GetScintillationCluster(k);
     }
  }

  double TimeSC = SecondSC->fTime;
  //if (TimeSC/1000 >  ED->fEventHeader.fTriggerOffset - 10 && TimeSC/1000 < ED->fEventHeader.fTriggerOffset + 10) {return false;} // the scintillation cluster must not be at t0

  bool foundCharge = false;
  for (int i = 0; i < ncl; i++) {
     double CL_Time = ED->GetChargeCluster(i)->fCollectionTime;
     if (CL_Time < TimeSC) {foundCharge = true; break;}
     //if (CL_Time/1000 > ED->fEventHeader.fTriggerOffset - 10 && CL_Time/1000 < ED->fEventHeader.fTriggerOffset + 10) {foundCharge = true; brek;}
  }

  //if (!foundCharge) {return false;} // one charge cluster must be at t0

  int k = 0;
  int i = 0;
  double dist = 0.0;
  double tclDiff = 0.0;
  bool found = false;

  double CL1_X = 0.0;
  double CL1_Y = 0.0;
  double CL1_Time = 0.0;

  double CL2_X = 0.0;
  double CL2_Y = 0.0;
  double CL2_Time = 0.0;

  while (k < ncl - 1) {
     CL1_X = ED->GetChargeCluster(k)->fX;
     CL1_Y = ED->GetChargeCluster(k)->fY;
     CL1_Time = ED->GetChargeCluster(k)->fCollectionTime;

     for (i = k+1; i < ncl; i++) {
        CL2_X = ED->GetChargeCluster(i)->fX;
        CL2_Y = ED->GetChargeCluster(i)->fY;
        CL2_Time = ED->GetChargeCluster(i)->fCollectionTime;

        dist = TMath::Sqrt((CL1_X - CL2_X)*(CL1_X - CL2_X) + (CL1_Y - CL2_Y)*(CL1_Y - CL2_Y));
        tclDiff = TMath::Abs(CL1_Time - CL2_Time) / 1000;

        if (CL1_Time > CL2_Time) {
           if (dist < R && tclDiff > tclDiffTresh) {found = true; break;}
        }
        else {
           if (dist < R && tclDiff > tclDiffTresh) {found = true; break;}
        }
     }

     if (found) {break;}
     k++;
  }

  if (!found) {return false;}

  EXOChargeCluster *FirstCC = 0; // this is the charge cluster of the beta decay
  EXOChargeCluster *SecondCC = 0; // this is the charge cluster of the alpha decay

  if (CL2_Time > CL1_Time) {
     FirstCC = ED->GetChargeCluster(k);
     SecondCC = ED->GetChargeCluster(i);
  }
  else {
     FirstCC = ED->GetChargeCluster(i);
     SecondCC = ED->GetChargeCluster(k);
  }

  if (FirstCC->fCollectionTime > TimeSC) {return false;} // with the time window cut the first charge cluster must be before the scintillation cluster
  if (SecondCC->fCollectionTime < TimeSC) {return false;} // with the time window cut the second charge cluster must be after the scintillation cluster

  fnr = ED->fRunNumber;
  fne = ED->fEventNumber;
  fdR = dist;
  fdt = SecondCC->fCollectionTime/1000 - FirstCC->fCollectionTime/1000;
  if (fdt < 110 || fdt > 980) {return false;} // apply time window cut (110us > dt < 980us)

  fercl1 = FirstCC->fRawEnergy;
  fercl2 = SecondCC->fRawEnergy;
  feccl1 = FirstCC->fCorrectedEnergy;
  feccl2 = SecondCC->fCorrectedEnergy;
  fepcl1 = FirstCC->fCorrectedEnergy * TMath::Exp(FirstCC->fDriftTime / elifetime);
  fepcl2 = SecondCC->fCorrectedEnergy * TMath::Exp(SecondCC->fDriftTime / elifetime);

  fcsc1 = 0;
  fcsc2 = (SecondSC->fCountsSumOnAPDPlaneOne + SecondSC->fCountsSumOnAPDPlaneTwo);
  if (fcsc2 < cscThresh) {return false;}

  // update purity corrected energy if values are set to zero
  if (FirstCC->fPurityCorrectedEnergy == 0) {FirstCC->fPurityCorrectedEnergy = fepcl1;}
  if (SecondCC->fPurityCorrectedEnergy == 0) {SecondCC->fPurityCorrectedEnergy = fepcl2;}

  fCL1 = -999;
  fCL2 = fercl2 / fcsc2;
  fpCL1 = -999;
  fpCL2 = fepcl2 / fcsc2;
  fPeakRatio = -999;
  fx1 = FirstCC->fX;
  fy1 = FirstCC->fY;
  fz1 = FirstCC->fZ;
  fx2 = SecondCC->fX;
  fy2 = SecondCC->fY;
  fz2 = SecondCC->fZ;
  fAlg = 2;

  fBetaCC = FirstCC;
  fAlphaCC = SecondCC;
  fBetaSC = FirstSC;
  fAlphaSC = SecondSC;

  return true;
}

// Require two scintillation cluster having a certain peak ratio and the second strong one having a certain absolute value.
// There must be at least one charge cluster after first scintillation cluster. If there is a charge cluster after second scintillation cluster
// that was assigned to first scintillation cluster, a search for a charge cluster matching in x-y plane and reassignment of scintillation cluster
// is performed.
bool EXOBackgroundAnalysis::ThirdAttempt(EXOEventData *ED)
{
  double R = 30; // radius within the position must match [mm]
  double tscDiffTresh = 5; // time difference that the two scintillation cluster must have [us]
  double PeakRatioThresh = 3; // the peak ratio must be greater than PeakRatioThresh
  double csc2Thresh = 10000; // the second scintillation cluster must have at least csc2Thresh photon counts

  int ncl = int(ED->GetNumChargeClusters());
  int nsc = int(ED->GetNumScintillationClusters());

  if (ncl < 1) {return false;} // require at least one charge cluster
  if (nsc < 2) {return false;} // require at least two scintillation cluster

  EXOScintillationCluster *FirstSC = 0; // this is the scintillation cluster of the beta decay
  EXOScintillationCluster *SecondSC = 0; // this is the scintillation cluster of the alpha decay

  if (nWaveFormPeaks < 2 && hasWaveforms) {return false;}
  if (nsc == 2) {
     double TimeSC1 = ED->GetScintillationCluster(0)->fTime;
     double TimeSC2 = ED->GetScintillationCluster(1)->fTime;

     if (TimeSC2 > TimeSC1) {
        FirstSC = ED->GetScintillationCluster(0);
        SecondSC = ED->GetScintillationCluster(1);
     }
     else {
        FirstSC = ED->GetScintillationCluster(1);
        SecondSC = ED->GetScintillationCluster(0);
     }
  }
  else {
     if (!hasWaveforms) {return false;}
     
     // if there are more than two scintillation cluster, compare with sum APD peaks
     if (nsc > 2 && nWaveFormPeaks > 2) {return false;} // there are definitely more than two scintillation cluster
     if (nsc > 2 && nWaveFormPeaks == 2) {

        // find first waveform peak in reconstructed cluster
        int m = 0;
        double TimeSC1 = 0.0;
        bool foundFirst = false;
        while (m < nsc) {
           TimeSC1 = ED->GetScintillationCluster(m)->fTime;
           if (WaveFormPeakX[0] > TMath::Nint(TimeSC1/1000) - 2 && WaveFormPeakX[0] < TMath::Nint(TimeSC1/1000) + 2) {foundFirst = true; break;}
           m++;
        }

        // find second waveform peak in reconstructed cluster
        int n = 0;
        double TimeSC2 = 0.0;
        bool foundSecond = false;
        while (n < nsc) {
           TimeSC2 = ED->GetScintillationCluster(n)->fTime;
           if (WaveFormPeakX[1] > TMath::Nint(TimeSC2/1000) - 2 && WaveFormPeakX[1] < TMath::Nint(TimeSC2/1000) + 2) {foundSecond = true; break;}
           n++;
        }

        if (foundFirst && foundSecond) {
           if (TimeSC2 > TimeSC1) {
              FirstSC = ED->GetScintillationCluster(m);
              SecondSC = ED->GetScintillationCluster(n);
           }
           else {
              FirstSC = ED->GetScintillationCluster(n);
              SecondSC = ED->GetScintillationCluster(m);
           }
        }
        else {return false;}
     }
  }

  // check if there are charges after second scintillation cluster and try to reassign them
  double dist = 0.0;
  double dist2 = 0.0;
  double newZ1 = 0.0;
  double newZ2 = 0.0;
  bool IsZeroZ1 = false;
  bool IsZeroZ2 = false;
  bool found = false;
  int k = 0;
  int i = 0;

  double CL1_X = 0.0;
  double CL1_Y = 0.0;
  double CL1_Z = 0.0;
  double CL1_Time = 0.0;

  double CL2_X = 0.0;
  double CL2_Y = 0.0;
  double CL2_Z = 0.0;
  double CL2_Time = 0.0;

  for (k = 0; k < ncl; k++) {
     double TimeCC = ED->GetChargeCluster(k)->fCollectionTime;
     if (TimeCC < SecondSC->fTime) {continue;} // look for charge cluster after second scintillation cluster

     CL1_X = ED->GetChargeCluster(k)->fX;
     CL1_Y = ED->GetChargeCluster(k)->fY;
     CL1_Z = ED->GetChargeCluster(k)->fZ;
     CL1_Time = ED->GetChargeCluster(k)->fCollectionTime;

     for (i = 0; i < ncl; i++) {
        if (i == k) {continue;}

        CL2_X = ED->GetChargeCluster(i)->fX;
        CL2_Y = ED->GetChargeCluster(i)->fY;
        CL2_Z = ED->GetChargeCluster(i)->fZ;
        CL2_Time = ED->GetChargeCluster(i)->fCollectionTime;

        dist = TMath::Sqrt((CL1_X - CL2_X)*(CL1_X - CL2_X) + (CL1_Y - CL2_Y)*(CL1_Y - CL2_Y));
        if (dist > R) {continue;}

        double dtcl1 = 0.0;
        double dtcl2 = 0.0;
        int uWireCH1 = 0;
        int uWireCH2 = 0;
        if (CL1_Time < CL2_Time) {
           dtcl1 = CL1_Time - FirstSC->fTime;
           dtcl2 = CL2_Time - SecondSC->fTime;

           // get first u-wire signal channel to determine detector half
           if (ED->GetChargeCluster(k)->GetNumUWireSignals() > 0 && ED->GetChargeCluster(i)->GetNumUWireSignals() > 0) {
              uWireCH1 = ED->GetChargeCluster(k)->GetUWireSignalAt(0)->fChannel;
              uWireCH2 = ED->GetChargeCluster(i)->GetUWireSignalAt(0)->fChannel;
           }
        }
        else {
           dtcl1 = CL2_Time - FirstSC->fTime;
           dtcl2 = CL1_Time - SecondSC->fTime;

           // get first u-wire signal channel to determine detector half
           if (ED->GetChargeCluster(k)->GetNumUWireSignals() > 0 && ED->GetChargeCluster(i)->GetNumUWireSignals() > 0) {
              uWireCH2 = ED->GetChargeCluster(k)->GetUWireSignalAt(0)->fChannel;
              uWireCH1 = ED->GetChargeCluster(i)->GetUWireSignalAt(0)->fChannel;
           }
        }

        // calculate new z position of reassigned cluster
        newZ1 = 198.4065 - dtcl1 * driftVelocity;
        newZ2 = 198.4065 - dtcl2 * driftVelocity;

        if (newZ1 < 0) {IsZeroZ1 = true;}
        if (newZ2 < 0) {IsZeroZ2 = true;}

        // assign detector half
        if (uWireCH1 >= 0 && uWireCH1 < 38) {newZ1 *= -1;}
        if (uWireCH2 >= 0 && uWireCH2 < 38) {newZ2 *= -1;}

        dist2 = TMath::Sqrt((CL1_X - CL2_X)*(CL1_X - CL2_X) + (CL1_Y - CL2_Y)*(CL1_Y - CL2_Y) + (newZ1 - newZ2)*(newZ1 - newZ2));
        if (dist2 < R) {found = true; break;}
     }

     if (found) {break;}
  }

  EXOChargeCluster *FirstCC = 0; // this is the charge cluster of the beta decay
  EXOChargeCluster *SecondCC = 0; // this is the charge cluster of the alpha decay

  fnr = ED->fRunNumber;
  fne = ED->fEventNumber;
  if (found) {fdR = dist2;}
  else {fdR = -999;}

  fdt = SecondSC->fTime/1000 - FirstSC->fTime/1000;
  //if (fdt < tscDiffTresh) {return false;}
  if (fdt < 110 || fdt > 980) {return false;} // apply time window cut (110us > dt < 980us)

  fcsc1 = (FirstSC->fCountsSumOnAPDPlaneOne + FirstSC->fCountsSumOnAPDPlaneTwo);
  fcsc2 = (SecondSC->fCountsSumOnAPDPlaneOne + SecondSC->fCountsSumOnAPDPlaneTwo);
  if (fcsc2 < csc2Thresh) {return false;}

  fPeakRatio = fcsc2 / fcsc1;
  //if (fPeakRatio < PeakRatioThresh) {return false;}

  if (found) {
     int alphaCharge = 0;
     if (CL2_Time > CL1_Time) {
        FirstCC = ED->GetChargeCluster(k);
        SecondCC = ED->GetChargeCluster(i);
        alphaCharge = i;
     }
     else {
        FirstCC = ED->GetChargeCluster(i);
        SecondCC = ED->GetChargeCluster(k);
        alphaCharge = k;
     }

     double tclDiff = SecondCC->fCollectionTime / 1000 - FirstCC->fCollectionTime / 1000;
     if (TMath::Abs(fdt - tclDiff) > 10) {return false;}

     fercl1 = FirstCC->fRawEnergy;
     fercl2 = SecondCC->fRawEnergy;
     feccl1 = FirstCC->fCorrectedEnergy;
     feccl2 = SecondCC->fCorrectedEnergy;
     fepcl1 = FirstCC->fCorrectedEnergy * TMath::Exp(FirstCC->fDriftTime / elifetime);
     fepcl2 = SecondCC->fCorrectedEnergy * TMath::Exp(SecondCC->fDriftTime / elifetime);

     // update purity corrected energy if values are set to zero
     if (FirstCC->fPurityCorrectedEnergy == 0) {FirstCC->fPurityCorrectedEnergy = fepcl1;}
     if (SecondCC->fPurityCorrectedEnergy == 0) {SecondCC->fPurityCorrectedEnergy = fepcl2;}

     double sumEnergy = 0.0;
     double sumPEnergy = 0.0;
     for (int i = 0; i < ncl; i++) {
        if (i != alphaCharge) {
           sumEnergy += ED->GetChargeCluster(i)->fCorrectedEnergy;
           sumPEnergy += ED->GetChargeCluster(i)->fCorrectedEnergy * TMath::Exp(ED->GetChargeCluster(i)->fDriftTime / elifetime);
        }
     }

     fCL1 = sumEnergy / fcsc1;
     fCL2 = fercl2 / fcsc2;
     fpCL1 = sumPEnergy / fcsc1;
     fpCL2 = fepcl2 / fcsc2;

     // set redetermined z positions to zero if they were negative
     if (IsZeroZ1) {newZ1 = 0.0;}
     if (IsZeroZ2) {newZ2 = 0.0;}

     fx1 = FirstCC->fX;
     fy1 = FirstCC->fY;
     fz1 = newZ1;
     fx2 = SecondCC->fX;
     fy2 = SecondCC->fY;
     fz2 = newZ2;

     // update new determined z position
     FirstCC->fZ = newZ1;
     SecondCC->fZ = newZ2;

     fAlg = 31;
  }
  else {
     double CalcTSC = 0.0;
     double TimeSC2 = SecondSC->fTime;
     for (int i = 0; i < ncl; i++) { // if we haven't found a position match all charge cluster must have first scintillation cluster as t0
        //double TimeAssocSC = ED->GetChargeCluster(i)->GetScintillationCluster()->fTime;

        EXOChargeCluster *charge_cluster = ED->GetChargeCluster(i);
        if (!charge_cluster) {continue;}

        EXOScintillationCluster *scint_cluster = charge_cluster->GetScintillationCluster();
        if (!scint_cluster) {continue;}

        double TimeAssocSC = scint_cluster->fTime;

        if (TimeAssocSC == TimeSC2) {return false;}
     }

     if (ncl == 1) {
        FirstCC = ED->GetChargeCluster(0);
        fercl1 = ED->GetChargeCluster(0)->fRawEnergy;
        feccl1 = ED->GetChargeCluster(0)->fCorrectedEnergy;
        fepcl1 = ED->GetChargeCluster(0)->fCorrectedEnergy * TMath::Exp(ED->GetChargeCluster(0)->fDriftTime / elifetime);
        fx1 = ED->GetChargeCluster(0)->fX;
        fy1 = ED->GetChargeCluster(0)->fY;
        fz1 = ED->GetChargeCluster(0)->fZ;

        // update purity corrected energy if value are set to zero
        if (FirstCC->fPurityCorrectedEnergy == 0) {FirstCC->fPurityCorrectedEnergy = fepcl1;}
     }
     else {
        fercl1 = -999;
        feccl1 = -999;
        fepcl1 = -999;

        fx1 = -999;
        fy1 = -999;
        fz1 = -999;
     }

     fercl2 = 0;
     feccl2 = 0;
     fepcl2 = 0;

     double sumEnergy = 0.0;
     double sumPEnergy = 0.0;
     for (int i = 0; i < ncl; i++) {
        sumEnergy += ED->GetChargeCluster(i)->fCorrectedEnergy;
        sumPEnergy += ED->GetChargeCluster(i)->fCorrectedEnergy * TMath::Exp(ED->GetChargeCluster(i)->fDriftTime / elifetime);
     }

     fCL1 = sumEnergy / fcsc1;
     fCL2 = 0;
     fpCL1 = sumPEnergy / fcsc1;
     fpCL2 = 0;

     fx2 = -999;
     fy2 = -999;
     fz2 = -999;

     fAlg = 3;
  }

  fBetaCC = FirstCC;
  fAlphaCC = SecondCC;
  fBetaSC = FirstSC;
  fAlphaSC = SecondSC;

  return true;
}

// Require one scintillation cluster above threshold and one charge cluster at t0. Scintillation cluster must be after t0. No charge
// after scintillation cluster
bool EXOBackgroundAnalysis::FourthAttempt(EXOEventData *ED)
{
  double tscDiffTresh = 10; // minimum time difference of scintillation cluster from t0 [us]
  double cscThresh = 10000; // the scintillation cluster must at least have cscThresh photon counts

  int ncl = int(ED->GetNumChargeClusters());
  int nsc = int(ED->GetNumScintillationClusters());

  if (ncl < 1) {return false;} // require at least one charge cluster
  if (nsc < 1) {return false;} // require one scintillation cluster

  EXOScintillationCluster *FirstSC = 0; // this is the scintillation cluster of the beta decay
  EXOScintillationCluster *SecondSC = 0; // this is the scintillation cluster of the alpha decay

  if (nWaveFormPeaks == 0 && hasWaveforms) {return false;} // check if the found scintillation cluster are not noise
  if (nsc == 1) {SecondSC = ED->GetScintillationCluster(0);}
  else {
     if (!hasWaveforms) {return false;}
     
     // if there are more than one scintillation cluster, compare with sum APD peaks
     if (nsc > 1 && nWaveFormPeaks > 1) {return false;} // there are definately more than one scintillation cluster
     if (nsc > 1 && nWaveFormPeaks == 1) {
        // find the one waveform peak in reconstructed cluster
        int k = 0;
        bool foundFirst = false;
        while (k < nsc) {
           double TimeSC = ED->GetScintillationCluster(k)->fTime;
           if (WaveFormPeakX[0] > TMath::Nint(TimeSC/1000) - 2 && WaveFormPeakX[0] < TMath::Nint(TimeSC/1000) + 2) {foundFirst = true; break;}
           k++;
        }
        SecondSC = ED->GetScintillationCluster(k);

        if (!foundFirst) {return false;}
     }
  }

  double TimeSC = SecondSC->fTime;
  //if (TimeSC/1000 >  ED->fEventHeader.fTriggerOffset - 10 && TimeSC/1000 < ED->fEventHeader.fTriggerOffset + 10) {return false;} // the scintillation cluster must not be at t0
  //if (TimeSC/1000 <= ED->fEventHeader.fTriggerOffset) {return false;} // scintillation cluster must be after t0

  // one charge cluster must be at t0
  bool found = false;
  for (int i = 0; i < ncl; i++) {
     double TimeCC = ED->GetChargeCluster(i)->fCollectionTime;
     if (TimeCC < TimeSC) {found = true;}
     if (TimeCC > TimeSC) {return false;}
     //if (TimeCC/1000 > ED->fEventHeader.fTriggerOffset - 10 && TimeCC/1000 < ED->fEventHeader.fTriggerOffset + 10) {found = true;}
  }

  if (!found) {return false;}

  EXOChargeCluster *FirstCC = 0; // this is the charge cluster of the beta decay
  EXOChargeCluster *SecondCC = 0; // this is the charge cluster of the alpha decay

  fnr = ED->fRunNumber;
  fne = ED->fEventNumber;
  fdR = -999;
  fdt = -999;

  if (ncl == 1) {
     FirstCC = ED->GetChargeCluster(0);
     fercl1 = ED->GetChargeCluster(0)->fRawEnergy;
     feccl1 = ED->GetChargeCluster(0)->fCorrectedEnergy;
     fepcl1 = ED->GetChargeCluster(0)->fCorrectedEnergy * TMath::Exp(ED->GetChargeCluster(0)->fDriftTime / elifetime);

     fx1 = ED->GetChargeCluster(0)->fX;
     fy1 = ED->GetChargeCluster(0)->fY;
     fz1 = ED->GetChargeCluster(0)->fZ;

     // update purity corrected energy if values are set to zero
     if (FirstCC->fPurityCorrectedEnergy == 0) {FirstCC->fPurityCorrectedEnergy = fepcl1;}
  }
  else {
     fercl1 = -999;
     feccl1 = -999;
     fepcl1 = -999;

     fx1 = -999;
     fy1 = -999;
     fz1 = -999;
  }

  fercl2 = 0;
  feccl2 = 0;
  fcsc1 = 0;

  fcsc2 = (SecondSC->fCountsSumOnAPDPlaneOne + SecondSC->fCountsSumOnAPDPlaneTwo);
  if (fcsc2 < cscThresh) {return false;}

  fCL1 = -999;
  fCL2 = 0;
  fpCL1 = -999;
  fpCL2 = 0;
  fPeakRatio = -999;
  fx2 = -999;
  fy2 = -999;
  fz2 = -999;
  fAlg = 4;

  fBetaCC = FirstCC;
  fAlphaCC = SecondCC;
  fBetaSC = FirstSC;
  fAlphaSC = SecondSC;

  return true;
}

bool EXOBackgroundAnalysis::FirstSubAttempt(EXOEventData *ED)
{
  double R = 10; // radius within the position must match [mm]
  double tclDiffTresh = 5;
  double cscThresh = 10000.0; // threshold for second scintillation cluster

  int ncl = int(ED->GetNumChargeClusters());
  int nsc = int(ED->GetNumScintillationClusters());

  if (ncl < 2) {return false;} // require at least two charge cluster
  if (nsc < 2) {return false;} // require at least two scintillation cluster

  EXOScintillationCluster *FirstSC = 0; // this is the scintillation cluster of the beta decay
  EXOScintillationCluster *SecondSC = 0; // this is the scintillation cluster of the alhpa decay

  if (nWaveFormPeaks < 2 && hasWaveforms) {return false;}
  if (nsc == 2) {
     double TimeSC1 = ED->GetScintillationCluster(0)->fTime;
     double TimeSC2 = ED->GetScintillationCluster(1)->fTime;

     if (TimeSC2 > TimeSC1) {
        FirstSC = ED->GetScintillationCluster(0);
        SecondSC = ED->GetScintillationCluster(1);
     }
     else {
        FirstSC = ED->GetScintillationCluster(1);
        SecondSC = ED->GetScintillationCluster(0);
     }
  }
  else {
     if (!hasWaveforms) {return false;}
     
     // if there are more than two scintillation cluster, compare with sum APD peaks
     if (nsc > 2 && nWaveFormPeaks > 2) {return false;} // there are definitely more than two scintillation cluster
     if (nsc > 2 && nWaveFormPeaks == 2) {

        // find first waveform peak in reconstructed cluster
        int m = 0;
        double TimeSC1 = 0.0;
        bool foundFirst = false;
        while (m < nsc) {
           TimeSC1 = ED->GetScintillationCluster(m)->fTime;
           if (WaveFormPeakX[0] > TMath::Nint(TimeSC1/1000) - 2 && WaveFormPeakX[0] < TMath::Nint(TimeSC1/1000) + 2) {foundFirst = true; break;}
           m++;
        }

        // find second waveform peak in reconstructed cluster
        int n = 0;
        double TimeSC2 = 0.0;
        bool foundSecond = false;
        while (n < nsc) {
           TimeSC2 = ED->GetScintillationCluster(n)->fTime;
           if (WaveFormPeakX[1] > TMath::Nint(TimeSC2/1000) - 2 && WaveFormPeakX[1] < TMath::Nint(TimeSC2/1000) + 2) {foundSecond = true; break;}
           n++;
        }

        if (foundFirst && foundSecond) {
           if (TimeSC2 > TimeSC1) {
              FirstSC = ED->GetScintillationCluster(m);
              SecondSC = ED->GetScintillationCluster(n);
           }
           else {
              FirstSC = ED->GetScintillationCluster(n);
              SecondSC = ED->GetScintillationCluster(m);
           }
        }
        else {return false;}
     }
  }

  int k = 0;
  int i = 0;
  double dist = 0.0;
  double dU = 0.0;
  double tclDiff = 0.0;
  bool found = false;

  double CL1_U = 0.0;
  double CL1_Z = 0.0;
  double CL1_Time = 0.0;

  double CL2_U = 0.0;
  double CL2_Z = 0.0;
  double CL2_Time = 0.0;

  double distMin = 1000000.0;

  if (FirstSC->GetNumChargeClusters() == 0) {return false;}
  if (!SecondSC) {return false;}
  if (SecondSC->GetNumChargeClusters() != 1) {return false;}

  EXOChargeCluster *FirstCC = 0; // this is the charge cluster of the beta decay
  EXOChargeCluster *SecondCC = SecondSC->GetChargeClusterAt(0); // this is the charge cluster of the alpha decay
  CL2_U = SecondCC->fU;
  CL2_Z = SecondCC->fZ;
  CL2_Time = SecondCC->fCollectionTime;

  int nclFirstCC = FirstSC->GetNumChargeClusters();
  for (int i = 0; i < nclFirstCC; i++) {
     CL1_U = FirstSC->GetChargeClusterAt(i)->fU;
     CL1_Z = FirstSC->GetChargeClusterAt(i)->fZ;
     if (TMath::Sqrt((CL1_U - CL2_U)*(CL1_U - CL2_U) + (CL1_Z - CL2_Z)*(CL1_Z - CL2_Z)) < distMin) {distMin = TMath::Sqrt((CL1_U - CL2_U)*(CL1_U - CL2_U) + (CL1_Z - CL2_Z)*(CL1_Z - CL2_Z)); FirstCC = FirstSC->GetChargeClusterAt(i);}
  }

/*  while (k < ncl - 1) {
     CL1_U = ED->GetChargeCluster(k)->fU;
     CL1_Z = ED->GetChargeCluster(k)->fZ;
     CL1_Time = ED->GetChargeCluster(k)->fCollectionTime;

     for (i = k+1; i < ncl; i++) {
        CL2_U = ED->GetChargeCluster(i)->fU;
        CL2_Z = ED->GetChargeCluster(i)->fZ;
        CL2_Time = ED->GetChargeCluster(i)->fCollectionTime;

        dist = TMath::Abs(CL1_Z - CL2_Z);
        dU = TMath::Abs(CL1_U - CL2_U);
        tclDiff = TMath::Abs(CL1_Time - CL2_Time) / 1000;
        if (dist < R && dU < 20 && tclDiff > tclDiffTresh) {found = true; break;}
     }

     if (found) {break;}
     k++;
  }

  if (!found) {return false;}

  EXOChargeCluster *FirstCC; // this is the charge cluster of the beta decay
  EXOChargeCluster *SecondCC; // this is the charge cluster of the alpha decay

  if (CL2_Time > CL1_Time) {
     FirstCC = ED->GetChargeCluster(k);
     SecondCC = ED->GetChargeCluster(i);
  }
  else {
     FirstCC = ED->GetChargeCluster(i);
     SecondCC = ED->GetChargeCluster(k);
  }*/

  fcsc1 = (FirstSC->fCountsSumOnAPDPlaneOne + FirstSC->fCountsSumOnAPDPlaneTwo);
  fcsc2 = (SecondSC->fCountsSumOnAPDPlaneOne + SecondSC->fCountsSumOnAPDPlaneTwo);

  if (fcsc2 < cscThresh) {return false;}

  //if (FirstSC->GetNumChargeClusters() == 1) {
     //FirstCC = FirstSC->GetChargeClusterAt(0);
     fercl1 = FirstCC->fRawEnergy;
     feccl1 = FirstCC->fCorrectedEnergy;
     fepcl1 = FirstCC->fCorrectedEnergy * TMath::Exp(FirstCC->fDriftTime / elifetime);
     fCL1 = feccl1/fcsc1;
     fpCL1 = fepcl1/fcsc1;
     fx1 = FirstCC->fX;
     fy1 = FirstCC->fY;
     fz1 = FirstCC->fZ;
/*  }
  else {
     fercl1 = -999;
     feccl1 = -999;
     fepcl1 = -999;
     fCL1 = -999;
     fpCL1 = -999;
     fx1 = -999;
     fy1 = -999;
     fz1 = -999;
  }*/
  //if (SecondSC->GetNumChargeClusters() == 1) {
     //SecondCC = SecondSC->GetChargeClusterAt(0);
     fercl2 = SecondCC->fRawEnergy;
     feccl2 = SecondCC->fCorrectedEnergy;
     fepcl2 = SecondCC->fCorrectedEnergy * TMath::Exp(SecondCC->fDriftTime / elifetime);
     fCL2 = feccl2/fcsc2;
     fpCL2 = fepcl2/fcsc2;
     fx2 = SecondCC->fX;
     fy2 = SecondCC->fY;
     fz2 = SecondCC->fZ;
/*  }
  else {
     fercl2 = -999;
     feccl2 = -999;
     fepcl2 = -999;
     fCL2 = -999;
     fpCL2 = -999;
     fx2 = -999;
     fy2 = -999;
     fz2 = -999;
  }*/

  fnr = ED->fRunNumber;
  fne = ED->fEventNumber;
  fdR = 0.0;
  fdt = SecondSC->fTime/1000 - FirstSC->fTime/1000;
  if (fdt < 110 || fdt > 980) {return false;} // apply time window cut (110us > dt < 980)

  fPeakRatio = fcsc2 / fcsc1;
  fAlg = 11;

  return true;
}

bool EXOBackgroundAnalysis::SecondSubAttempt(EXOEventData *ED)
{
  double cscThresh = 10000; // the scintillation cluster must at least have cscThresh photon counts
  double tclDiffTresh = 5;

  int ncl = int(ED->GetNumChargeClusters());
  int nsc = int(ED->GetNumScintillationClusters());

  if (ncl < 2) {return false;} // require at least two charge cluster
  if (nsc == 0) {return false;} // require one scintillation cluster

  EXOScintillationCluster *FirstSC = 0; // this is the scintillation cluster of the beta decay
  EXOScintillationCluster *SecondSC = 0; // this is the scintillation cluster of the alpha decay

  if (nWaveFormPeaks == 0 && hasWaveforms) {return false;}
  if (nsc == 1) {SecondSC = ED->GetScintillationCluster(0);}
  else {
     if (!hasWaveforms) {return false;}
     
     // if there are more than one scintillation cluster, compare with sum APD peaks
     if (nsc > 1 && nWaveFormPeaks > 1) {return false;} // there are definitely more than one scintillation cluster
     if (nsc > 1 && nWaveFormPeaks == 1) { // find the one waveform peak in reconstructed cluster
        int k = 0;
        bool foundFirst = false;
        while (k < nsc) {
           double TimeSC = ED->GetScintillationCluster(k)->fTime;
           if (WaveFormPeakX[0] > TMath::Nint(TimeSC/1000) - 2 && WaveFormPeakX[0] < TMath::Nint(TimeSC/1000) + 2) {foundFirst = true; break;}
           k++;
        }

        if (!foundFirst) {return false;}
        SecondSC = ED->GetScintillationCluster(k);
     }
  }

  double TimeSC = SecondSC->fTime;
  //if (TimeSC/1000 >  ED->fEventHeader.fTriggerOffset - 10 && TimeSC/1000 < ED->fEventHeader.fTriggerOffset + 10) {return false;} // the scintillation cluster must not be at t0

  bool foundCharge = false;
  for (int i = 0; i < ncl; i++) {
     double CL_Time = ED->GetChargeCluster(i)->fCollectionTime;
     if (CL_Time < TimeSC) {foundCharge = true; break;}
     //if (CL_Time/1000 > ED->fEventHeader.fTriggerOffset - 10 && CL_Time/1000 < ED->fEventHeader.fTriggerOffset + 10) {foundCharge = true; break;}
  }

  if (!foundCharge) {return false;} // one charge cluster must be at t0

  int k = 0;
  int i = 0;
  double dist = 0.0;
  double dU = 0.0;
  double tclDiff = 0.0;
  bool found = false;

  double CL1_U = 0.0;
  double CL1_Time = 0.0;

  double CL2_U = 0.0;
  double CL2_Time = 0.0;

  while (k < ncl - 1) {
     CL1_U = ED->GetChargeCluster(k)->fU;
     CL1_Time = ED->GetChargeCluster(k)->fCollectionTime;

     for (i = k+1; i < ncl; i++) {
        CL2_U = ED->GetChargeCluster(i)->fU;
        CL2_Time = ED->GetChargeCluster(i)->fCollectionTime;

        dU = TMath::Abs(CL1_U - CL2_U);
        tclDiff = TMath::Abs(CL1_Time - CL2_Time) / 1000;
        if (dU < 20 && tclDiff > tclDiffTresh) {found = true; break;}
     }

     if (found) {break;}
     k++;
  }

  if (!found) {return false;}

  EXOChargeCluster *FirstCC = 0; // this is the charge cluster of the beta decay
  EXOChargeCluster *SecondCC = 0; // this is the charge cluster of the alpha decay

  if (CL2_Time > CL1_Time) {
     FirstCC = ED->GetChargeCluster(k);
     SecondCC = ED->GetChargeCluster(i);
  }
  else {
     FirstCC = ED->GetChargeCluster(i);
     SecondCC = ED->GetChargeCluster(k);
  }

  fcsc1 = 0;
  fcsc2 = (SecondSC->fCountsSumOnAPDPlaneOne + SecondSC->fCountsSumOnAPDPlaneTwo);
  if (fcsc2 < cscThresh) {return false;}

  fnr = ED->fRunNumber;
  fne = ED->fEventNumber;
  fdR = -999;
  fdt = SecondCC->fCollectionTime/1000.0 - FirstCC->fCollectionTime/1000.0;
  if (fdt < 110 || fdt > 980) {return false;} // apply time window cut (110us > dt < 980us)

  fercl1 = FirstCC->fRawEnergy;
  feccl1 = FirstCC->fCorrectedEnergy;
  fepcl1 = FirstCC->fCorrectedEnergy * TMath::Exp(FirstCC->fDriftTime / elifetime);
  fCL1 = -999;
  fpCL1 = -999;
  fx1 = FirstCC->fX;
  fy1 = FirstCC->fY;
  fz1 = FirstCC->fZ;

  //if (SecondSC->GetNumChargeClusters() != 1) {return false;}

     //SecondCC = SecondSC->GetChargeClusterAt(0);
     fercl2 = SecondCC->fRawEnergy;
     feccl2 = SecondCC->fCorrectedEnergy;
     fepcl2 = SecondCC->fCorrectedEnergy * TMath::Exp(SecondCC->fDriftTime / elifetime);
     fCL2 = fercl2 / fcsc2;
     fpCL2 = fepcl2 / fcsc2;
     fx2 = SecondCC->fX;
     fy2 = SecondCC->fY;
     fz2 = SecondCC->fZ;
/*  }
  else {
     fercl2 = -999;
     feccl2 = -999;
     fepcl2 = -999;
     fCL2 = -999;
     fpCL2 = -999;
     fx2 = -999;
     fy2 = -999;
     fz2 = -999;
  }*/

  fPeakRatio = -999;
  fAlg = 21;

  return true;
}
