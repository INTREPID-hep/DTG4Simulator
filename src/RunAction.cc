#include "RunAction.hh"

#include "G4Run.hh"
#include "G4AnalysisManager.hh"
#include "G4GenericMessenger.hh"

#include "DTSimConstants.hh"

namespace DTSim
{
RunAction::RunAction()
 : G4UserRunAction(),
   fMessenger(nullptr),
   fOutputFileName(DTSim::kOutputFileName)
{
    auto analysisManager = G4AnalysisManager::Instance();
    
    // Default settings
    analysisManager->SetNtupleDirectoryName("DTG4SimNTuple");
    analysisManager->SetNtupleMerging(true);

    // Define commands
    fMessenger = new G4GenericMessenger(this, "/DTSim/run/", "Run control");
    fMessenger->DeclareProperty("setOutputFileName", fOutputFileName, "Set output file name (without extension)");
    
    // Create ntuple for DTG4Sim hits - one row per event with vector branches
    analysisManager->CreateNtuple("DTG4Tree", "DTG4Tree");
    analysisManager->CreateNtupleIColumn("event_eventNumber");
    analysisManager->CreateNtupleIColumn("simHit_nSimHits");
    analysisManager->CreateNtupleIColumn("simHit_PDG", fHit_PDG);
    analysisManager->CreateNtupleIColumn("simHit_q", fHit_Charge);
    analysisManager->CreateNtupleIColumn("simHit_wheel", fHit_Wheel);
    analysisManager->CreateNtupleIColumn("simHit_sector", fHit_Sector);
    analysisManager->CreateNtupleIColumn("simHit_station", fHit_Station);
    analysisManager->CreateNtupleIColumn("simHit_superlayer", fHit_SuperLayer);
    analysisManager->CreateNtupleIColumn("simHit_layer", fHit_Layer);
    analysisManager->CreateNtupleIColumn("simHit_cell", fHit_Wire);
    analysisManager->CreateNtupleDColumn("simHit_xlocal", fHit_XLocal);
    analysisManager->CreateNtupleDColumn("simHit_ylocal", fHit_YLocal);
    analysisManager->CreateNtupleDColumn("simHit_zlocal", fHit_ZLocal);
    analysisManager->CreateNtupleDColumn("simHit_time", fHit_Time);
    analysisManager->CreateNtupleDColumn("simHit_edep", fHit_Edep);
    analysisManager->CreateNtupleIColumn("simHit_process_type", fHit_ProcessType);
    analysisManager->CreateNtupleIColumn("simHit_trackId", fHit_TrackID);
    analysisManager->CreateNtupleIColumn("simHit_parentId", fHit_ParentID);
    analysisManager->CreateNtupleDColumn("simHit_trackLength", fHit_TrackLength);
    analysisManager->CreateNtupleDColumn("simHit_vertexKineticEnergy", fHit_VertexKineticEnergy);
    analysisManager->CreateNtupleDColumn("simHit_vertexPosX", fHit_VertexPosX);
    analysisManager->CreateNtupleDColumn("simHit_vertexPosY", fHit_VertexPosY);
    analysisManager->CreateNtupleDColumn("simHit_vertexPosZ", fHit_VertexPosZ);
    // Digi columns
    analysisManager->CreateNtupleIColumn("digi_nDigis");
    analysisManager->CreateNtupleIColumn("digi_wheel", fDigi_Wheel);
    analysisManager->CreateNtupleIColumn("digi_sector", fDigi_Sector);
    analysisManager->CreateNtupleIColumn("digi_station", fDigi_Station);
    analysisManager->CreateNtupleIColumn("digi_superlayer", fDigi_SuperLayer);
    analysisManager->CreateNtupleIColumn("digi_layer", fDigi_Layer);
    analysisManager->CreateNtupleIColumn("digi_cell", fDigi_Wire);
    analysisManager->CreateNtupleIColumn("digi_TDC", fDigi_TDC);
    analysisManager->CreateNtupleIColumn("digi_trackId", fDigi_TrackID);
    // Generator columns
    analysisManager->CreateNtupleIColumn("gen_nGenParts");
    analysisManager->CreateNtupleIColumn("gen_pdgId", fGen_PDG);
    analysisManager->CreateNtupleIColumn("gen_charge", fGen_Charge);
    analysisManager->CreateNtupleDColumn("gen_pt", fGen_Pt);
    analysisManager->CreateNtupleDColumn("gen_eta", fGen_Eta);
    analysisManager->CreateNtupleDColumn("gen_phi", fGen_Phi);
    // Segment columns
    analysisManager->CreateNtupleIColumn("seg_nSegments");
    analysisManager->CreateNtupleIColumn("seg_wheel", fSeg_Wheel);
    analysisManager->CreateNtupleIColumn("seg_sector", fSeg_Sector);
    analysisManager->CreateNtupleIColumn("seg_station", fSeg_Station);
    analysisManager->CreateNtupleDColumn("seg_localPosX", fSeg_LocalPosX);
    analysisManager->CreateNtupleDColumn("seg_localPosY", fSeg_LocalPosY);
    analysisManager->CreateNtupleDColumn("seg_localPosZ", fSeg_LocalPosZ);
    analysisManager->CreateNtupleDColumn("seg_localDirX", fSeg_LocalDirX);
    analysisManager->CreateNtupleDColumn("seg_localDirY", fSeg_LocalDirY);
    analysisManager->CreateNtupleDColumn("seg_localDirZ", fSeg_LocalDirZ);
    analysisManager->CreateNtupleDColumn("seg_globalPosX", fSeg_GlobalPosX);
    analysisManager->CreateNtupleDColumn("seg_globalPosY", fSeg_GlobalPosY);
    analysisManager->CreateNtupleDColumn("seg_globalPosZ", fSeg_GlobalPosZ);
    analysisManager->CreateNtupleDColumn("seg_globalDirX", fSeg_GlobalDirX);
    analysisManager->CreateNtupleDColumn("seg_globalDirY", fSeg_GlobalDirY);
    analysisManager->CreateNtupleDColumn("seg_globalDirZ", fSeg_GlobalDirZ);
    analysisManager->FinishNtuple();
}

RunAction::~RunAction()
{
    delete fMessenger;
}

void RunAction::BeginOfRunAction(const G4Run* run)
{
    auto analysisManager = G4AnalysisManager::Instance();
    
    G4int runID = run->GetRunID();
    G4cout << "### Run " << runID << " start." << G4endl;
    
    
    // Open file first
    G4String fileName = fOutputFileName;
    // If using default name, append run ID to avoid overwriting in sequential runs
    if (fileName == DTSim::kOutputFileName) {
        fileName += "_" + std::to_string(runID);
    }
    
    analysisManager->OpenFile(fileName + ".root");
}

void RunAction::EndOfRunAction(const G4Run* run)
{
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();

    G4int runID = run->GetRunID();
    G4cout << "### Run " << runID << " end." << G4endl;
}

}