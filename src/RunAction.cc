#include "RunAction.hh"
#include "G4Run.hh"
#include "G4AnalysisManager.hh"

namespace DTSim
{
RunAction::RunAction()
 : G4UserRunAction()
{
    auto analysisManager = G4AnalysisManager::Instance();
    
    // Default settings
    analysisManager->SetNtupleDirectoryName("DTG4SimNTuple");
    analysisManager->SetNtupleMerging(true);
    
    // Create ntuple for DTG4Sim hits - one row per event with vector branches
    analysisManager->CreateNtuple("DTG4Tree", "DTG4Tree");
    analysisManager->CreateNtupleIColumn("event_eventNumber");
    analysisManager->CreateNtupleIColumn("g4dtSimHit_nSimHits");
    analysisManager->CreateNtupleIColumn("g4dtSimHit_PDG", fHit_PDG);
    analysisManager->CreateNtupleIColumn("g4dtSimHit_q", fHit_Charge);
    analysisManager->CreateNtupleIColumn("g4dtSimHit_wheel", fHit_Wheel);
    analysisManager->CreateNtupleIColumn("g4dtSimHit_sector", fHit_Sector);
    analysisManager->CreateNtupleIColumn("g4dtSimHit_station", fHit_Station);
    analysisManager->CreateNtupleIColumn("g4dtSimHit_superlayer", fHit_SuperLayer);
    analysisManager->CreateNtupleIColumn("g4dtSimHit_layer", fHit_Layer);
    analysisManager->CreateNtupleIColumn("g4dtSimHit_cell", fHit_Wire);
    analysisManager->CreateNtupleDColumn("g4dtSimHit_xlocal", fHit_XLocal);
    analysisManager->CreateNtupleDColumn("g4dtSimHit_ylocal", fHit_YLocal);
    analysisManager->CreateNtupleDColumn("g4dtSimHit_zlocal", fHit_ZLocal);
    analysisManager->CreateNtupleDColumn("g4dtSimHit_time", fHit_Time);
    analysisManager->CreateNtupleDColumn("g4dtSimHit_edep", fHit_Edep);
    analysisManager->CreateNtupleIColumn("g4dtSimHit_process_type", fHit_ProcessType);
    analysisManager->FinishNtuple();
}

void RunAction::BeginOfRunAction(const G4Run* run)
{
    auto analysisManager = G4AnalysisManager::Instance();
    
    G4int runID = run->GetRunID();
    G4cout << "### Run " << runID << " start." << G4endl;
    
    
    // Open file first
    analysisManager->OpenFile("DTG4Simulation_" + std::to_string(runID) + ".root");
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