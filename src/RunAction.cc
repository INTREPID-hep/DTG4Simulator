#include "RunAction.hh"
#include "G4Run.hh"
#include "G4AnalysisManager.hh"

namespace DTSim
{
RunAction::RunAction()
 : G4UserRunAction()
{
    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
    // Default settings
    analysisManager->SetNtupleDirectoryName("DTG4SimNTuple");
    analysisManager->SetNtupleMerging(true);
    // Create ntuple for DTG4Sim hits
    analysisManager->CreateNtuple("DTG4Tree", "DTG4Tree");
    analysisManager->CreateNtupleIColumn("g4dtsimHit_eventNumber");
    analysisManager->CreateNtupleIColumn("g4dtsimHit_PDG");
    analysisManager->CreateNtupleIColumn("g4dtsimHit_q");
    analysisManager->CreateNtupleIColumn("g4dtsimHit_wheel");
    analysisManager->CreateNtupleIColumn("g4dtsimHit_sector");
    analysisManager->CreateNtupleIColumn("g4dtsimHit_station");
    analysisManager->CreateNtupleIColumn("g4dtsimHit_superlayer");
    analysisManager->CreateNtupleIColumn("g4dtsimHit_layer");
    analysisManager->CreateNtupleIColumn("g4dtsimHit_cell");
    analysisManager->CreateNtupleDColumn("g4dtsimHit_xlocal");
    analysisManager->CreateNtupleDColumn("g4dtsimHit_ylocal");
    analysisManager->CreateNtupleDColumn("g4dtsimHit_zlocal");
    analysisManager->CreateNtupleDColumn("g4dtsimHit_timewithdrift");
    analysisManager->CreateNtupleDColumn("g4dtsimHit_edep");  // ← Add this
    analysisManager->FinishNtuple();
}
RunAction::~RunAction()
{
}

void RunAction::BeginOfRunAction(const G4Run* run)
{
    auto analysisManager = G4AnalysisManager::Instance();

    G4int runID = run->GetRunID();
    G4cout << "### Run " << runID << " start." << G4endl;
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