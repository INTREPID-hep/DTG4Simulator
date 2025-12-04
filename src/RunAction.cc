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
    // Digi columns
    analysisManager->CreateNtupleIColumn("digi_nDigis");
    analysisManager->CreateNtupleIColumn("digi_wheel", fDigi_Wheel);
    analysisManager->CreateNtupleIColumn("digi_sector", fDigi_Sector);
    analysisManager->CreateNtupleIColumn("digi_station", fDigi_Station);
    analysisManager->CreateNtupleIColumn("digi_superlayer", fDigi_SuperLayer);
    analysisManager->CreateNtupleIColumn("digi_layer", fDigi_Layer);
    analysisManager->CreateNtupleIColumn("digi_cell", fDigi_Wire);
    analysisManager->CreateNtupleIColumn("digi_TDC", fDigi_TDC);
    // Generator columns
    analysisManager->CreateNtupleIColumn("gen_nGenParts");
    analysisManager->CreateNtupleIColumn("gen_pdgId", fGen_PDG);
    analysisManager->CreateNtupleIColumn("gen_charge", fGen_Charge);
    analysisManager->CreateNtupleDColumn("gen_pt", fGen_Pt);
    analysisManager->CreateNtupleDColumn("gen_eta", fGen_Eta);
    analysisManager->CreateNtupleDColumn("gen_phi", fGen_Phi);
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