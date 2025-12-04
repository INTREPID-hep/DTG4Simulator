#include "EventAction.hh"
#include "RunAction.hh"
#include "DriftCellHit.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4AnalysisManager.hh"

namespace DTSim
{

EventAction::EventAction(RunAction* runAction)
 : G4UserEventAction(), fRunAction(runAction)
{
}

void EventAction::BeginOfEventAction(const G4Event* event)
{
    // Clear vectors at the beginning of each event
    clearVectors();
}

void EventAction::EndOfEventAction(const G4Event* event)
{
    // Fill hits ntuple
    FillHitsNtuple(event);
}

void EventAction::clearVectors()
{
    fRunAction->fHit_PDG.clear();
    fRunAction->fHit_Charge.clear();
    fRunAction->fHit_Wheel.clear();
    fRunAction->fHit_Sector.clear();
    fRunAction->fHit_Station.clear();
    fRunAction->fHit_SuperLayer.clear();
    fRunAction->fHit_Layer.clear();
    fRunAction->fHit_Wire.clear();
    fRunAction->fHit_XLocal.clear();
    fRunAction->fHit_YLocal.clear();
    fRunAction->fHit_ZLocal.clear();
    fRunAction->fHit_Time.clear();
    fRunAction->fHit_Edep.clear();
    fRunAction->fHit_ProcessType.clear();
}

void EventAction::FillHitsNtuple(const G4Event* event)
{
    // Get hits collection from the event
    G4HCofThisEvent* hce = event->GetHCofThisEvent();
    if (!hce) return;

    // Get the collection ID for DriftCellHitsCollection
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    G4int hcID = sdManager->GetCollectionID("DriftCellHitsCollection");
    if (hcID < 0) return;

    // Retrieve the hits collection
    DriftCellHitsCollection* hitsCollection = 
        static_cast<DriftCellHitsCollection*>(hce->GetHC(hcID));
    
    if (!hitsCollection) return;

    G4int nHits = hitsCollection->entries();
    
    G4cout << "\n=== EventAction: " << nHits << " hits collected in event " 
           << event->GetEventID() << " ===" << G4endl;
    
    // Loop over all hits and fill RunAction's vectors
    for (G4int i = 0; i < nHits; i++) {
        DriftCellHit* hit = (*hitsCollection)[i];
        CellID cellID = hit->GetCellID();
        G4ThreeVector localPos = hit->GetLocalPos();
        
        // Fill vectors with push_back - no need to call FillNtuple for vectors!
        fRunAction->fHit_PDG.push_back(hit->GetPDG());
        fRunAction->fHit_Charge.push_back(hit->GetCharge());
        fRunAction->fHit_Wheel.push_back(cellID.wheel);
        fRunAction->fHit_Sector.push_back(cellID.sector);
        fRunAction->fHit_Station.push_back(cellID.station);
        fRunAction->fHit_SuperLayer.push_back(cellID.superLayer);
        fRunAction->fHit_Layer.push_back(cellID.layer);
        fRunAction->fHit_Wire.push_back(cellID.wire);
        fRunAction->fHit_XLocal.push_back(localPos.x());
        fRunAction->fHit_YLocal.push_back(localPos.y());
        fRunAction->fHit_ZLocal.push_back(localPos.z());
        fRunAction->fHit_Time.push_back(hit->GetTimeDrift());
        fRunAction->fHit_Edep.push_back(hit->GetEnergyDeposit());
        fRunAction->fHit_ProcessType.push_back(hit->GetProcessType());
    }
    
    // Fill only scalar columns (event number and nHits)
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleIColumn(0, event->GetEventID());
    analysisManager->FillNtupleIColumn(1, nHits);
    
    // Add one row per event - vectors are automatically saved!
    analysisManager->AddNtupleRow(0);
}

}
