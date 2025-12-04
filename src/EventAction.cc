#include "EventAction.hh"
#include "RunAction.hh"
#include "DriftCellHit.hh"
#include "DriftCellDigi.hh"
#include "DriftCellDigitizer.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4DCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4DigiManager.hh"
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
    // Trigger digitization
    G4DigiManager* digiMan = G4DigiManager::GetDMpointer();
    DriftCellDigitizer* digitizer = 
    static_cast<DriftCellDigitizer*>(digiMan->FindDigitizerModule("DriftCellDigitizer"));
    
    if (digitizer) {
        digitizer->Digitize();
    }

    //get analysis manager
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    
    // Fill event ID
    analysisManager->FillNtupleIColumn(0, 0, event->GetEventID());

    // Fill generator, hits and digis branches
    fillGenBranches(event, analysisManager);
    fillHitBranches(event, analysisManager);
    fillDigiBranches(event, analysisManager);
    
    // Write ntuple row
    analysisManager->AddNtupleRow(0);
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
    
    // Clear digi vectors
    fRunAction->fDigi_Wheel.clear();
    fRunAction->fDigi_Sector.clear();
    fRunAction->fDigi_Station.clear();
    fRunAction->fDigi_SuperLayer.clear();
    fRunAction->fDigi_Layer.clear();
    fRunAction->fDigi_Wire.clear();
    fRunAction->fDigi_TDC.clear();
    
    // Clear gen vectors
    fRunAction->fGen_PDG.clear();
    fRunAction->fGen_Charge.clear();
    fRunAction->fGen_Pt.clear();
    fRunAction->fGen_Eta.clear();
    fRunAction->fGen_Phi.clear();
}

void EventAction::fillHitBranches(const G4Event* event, G4AnalysisManager* analysisManager)
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
    
    G4cout << "=== EventAction: " << nHits << " hits collected in event " 
           << event->GetEventID() << " ===" << G4endl;
    
    // Fill number of hits
    analysisManager->FillNtupleIColumn(0, 1, nHits);

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
}

void EventAction::fillDigiBranches(const G4Event* event, G4AnalysisManager* analysisManager)
{
    // Get digi collection from event
    G4DCofThisEvent* dce = event->GetDCofThisEvent();
    if (!dce) {
        // No digis - fill zero
        analysisManager->FillNtupleIColumn(0, 16, 0);
        return;
    }
    
    G4DigiManager* digiMan = G4DigiManager::GetDMpointer();
    G4int dcID = digiMan->GetDigiCollectionID("DriftCellDigiCollection");
    if (dcID < 0) {
        analysisManager->FillNtupleIColumn(0, 16, 0);
        return;
    }
    
    DriftCellDigiCollection* digiCollection = 
        static_cast<DriftCellDigiCollection*>(dce->GetDC(dcID));
    if (!digiCollection) {
        analysisManager->FillNtupleIColumn(0, 16, 0);
        return;
    }
    
    G4int nDigis = digiCollection->entries();
    
    G4cout << "=== EventAction: " << nDigis << " digis created in event " 
           << event->GetEventID() << " ===" << G4endl;
    
    // Fill scalar column for nDigis (column 16)
    analysisManager->FillNtupleIColumn(0, 16, nDigis);
    
    // Loop over all digis and fill RunAction's vectors
    for (G4int i = 0; i < nDigis; i++) {
        DriftCellDigi* digi = (*digiCollection)[i];
        CellID cellID = digi->GetCellID();
        
        // Fill vectors with push_back
        fRunAction->fDigi_Wheel.push_back(cellID.wheel);
        fRunAction->fDigi_Sector.push_back(cellID.sector);
        fRunAction->fDigi_Station.push_back(cellID.station);
        fRunAction->fDigi_SuperLayer.push_back(cellID.superLayer);
        fRunAction->fDigi_Layer.push_back(cellID.layer);
        fRunAction->fDigi_Wire.push_back(cellID.wire);
        fRunAction->fDigi_TDC.push_back(digi->GetTDC());
    }
}

void EventAction::fillGenBranches(const G4Event* event, G4AnalysisManager* analysisManager)
{
    G4int nPrimaries = 0;
    
    // Loop over all primary vertices
    G4int nVertices = event->GetNumberOfPrimaryVertex();
    for (G4int iVtx = 0; iVtx < nVertices; iVtx++) {
        G4PrimaryVertex* vertex = event->GetPrimaryVertex(iVtx);
        if (!vertex) continue;
        
        // Loop over all primaries in this vertex
        G4PrimaryParticle* primary = vertex->GetPrimary();
        while (primary) {
            nPrimaries++;
            
            // Get momentum components
            G4double px = primary->GetPx();
            G4double py = primary->GetPy();
            G4double pz = primary->GetPz();
            
            // Calculate pt, eta, phi
            G4double pt = std::sqrt(px*px + py*py);
            G4double p = std::sqrt(px*px + py*py + pz*pz);
            G4double eta = 0.5 * std::log((p + pz) / (p - pz));
            G4double phi = std::atan2(py, px);
            
            // Fill vectors
            fRunAction->fGen_PDG.push_back(primary->GetPDGcode());
            fRunAction->fGen_Charge.push_back(primary->GetCharge());
            fRunAction->fGen_Pt.push_back(pt);
            fRunAction->fGen_Eta.push_back(eta);
            fRunAction->fGen_Phi.push_back(phi);
            
            // Move to next primary in this vertex
            primary = primary->GetNext();
        }
    }
    
    // Fill scalar column for number of primaries (column 24)
    analysisManager->FillNtupleIColumn(0, 24, nPrimaries);
    
    G4cout << "=== EventAction: " << nPrimaries << " primary particles generated in event " 
           << event->GetEventID() << " ===" << G4endl;
}

}
