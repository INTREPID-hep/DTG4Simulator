#include "EventAction.hh"

#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4DCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4DigiManager.hh"
#include "G4SystemOfUnits.hh"

#include "RunAction.hh"
#include "DriftCellHit.hh"
#include "DriftCellDigi.hh"
#include "DriftCellDigitizer.hh"
#include "DTSegment.hh"

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
    fillSegmentBranches(event, analysisManager);
    
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
    fRunAction->fHit_TrackID.clear();
    fRunAction->fHit_ParentID.clear();
    fRunAction->fHit_TrackLength.clear();
    fRunAction->fHit_VertexKineticEnergy.clear();
    fRunAction->fHit_VertexPosX.clear();
    fRunAction->fHit_VertexPosY.clear();
    fRunAction->fHit_VertexPosZ.clear();
    
    // Clear digi vectors
    fRunAction->fDigi_Wheel.clear();
    fRunAction->fDigi_Sector.clear();
    fRunAction->fDigi_Station.clear();
    fRunAction->fDigi_SuperLayer.clear();
    fRunAction->fDigi_Layer.clear();
    fRunAction->fDigi_Wire.clear();
    fRunAction->fDigi_TDC.clear();
    fRunAction->fDigi_TrackID.clear();
    
    // Clear gen vectors
    fRunAction->fGen_PDG.clear();
    fRunAction->fGen_Charge.clear();
    fRunAction->fGen_Pt.clear();
    fRunAction->fGen_Eta.clear();
    fRunAction->fGen_Phi.clear();
    fRunAction->fGen_RadEnergy.clear();
    fRunAction->fGen_nSecondaries.clear();
    
    // Clear segment vectors
    fRunAction->fSeg_Wheel.clear();
    fRunAction->fSeg_Sector.clear();
    fRunAction->fSeg_Station.clear();
    fRunAction->fSeg_LocalPosX.clear();
    fRunAction->fSeg_LocalPosY.clear();
    fRunAction->fSeg_LocalPosZ.clear();
    fRunAction->fSeg_LocalDirX.clear();
    fRunAction->fSeg_LocalDirY.clear();
    fRunAction->fSeg_LocalDirZ.clear();
    fRunAction->fSeg_GlobalPosX.clear();
    fRunAction->fSeg_GlobalDirY.clear();
    fRunAction->fSeg_GlobalDirZ.clear();

    // Reset maps
    fRadiatedEnergyMap.clear();
    fSecondaryCountMap.clear();
}

void EventAction::fillHitBranches(const G4Event* event, G4AnalysisManager* analysisManager)
{
    // Get hits collection from the event
    G4HCofThisEvent* hce = event->GetHCofThisEvent();
    if (!hce) {
        analysisManager->FillNtupleIColumn(0, 1, 0);
        return;
    }

    // Get the collection ID for DriftCellHitsCollection
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    G4int hcID = sdManager->GetCollectionID("DriftCellHitsCollection");
    if (hcID < 0) {
        // SD is disabled or not registered
        if (event->GetEventID() == 0) {
            G4cout << "EventAction: DriftCellHitsCollection not found (DriftCell SD may be disabled)" << G4endl;
        }
        analysisManager->FillNtupleIColumn(0, 1, 0);
        return;
    }

    // Retrieve the hits collection
    DriftCellHitsCollection* hitsCollection = 
        static_cast<DriftCellHitsCollection*>(hce->GetHC(hcID));
    
    if (!hitsCollection) {
        analysisManager->FillNtupleIColumn(0, 1, 0);
        return;
    }

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
        fRunAction->fHit_XLocal.push_back(localPos.x()/mm);
        fRunAction->fHit_YLocal.push_back(localPos.y()/mm);
        fRunAction->fHit_ZLocal.push_back(localPos.z()/mm);
        fRunAction->fHit_Time.push_back(hit->GetTimeDrift()/ns);
        fRunAction->fHit_Edep.push_back(hit->GetEnergyDeposit()/keV);
        fRunAction->fHit_ProcessType.push_back(hit->GetProcessType());
        fRunAction->fHit_TrackID.push_back(hit->GetTrackID());
        fRunAction->fHit_ParentID.push_back(hit->GetParentID());
        fRunAction->fHit_TrackLength.push_back(hit->GetTrackLength()/mm);
        fRunAction->fHit_VertexKineticEnergy.push_back(hit->GetVertexKineticEnergy()/keV);
        fRunAction->fHit_VertexPosX.push_back(hit->GetVertexPos().x()/mm);
        fRunAction->fHit_VertexPosY.push_back(hit->GetVertexPos().y()/mm);
        fRunAction->fHit_VertexPosZ.push_back(hit->GetVertexPos().z()/mm);
    }    
}

void EventAction::fillDigiBranches(const G4Event* event, G4AnalysisManager* analysisManager)
{
    // Get digi collection from event
    G4DCofThisEvent* dce = event->GetDCofThisEvent();
    if (!dce) {
        // No digis - fill zero
        analysisManager->FillNtupleIColumn(0, 23, 0);
        return;
    }
    
    G4DigiManager* digiMan = G4DigiManager::GetDMpointer();
    G4int dcID = digiMan->GetDigiCollectionID("DriftCellDigiCollection");
    if (dcID < 0) {
        if (event->GetEventID() == 0) {
            G4cout << "EventAction: DriftCellDigiCollection not found (digitizer may not have run)" << G4endl;
        }
        analysisManager->FillNtupleIColumn(0, 23, 0);
        return;
    }
    
    DriftCellDigiCollection* digiCollection = 
        static_cast<DriftCellDigiCollection*>(dce->GetDC(dcID));
    if (!digiCollection) {
        analysisManager->FillNtupleIColumn(0, 23, 0);
        return;
    }
    
    G4int nDigis = digiCollection->entries();
    
    G4cout << "=== EventAction: " << nDigis << " digis created in event " 
           << event->GetEventID() << " ===" << G4endl;
    
    // Fill scalar column for nDigis
    analysisManager->FillNtupleIColumn(0, 23, nDigis);
    
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
        fRunAction->fDigi_TrackID.push_back(digi->GetTrackID());
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
            fRunAction->fGen_Pt.push_back(pt/GeV);
            fRunAction->fGen_Eta.push_back(eta);
            fRunAction->fGen_Phi.push_back(phi);
            
            // Fill radiation info for this primary
            G4int trackID = primary->GetTrackID();
            fRunAction->fGen_RadEnergy.push_back(fRadiatedEnergyMap[trackID]/keV);
            fRunAction->fGen_nSecondaries.push_back(fSecondaryCountMap[trackID]);
            
            // Move to next primary in this vertex
            primary = primary->GetNext();
        }
    }
    
    // Fill scalar column for number of primaries
    analysisManager->FillNtupleIColumn(0, 32, nPrimaries);
    
    G4cout << "=== EventAction: " << nPrimaries << " primary particles generated in event " 
           << event->GetEventID() << " ===" << G4endl;
}


void EventAction::fillSegmentBranches(const G4Event* event, G4AnalysisManager* analysisManager)
{
    // Get hits collection from the event
    G4HCofThisEvent* hce = event->GetHCofThisEvent();
    if (!hce) {
        // No hits collection - fill zero segments
        analysisManager->FillNtupleIColumn(0, 40, 0);
        return;
    }

    // Get the collection ID for DTSegmentCollection
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    G4int hcID = sdManager->GetCollectionID("DTSegmentCollection");
    if (hcID < 0) {
        // SD is disabled or not registered
        if (event->GetEventID() == 0) {
            G4cout << "EventAction: DTSegmentCollection not found (Station SD may be disabled)" << G4endl;
        }
        analysisManager->FillNtupleIColumn(0, 40, 0);
        return;
    }

    // Retrieve the segment collection
    DTSegmentCollection* segmentCollection = 
        static_cast<DTSegmentCollection*>(hce->GetHC(hcID));
    
    if (!segmentCollection) {
        analysisManager->FillNtupleIColumn(0, 40, 0);
        return;
    }

    G4int nSegments = segmentCollection->entries();
    
    G4cout << "=== EventAction: " << nSegments << " segments collected in event " 
           << event->GetEventID() << " ===" << G4endl;
    
    // Fill scalar column for number of segments
    analysisManager->FillNtupleIColumn(0, 40, nSegments);

    // Loop over all segments and fill RunAction's vectors
    for (G4int i = 0; i < nSegments; i++) {
        DTSegment* segment = (*segmentCollection)[i];
        StationID stationID = segment->GetStationID();
        G4ThreeVector localPos = segment->GetLocalPos();
        G4ThreeVector localDir = segment->GetLocalDir();
        G4ThreeVector globalPos = segment->GetGlobalPos();
        G4ThreeVector globalDir = segment->GetGlobalDir();
        
        // Fill vectors with push_back
        fRunAction->fSeg_Wheel.push_back(stationID.wheel);
        fRunAction->fSeg_Sector.push_back(stationID.sector);
        fRunAction->fSeg_Station.push_back(stationID.station);
        fRunAction->fSeg_LocalPosX.push_back(localPos.x()/mm);
        fRunAction->fSeg_LocalPosY.push_back(localPos.y()/mm);
        fRunAction->fSeg_LocalPosZ.push_back(localPos.z()/mm);
        fRunAction->fSeg_LocalDirX.push_back(localDir.x());
        fRunAction->fSeg_LocalDirY.push_back(localDir.y());
        fRunAction->fSeg_LocalDirZ.push_back(localDir.z());
        fRunAction->fSeg_GlobalPosX.push_back(globalPos.x()/mm);
        fRunAction->fSeg_GlobalPosY.push_back(globalPos.y()/mm);
        fRunAction->fSeg_GlobalPosZ.push_back(globalPos.z()/mm);
        fRunAction->fSeg_GlobalDirX.push_back(globalDir.x());
        fRunAction->fSeg_GlobalDirY.push_back(globalDir.y());
        fRunAction->fSeg_GlobalDirZ.push_back(globalDir.z());
    }    
}

}
