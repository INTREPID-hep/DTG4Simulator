#include "DriftCellSD.hh"
#include "DTSimUtils.hh"

#include "G4RunManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include <sstream>

#include "DTSimConstants.hh"

namespace DTSim
{

DriftCellSD::DriftCellSD(const G4String& name)
  : G4VSensitiveDetector(name)
{
    collectionName.insert("DriftCellHitsCollection");
}

void DriftCellSD::Initialize(G4HCofThisEvent* hce)
{
    // Create hits collection
    fHitsCollection = new DriftCellHitsCollection(SensitiveDetectorName, collectionName[0]);

    // Add this collection in hce
    G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    hce->AddHitsCollection(hcID, fHitsCollection);
}

G4bool DriftCellSD::ProcessHits(G4Step* step, G4TouchableHistory* history)
{
    // Don't create hit if neutral particle or energy below threshold
    auto charge = step->GetTrack()->GetDefinition()->GetPDGCharge();
    auto edep = step->GetTotalEnergyDeposit();

    if (charge == 0 || edep < DTSim::kMinEnergyDeposit) return true;

    auto preStepPoint = step->GetPreStepPoint();
    auto postStepPoint = step->GetPostStepPoint();

    // Use PreStepPoint to get the volume where the step started (the sensitive volume)
    auto touchable = preStepPoint->GetTouchable();
    auto physVol = touchable->GetVolume();
    auto copyNo = physVol->GetCopyNo();
    auto volumeName = physVol->GetLogicalVolume()->GetName();
    // Decode complete cell identification in one call
    CellID cellID = DecodeCellID(volumeName, copyNo);
    // Check if decoding was successful
    if (!cellID.isValid()) {
        G4cerr << "Warning: Failed to decode cell ID for volume " << volumeName 
               << " copy " << copyNo << G4endl;
        return false;
    }

    // Calculate positions in different coordinate frames
    G4ThreeVector worldPos, cellStationPos, cellLocalPos;
    // middle point of the step
    worldPos = (preStepPoint->GetPosition() + postStepPoint->GetPosition()) * 0.5;
    
    // Transform to Station frame (one level up from DriftCell)
    cellStationPos = touchable->GetHistory()->GetTransform(1).TransformPoint(worldPos);
    
    // Transform to DriftCell's local frame (accounts for cell rotation)
    cellLocalPos = touchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);
    
    // Calculate time with drift
    G4double midTime = (preStepPoint->GetGlobalTime() + postStepPoint->GetGlobalTime()) * 0.5;
    // Use radial distance from wire (cell center)
    G4double radialDistance = sqrt(cellLocalPos.x()*cellLocalPos.x() + cellLocalPos.z()*cellLocalPos.z());

    auto timeDrift = midTime + radialDistance / kDriftVelocity;
        
    // Get other particle properties
    auto pdgID  = step->GetTrack()->GetDefinition()->GetPDGEncoding();
    const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();
    G4int processType = process ? process->GetProcessType() : -999;
    // Get current event ID
    G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

    // Create a new hit and fill it
    DriftCellHit* hit = new DriftCellHit();

    hit->SetEventID(evt);
    hit->SetPDG(pdgID);
    hit->SetCharge(charge);
    hit->SetProcessType(processType);
    hit->SetCellID(cellID);
    hit->SetLocalPos(cellStationPos);  // Position in Station frame
    hit->SetGlobalPos(worldPos);
    // Get cell center position in world coordinates
    G4ThreeVector cellCenter = touchable->GetHistory()->GetTopTransform().Inverse().TransformPoint(G4ThreeVector(0,0,0));
    hit->SetCellCenterPos(cellCenter);
    hit->SetTimeDrift(timeDrift);
    hit->SetEnergyDeposit(edep);

    // Add hit to collection
    fHitsCollection->insert(hit);

    return true;
}

CellID DriftCellSD::DecodeCellID(const G4String& volumeName, G4int copyNo) const
{
    CellID cellID;
    
    // Volume name format: DriftCell_W{wheel}_Sec{sector}_St{station}_SL{sl}_{encoding}
    // Example: DriftCell_W-1_Sec1_St2_SL1_60605959

    // Extract common station identification using utility function
    cellID.wheel = ExtractIntAfterToken(volumeName, "_W");
    cellID.sector = ExtractIntAfterToken(volumeName, "_Sec");
    cellID.station = ExtractIntAfterToken(volumeName, "_St");
    cellID.superLayer = ExtractIntAfterToken(volumeName, "_SL");
    
    // Check if basic parsing succeeded
    if (cellID.wheel == -999 || cellID.sector == -999 || 
        cellID.station == -999 || cellID.superLayer == -999) {
        G4cerr << "Warning: Could not parse basic cell ID from " << volumeName << G4endl;
        return cellID;
    }
    
    // Extract cells-per-layer encoding from volume name
    size_t lastUnderscore = volumeName.rfind('_');
    if (lastUnderscore == std::string::npos) {
        G4cerr << "Warning: Could not find encoding in " << volumeName << G4endl;
        return cellID;
    }
    
    G4String encoding = volumeName.substr(lastUnderscore + 1);
    
    // Parse encoding: each 2 characters = one layer count
    std::vector<G4int> cellsPerLayer;
    for (size_t i = 0; i < encoding.length(); i += 2) {
        if (i + 1 < encoding.length()) {
            G4String countStr = encoding.substr(i, 2);
            G4int count = std::stoi(countStr);
            cellsPerLayer.push_back(count);
        }
    }
    
    if (cellsPerLayer.empty()) {
        G4cerr << "Warning: Empty layer structure from encoding " << encoding << G4endl;
        return cellID;
    }
    
    // Decode layer and wire from copy number
    G4int cumulativeCells = 0;
    for (size_t layerIdx = 0; layerIdx < cellsPerLayer.size(); ++layerIdx) {
        G4int layerStartCopy = cumulativeCells + 1;
        cumulativeCells += cellsPerLayer[layerIdx];
        
        if (copyNo <= cumulativeCells) {
            cellID.layer = layerIdx + 1;
            cellID.wire = copyNo - layerStartCopy + 1;
            return cellID;
        }
    }
    
    G4cerr << "Warning: CopyNo " << copyNo << " out of range for " << volumeName << G4endl;
    return cellID;
}

void DriftCellSD::EndOfEvent(G4HCofThisEvent* hce)
{
    // Optional: print summary or debug information
    G4int nHits = fHitsCollection->entries();
    if (nHits > 0) {
        G4cout << "DriftCellSD: " << nHits << " hits produced" << G4endl;
    }
}

}