#include "DriftCellSD.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include <sstream>

namespace DTSim
{

DriftCellSD::DriftCellSD(const G4String& name, const G4String& hitsCollectionName)
  : G4VSensitiveDetector(name)
{
    collectionName.insert(hitsCollectionName);
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

    size_t pos = 0;
    
    // Extract wheel
    pos = volumeName.find("_W");
    if (pos != std::string::npos) {
        pos += 2;  // Skip "_W"
        size_t endPos = volumeName.find("_", pos);
        G4String wheelStr = volumeName.substr(pos, endPos - pos);
        cellID.wheel = std::stoi(wheelStr);
        pos = endPos;
    } else {
        G4cerr << "Warning: Could not parse wheel from " << volumeName << G4endl;
        return cellID;
    }
    
    // Extract sector
    pos = volumeName.find("_Sec", pos);
    if (pos != std::string::npos) {
        pos += 4;  // Skip "_Sec"
        size_t endPos = volumeName.find("_", pos);
        G4String sectorStr = volumeName.substr(pos, endPos - pos);
        cellID.sector = std::stoi(sectorStr);
        pos = endPos;
    } else {
        G4cerr << "Warning: Could not parse sector from " << volumeName << G4endl;
        return cellID;
    }
    
    // Extract station
    pos = volumeName.find("_St", pos);
    if (pos != std::string::npos) {
        pos += 3;  // Skip "_St"
        size_t endPos = volumeName.find("_", pos);
        G4String stationStr = volumeName.substr(pos, endPos - pos);
        cellID.station = std::stoi(stationStr);
        pos = endPos;
    } else {
        G4cerr << "Warning: Could not parse station from " << volumeName << G4endl;
        return cellID;
    }
    
    // Extract superlayer
    pos = volumeName.find("_SL", pos);
    if (pos != std::string::npos) {
        pos += 3;  // Skip "_SL"
        size_t endPos = volumeName.find("_", pos);
        G4String slStr = volumeName.substr(pos, endPos - pos);
        cellID.superLayer = std::stoi(slStr);
        pos = endPos;
    } else {
        G4cerr << "Warning: Could not parse superlayer from " << volumeName << G4endl;
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
    G4int nHits = fHitsCollection->entries();
    
    if (nHits > 0) {
        G4cout << "\n=== DriftCellSD: " << nHits << " hits collected in this event ===" << G4endl;
        
        G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
        
        for (G4int i = 0; i < nHits; i++) {
            DriftCellHit* hit = (*fHitsCollection)[i];
            CellID cellID = hit->GetCellID();
            G4ThreeVector localPos = hit->GetLocalPos();
            
            // Fill ntuple
            analysisManager->FillNtupleIColumn(0, 0, hit->GetEventID());
            analysisManager->FillNtupleIColumn(0, 1, hit->GetPDG());
            analysisManager->FillNtupleIColumn(0, 2, hit->GetCharge());
            analysisManager->FillNtupleIColumn(0, 3, cellID.wheel);
            analysisManager->FillNtupleIColumn(0, 4, cellID.sector);
            analysisManager->FillNtupleIColumn(0, 5, cellID.station);
            analysisManager->FillNtupleIColumn(0, 6, cellID.superLayer);
            analysisManager->FillNtupleIColumn(0, 7, cellID.layer);
            analysisManager->FillNtupleIColumn(0, 8, cellID.wire);
            analysisManager->FillNtupleDColumn(0, 9, localPos.x());
            analysisManager->FillNtupleDColumn(0, 10, localPos.y());
            analysisManager->FillNtupleDColumn(0, 11, localPos.z());
            analysisManager->FillNtupleDColumn(0, 12, hit->GetTimeDrift());
            analysisManager->FillNtupleDColumn(0, 13, hit->GetEnergyDeposit());
            analysisManager->FillNtupleIColumn(0, 14, hit->GetProcessType());
            
            analysisManager->AddNtupleRow(0);
        }
    }
}

}