#include "DriftCellSD.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include <sstream>

namespace DTSim
{

DriftCellSD::DriftCellSD(const G4String& name)
  : G4VSensitiveDetector(name)
{
}

DriftCellSD::~DriftCellSD()
{
}

void DriftCellSD::Initialize(G4HCofThisEvent* hce)
{
    // Implement initialization logic here
}

void DriftCellSD::EndOfEvent(G4HCofThisEvent* hce)
{
    // Implement end-of-event logic here
}

G4bool DriftCellSD::ProcessHits(G4Step* step, G4TouchableHistory* history)
{
    auto charge = step->GetTrack()->GetDefinition()->GetPDGCharge();
    auto pdgID  = step->GetTrack()->GetDefinition()->GetPDGEncoding();

    if (charge==0.) return true;

    auto preStepPoint = step->GetPreStepPoint();

    auto touchable = step->GetPreStepPoint()->GetTouchable();
    auto physVol = touchable->GetVolume();
    auto copyNo = physVol->GetCopyNo();
    auto volumeName = physVol->GetLogicalVolume()->GetName();

    auto worldPos = preStepPoint->GetPosition();
    auto localPos = touchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);

    // Get position with respect to the center of the cell:
    G4ThreeVector origin(0., 0., 0.);
    auto detectorCenter = touchable->GetHistory()->GetTopTransform().Inverse().TransformPoint(origin);
    auto hitOffset = worldPos - detectorCenter;

    auto timeDrift = GetTimeWithDrift(preStepPoint->GetGlobalTime(), hitOffset);
    
    // Decode complete cell identification in one call
    CellID cellID = DecodeCellID(volumeName, copyNo);
    
    // Check if decoding was successful
    if (!cellID.isValid()) {
        G4cerr << "Warning: Failed to decode cell ID for volume " << volumeName 
               << " copy " << copyNo << G4endl;
        return false;
    }
    
    G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

    // Pretty print using the << operator
    G4cout << "DriftCellSD Hit in Event " << evt << ": "
            << " PDG=" << pdgID
            << " q=" << charge
            << " CellID=" << cellID  // ← Uses the overloaded << operator!
            << " LocalPos=" << localPos
            << " TimeWithDrift=" << timeDrift
            << G4endl;

    // Fill ntuple with all information
    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
    
    analysisManager->FillNtupleIColumn(0, 0, evt);
    analysisManager->FillNtupleIColumn(0, 1, pdgID);
    analysisManager->FillNtupleIColumn(0, 2, charge);
    analysisManager->FillNtupleIColumn(0, 3, cellID.wheel);
    analysisManager->FillNtupleIColumn(0, 4, cellID.sector);
    analysisManager->FillNtupleIColumn(0, 5, cellID.station);
    analysisManager->FillNtupleIColumn(0, 6, cellID.superLayer);
    analysisManager->FillNtupleIColumn(0, 7, cellID.layer);
    analysisManager->FillNtupleIColumn(0, 8, cellID.wire);
    analysisManager->FillNtupleDColumn(0, 9, localPos.x());
    analysisManager->FillNtupleDColumn(0, 10, localPos.y());
    analysisManager->FillNtupleDColumn(0, 11, localPos.z());
    analysisManager->FillNtupleDColumn(0, 12, timeDrift);
    
    analysisManager->AddNtupleRow(0);

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
    // Find the last underscore to get the encoding part
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
            cellID.layer = layerIdx + 1;  // Layer IDs start from 1
            cellID.wire = copyNo - layerStartCopy + 1;  // Wire ID within layer (1-indexed)
            return cellID;  // Success!
        }
    }
    
    // If we reach here, copyNo was out of range
    G4cerr << "Warning: CopyNo " << copyNo << " out of range for " << volumeName << G4endl;
    return cellID;
}

G4double DriftCellSD::GetTimeWithDrift(G4double time, G4ThreeVector distance) const
{ 
    return time + abs(distance.x()) / kDriftVelocity;
}

}