#include "DriftCellSD.hh"
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

    // Get position with respect of the center of the cell:
    G4ThreeVector origin(0., 0., 0.);
    auto detectorCenter = touchable->GetHistory()->GetTopTransform().Inverse().TransformPoint(origin);
    auto hitOffset = worldPos - detectorCenter;

    auto timeDrift = GetTimeWithDrift(preStepPoint->GetGlobalTime(), hitOffset);
    
    // Extract layer structure from volume Name
    ExtractLayerStructure(volumeName);
    
    // Decode both IDs
    G4int layerID = -1;
    G4int cellID = -1;
    DecodeCopyNumber(copyNo, layerID, cellID);

    G4cout << "DriftCellSD Hit: "
            << " PDG=" << pdgID
            << " q=" << charge
            << " CopyNo=" << copyNo
            << " LayerID=" << layerID
            << " CellID=" << cellID
            << " LocalPos=" << localPos
            << " TimeWithDrift=" << timeDrift
            << G4endl;

    return true;
}

void DriftCellSD::ExtractLayerStructure(const G4String& volumeName) const
{
    // Check if we already have the structure for this volume 
    if (fVolumeName == volumeName && !fCellsPerLayer.empty()) {
        return; // Already cached
    }
    
    // Volume name format: DriftCell_W{wheel}_Sec{sector}_St{station}_SL{sl}_{encoding}
    // Encoding format: 2 digits per layer (e.g., "60605959" = [60, 60, 59, 59])
    
    fCellsPerLayer.clear();
    fVolumeName = volumeName;
    
    // Find the last underscore to get the encoding part
    size_t lastUnderscore = volumeName.rfind('_');
    if (lastUnderscore == std::string::npos) {
        G4cerr << "Warning: Could not parse volume name: " << volumeName << G4endl;
        return;
    }
    
    G4String encoding = volumeName.substr(lastUnderscore + 1);
    
    // Parse encoding: each 2 characters = one layer count
    for (size_t i = 0; i < encoding.length(); i += 2) {
        if (i + 1 < encoding.length()) {
            G4String countStr = encoding.substr(i, 2);
            G4int count = std::stoi(countStr);
            fCellsPerLayer.push_back(count);
        }
    }
}

void DriftCellSD::DecodeCopyNumber(G4int copyNo, G4int& layerID, G4int& cellID) const
{
    // Initialize to invalid values
    layerID = -1;
    cellID = -1;
    
    if (fCellsPerLayer.empty()) {
        G4cerr << "Warning: Empty cached layer structure" << G4endl;
        return;
    }
    
    // Find which layer this copyNo belongs to and the cell within that layer
    G4int cumulativeCells = 0;
    for (size_t layerIdx = 0; layerIdx < fCellsPerLayer.size(); ++layerIdx) {
        G4int layerStartCopy = cumulativeCells + 1;
        cumulativeCells += fCellsPerLayer[layerIdx];
        
        if (copyNo <= cumulativeCells) {
            layerID = layerIdx + 1;  // Layer IDs start from 1
            cellID = copyNo - layerStartCopy + 1;  // Cell ID within layer (1-indexed)
            return;  // Found it, exit early
        }
    }
    
    // If we reach here, copyNo was out of range
    G4cerr << "Warning: CopyNo " << copyNo << " out of range" << G4endl;
}

G4double DriftCellSD::GetTimeWithDrift(G4double time, G4ThreeVector distance) const
{ 
    return time + abs(distance.x()) / kDriftVelocity; //This is in the correct units  (ns)
}

}