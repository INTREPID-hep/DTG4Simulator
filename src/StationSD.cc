#include "DTSimLogger.hh"
#include "StationSD.hh"
#include "DTSimUtils.hh"

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

namespace DTSim
{

StationSD::StationSD(const G4String& name)
  : G4VSensitiveDetector(name)
{
    collectionName.insert("DTSegmentCollection");
}

void StationSD::Initialize(G4HCofThisEvent* hce)
{
    // Create segment collection
    fSegmentCollection = new DTSegmentCollection(SensitiveDetectorName, collectionName[0]);

    // Add this collection in hce
    G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    hce->AddHitsCollection(hcID, fSegmentCollection);

    // Clear the entry map at the beginning of each event
    fEntryStepMap.clear();
}

G4bool StationSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    // Only process muons (PDG = ±13)
    G4int pdg = step->GetTrack()->GetDefinition()->GetPDGEncoding();
    if (std::abs(pdg) != 13) return true;

    G4int trackID = step->GetTrack()->GetTrackID();
    auto preStepPoint = step->GetPreStepPoint();
    auto postStepPoint = step->GetPostStepPoint();
    
    // Check if this is the first time we see this track (entry into station)
    auto it = fEntryStepMap.find(trackID);
    G4bool hasEntry = (it != fEntryStepMap.end());
    
    if (!hasEntry) {
        // First hit in this station - store entry position (copy the value, not pointer!)
        G4ThreeVector entryPos = preStepPoint->GetPosition();
        fEntryStepMap[trackID] = entryPos;
    }
    
    // Check if particle is exiting the station
    G4VPhysicalVolume* postVolume = postStepPoint->GetPhysicalVolume();
    G4bool isExiting = false;
    
    if (postVolume) {
        G4String postName = postVolume->GetLogicalVolume()->GetName();
        // Exiting if next volume is world or yoke (not a station)
        isExiting = G4StrUtil::contains(postName, "world") || G4StrUtil::contains(postName, "Yoke");
    } else {
        // postVolume is null - particle killed or leaving world
        isExiting = true;
    }
    
    if (isExiting && hasEntry) {
        // Particle is exiting and we have stored entry position
        G4ThreeVector entryPosGlobal = it->second;
        G4ThreeVector exitPosGlobal = postStepPoint->GetPosition();
        
        // Get station ID from the current (pre) step volume
        auto touchable = preStepPoint->GetTouchable();
        auto physVol = touchable->GetVolume();
        G4String volumeName = physVol->GetLogicalVolume()->GetName();
        StationID stationID = DecodeStationID(volumeName);
        
         if (!stationID.isValid()) {
             LogWarn("StationSD") << "Failed to decode station ID from " 
                 << volumeName << G4endl;
             fEntryStepMap.erase(it);
             return false;
        }
        
        // Calculate segment properties in global coordinates
        G4ThreeVector posGlobal = (entryPosGlobal + exitPosGlobal) * 0.5;
        G4ThreeVector dirGlobal = (exitPosGlobal - entryPosGlobal).unit();
        
        // Transform to local station coordinates
        auto transform = touchable->GetHistory()->GetTopTransform();
        G4ThreeVector posLocal = transform.TransformPoint(posGlobal);
        G4ThreeVector dirLocal = transform.TransformAxis(dirGlobal);
        
        // Create and fill the segment
        DTSegment* segment = new DTSegment();
        segment->SetStationID(stationID);
        segment->SetLocalPos(posLocal);
        segment->SetLocalDir(dirLocal);
        segment->SetGlobalPos(posGlobal);
        segment->SetGlobalDir(dirGlobal);
        segment->SetEntryPos(entryPosGlobal);
        segment->SetExitPos(exitPosGlobal);
        
        // Add segment to collection
        fSegmentCollection->insert(segment);

        LogDebug("StationSD") << *segment << G4endl;
        
        // Clean up entry from map
        fEntryStepMap.erase(it);
    }

    return true;
}

StationID StationSD::DecodeStationID(const G4String& volumeName) const
{
    // Parse volume name pattern: "Station_W0_Sec1_St1"
    StationID stationID;
    
    stationID.wheel = ExtractIntAfterToken(volumeName, "_W");
    stationID.sector = ExtractIntAfterToken(volumeName, "_Sec");
    stationID.station = ExtractIntAfterToken(volumeName, "_St");
    
    if (!stationID.isValid()) {
        LogWarn("StationSD") << "Failed to decode station ID from " 
               << volumeName << G4endl;
    }
    
    return stationID;
}

void StationSD::EndOfEvent(G4HCofThisEvent*)
{
    // Clear the entry map (handles cases where particles entered but didn't exit)
    if (!fEntryStepMap.empty()) {
        LogDebug("StationSD") << "StationSD: Clearing " << fEntryStepMap.size() 
               << " unmatched entry points at end of event" << G4endl;
        fEntryStepMap.clear();
    }
    // Optional: print summary
    G4int nSegments = fSegmentCollection->entries();
    if (nSegments > 0) {
        LogDebug("StationSD") << "StationSD: Collected " << nSegments << " segments in this event" << G4endl;
    }
}

}
