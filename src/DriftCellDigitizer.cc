#include "DriftCellDigitizer.hh"

#include "G4DigiManager.hh"
#include "G4RunManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "DriftCellHit.hh"
#include "DTSimConstants.hh"

#include <unordered_map>
#include <vector>

namespace DTSim
{

// Helper for sorting hits by time
static bool CompareHitTimes(const DriftCellHit* a, const DriftCellHit* b) {
    return a->GetTimeDrift() < b->GetTimeDrift();
}

DriftCellDigitizer::DriftCellDigitizer(G4String name)
 : G4VDigitizerModule(name)
{
  collectionName.push_back("DriftCellDigiCollection");
}

DriftCellDigitizer::~DriftCellDigitizer()
{}

void DriftCellDigitizer::Digitize()
{
  // -------------------- Initialization --------------------
  // Create digi collection using the registered collection name
  fDigiCollection = new DriftCellDigiCollection(moduleName, collectionName[0]);

  // Get the DigiManager
  G4DigiManager* digiMan = G4DigiManager::GetDMpointer();

  // Get current event
  const G4Event* currentEvent = G4RunManager::GetRunManager()->GetCurrentEvent();
  if (!currentEvent) {
    G4cerr << "DriftCellDigitizer: No current event" << G4endl;
    StoreDigiCollection(fDigiCollection);
    return;
  }

  // Get hits collection ID (only once)
  if (fHCID < 0) {
    fHCID = digiMan->GetHitsCollectionID("DriftCellHitsCollection");
  }

  if (fHCID < 0) {
    // This is expected when DriftCell SD is disabled - not an error
    if (currentEvent->GetEventID() == 0) {
      G4cout << "DriftCellDigitizer: DriftCellHitsCollection not found (SD may be disabled)" << G4endl;
    }
    StoreDigiCollection(fDigiCollection);
    return;
  }

  // Get hits collection from event
  G4HCofThisEvent* hce = currentEvent->GetHCofThisEvent();
  if (!hce) {
    StoreDigiCollection(fDigiCollection);
    return;
  }

  DriftCellHitsCollection* hitsCollection = 
    static_cast<DriftCellHitsCollection*>(hce->GetHC(fHCID));
  
  if (!hitsCollection) {
    StoreDigiCollection(fDigiCollection);
    return;
  }
  // -------------------- Digitization --------------------
  G4int nHits = hitsCollection->entries();
  
  // Map to group hits by CellID
  // Using unordered_map (hash map) for O(1) access
  std::unordered_map<CellID, std::vector<DriftCellHit*>> cellHitsMap;

  // 1. Group hits by cell
  for (G4int i = 0; i < nHits; i++) {
    DriftCellHit* hit = (*hitsCollection)[i];
    cellHitsMap[hit->GetCellID()].push_back(hit);
  }

  // 2. Process each cell
  for (auto& entry : cellHitsMap) {
    std::vector<DriftCellHit*>& hits = entry.second;

    // Sort hits by time
    std::sort(hits.begin(), hits.end(), CompareHitTimes);

    // Take the first hit (earliest time)
    // TODO: Implement dead time logic here (process hits[1..N])
    DriftCellHit* hit = hits[0];
    
    // Apply detection efficiency
    if (G4UniformRand() > kEfficiency) {
      continue;  // Hit not detected
    }
    
    // Get drift time from hit
    G4double driftTime = hit->GetTimeDrift();
    
    // Apply time resolution smearing (Gaussian)
    G4double smearedTime = driftTime + G4RandGauss::shoot(0., kTimeResolution);
    
    // Ensure non-negative time
    if (smearedTime < 0.) smearedTime = 0.;
    
    // Convert to TDC counts
    G4int tdc = static_cast<G4int>(smearedTime / kTDCResolution);
    
    // Create digi
    DriftCellDigi* digi = new DriftCellDigi();
    digi->SetEventID(hit->GetEventID());
    digi->SetCellID(hit->GetCellID());
    digi->SetTDC(tdc);
    digi->SetGlobalPos(hit->GetCellCenterPos());  // Use cell center for visualization
    
    // Add digi to collection
    fDigiCollection->insert(digi);
  }

  // Store the digi collection
  StoreDigiCollection(fDigiCollection);
  
  if (fDigiCollection->entries() > 0) {
    G4cout << "DriftCellDigitizer: Created " << fDigiCollection->entries() 
           << " digis" << G4endl;
  }
}

}
