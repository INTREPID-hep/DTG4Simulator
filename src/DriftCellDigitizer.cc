#include "DTSimLogger.hh"
#include "DriftCellDigitizer.hh"

#include "G4DigiManager.hh"
#include "G4RunManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"
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
 : G4VDigitizerModule(name),
   fMessenger(nullptr),
   fEfficiency(DTSim::kEfficiency),
   fTimeResolution(DTSim::kTimeResolution),
   fTDCResolution(DTSim::kTDCResolution)
{
  collectionName.push_back("DriftCellDigiCollection");
  DefineCommands();
}

DriftCellDigitizer::~DriftCellDigitizer()
{
  delete fMessenger;
}

void DriftCellDigitizer::DefineCommands()
{
    fMessenger = new G4GenericMessenger(this, "/DTSim/digitizer/",
                                        "Digitization response parameters");
    
    fMessenger->DeclareProperty("setEfficiency", fEfficiency,
                                "Set detection efficiency (0.0 to 1.0)");
    
    fMessenger->DeclarePropertyWithUnit("setTimeResolution", "ns",
                                        fTimeResolution,
                                        "Set time resolution (sigma)");
    
    fMessenger->DeclarePropertyWithUnit("setTDCResolution", "ns",
                                        fTDCResolution,
                                        "Set TDC bin resolution");
}

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
    LogError("DriftCellDigitizer") << "No current event" << G4endl;
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
      LogDebug("DriftCellDigitizer") << "DriftCellHitsCollection not found (SD may be disabled)" << G4endl;
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
  LogDebug("DriftCellDigitizer") << "Processing " << nHits << " hits for digitization" << G4endl;
  
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
    LogDebug("DriftCellDigitizer") << "Digitizing cell: " << entry.first << ", hits: " << entry.second.size() << G4endl;
    std::vector<DriftCellHit*>& hits = entry.second;

    // Sort hits by time
    std::sort(hits.begin(), hits.end(), CompareHitTimes);

    // Take the first hit (earliest time)
    // TODO: Implement dead time logic here (process hits[1..N])
    DriftCellHit* hit = hits[0];
    
    // Apply detection efficiency
    if (G4UniformRand() > fEfficiency) {
      continue;  // Hit not detected
    }
    
    // Get drift time from hit
    G4double driftTime = hit->GetTimeDrift();
    
    // Apply time resolution smearing (Gaussian)
    G4double smearedTime = driftTime + G4RandGauss::shoot(0., fTimeResolution);
    
    // Ensure non-negative time
    if (smearedTime < 0.) smearedTime = 0.;
    
    // Convert to TDC counts
    G4int tdc = static_cast<G4int>(smearedTime / fTDCResolution);
    
    // Create digi
    DriftCellDigi* digi = new DriftCellDigi();
    digi->SetEventID(hit->GetEventID());
    digi->SetCellID(hit->GetCellID());
    digi->SetTDC(tdc);
    digi->SetGlobalPos(hit->GetCellCenterPos());  // Use cell center for visualization
    digi->SetTrackID(hit->GetTrackID()); // Link to the particle that caused the hit
    
    // Add digi to collection
    fDigiCollection->insert(digi);
  }

  // Store the digi collection
  StoreDigiCollection(fDigiCollection);
  
  if (fDigiCollection->entries() > 0) {
    LogDebug("DriftCellDigitizer") << "Created " << fDigiCollection->entries() 
           << " digis" << G4endl;
  }
}

}
