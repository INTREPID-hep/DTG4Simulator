#ifndef DTSimDriftCellSD_hh
#define DTSimDriftCellSD_hh 1

#include "G4VSensitiveDetector.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include <vector>

namespace DTSim
{

class DriftCellSD : public G4VSensitiveDetector
{
  public:
      DriftCellSD(const G4String& name);
      ~DriftCellSD() override;

      void Initialize(G4HCofThisEvent* hce) override;
      void EndOfEvent(G4HCofThisEvent* hce) override;

      G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;

  private:
      G4double GetTimeWithDrift(G4double time, G4ThreeVector distance) const;
      
      // Helper to extract and cells-per-layer info from volume name
      void ExtractLayerStructure(const G4String& volumeName) const;
      
      // Decode both layer and cell ID
      void DecodeCopyNumber(G4int copyNo, G4int& layerID, G4int& cellID) const;

      // Cache for layer structure (mutable to allow caching in const methods)
      mutable G4String fVolumeName;
      mutable std::vector<G4int> fCellsPerLayer;

      // constants
      const G4double kDriftVelocity = 54.0*micrometer/nanosecond;
};

}

#endif