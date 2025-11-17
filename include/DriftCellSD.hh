#ifndef DTSimDriftCellSD_hh
#define DTSimDriftCellSD_hh 1

#include "G4VSensitiveDetector.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "DTSimTypes.hh"
#include "DTSimConstants.hh"
#include "DriftCellHit.hh"
#include <vector>

namespace DTSim
{

class DriftCellSD : public G4VSensitiveDetector
{
  public:
      DriftCellSD(const G4String& name, const G4String& hitsCollectionName);
      ~DriftCellSD() override = default;

      void Initialize(G4HCofThisEvent* hce) override;
      void EndOfEvent(G4HCofThisEvent* hce) override;

      G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;

  private:
      G4double GetTimeWithDrift(G4double time, G4ThreeVector distance) const;
      
      CellID DecodeCellID(const G4String& volumeName, G4int copyNo) const;

      DriftCellHitsCollection* fHitsCollection = nullptr;
};

}

#endif