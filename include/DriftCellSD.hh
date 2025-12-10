#ifndef DTSimDriftCellSD_hh
#define DTSimDriftCellSD_hh 1

#include "G4VSensitiveDetector.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

#include "DTSimTypes.hh"
#include "DriftCellHit.hh"

namespace DTSim
{

class DriftCellSD : public G4VSensitiveDetector
{
  public:
      DriftCellSD(const G4String& name);
      ~DriftCellSD() override = default;

      void Initialize(G4HCofThisEvent* hce) override;
      void EndOfEvent(G4HCofThisEvent* hce) override;

      G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;

  private:
      CellID DecodeCellID(const G4String& volumeName, G4int copyNo) const;
      bool PassHitCriteria(const G4Step* step) const;
      void ApplyElectrostaticConfinement(G4Step* step);
      void EmulateWallCrossing(G4Step* step);

      DriftCellHitsCollection* fHitsCollection = nullptr;
};

}

#endif