#ifndef DTSimDriftCellSD_hh
#define DTSimDriftCellSD_hh 1

#include "G4VSensitiveDetector.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;
class G4GenericMessenger;

#include "DTSimTypes.hh"
#include "DriftCellHit.hh"

namespace DTSim
{

class RunAction;

class DriftCellSD : public G4VSensitiveDetector
{
  public:
      DriftCellSD(const G4String& name);
      ~DriftCellSD() override;

      void Initialize(G4HCofThisEvent* hce) override;
      void EndOfEvent(G4HCofThisEvent* hce) override;

      G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;

  private:
      void DefineCommands();
      
      CellID DecodeCellID(const G4String& volumeName, G4int copyNo) const;
      bool PassHitCriteria(const G4Step* step) const;
      G4bool ApplyElectrostaticConfinement(G4Step* step);
      G4bool EmulateWallCrossing(G4Step* step);
      G4bool ApplyWireCut(G4Step* step);

      const RunAction* fRunAction;

      DriftCellHitsCollection* fHitsCollection = nullptr;
      G4GenericMessenger* fMessenger;
      
      // Configurable physics parameters
      G4double fDriftVelocity;
      G4double fMinEnergyDeposit;
      G4double fCellBarrierEnergy;
      G4double fWallEnergyLoss;
      G4double fWireCutRadius;
      
      // Physics model enable/disable flags
      G4bool fEnableElectrostaticConfinement;
      G4bool fEnableWallCrossing;
      G4bool fEnableWireCut;
};

}

#endif