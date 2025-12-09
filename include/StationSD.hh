#ifndef DTSimStationSD_hh
#define DTSimStationSD_hh 1

#include "G4VSensitiveDetector.hh"
#include <map>

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;
class G4StepPoint;

#include "DTSimTypes.hh"
#include "DTSegment.hh"

namespace DTSim
{

class StationSD : public G4VSensitiveDetector
{
  public:
      StationSD(const G4String& name);
      ~StationSD() override = default;

      void Initialize(G4HCofThisEvent* hce) override;
      void EndOfEvent(G4HCofThisEvent* hce) override;

      G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;

  private:
      StationID DecodeStationID(const G4String& volumeName) const;

      DTSegmentCollection* fSegmentCollection = nullptr;
      std::map<G4int, G4ThreeVector> fEntryStepMap;  // TrackID -> entry position
};

}

#endif
