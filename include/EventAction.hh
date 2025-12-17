#ifndef DTSimEventAction_h
#define DTSimEventAction_h 1

#include "G4UserEventAction.hh"
#include "G4AnalysisManager.hh"

#include <map>

class G4Event;

namespace DTSim
{

class RunAction;

class EventAction : public G4UserEventAction
{
  public:
    EventAction(RunAction* runAction);
    ~EventAction() = default;

    void BeginOfEventAction(const G4Event* event) override;
    void EndOfEventAction(const G4Event* event) override;

    // Methods to accumulate radiation info
    void AddRadiatedEnergy(G4int trackID, G4double energy) { fRadiatedEnergyMap[trackID] += energy; }
    void AddSecondaryCount(G4int trackID, G4int count) { fSecondaryCountMap[trackID] += count; }

  private:
    RunAction* fRunAction;
    
    // Radiation counters mapped by TrackID
    std::map<G4int, G4double> fRadiatedEnergyMap;
    std::map<G4int, G4int> fSecondaryCountMap;
    
    // Helper functions for filling ntuples
    void clearVectors();
    void fillGenBranches(const G4Event* event, G4AnalysisManager* analysisManager);
    void fillHitBranches(const G4Event* event, G4AnalysisManager* analysisManager);
    void fillDigiBranches(const G4Event* event, G4AnalysisManager* analysisManager);
    void fillSegmentBranches(const G4Event* event, G4AnalysisManager* analysisManager);
};

}

#endif
