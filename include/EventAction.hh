#ifndef DTSimEventAction_h
#define DTSimEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

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

  private:
    RunAction* fRunAction;
    
    // Helper functions for filling ntuples
    void clearVectors();
    void FillHitsNtuple(const G4Event* event);
};

}

#endif
