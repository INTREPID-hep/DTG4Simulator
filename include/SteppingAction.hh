#ifndef SteppingAction_hh
#define SteppingAction_hh 1

#include "G4UserSteppingAction.hh"

namespace DTSim
{

class EventAction;
class RunAction;

class SteppingAction : public G4UserSteppingAction
{
  public:
    SteppingAction();
    ~SteppingAction() = default;

    void UserSteppingAction(const G4Step*) override;

  private:
    const EventAction* fEventAction;
    const RunAction* fRunAction;
};

}

#endif
