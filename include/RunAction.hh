#ifndef RunAction_hh
#define RunAction_hh

#include "G4UserRunAction.hh"
#include "G4Run.hh"

namespace DTSim
{

class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    ~RunAction() override;

    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;
};

}

#endif