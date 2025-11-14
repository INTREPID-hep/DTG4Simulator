#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"

namespace DTSim
{

void ActionInitialization::BuildForMaster() const
{
  DTSim::RunAction* runAction = new DTSim::RunAction;
  SetUserAction(runAction);

}

void ActionInitialization::Build() const
{
  DTSim::PrimaryGeneratorAction* generator = new DTSim::PrimaryGeneratorAction;
  SetUserAction(generator);

  DTSim::RunAction* runAction = new DTSim::RunAction;
  SetUserAction(runAction);
}

}
