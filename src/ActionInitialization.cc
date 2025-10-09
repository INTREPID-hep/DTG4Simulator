#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"

namespace DTSim
{

void ActionInitialization::BuildForMaster() const
{
}

void ActionInitialization::Build() const
{
  DTSim::PrimaryGeneratorAction* generator = new DTSim::PrimaryGeneratorAction;
  SetUserAction(generator);
}

}
