#include "SteppingAction.hh"
#include "EventAction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"

namespace DTSim
{

SteppingAction::SteppingAction(EventAction* eventAction)
 : G4UserSteppingAction(),
   fEventAction(eventAction)
{}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    G4Track* track = step->GetTrack();
    
    // We only care about primary particles (ParentID = 0)
    if (track->GetParentID() != 0) return;

    // Check for secondaries produced in this step
    const std::vector<const G4Track*>* secondaries = step->GetSecondaryInCurrentStep();

    if (!secondaries || secondaries->empty()) return;

    G4double radiatedEnergy = 0.0;
    G4int nElectrons = 0;

    for (const auto* secondary : *secondaries) {
        // Check if secondary is an electron or positron (PDG +/- 11)
        if (std::abs(secondary->GetDefinition()->GetPDGEncoding()) == 11) {
            radiatedEnergy += secondary->GetKineticEnergy();
            nElectrons++;
        }
    }

    if (nElectrons > 0) {
        fEventAction->AddRadiatedEnergy(track->GetTrackID(), radiatedEnergy);
        fEventAction->AddSecondaryCount(track->GetTrackID(), nElectrons);
    }
}

}
