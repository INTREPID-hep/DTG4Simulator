#ifndef DTSimPrimaryGeneratorAction_h
#define DTSimPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4Types.hh"

class G4ParticleGun;
class G4GenericMessenger;
class G4Event;
class G4ParticleDefinition;

#include "DTSimConstants.hh"

namespace DTSim
{

/// Primary generator
///
/// A single particle is generated.
/// User can select
/// - the initial momentum and angle
/// - the momentum and angle spreads
/// - random selection of a particle type from proton, kaon+, pi+, muon+, e+


class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event*) override;

    void SetMomentum(G4double val) { fMomentum = val; }
    G4double GetMomentum() const { return fMomentum; }

    void SetSigmaMomentum(G4double val) { fSigmaMomentum = val; }
    G4double GetSigmaMomentum() const { return fSigmaMomentum; }

    void SetSigmaAngle(G4double val) { fSigmaAngle = val; }
    G4double GetSigmaAngle() const { return fSigmaAngle; }

    void SetRandomize(G4bool val) { fRandomizePrimary = val; }
    G4bool GetRandomize() const { return fRandomizePrimary; }

  private:
    void DefineCommands();

    G4ParticleGun* fParticleGun = nullptr;
    G4GenericMessenger* fMessenger = nullptr;
    G4ParticleDefinition* fPositron = nullptr;
    G4ParticleDefinition* fMuon = nullptr;
    G4ParticleDefinition* fPion = nullptr;
    G4ParticleDefinition* fKaon = nullptr;
    G4ParticleDefinition* fProton = nullptr;
    G4double fMomentum = DTSim::kDefaultMomentum;
    G4double fSigmaMomentum = DTSim::kDefaultSigmaMomentum;
    G4double fSigmaAngle = DTSim::kDefaultSigmaAngle;
    G4bool fRandomizePrimary = true;
};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
