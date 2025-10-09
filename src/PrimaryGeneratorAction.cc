#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

namespace DTSim
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  G4int nofParticles = 1;
  fParticleGun  = new G4ParticleGun(nofParticles);

  auto particleTable = G4ParticleTable::GetParticleTable();
  fPositron = particleTable->FindParticle("e+");
  fMuon = particleTable->FindParticle("mu+");
  fPion = particleTable->FindParticle("pi+");
  fKaon = particleTable->FindParticle("kaon+");
  fProton = particleTable->FindParticle("proton");

  // default particle kinematics
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,-5.*m)); // <---change this latter
  fParticleGun->SetParticleDefinition(fMuon);

  // define commands for this class
  DefineCommands();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fMessenger;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  G4ParticleDefinition* particle;
  if (fRandomizePrimary) { 
    // if randomizing, select a particle type at random
    auto i = (int)(5. * G4UniformRand());
    switch(i) {
      case 0:
          particle = fPositron;
          break;
      case 1:
          particle = fMuon;
          break;
      case 2:
          particle = fPion;
          break;
      case 3:
          particle = fKaon;
          break;
      default:
          particle = fProton;
          break;
    }
    fParticleGun->SetParticleDefinition(particle);
  }
  else {
    // if not randomizing, use the particle selected with /gun/particle command or muon by default
    particle = fParticleGun->GetParticleDefinition();
  }

  auto pp = fMomentum + (G4UniformRand()-0.5)*fSigmaMomentum;
  auto mass = particle->GetPDGMass();
  auto ekin = std::sqrt(pp*pp+mass*mass)-mass;
  fParticleGun->SetParticleEnergy(ekin);

  auto angle = (G4UniformRand()-0.5)*fSigmaAngle;
  fParticleGun->SetParticleMomentumDirection(
                  G4ThreeVector(std::sin(angle),0.,std::cos(angle)));

  fParticleGun->GeneratePrimaryVertex(event);
}

void PrimaryGeneratorAction::DefineCommands()
{
  // Define /DTSim/generator command directory using generic messenger class
  fMessenger
    = new G4GenericMessenger(this,
                             "/DTSim/generator/",
                             "Primary generator control");

  // momentum command
  auto& momentumCmd
    = fMessenger->DeclarePropertyWithUnit("momentum", "GeV", fMomentum,
        "Mean momentum of primaries.");
  momentumCmd.SetParameterName("p", true);
  momentumCmd.SetRange("p>=0.");
  momentumCmd.SetDefaultValue("1000.");

  // sigmaMomentum command
  auto& sigmaMomentumCmd
    = fMessenger->DeclarePropertyWithUnit("sigmaMomentum",
        "MeV", fSigmaMomentum, "Sigma momentum of primaries.");
  sigmaMomentumCmd.SetParameterName("sp", true);
  sigmaMomentumCmd.SetRange("sp>=0.");
  sigmaMomentumCmd.SetDefaultValue("50.");

  // sigmaAngle command
  auto& sigmaAngleCmd
    = fMessenger->DeclarePropertyWithUnit("sigmaAngle", "deg", fSigmaAngle,
        "Sigma angle divergence of primaries.");
  sigmaAngleCmd.SetParameterName("t", true);
  sigmaAngleCmd.SetRange("t>=0.");
  sigmaAngleCmd.SetDefaultValue("2.");

  // randomizePrimary command
  auto& randomCmd
    = fMessenger->DeclareProperty("randomizePrimary", fRandomizePrimary);
  G4String guidance
    = "Boolean flag for randomizing primary particle types.\n";
  guidance
    += "In case this flag is false, you can select the primary particle\n";
  guidance += "  with /gun/particle command.";
  randomCmd.SetGuidance(guidance);
  randomCmd.SetParameterName("flg", true);
  randomCmd.SetDefaultValue("true");
}

}
