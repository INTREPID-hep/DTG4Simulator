#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "DTSimLogger.hh"

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
    LogDebug("PrimaryGeneratorAction") << "Randomizing primary particle type" << G4endl;
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

  // Set direction using spherical coordinates (theta, phi)
  auto theta = fTheta + (G4UniformRand()-0.5)*fSigmaTheta;
  auto phi = fPhi + (G4UniformRand()-0.5)*fSigmaPhi;
  auto sinTheta = std::sin(theta);
  auto cosTheta = std::cos(theta);
  auto sinPhi = std::sin(phi);
  auto cosPhi = std::cos(phi);
  fParticleGun->SetParticleMomentumDirection(
                  G4ThreeVector(sinTheta*cosPhi, sinTheta*sinPhi, cosTheta));

  // Set position with gaussian spread
  auto x = fPosition.x() + G4RandGauss::shoot(0., fSigmaPosition.x());
  auto y = fPosition.y() + G4RandGauss::shoot(0., fSigmaPosition.y());
  auto z = fPosition.z() + G4RandGauss::shoot(0., fSigmaPosition.z());
  fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));

  fParticleGun->GeneratePrimaryVertex(event);
   // Log selected particle and kinematics
  LogDebug("PrimaryGeneratorAction") << "Generating primary: "
      << particle->GetParticleName() << ", momentum = " << fMomentum/GeV << " GeV, "
      << "theta = " << theta/deg << " deg, phi = " << phi/deg << " deg, position = ("
      << x << ", " << y << ", " << z << ") m" << G4endl;
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

  // theta command (polar angle)
  auto& thetaCmd
    = fMessenger->DeclarePropertyWithUnit("theta", "deg", fTheta,
        "Polar angle (0=+Z axis, 90=XY plane, 180=-Z axis).");
  thetaCmd.SetParameterName("theta", true);
  thetaCmd.SetRange("theta>=0. && theta<=180.");
  thetaCmd.SetDefaultValue("90.");

  // phi command (azimuthal angle)
  auto& phiCmd
    = fMessenger->DeclarePropertyWithUnit("phi", "deg", fPhi,
        "Azimuthal angle (rotation around Z axis, 0=+X, 90=+Y).");
  phiCmd.SetParameterName("phi", true);
  phiCmd.SetRange("phi>=-180. && phi<=360.");
  phiCmd.SetDefaultValue("0.");

  // sigmaTheta command
  auto& sigmaThetaCmd
    = fMessenger->DeclarePropertyWithUnit("sigmaTheta", "deg", fSigmaTheta,
        "Sigma spread of polar angle.");
  sigmaThetaCmd.SetParameterName("st", true);
  sigmaThetaCmd.SetRange("st>=0.");
  sigmaThetaCmd.SetDefaultValue("2.");

  // sigmaPhi command
  auto& sigmaPhiCmd
    = fMessenger->DeclarePropertyWithUnit("sigmaPhi", "deg", fSigmaPhi,
        "Sigma spread of azimuthal angle.");
  sigmaPhiCmd.SetParameterName("sp", true);
  sigmaPhiCmd.SetRange("sp>=0.");
  sigmaPhiCmd.SetDefaultValue("2.");

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

  // position command
  auto& positionCmd
    = fMessenger->DeclarePropertyWithUnit("position", "m", fPosition,
        "Mean position of primaries (x, y, z).");
  positionCmd.SetParameterName("x", "y", "z", true);
  positionCmd.SetDefaultValue("0. 0. -2.5");

  // sigmaPosition command
  auto& sigmaPositionCmd
    = fMessenger->DeclarePropertyWithUnit("sigmaPosition", "cm", fSigmaPosition,
        "Sigma spread of primary position (x, y, z).");
  sigmaPositionCmd.SetParameterName("sx", "sy", "sz", true);
  sigmaPositionCmd.SetDefaultValue("0. 0. 0.");
}

}
