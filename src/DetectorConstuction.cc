#include "DetectorConstruction.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4tgbVolumeMgr.hh"
#include "G4tgrVolume.hh"
#include "G4tgrMessenger.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4String.hh"
#include "G4SDManager.hh"
#include "DriftCellSD.hh"
#include "MagneticField.hh"
#include "G4FieldBuilder.hh"
#include "DTSimConstants.hh"

namespace DTSim
{

DetectorConstruction::DetectorConstruction()
{
    fMessenger = new G4tgrMessenger;
}

DetectorConstruction::~DetectorConstruction()
{
    delete fMessenger;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    fMessenger->SetVerboseLevel(1);
    G4tgbVolumeMgr* volMgr = G4tgbVolumeMgr::GetInstance();

    // Read geometry from text file
    G4String geometryFile = "geometry/geometry_concentrator.tg";
    volMgr->AddTextFile(geometryFile);
    
    // Construct the geometry
    const G4VPhysicalVolume* worldPhys = volMgr->ReadAndConstructDetector();

    // Find all DriftCell logical volumes: Necessary for setting sensitive volumes
    auto logVolStore = G4LogicalVolumeStore::GetInstance();
    driftCellsLogicals.clear();
    
    for (auto* logVol : *logVolStore) {
        G4String name = logVol->GetName();
        if (G4StrUtil::contains(name, "DriftCell")) {
            driftCellsLogicals.push_back(logVol);
        }
        else if (G4StrUtil::contains(name, "Yoke")) {
            yokeLogicals.push_back(logVol);
        }
    }
    
    return const_cast<G4VPhysicalVolume*>(worldPhys);
}

void DetectorConstruction::ConstructSDandField()
{
    DTSim::DriftCellSD* driftCellSD = new DTSim::DriftCellSD("/DriftCellSD", "DriftCellHitsCollection");
    auto sdManager = G4SDManager::GetSDMpointer();
    sdManager->AddNewDetector(driftCellSD);
    
    for (auto* logVol : driftCellsLogicals) {
        logVol->SetSensitiveDetector(driftCellSD);
    }

    // Magnetic field
    G4FieldBuilder* fieldBuilder = G4FieldBuilder::Instance();

    // Define global magnetic field
    DTSim::MagneticField* worldMagField = new DTSim::MagneticField();
    fieldBuilder->SetGlobalField(worldMagField);
    // Define yoke magnetic field
    G4MagneticField* yokeMagField = new G4UniformMagField(G4ThreeVector(0., 0., DTSim::kYokeMagneticField));;
    for (auto* logVol : yokeLogicals) {
        fieldBuilder->CreateFieldParameters(logVol->GetName());
        fieldBuilder->SetLocalField(yokeMagField, logVol);
    }

    // Construct all Geant4 field objects
    fieldBuilder->ConstructFieldSetup();

}

} // namespace DTSim
