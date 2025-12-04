#include "DetectorConstruction.hh"
#include "CommandLineParser.hh"
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
#include "DTSimConstants.hh"
#include "G4FieldBuilder.hh"
#include "G4UniformMagField.hh"

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
    fDriftCellsLogicals.clear();
    fYokeLogicals.clear();

    for (auto* logVol : *logVolStore) {
        G4String name = logVol->GetName();
        if (G4StrUtil::contains(name, "DriftCell")) {
            fDriftCellsLogicals.push_back(logVol);
        }
        else if (G4StrUtil::contains(name, "Yoke")) {
            fYokeLogicals.push_back(logVol);
        }
    }
    
    return const_cast<G4VPhysicalVolume*>(worldPhys);
}

void DetectorConstruction::ConstructSDandField()
{
    // ========== Setup Sensitive Detectors ==========
    if (!fDriftCellsLogicals.empty()) {
        DTSim::DriftCellSD* driftCellSD = new DTSim::DriftCellSD("/DriftCellSD");
        auto sdManager = G4SDManager::GetSDMpointer();
        sdManager->AddNewDetector(driftCellSD);
        
        
        for (auto* logVol : fDriftCellsLogicals) {
            logVol->SetSensitiveDetector(driftCellSD);
        }
    }

    // ========== Setup Magnetic Fields ==========
    // Check if magnetic field option was enabled via command line
    auto* parser = DTSim::CommandLineParser::Instance();
    G4bool enableMagneticField = parser->HasFlag("-B");
    
    if (!enableMagneticField) {
        G4cout << "Magnetic field disabled (use -B flag to enable)" << G4endl;
        return;
    }
    
    G4cout << "Building magnetic field configuration..." << G4endl;
    
    //    This automatically creates UI commands under /field/
    auto fieldBuilder = G4FieldBuilder::Instance();
    
    // world magnetic field
    // DTSim::MagneticField* worldMagField = new DTSim::MagneticField(); currently unused
    G4MagneticField* worldMagField = new G4UniformMagField(
        G4ThreeVector(0., 0., DTSim::kOutMagneticField)
    );
    fieldBuilder->SetGlobalField(worldMagField);
    
    if (!fYokeLogicals.empty()) {
        // Create uniform field for yoke
        G4MagneticField* yokeMagField = new G4UniformMagField(
            G4ThreeVector(0., 0., DTSim::kYokeMagneticField)
        );
        // This creates UI commands under /field/Yokes_...
        fieldBuilder->SetLocalField(yokeMagField, fYokeLogicals[0]);
        
        // Share the same FieldManager with all other yoke volumes
        G4FieldManager* yokeFieldMgr = fYokeLogicals[0]->GetFieldManager();
        for (size_t i = 1; i < fYokeLogicals.size(); ++i) {
            fYokeLogicals[i]->SetFieldManager(yokeFieldMgr, false);
        }
    }
    
    fieldBuilder->ConstructFieldSetup();
}

} // namespace DTSim
