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
#include "G4LogicalVolumeStore.hh"
#include "G4String.hh"
#include "G4SDManager.hh"
#include "G4FieldBuilder.hh"
#include "G4UniformMagField.hh"
#include "G4GenericMessenger.hh"
#include "G4Region.hh"
#include "G4UserLimits.hh"

#include "DriftCellSD.hh"
#include "StationSD.hh"
#include "DTSimConstants.hh"
#include "DTSimLogger.hh"

namespace DTSim
{

DetectorConstruction::DetectorConstruction()
  : fDetMessenger(nullptr),
    fUseBField(true),
    fGlobalField(DTSim::kGlobalMagneticField),
    fYokeField(DTSim::kYokeMagneticField),
    fGeometryFileName(DTSim::kGeometryFileName),
    fEnableDriftSD(true),
    fEnableStationSD(true),
    // Station limits
    fStationStepMax(DBL_MAX),
    fStationTrakMax(DBL_MAX),
    fStationTimeMax(DBL_MAX),
    fStationEkinMin(0.),
    fStationRangMin(0.),
    // Yoke limits
    fYokeStepMax(DBL_MAX),
    fYokeTrakMax(DBL_MAX),
    fYokeTimeMax(DBL_MAX),
    fYokeEkinMin(0.),
    fYokeRangMin(0.)
{
    DefineCommands();
}

DetectorConstruction::~DetectorConstruction()
{
    delete fDetMessenger;
    delete fStationLimits;
    delete fYokeLimits;
}

void DetectorConstruction::DefineCommands()
{
    fDetMessenger = new G4GenericMessenger(this, "/DTSim/detector/", 
                                            "Detector construction control");
    
    // Magnetic field commands
    fDetMessenger->DeclareProperty("useBField", fUseBField,
                                   "Enable/disable magnetic field");
    
    fDetMessenger->DeclarePropertyWithUnit("BField/setGlobal", "tesla", 
                                           fGlobalField,
                                           "Set global magnetic field vector (x,y,z)");
    
    fDetMessenger->DeclarePropertyWithUnit("BField/setYoke", "tesla", 
                                           fYokeField,
                                           "Set yoke magnetic field vector (x,y,z)");
    
    // Geometry file command
    fDetMessenger->DeclareProperty("setGeometryFile", fGeometryFileName,
                                   "Set path to geometry text file");
    
    // Sensitive detector commands
    fDetMessenger->DeclareProperty("enableDriftSD", fEnableDriftSD,
                                   "Enable/disable Drift Cell sensitive detector");
    
    fDetMessenger->DeclareProperty("enableStationSD", fEnableStationSD,
                                   "Enable/disable Station sensitive detector");

    // Station User Limits
    fDetMessenger->DeclarePropertyWithUnit("setStationStepMax", "mm", 
                                         fStationStepMax, 
                                         "Set max step limit for StationRegion");

    fDetMessenger->DeclarePropertyWithUnit("setStationTrakMax", "mm", 
                                         fStationTrakMax, 
                                         "Set max track length for StationRegion");

    fDetMessenger->DeclarePropertyWithUnit("setStationTimeMax", "ns", 
                                         fStationTimeMax, 
                                         "Set max time for StationRegion");

    fDetMessenger->DeclarePropertyWithUnit("setStationEkinMin", "MeV", 
                                         fStationEkinMin, 
                                         "Set min kinetic energy for StationRegion");

    fDetMessenger->DeclarePropertyWithUnit("setStationRangMin", "mm", 
                                         fStationRangMin, 
                                         "Set min range for StationRegion");

    // Yoke User Limits
    fDetMessenger->DeclarePropertyWithUnit("setYokeStepMax", "mm", 
                                         fYokeStepMax, 
                                         "Set max step limit for YokeRegion");

    fDetMessenger->DeclarePropertyWithUnit("setYokeTrakMax", "mm", 
                                         fYokeTrakMax, 
                                         "Set max track length for YokeRegion");

    fDetMessenger->DeclarePropertyWithUnit("setYokeTimeMax", "ns", 
                                         fYokeTimeMax, 
                                         "Set max time for YokeRegion");

    fDetMessenger->DeclarePropertyWithUnit("setYokeEkinMin", "MeV", 
                                         fYokeEkinMin, 
                                         "Set min kinetic energy for YokeRegion");

    fDetMessenger->DeclarePropertyWithUnit("setYokeRangMin", "mm", 
                                         fYokeRangMin, 
                                         "Set min range for YokeRegion");
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    // NOTE: Configuration parameters (enableDriftSD, BField, etc.) should be set
    // BEFORE calling /run/initialize. Runtime geometry reinitialization is not supported.
    
    G4tgbVolumeMgr* volMgr = G4tgbVolumeMgr::GetInstance();

    // Read geometry from text file (now configurable)
    volMgr->AddTextFile(fGeometryFileName);
    LogInfo("DetectorConstruction") << "Loading geometry from: " << fGeometryFileName << G4endl;
    
    // Construct the geometry
    const G4VPhysicalVolume* worldPhys = volMgr->ReadAndConstructDetector();

    // Find all DriftCell logical volumes: Necessary for setting sensitive volumes
    auto logVolStore = G4LogicalVolumeStore::GetInstance();
    fDriftCellsLogicals.clear();
    fYokeLogicals.clear();
    fStationLogicals.clear();

    // Create a region for DriftCells to allow setting specific cuts
    G4Region* yokeRegion = new G4Region("YokeRegion");
    fYokeLimits = new G4UserLimits(fYokeStepMax,
                                   fYokeTrakMax,
                                   fYokeTimeMax,
                                   fYokeEkinMin,
                                   fYokeRangMin);
    yokeRegion->SetUserLimits(fYokeLimits);

    G4Region* stationRegion = new G4Region("StationRegion");
    fStationLimits = new G4UserLimits(fStationStepMax,
                                      fStationTrakMax,
                                      fStationTimeMax,
                                      fStationEkinMin,
                                      fStationRangMin);
    stationRegion->SetUserLimits(fStationLimits);

    for (auto* logVol : *logVolStore) {
        G4String name = logVol->GetName();
        if (G4StrUtil::contains(name, "DriftCell")) {
            fDriftCellsLogicals.push_back(logVol);
        }
        else if (G4StrUtil::contains(name, "Yoke")) {
            fYokeLogicals.push_back(logVol);
            yokeRegion->AddRootLogicalVolume(logVol);
        }
        else if (G4StrUtil::contains(name, "Station")) {
            fStationLogicals.push_back(logVol);
            stationRegion->AddRootLogicalVolume(logVol);
        }
    }
    
    return const_cast<G4VPhysicalVolume*>(worldPhys);
}

void DetectorConstruction::ConstructSDandField()
{
    // ========== Setup Sensitive Detectors ==========
    auto sdManager = G4SDManager::GetSDMpointer();
    
    // DriftCell sensitive detector (conditionally enabled)
    if (!fEnableDriftSD) {
        LogInfo("DetectorConstruction") << "DriftCell SD disabled by configuration" << G4endl;
    } else if (!fDriftCellsLogicals.empty()) {
        LogInfo("DetectorConstruction") << "DriftCell SD enabled: " << fDriftCellsLogicals.size() << " volumes" << G4endl;
        DTSim::DriftCellSD* driftCellSD = new DTSim::DriftCellSD("/DriftCellSD");
        sdManager->AddNewDetector(driftCellSD);
        
        for (auto* logVol : fDriftCellsLogicals) {
            logVol->SetSensitiveDetector(driftCellSD);
        }
    } else {
        LogWarn("DetectorConstruction") << "No DriftCell logical volumes found!" << G4endl;
    }

    // Station sensitive detector (conditionally enabled)
    if (!fEnableStationSD) {
        LogInfo("DetectorConstruction") << "Station SD disabled by configuration" << G4endl;
    } else if (!fStationLogicals.empty()) {
        LogInfo("DetectorConstruction") << "Station SD enabled: " << fStationLogicals.size() << " volumes" << G4endl;
        DTSim::StationSD* stationSD = new DTSim::StationSD("/StationSD");
        sdManager->AddNewDetector(stationSD);
        for (auto* logVol : fStationLogicals) {
            logVol->SetSensitiveDetector(stationSD);
        }
    } else {
        LogWarn("DetectorConstruction") << "No Station logical volumes found!" << G4endl;
    }

    // ========== Setup Magnetic Fields ==========
    if (!fUseBField) {
        LogInfo("DetectorConstruction") << "Magnetic field disabled by configuration" << G4endl;
        return;
    }
    LogDebug("DetectorConstruction") << "Building magnetic field configuration..." << G4endl;
    auto fieldBuilder = G4FieldBuilder::Instance(); // This automatically creates UI commands under /field/
    
    // world magnetic field (use configured value)
    G4MagneticField* worldMagField = new G4UniformMagField(fGlobalField);
    fieldBuilder->SetGlobalField(worldMagField);
    LogDebug("DetectorConstruction") << "  Global field: (" << fGlobalField.x()/tesla << ", " 
        << fGlobalField.y()/tesla << ", " << fGlobalField.z()/tesla << ") T" << G4endl;
    
    LogDebug("DetectorConstruction") << "  Yoke field: (" << fYokeField.x()/tesla << ", " 
        << fYokeField.y()/tesla << ", " << fYokeField.z()/tesla << ") T" << G4endl;
    if (!fYokeLogicals.empty()) {
        // Create uniform field for yoke (use configured value)
        G4MagneticField* yokeMagField = new G4UniformMagField(fYokeField);
        // This creates UI commands under /field/Yokes_...
        fieldBuilder->SetLocalField(yokeMagField, fYokeLogicals[0]);
        
        // Share the same FieldManager with all other yoke volumes
        G4FieldManager* yokeFieldMgr = fYokeLogicals[0]->GetFieldManager();
        for (size_t i = 1; i < fYokeLogicals.size(); ++i) {
            fYokeLogicals[i]->SetFieldManager(yokeFieldMgr, false);
        }
    } else {
        LogWarn("DetectorConstruction") << "No Yoke logical volumes found!" << G4endl;
    }
    
    fieldBuilder->ConstructFieldSetup();
}

} // namespace DTSim
