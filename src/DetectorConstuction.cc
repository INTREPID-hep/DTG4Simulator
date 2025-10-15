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
#include "G4tgrMessenger.hh"
#include "G4LogicalVolumeStore.hh"


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
    // Read materials from text file
    G4String materialsFile = "geometry/materials.tg";
    fMessenger->SetVerboseLevel(1);
    
    G4tgbVolumeMgr* volMgr = G4tgbVolumeMgr::GetInstance();
    volMgr->AddTextFile(materialsFile);
    
    // Read geometry from text file
    G4String geometryFile = "geometry/dt_geometry.txt";
    volMgr->AddTextFile(geometryFile);
    
    // Construct the geometry
    const G4VPhysicalVolume* worldPhys = volMgr->ReadAndConstructDetector();
    
    // Set visualization attributes
    SetVisualizationAttributes();
    
    return const_cast<G4VPhysicalVolume*>(worldPhys);
}


void DetectorConstruction::SetVisualizationAttributes()
{
    // G4Colour  white   ()              ;  // white
    // G4Colour  white   (1.0, 1.0, 1.0) ;  // white
    // G4Colour  gray    (0.5, 0.5, 0.5) ;  // gray
    // G4Colour  black   (0.0, 0.0, 0.0) ;  // black
    // G4Colour  red     (1.0, 0.0, 0.0) ;  // red
    // G4Colour  green   (0.0, 1.0, 0.0) ;  // green
    // G4Colour  blue    (0.0, 0.0, 1.0) ;  // blue
    // G4Colour  cyan    (0.0, 1.0, 1.0) ;  // cyan
    // G4Colour  magenta (1.0, 0.0, 1.0) ;  // magenta
    // G4Colour  yellow  (1.0, 1.0, 0.0) ;  // yellow

    // Get logical volume store
    auto logVolStore = G4LogicalVolumeStore::GetInstance();
    
    // World invisible
    G4LogicalVolume* worldLV = logVolStore->GetVolume("world");
    if (worldLV) {
        worldLV->SetVisAttributes(G4Colour::Red());
    }
    
    // DT Frame - semi-transparent blue
    G4LogicalVolume* frameLV = logVolStore->GetVolume("DTFrame");
    if (frameLV) {
        auto frameVis = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.5));
        frameVis->SetForceSolid(false);
        frameLV->SetVisAttributes(frameVis);
    }
    
    // Drift cells - yellow
    G4LogicalVolume* cellLV = logVolStore->GetVolume("DriftCell");
    if (cellLV) {
        auto cellVis = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0, 0.7));
        cellVis->SetForceSolid(true);
        cellLV->SetVisAttributes(cellVis);
    }

    // Axis markers - white and solid
    G4LogicalVolume* axisXLV = logVolStore->GetVolume("AxisX");
    if (axisXLV) {
        auto axisXVis = new G4VisAttributes(G4Colour::White());
        axisXVis->SetForceSolid(true);
        axisXLV->SetVisAttributes(axisXVis);
    }
    
    G4LogicalVolume* axisYLV = logVolStore->GetVolume("AxisY");
    if (axisYLV) {
        auto axisYVis = new G4VisAttributes(G4Colour::White());
        axisYVis->SetForceSolid(true);
        axisYLV->SetVisAttributes(axisYVis);
    }
    
    G4LogicalVolume* axisZLV = logVolStore->GetVolume("AxisZ");
    if (axisZLV) {
        auto axisZVis = new G4VisAttributes(G4Colour::White());
        axisZVis->SetForceSolid(true);
        axisZLV->SetVisAttributes(axisZVis);
    }
    
}

}
