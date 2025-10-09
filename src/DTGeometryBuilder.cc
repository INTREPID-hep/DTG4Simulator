//
// DT Geometry Builder Implementation
// Uses data from auto-generated DTGeometryData.hh
//

#include "DTGeometryBuilder.hh"
#include "DTGeometryData.hh"  // Auto-generated geometry data
#include "DTCellParameterisation.hh"  // Cell parameterisation
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include <sstream>

namespace DTSim {

DTGeometryBuilder::DTGeometryBuilder() {
}

DTGeometryBuilder::~DTGeometryBuilder() {}

G4LogicalVolume* DTGeometryBuilder::BuildStation(int wheel, int sector, int station,
                                                  G4LogicalVolume* motherVolume,
                                                  bool checkOverlaps) {
    // Get geometry data
    const StationGeometry* stGeom = GetStationGeometry(wheel, sector, station);
    if (!stGeom) {
        G4cerr << "ERROR: No geometry data for Wh=" << wheel 
               << " Sec=" << sector << " St=" << station << G4endl;
        return nullptr;
    }
    
    // Get gas mixture material
    auto gasMixture = G4Material::GetMaterial("GasMixture");
    if (!gasMixture) {
        G4cerr << "ERROR: GasMixture material not found!" << G4endl;
        return nullptr;
    }
    
    // Create station volume
    std::ostringstream stName;
    stName << "DTStation_Wh" << wheel << "_Sec" << sector << "_St" << station;
    
    auto stationSolid = new G4Box(stName.str() + "_solid",
        stGeom->width * cm / 2.0,
        stGeom->height * cm / 2.0,
        stGeom->length * cm / 2.0);
    
    auto stationLogical = new G4LogicalVolume(stationSolid, gasMixture,
        stName.str() + "_logical");
    
    // Place station in mother volume
    // Determine face orientation based on wheel and sector
    // int face_orientation;
    // if (wheel < 0) {
    //     face_orientation = -1;
    // } else if (wheel > 0) {
    //     face_orientation = 1;
    // } else if (sector == 1 || sector == 4 || sector == 5 || sector == 8 || sector == 9 || sector == 12 || sector == 13) {
    //     face_orientation = -1;
    // } else {
    //     face_orientation = 1;
    // }
    
    // Create rotation matrix (needs to be a pointer with 'new')
    // G4RotationMatrix* rotm = new G4RotationMatrix();
    // rotm->rotateX(-1 * face_orientation * 90 * deg);

    new G4PVPlacement(nullptr,
        G4ThreeVector(stGeom->centerX * cm, stGeom->centerY * cm, stGeom->centerZ * cm),
        stationLogical, stName.str() + "_physical",
        motherVolume, false, 0, checkOverlaps);
    
    // Build superlayers
    for (const auto& slGeom : stGeom->superlayers) {
        BuildSuperLayer(slGeom, stationLogical, checkOverlaps);
        // break; // For now, build only the first superlayer for testing
    }
    
    G4cout << "Built DT Station: Wh=" << wheel << " Sec=" << sector 
           << " St=" << station << " with " << stGeom->superlayers.size() 
           << " superlayers" << G4endl;
    
    return stationLogical;
}

G4LogicalVolume* DTGeometryBuilder::BuildSuperLayer(const SuperLayerGeometry& slGeom,
                                                     G4LogicalVolume* motherVolume,
                                                     bool checkOverlaps) {
    auto gasMixture = G4Material::GetMaterial("GasMixture");
    
    // Create superlayer volume
    std::ostringstream slName;
    slName << "SuperLayer" << slGeom.slNumber;
    
    auto slSolid = new G4Box(slName.str() + "_solid",
        slGeom.width * cm / 2.0,
        slGeom.height * cm / 2.0,
        slGeom.length * cm / 2.0);
    
    auto slLogical = new G4LogicalVolume(slSolid, gasMixture,
        slName.str() + "_logical");
    
    // Apply rotation if needed (SL2 is rotated 90 degrees)
    G4RotationMatrix* rotation = nullptr;
    // if (slGeom.rotation != 0.0) {
    //     rotation = new G4RotationMatrix();
    //     rotation->rotateZ(slGeom.rotation * deg);
    // }
    
    // Place superlayer
    new G4PVPlacement(rotation,
        G4ThreeVector(slGeom.centerX * cm, slGeom.centerZ * cm, slGeom.centerY * cm),
        slLogical, slName.str() + "_physical",
        motherVolume, false, slGeom.slNumber, checkOverlaps);
    
    // Build cells using parameterization (efficient!)
    BuildCellsParameterised(slGeom, slLogical, checkOverlaps);
    
    // Set visualization
    G4VisAttributes blue(G4Colour::Blue());
    slLogical->SetVisAttributes(blue);
    
    return slLogical;
}

void DTGeometryBuilder::BuildCellsParameterised(const SuperLayerGeometry& slGeom,
                                                G4LogicalVolume* slLogical,
                                                bool checkOverlaps) {
    auto gasMixture = G4Material::GetMaterial("GasMixture");
    
    if (slGeom.cells.empty()) {
        G4cerr << "ERROR: No cells found in SuperLayer " << slGeom.slNumber << G4endl;
        return;
    }
    
    // Get cell dimensions from first cell (all cells have same dimensions)
    const auto& firstCell = slGeom.cells[0];
    
    // Create ONE logical volume for all cells (memory efficient!)
    std::ostringstream cellName;
    cellName << "Cell_SL" << slGeom.slNumber;
    
    auto cellSolid = new G4Box(cellName.str() + "_solid",
        firstCell.width * cm / 2.0,
        firstCell.height * cm / 2.0,
        firstCell.length * cm / 2.0);
    
    auto cellLogical = new G4LogicalVolume(cellSolid, gasMixture,
        cellName.str() + "_logical");
    
    // Create parameterisation with mplDTs geometry
    auto cellParam = new DTCellParameterisation(slGeom);
    G4int totalCells = cellParam->GetNumberOfCells();
    
    // Place all cells using ONE logical volume (replicated with different positions)
    new G4PVParameterised(cellName.str() + "_physical",
        cellLogical,      // The ONE logical volume
        slLogical,        // Mother volume
        kUndefined,       // Axis (not used for custom parameterisation)
        totalCells,       // Number of copies
        cellParam,        // Parameterisation object
        checkOverlaps);
    
    G4cout << "Built " << totalCells << " cells in SuperLayer " << slGeom.slNumber 
           << " using parameterisation (1 logical volume reused)" << G4endl;
}

const StationGeometry* DTGeometryBuilder::GetStationGeometry(int wheel, int sector, int station) const {
    return DTGeometryData::Instance()->GetStationGeometry(wheel, sector, station);
}

std::vector<std::tuple<int,int,int>> DTGeometryBuilder::GetAvailableStations() const {
    std::vector<std::tuple<int,int,int>> stations;
    const auto& allStations = DTGeometryData::Instance()->GetAllStations();
    for (const auto& entry : allStations) {
        stations.push_back(entry.first);
    }
    return stations;
}

} // namespace DTSim
