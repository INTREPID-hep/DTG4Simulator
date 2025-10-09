//
// Auto-generated from mplDTs Python package
// DO NOT EDIT - Regenerate using generate_geometry_builder.py
//

#ifndef DTGeometryBuilder_h
#define DTGeometryBuilder_h 1

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include <map>
#include <vector>

namespace DTSim {

/// Structure to hold cell geometry data
struct CellGeometry {
    int cellNumber;
    int layerNumber;
    double x, y, z;
    double width, height, length;
};

/// Structure to hold superlayer geometry data
struct SuperLayerGeometry {
    int slNumber;
    double centerX, centerY, centerZ;
    double width, height, length;
    double rotation;                   // degrees
    std::vector<CellGeometry> cells; 
};

/// Structure to hold complete station geometry data
struct StationGeometry {
    int wheel, sector, station;
    double centerX, centerY, centerZ;
    double width, height, length;
    std::vector<SuperLayerGeometry> superlayers;
};

class DTGeometryBuilder {
public:
    DTGeometryBuilder();
    ~DTGeometryBuilder();
    
    /// Build a complete DT station
    /// @param wheel: Wheel number (-2 to +2)
    /// @param sector: Sector number (1-12 for MB1/MB2/MB3, 1-14 for MB4)
    /// @param station: Station number (1-4)
    /// @param motherVolume: Mother logical volume to place station in
    /// @param checkOverlaps: Whether to check for volume overlaps
    /// @return Station logical volume
    G4LogicalVolume* BuildStation(int wheel, int sector, int station,
                                   G4LogicalVolume* motherVolume,
                                   bool checkOverlaps = true);
    
    /// Get station geometry data
    const StationGeometry* GetStationGeometry(int wheel, int sector, int station) const;
    
    /// Get all available station configurations
    std::vector<std::tuple<int,int,int>> GetAvailableStations() const;

private:
    /// Build a single superlayer
    G4LogicalVolume* BuildSuperLayer(const SuperLayerGeometry& slGeom,
                                      G4LogicalVolume* motherVolume,
                                      bool checkOverlaps);

    /// Build cells using parameterisation (efficient - recommended)
    void BuildCellsParameterised(const SuperLayerGeometry& slGeom,
                                 G4LogicalVolume* slLogical,
                                 bool checkOverlaps);
};

} // namespace DTSim

#endif
