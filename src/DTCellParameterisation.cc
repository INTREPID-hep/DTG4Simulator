//
// DT Cell Parameterisation Implementation
//

#include "DTCellParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include <iostream>

namespace DTSim {

DTCellParameterisation::DTCellParameterisation(const SuperLayerGeometry& slGeom)
    : G4VPVParameterisation(),
      fSLGeometry(slGeom),
      fTotalCells(slGeom.cells.size())  // Simple: just count cells directly!
{
    G4cout << "DTCellParameterisation created for SuperLayer " << slGeom.slNumber 
           << " with " << fTotalCells << " total cells" << G4endl;
}

DTCellParameterisation::~DTCellParameterisation()
{
}

void DTCellParameterisation::ComputeTransformation(
    const G4int copyNo, G4VPhysicalVolume* physVol) const
{
    // Simple direct indexing: copyNo is the cell index!
    if (copyNo < 0 || copyNo >= fTotalCells) {
        G4cerr << "ERROR: Invalid copy number " << copyNo 
               << " (total cells: " << fTotalCells << ")" << G4endl;
        return;
    }
    
    // Get cell geometry directly from the flat array
    const auto& cellGeom = fSLGeometry.cells[copyNo];
    
    // Set cell position (convert from cm to Geant4 units)
    // Local DT framework: x = x, y = z, z = -y
    G4ThreeVector position(
        cellGeom.x * cm,
        cellGeom.z * cm,
        cellGeom.y * cm
    );
    
    physVol->SetTranslation(position);
}

} // namespace DTSim
