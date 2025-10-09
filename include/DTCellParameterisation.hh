//
// DT Cell Parameterisation using mplDTs geometry
//

#ifndef DTCellParameterisation_h
#define DTCellParameterisation_h 1

#include "DTGeometryBuilder.hh"
#include "G4VPVParameterisation.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"

namespace DTSim {

/// Parameterisation class for DT cells
/// Uses geometry data from mplDTs to position cells efficiently
/// Simple indexing: copyNo directly corresponds to cell index
class DTCellParameterisation : public G4VPVParameterisation
{
public:
    DTCellParameterisation(const SuperLayerGeometry& slGeom);
    ~DTCellParameterisation() override;
    
    /// Compute transformation for each cell copy
    void ComputeTransformation(const G4int copyNo, 
                              G4VPhysicalVolume* physVol) const override;
    
    /// Get total number of cells in this SuperLayer
    G4int GetNumberOfCells() const { return fTotalCells; }

private:
    const SuperLayerGeometry& fSLGeometry;
    G4int fTotalCells;
};

} // namespace DTSim

#endif
