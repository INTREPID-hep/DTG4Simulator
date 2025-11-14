#ifndef DTSimDriftCellSD_hh
#define DTSimDriftCellSD_hh 1

#include "G4VSensitiveDetector.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include <vector>
#include <ostream>

namespace DTSim
{

// Struct to hold complete cell identification
struct CellID {
    G4int wheel;        // Valid range: -2 to +2
    G4int sector;       // Valid range: 1 to 12 (or 14)
    G4int station;      // Valid range: 1 to 4
    G4int superLayer;   // Valid range: 1 to 3
    G4int layer;        // Valid range: 1 to 4
    G4int wire;         // Valid range: 1 to ~60 (depends on geometry)
    
    // Constructor - initialize to clearly invalid values
    CellID() : wheel(-999), sector(-999), station(-999), 
               superLayer(-999), layer(-999), wire(-999) {}
    
    // Check if valid (all values must be in their valid ranges)
    bool isValid() const {
        return wheel >= -2 && wheel <= 2 &&
               sector >= 1 && sector <= 14 &&
               station >= 1 && station <= 4 &&
               superLayer >= 1 && superLayer <= 3 &&
               layer >= 1 && layer <= 4 &&
               wire >= 0;  // wire can vary, just check positive
    }
    
    // Output operator for easy printing
    friend std::ostream& operator<<(std::ostream& os, const CellID& id) {
        os << "W" << id.wheel 
           << "_Sec" << id.sector 
           << "_St" << id.station 
           << "_SL" << id.superLayer 
           << "_L" << id.layer 
           << "_Wire" << id.wire;
        return os;
    }
};

class DriftCellSD : public G4VSensitiveDetector
{
  public:
      DriftCellSD(const G4String& name);
      ~DriftCellSD() override;

      void Initialize(G4HCofThisEvent* hce) override;
      void EndOfEvent(G4HCofThisEvent* hce) override;

      G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;

  private:
      G4double GetTimeWithDrift(G4double time, G4ThreeVector distance) const;
      
      // Decode complete cell ID from volume name and copy number
      CellID DecodeCellID(const G4String& volumeName, G4int copyNo) const;

      // constants
      const G4double kDriftVelocity = 54.0*micrometer/nanosecond;
};

}

#endif