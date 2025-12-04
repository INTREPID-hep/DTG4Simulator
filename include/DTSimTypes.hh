#ifndef DTSimTypes_hh
#define DTSimTypes_hh 1

#include "G4Types.hh"
#include "DTSimConstants.hh"  // ← Include constants for validation
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
    
    // Check if valid using constants
    bool isValid() const {
        return wheel >= kMinWheel && wheel <= kMaxWheel &&
               sector >= kMinSector && sector <= kMaxSector &&
               station >= kMinStation && station <= kMaxStation &&
               superLayer >= kMinSuperLayer && superLayer <= kMaxSuperLayer &&
               layer >= kMinLayer && layer <= kMaxLayer &&
               wire >= 0;  // wire can vary, just check positive
    }
    
    // Equality operator
    bool operator==(const CellID& other) const {
        return wheel == other.wheel &&
               sector == other.sector &&
               station == other.station &&
               superLayer == other.superLayer &&
               layer == other.layer &&
               wire == other.wire;
    }
    
    bool operator!=(const CellID& other) const {
        return !(*this == other);
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

}

#endif