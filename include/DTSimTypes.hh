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

    // Unique hash for using CellID in unordered_map
    // Packing strategy: W(3bit) Sec(4bit) St(3bit) SL(2bit) L(3bit) Wire(8bit)
    // Or simple decimal packing since ranges are small
    std::size_t getHash() const {
        // Wheel: -2..2 -> 0..4 (add 2)
        // Sector: 1..14
        // Station: 1..4
        // SL: 1..3
        // Layer: 1..4
        // Wire: 1..100
        
        std::size_t hash = 0;
        hash += (wheel + 2);
        hash = hash * 100 + sector;
        hash = hash * 10 + station;
        hash = hash * 10 + superLayer;
        hash = hash * 10 + layer;
        hash = hash * 1000 + wire; // Allow up to 999 wires
        return hash;
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

// Struct to hold station identification (for truth segments)
struct StationID {
    G4int wheel;        // Valid range: -2 to +2
    G4int sector;       // Valid range: 1 to 12 (or 14)
    G4int station;      // Valid range: 1 to 4
    
    // Constructor - initialize to clearly invalid values
    StationID() : wheel(-999), sector(-999), station(-999) {}
    
    StationID(G4int w, G4int sec, G4int st) 
        : wheel(w), sector(sec), station(st) {}
    
    // Check if valid using constants
    bool isValid() const {
        return wheel >= kMinWheel && wheel <= kMaxWheel &&
               sector >= kMinSector && sector <= kMaxSector &&
               station >= kMinStation && station <= kMaxStation;
    }
    
    // Equality operator
    bool operator==(const StationID& other) const {
        return wheel == other.wheel &&
               sector == other.sector &&
               station == other.station;
    }
    
    bool operator!=(const StationID& other) const {
        return !(*this == other);
    }

    // Unique hash for using StationID in unordered_map
    std::size_t getHash() const {
        std::size_t hash = 0;
        hash += (wheel + 2);
        hash = hash * 100 + sector;
        hash = hash * 10 + station;
        return hash;
    }
    
    // Output operator for easy printing
    friend std::ostream& operator<<(std::ostream& os, const StationID& id) {
        os << "W" << id.wheel 
           << "_Sec" << id.sector 
           << "_St" << id.station;
        return os;
    }
};

}

// Specialization of std::hash for CellID and StationID
namespace std {
    template <>
    struct hash<DTSim::CellID> {
        std::size_t operator()(const DTSim::CellID& id) const {
            return id.getHash();
        }
    };

    template <>
    struct hash<DTSim::StationID> {
        std::size_t operator()(const DTSim::StationID& id) const {
            return id.getHash();
        }
    };
}

#endif