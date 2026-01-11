# Future Plan: Flexible User Limits Implementation

## Overview
The current implementation uses hardcoded member variables and UI commands for specific regions (Station, Yoke). This plan proposes a more flexible "Smart Parser" approach that allows setting limits for **any** region dynamically without modifying the C++ code.

## Proposed UI Syntax
Instead of specific commands like `/DTSim/detector/setStationStepLimit`, we will use generic commands that take the region name as an argument.

```plaintext
# Syntax: /DTSim/limits/setStepLimit <RegionName> <Value> <Unit>

/DTSim/limits/setStepLimit StationRegion 6.5 mm
/DTSim/limits/setStepLimit YokeRegion 10.0 cm
/DTSim/limits/setTrackLimit DriftCellRegion 50.0 cm
/DTSim/limits/setMinEkine StationRegion 100.0 keV
```

## C++ Implementation Strategy

### 1. Data Structure
Replace individual `G4UserLimits*` pointers with a map.

```cpp
// DetectorConstruction.hh
#include <map>

class DetectorConstruction : ... {
private:
    // Map region names to their limit objects
    std::map<G4String, G4UserLimits*> fLimitMap;
    
    // Helper to manage the map
    G4UserLimits* GetOrCreateLimits(const G4String& regionName);
    
    // ...
};
```

### 2. Command Handling
Define a single method for each limit type that parses the string argument.

```cpp
// DetectorConstruction.cc

void DetectorConstruction::SetStepLimit(G4String command)
{
    // 1. Parse the command string
    std::stringstream ss(command);
    G4String regionName;
    G4double value;
    G4String unitName;
    
    ss >> regionName >> value >> unitName;
    
    // 2. Convert unit
    G4double limit = value * G4UnitDefinition::GetValueOf(unitName);

    // 3. Get or Create UserLimits for this region
    G4UserLimits* limits = GetOrCreateLimits(regionName);
    
    // 4. Apply
    limits->SetMaxAllowedStep(limit);
    
    // 5. If geometry is already built, ensure the region has this limit attached
    // (This handles the case where limits are set AFTER initialization)
    G4Region* region = G4RegionStore::GetInstance()->GetRegion(regionName, false);
    if (region) {
        region->SetUserLimits(limits);
    }
}
```

### 3. Helper Method
```cpp
G4UserLimits* DetectorConstruction::GetOrCreateLimits(const G4String& regionName)
{
    if (fLimitMap.find(regionName) == fLimitMap.end()) {
        fLimitMap[regionName] = new G4UserLimits();
    }
    return fLimitMap[regionName];
}
```

### 4. Cleanup
In the destructor, iterate through the map and delete all `G4UserLimits` objects.

```cpp
DetectorConstruction::~DetectorConstruction() {
    for (auto const& [name, limit] : fLimitMap) {
        delete limit;
    }
}
```
