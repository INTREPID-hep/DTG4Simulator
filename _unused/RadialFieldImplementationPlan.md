# Radial Electric Field Implementation Plan

To implement a radial electric field (e.g., 2 kV/cm) pointing to the center of each `DriftCell`, we need a custom field class that handles the coordinate transformation from Global (Geant4 standard) to Local (where the radial direction is defined).

## 1. The "Smart" Field Class

Create a new header (e.g., `include/DriftCellField.hh`) or include this class definition directly. This class inherits from `G4ElectroMagneticField` to support both the existing B-field and the new E-field.

```cpp
#include "G4ElectroMagneticField.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4SystemOfUnits.hh"

class DriftCellField : public G4ElectroMagneticField {
public:
    // Constructor takes the global B-field (to preserve it) and the E-field magnitude
    DriftCellField(G4ThreeVector globalB, G4double eMagnitude) 
        : fBField(globalB), fEMagnitude(eMagnitude) {}

    virtual void GetFieldValue(const G4double Point[4], G4double *Bfield) const override {
        // 1. Set the Magnetic Field (Bfield[0..2]) from the global value
        Bfield[0] = fBField.x();
        Bfield[1] = fBField.y();
        Bfield[2] = fBField.z();

        // 2. Get Local Coordinates to calculate Radial E-field
        // We use the Navigator to find where this point is in the current volume's frame
        auto* navigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
        G4ThreeVector globalPos(Point[0], Point[1], Point[2]);
        
        // Transform Global Point -> Local Point
        G4ThreeVector localPos = navigator->GetGlobalToLocalTransform().TransformPoint(globalPos);

        // 3. Calculate Radial Electric Field
        // Assuming the wire is along the Z-axis in the local DriftCell coordinates
        G4double r = std::sqrt(localPos.x()*localPos.x() + localPos.y()*localPos.y());
        
        if (r > 0.0) {
            // "Radially pointing to the center" means direction is opposite to position vector
            // Magnitude: 2 kV/cm (constant)
            G4double Ex_local = -fEMagnitude * (localPos.x() / r);
            G4double Ey_local = -fEMagnitude * (localPos.y() / r);
            G4double Ez_local = 0.0;

            G4ThreeVector localE(Ex_local, Ey_local, Ez_local);

            // 4. Transform E-field Vector back to Global
            // Note: Vectors transform differently than Points (use Inverse().TransformAxis)
            G4ThreeVector globalE = navigator->GetGlobalToLocalTransform().Inverse().TransformAxis(localE);

            // Set Electric Field (Bfield[3..5])
            Bfield[3] = globalE.x();
            Bfield[4] = globalE.y();
            Bfield[5] = globalE.z();
        } else {
            Bfield[3] = 0.0; Bfield[4] = 0.0; Bfield[5] = 0.0;
        }
    }

    // CRITICAL: Tell Geant4 this field does work (changes energy)
    virtual G4bool DoesFieldChangeEnergy() const override { return true; }

private:
    G4ThreeVector fBField;
    G4double fEMagnitude;
};
```

## 2. Integration in DetectorConstruction

Update `DetectorConstruction::ConstructSDandField` to instantiate and assign this field to the DriftCell logical volumes.

```cpp
void DetectorConstruction::ConstructSDandField()
{
    // ... (Existing SD setup) ...

    // ... (Existing Global Field setup) ...
    // Ensure fGlobalField is set before creating the drift field

    // --- Setup Drift Cell Electric Field ---
    if (!fDriftCellsLogicals.empty()) {
        G4cout << "Creating Radial Electric Field for Drift Cells..." << G4endl;
        
        // Create the field: Superimpose Global B + Radial E (2 kV/cm)
        // Note: Ensure 'kilovolt' and 'cm' are available (G4SystemOfUnits.hh)
        auto* driftField = new DriftCellField(fGlobalField, 2.0*kilovolt/cm);
        
        // Attach this field to ALL DriftCell logical volumes
        // The 'SetLocalField' helper handles the FieldManager creation for you.
        for (auto* logVol : fDriftCellsLogicals) {
            G4FieldBuilder::Instance()->SetLocalField(driftField, logVol);
        }
    }
    
    // Construct the rest of the field setup
    G4FieldBuilder::Instance()->ConstructFieldSetup();
}
```

## 3. Key Implementation Details

*   **Superposition**: The `DriftCellField` must manually re-apply the global magnetic field because assigning a local field manager overrides the global one for that volume.
*   **Navigator**: Using `GetNavigatorForTracking()` allows a single field object to serve all placed copies of the DriftCell. It calculates the field based on the *current* volume's local coordinate system.
*   **Equation of Motion**: By inheriting from `G4ElectroMagneticField` and returning `true` for `DoesFieldChangeEnergy()`, Geant4 should automatically select the correct equation of motion (e.g., `G4EqMagElectricField`) to simulate the particle's trajectory under both electric and magnetic forces.
