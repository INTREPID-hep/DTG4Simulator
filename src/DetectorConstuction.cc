#include "DetectorConstruction.hh"
#include "DTGeometryBuilder.hh"  // Auto-generated geometry builder

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4Box.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"

#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"


namespace DTSim
{

DetectorConstruction::DetectorConstruction() {}

DetectorConstruction::~DetectorConstruction(){}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    // Construct materials
    ConstructMaterials();

    // Construct world volume
    auto worldPhysical = ConstructWorld();
    
    // Build DT Station using mplDTs geometry
    DTGeometryBuilder builder;
    
    // Build MB2 station: Wheel=-1, Sector=1, Station=2
    fDTStationLogical = builder.BuildStation(-1, 1, 2, fWorldLogical, false);
    
    if (!fDTStationLogical) {
        G4cerr << "ERROR: Failed to build DT Station!" << G4endl;
        return worldPhysical;
    }
    // Optionally build more stations:
    // builder.BuildStation(-1, 1, 1, fWorldLogical, true);  // MB1
    // builder.BuildStation(-1, 1, 3, fWorldLogical, true);  // MB3
    
    // Set visualization attributes
    SetVisualizationAttributes();

    return worldPhysical;
}

G4Material* DetectorConstruction::GetMaterial(const G4String& name)
{
    auto material = G4Material::GetMaterial(name);
    if (!material) {
        G4ExceptionDescription msg;
        msg << "Material " << name << " not found!";
        G4Exception("DetectorConstruction::GetMaterial()",
                    "MyCode0001", FatalException, msg);
    }
    return material;
}

void DetectorConstruction::ConstructMaterials()
{
  auto nistManager = G4NistManager::Instance();

  // Air 
  auto air = nistManager->FindOrBuildMaterial("G4_AIR");
  G4double air_density = air->GetDensity();

  // Gas mixture
  auto gas_mixture = new G4Material("GasMixture", air_density, 2);
  gas_mixture->AddMaterial(nistManager->FindOrBuildMaterial("G4_Ar"), 85*perCent);
  gas_mixture->AddMaterial(nistManager->FindOrBuildMaterial("G4_CARBON_DIOXIDE"), 15*perCent);  
  // Aluminium - Honeycomb (GAP)
  nistManager->FindOrBuildMaterial("G4_Al");
  
  // iron Yoke
  nistManager->FindOrBuildMaterial("G4_Fe");
  
  // Vacuum "Galactic"
  nistManager->FindOrBuildMaterial("G4_Galactic");

  // Vacuum "Air with low density"
  G4double density = 1.0e-5*air_density;
  nistManager
    ->BuildMaterialWithNewDensity("Air_lowDensity", "G4_AIR", density);

  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

G4VPhysicalVolume* DetectorConstruction::ConstructWorld()
{
    // Get materials when needed
    auto defaultMaterial = GetMaterial("G4_Galactic");
    
    // CMS cavern 8m x 8m x 15m
    auto worldSolid = new G4Box("worldBox", 0.5 * 8*m, 0.5 * 8*m, 0.5 * 15*m);
    fWorldLogical = new G4LogicalVolume(worldSolid, defaultMaterial, "worldLogical");
    auto worldPhysical = new G4PVPlacement(
      nullptr, G4ThreeVector(), fWorldLogical, "worldPhysical", nullptr, false, 0, true
    );
    
  return worldPhysical;
}

void DetectorConstruction::SetVisualizationAttributes()
{
  G4VisAttributes invisible(G4VisAttributes::GetInvisible());
  G4VisAttributes blue(G4Colour::Blue());
  G4VisAttributes green(G4Colour::Green());
  G4VisAttributes red(G4Colour::Red());

  fWorldLogical->SetVisAttributes(red);
  fDTStationLogical->SetVisAttributes(green);
}

}
