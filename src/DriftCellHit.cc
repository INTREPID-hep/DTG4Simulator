#include "DriftCellHit.hh"
#include "DTSimTypes.hh"
#include "RunAction.hh"

#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4RunManager.hh"
#include <iomanip>


namespace DTSim
{

G4ThreadLocal G4Allocator<DriftCellHit>* DriftCellHitAllocator = nullptr;

DriftCellHit::DriftCellHit()
 : G4VHit(),
   fPDG(-999),
   fCharge(-999),
   fProcessType(-999),
   fCellID(DTSim::CellID()),
   fLocalPos(G4ThreeVector()),
   fGlobalPos(G4ThreeVector()),
   fTimeDrift(-1),
   fEnergyDeposit(-1),
   fTrackID(-1),
   fParentID(-1),
   fTrackLength(-1),
   fVertexKineticEnergy(-1),
   fVertexPos(G4ThreeVector())
{
    auto* runManager = G4RunManager::GetRunManager();
    fRunAction = static_cast<const RunAction*>(runManager->GetUserRunAction());
}

G4bool DriftCellHit::operator==(const DriftCellHit& right) const
{
    // Two hits are equal if they match in all key physics properties
    return (fCellID == right.fCellID &&
            fPDG == right.fPDG &&
            fCharge == right.fCharge &&
            std::abs(fTimeDrift - right.fTimeDrift) < .1);  // Within .1 ns
}

void DriftCellHit::Draw()
{
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if (pVVisManager) {
        // Draw a small circle at the local position
        G4Circle circle(fGlobalPos);
        circle.SetScreenSize(4.);
        circle.SetFillStyle(G4Circle::filled);
        G4VisAttributes attribs(G4Colour::Yellow());
        circle.SetVisAttributes(attribs);
        pVVisManager->Draw(circle);
    }
}

void DriftCellHit::Print()
{
    G4cout << "DriftCellHit: "
           << " PDG=" << fPDG
           << " q=" << fCharge
           << " CellID=" << fCellID
           << " LocalPos=" << fLocalPos/cm << " cm"
           << " TimeDrift=" << G4BestUnit(fTimeDrift, "Time")
           << " GlobalPos=" << fGlobalPos/cm << " cm";
    if (fRunAction->IsExtendedActive()) {
        G4cout << " ProcessType=" << fProcessType
               << " Edep=" << G4BestUnit(fEnergyDeposit, "Energy");    
    }
    G4cout << G4endl;
}

void DriftCellHit::Print(std::ostream& os) const {
    os << "DriftCellHit: "
       << " PDG=" << fPDG
       << " q=" << fCharge
       << " CellID=" << fCellID
       << " TimeDrift=" << G4BestUnit(fTimeDrift, "Time")
       << " LocalPos=" << fLocalPos/cm << " cm"
       << " GlobalPos=" << fGlobalPos/cm << " cm";
    if (fRunAction->IsExtendedActive()) {
        os << " ProcessType=" << fProcessType
           << " Edep=" << G4BestUnit(fEnergyDeposit, "Energy");    
    }
}

std::ostream& operator<<(std::ostream& os, const DriftCellHit& hit) {
    hit.Print(os);
    return os;
}

}