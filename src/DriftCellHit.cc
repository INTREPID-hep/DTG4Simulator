#include "DriftCellHit.hh"

#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include <iomanip>

#include "DTSimTypes.hh"

namespace DTSim
{

G4ThreadLocal G4Allocator<DriftCellHit>* DriftCellHitAllocator = nullptr;

DriftCellHit::DriftCellHit()
 : G4VHit(),
   fEventID(-1),
   fPDG(-999),
   fCharge(-999),
   fProcessType(-999),
   fCellID(DTSim::CellID()),
   fLocalPos(G4ThreeVector()),
   fGlobalPos(G4ThreeVector()),
   fTimeDrift(-1),
   fEnergyDeposit(-1)
{}

G4bool DriftCellHit::operator==(const DriftCellHit& right) const
{
    // Two hits are equal if they match in all key physics properties
    return (fEventID == right.fEventID && 
            fCellID == right.fCellID &&
            fPDG == right.fPDG &&
            fCharge == right.fCharge &&
            fProcessType == right.fProcessType &&
            std::abs(fTimeDrift - right.fTimeDrift) < .1 &&  // Within .1 ns
            std::abs(fEnergyDeposit - right.fEnergyDeposit) < 1e-6);  // Within 1 eV
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
    G4cout << "DriftCellHit: Event " << fEventID
           << " PDG=" << fPDG
           << " ProcessType=" << fProcessType
           << " q=" << fCharge
           << " Edep=" << G4BestUnit(fEnergyDeposit, "Energy")
           << " CellID=" << fCellID
           << " LocalPos=" << fLocalPos/cm << " cm"
           << " GlobalPos=" << fGlobalPos/cm << " cm"
           << " TimeDrift=" << G4BestUnit(fTimeDrift, "Time")
           << G4endl;
}

}