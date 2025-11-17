#include "DriftCellHit.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

namespace DTSim
{

G4ThreadLocal G4Allocator<DriftCellHit>* DriftCellHitAllocator = nullptr;

G4bool DriftCellHit::operator==(const DriftCellHit& right) const
{
    return (this == &right) ? true : false;
}

void DriftCellHit::Draw()
{
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if (pVVisManager) {
        // Draw a small circle at the local position
        G4Circle circle(fLocalPos);
        circle.SetScreenSize(4.);
        circle.SetFillStyle(G4Circle::filled);
        G4VisAttributes attribs(G4Colour::Red());
        circle.SetVisAttributes(attribs);
        pVVisManager->Draw(circle);
    }
}

void DriftCellHit::Print()
{
    G4cout << "DriftCellHit: Event " << fEventID
           << " PDG=" << fPDG
           << " q=" << fCharge
           << " CellID=" << fCellID
           << " LocalPos=" << fLocalPos
           << " TimeDrift=" << G4BestUnit(fTimeDrift, "Time")
           << " Edep=" << G4BestUnit(fEnergyDeposit, "Energy") 
           << G4endl;
}

}