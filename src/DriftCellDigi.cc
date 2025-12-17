#include "DriftCellDigi.hh"

#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

namespace DTSim
{

G4ThreadLocal G4Allocator<DriftCellDigi>* DriftCellDigiAllocator = nullptr;

DriftCellDigi::DriftCellDigi()
 : G4VDigi(), fEventID(-1), fTDC(-1), fGlobalPos(G4ThreeVector()), fTrackID(-1)
{}

DriftCellDigi::DriftCellDigi(const DriftCellDigi& right)
 : G4VDigi()
{
  fEventID = right.fEventID;
  fCellID = right.fCellID;
  fTDC = right.fTDC;
  fGlobalPos = right.fGlobalPos;
  fTrackID = right.fTrackID;
}

DriftCellDigi::~DriftCellDigi()
{}

DriftCellDigi& DriftCellDigi::operator=(const DriftCellDigi& right)
{
  fEventID = right.fEventID;
  fCellID = right.fCellID;
  fTDC = right.fTDC;
  fGlobalPos = right.fGlobalPos;
  fTrackID = right.fTrackID;
  return *this;
}

G4bool DriftCellDigi::operator==(const DriftCellDigi& right) const
{
  return (fEventID == right.fEventID && fCellID == right.fCellID && fTDC == right.fTDC);
}

void DriftCellDigi::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if (pVVisManager) {
    // Draw a filled square at the digi position
    G4Square square(fGlobalPos);
    square.SetScreenSize(6.);  // Slightly larger than hit circles
    square.SetFillStyle(G4Square::filled);
    G4VisAttributes attribs(G4Colour::Red());  // Red for digis
    square.SetVisAttributes(attribs);
    pVVisManager->Draw(square);
  }
}

void DriftCellDigi::Print()
{
  G4cout << "DriftCellDigi - Event: " << fEventID
         << " CellID: " << fCellID 
         << " TDC: " << fTDC
         << " GlobalPos: " << fGlobalPos/cm << " cm"
         << G4endl;
}

}
