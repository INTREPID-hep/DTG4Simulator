#include "DriftCellDigi.hh"
#include "RunAction.hh"

#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4RunManager.hh"

namespace DTSim
{

G4ThreadLocal G4Allocator<DriftCellDigi>* DriftCellDigiAllocator = nullptr;

DriftCellDigi::DriftCellDigi()
 : G4VDigi(), fTDC(-1), fGlobalPos(G4ThreeVector()), fParentPDG(-1), fTrackID(-1)
{
    auto* runManager = G4RunManager::GetRunManager();
    fRunAction = static_cast<const RunAction*>(runManager->GetUserRunAction());
}

DriftCellDigi::DriftCellDigi(const DriftCellDigi& right)
 : G4VDigi()
{
  fCellID = right.fCellID;
  fTDC = right.fTDC;
  fGlobalPos = right.fGlobalPos;
  fParentPDG = right.fParentPDG;
  fTrackID = right.fTrackID;
}

DriftCellDigi::~DriftCellDigi()
{}

DriftCellDigi& DriftCellDigi::operator=(const DriftCellDigi& right)
{
  fCellID = right.fCellID;
  fTDC = right.fTDC;
  fGlobalPos = right.fGlobalPos;
  fParentPDG = right.fParentPDG;
  fTrackID = right.fTrackID;
  return *this;
}

G4bool DriftCellDigi::operator==(const DriftCellDigi& right) const
{
  return (fCellID == right.fCellID && fTDC == right.fTDC);
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
  G4cout << "DriftCellDigi: "
         << " CellID: " << fCellID 
         << " TDC: " << fTDC
         << " GlobalPos: " << fGlobalPos/cm << " cm"
         << " ParentPDG: " << fParentPDG;
  if (fRunAction->IsExtendedActive())
         G4cout << " TrackID: " << fTrackID;
  G4cout << G4endl;
}

void DriftCellDigi::Print(std::ostream& os) const {
  os << "DriftCellDigi: "
     << " CellID: " << fCellID
     << " TDC: " << fTDC
     << " GlobalPos: " << fGlobalPos/cm << " cm"
      << " ParentPDG: " << fParentPDG;
  if (fRunAction->IsExtendedActive()) os << " TrackID: " << fTrackID;
}

std::ostream& operator<<(std::ostream& os, const DriftCellDigi& digi) {
  digi.Print(os);
  return os;
}

}
