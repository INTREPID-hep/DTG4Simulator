#include "DTSegment.hh"

#include "G4Circle.hh"
#include "G4Polyline.hh"
#include "G4Colour.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include <iomanip>

namespace DTSim
{

G4ThreadLocal G4Allocator<DTSegment>* DTSegmentAllocator = nullptr;

DTSegment::DTSegment()
 : G4VHit(),
   fStationID(StationID()),
   fLocalPos(G4ThreeVector()),
   fLocalDir(G4ThreeVector()),
   fGlobalPos(G4ThreeVector()),
   fGlobalDir(G4ThreeVector()),
   fEntryPos(G4ThreeVector()),
   fExitPos(G4ThreeVector())
{}

G4bool DTSegment::operator==(const DTSegment& right) const
{
    // Two segments are equal if they match in station ID and position
    return (fStationID == right.fStationID &&
            (fGlobalPos - right.fGlobalPos).mag() < 1e-6);  // Within 1 nm
}

void DTSegment::Draw()
{
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if (!pVVisManager) return;
    // Draw the trajectory line from entry to exit
    G4Polyline line;
    line.push_back(fEntryPos);
    line.push_back(fExitPos);
    G4VisAttributes lineAttribs(G4Colour::Magenta());
    lineAttribs.SetLineWidth(10.);
    line.SetVisAttributes(lineAttribs);  // Pass pointer, not copy
    pVVisManager->Draw(line);

}

void DTSegment::Print()
{
    G4cout << "DTSegment: " << fStationID
           << " LocalPos=" << fLocalPos/cm << " cm"
           << " LocalDir=" << fLocalDir
           << " GlobalPos=" << fGlobalPos/cm << " cm"
           << " GlobalDir=" << fGlobalDir
           << G4endl;
}

}
