#ifndef DTSimDTSegment_hh
#define DTSimDTSegment_hh 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Point3D.hh"

#include "DTSimTypes.hh"

namespace DTSim
{

class DTSegment : public G4VHit
{
  public:
      DTSegment();
      DTSegment(const DTSegment&) = default;
      ~DTSegment() override = default;

      // Operators
      DTSegment& operator=(const DTSegment&) = default;
      G4bool operator==(const DTSegment&) const;

      inline void* operator new(size_t);
      inline void  operator delete(void*);

      // methods from base class
      void Draw() override;
      void Print() override;
       void Print(std::ostream& os) const;

      // Setters
      void SetStationID(const StationID& id)      { fStationID = id; }
      void SetLocalPos(const G4ThreeVector& pos)  { fLocalPos = pos; }
      void SetLocalDir(const G4ThreeVector& dir)  { fLocalDir = dir; }
      void SetGlobalPos(const G4ThreeVector& pos) { fGlobalPos = pos; }
      void SetGlobalDir(const G4ThreeVector& dir) { fGlobalDir = dir; }
      void SetEntryPos(const G4ThreeVector& pos)  { fEntryPos = pos; }
      void SetExitPos(const G4ThreeVector& pos)   { fExitPos = pos; }

      // Getters
      StationID GetStationID() const        { return fStationID; }
      G4int GetWheel() const                { return fStationID.wheel; }
      G4int GetSector() const               { return fStationID.sector; }
      G4int GetStation() const              { return fStationID.station; }
      G4ThreeVector GetLocalPos() const     { return fLocalPos; }
      G4ThreeVector GetLocalDir() const     { return fLocalDir; }
      G4ThreeVector GetGlobalPos() const    { return fGlobalPos; }
      G4ThreeVector GetGlobalDir() const    { return fGlobalDir; }
      G4ThreeVector GetEntryPos() const     { return fEntryPos; }
      G4ThreeVector GetExitPos() const      { return fExitPos; }

  private:
      StationID fStationID;
      G4ThreeVector fLocalPos;
      G4ThreeVector fLocalDir;
      G4ThreeVector fGlobalPos;
      G4ThreeVector fGlobalDir;
      G4Point3D fEntryPos;  // Entry position (global) for visualization
      G4Point3D fExitPos;   // Exit position (global) for visualization
};

// Define hits collection type
typedef G4THitsCollection<DTSegment> DTSegmentCollection;

// Memory allocation
extern G4ThreadLocal G4Allocator<DTSegment>* DTSegmentAllocator;

// Global streaming operator
std::ostream& operator<<(std::ostream& os, const DTSegment& seg);

inline void* DTSegment::operator new(size_t)
{
    if (!DTSegmentAllocator) {
        DTSegmentAllocator = new G4Allocator<DTSegment>;
    }
    return (void*)DTSegmentAllocator->MallocSingle();
}

inline void DTSegment::operator delete(void* hit)
{
    DTSegmentAllocator->FreeSingle((DTSegment*)hit);
}

}

#endif
