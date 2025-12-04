#ifndef DTSimDriftCellDigi_hh
#define DTSimDriftCellDigi_hh 1

#include "G4VDigi.hh"
#include "G4TDigiCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

#include "DTSimTypes.hh"

namespace DTSim
{

class DriftCellDigi : public G4VDigi
{
  public:
    DriftCellDigi();
    DriftCellDigi(const DriftCellDigi&);
    ~DriftCellDigi() override;

    DriftCellDigi& operator=(const DriftCellDigi&);
    G4bool operator==(const DriftCellDigi&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    void Draw() override;
    void Print() override;

    // Setters
    void SetEventID(G4int id) { fEventID = id; }
    void SetCellID(const CellID& id) { fCellID = id; }
    void SetTDC(G4int tdc) { fTDC = tdc; }
    void SetGlobalPos(const G4ThreeVector& pos) { fGlobalPos = pos; }
    
    // Getters
    G4int GetEventID() const { return fEventID; }
    CellID GetCellID() const { return fCellID; }
    G4int GetTDC() const { return fTDC; }
    G4ThreeVector GetGlobalPos() const { return fGlobalPos; }

  private:
    G4int fEventID;
    CellID fCellID;
    G4int fTDC;              // Time-to-Digital Converter value (in TDC counts)
    G4ThreeVector fGlobalPos; // Global position of the hit for visualization
};

// Define digi collection type
typedef G4TDigiCollection<DriftCellDigi> DriftCellDigiCollection;

// Memory allocation
extern G4ThreadLocal G4Allocator<DriftCellDigi>* DriftCellDigiAllocator;

inline void* DriftCellDigi::operator new(size_t)
{
  if (!DriftCellDigiAllocator)
    DriftCellDigiAllocator = new G4Allocator<DriftCellDigi>;
  return (void*)DriftCellDigiAllocator->MallocSingle();
}

inline void DriftCellDigi::operator delete(void* digi)
{
  DriftCellDigiAllocator->FreeSingle((DriftCellDigi*)digi);
}

}

#endif
