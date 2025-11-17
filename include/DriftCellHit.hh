#ifndef DTSimDriftCellHit_hh
#define DTSimDriftCellHit_hh 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "DTSimTypes.hh"

namespace DTSim
{

class DriftCellHit : public G4VHit
{
  public:
      DriftCellHit() = default;
      DriftCellHit(const DriftCellHit&) = default;
      ~DriftCellHit() override = default;           

      // Operators
      DriftCellHit& operator=(const DriftCellHit&) = default;
      G4bool operator==(const DriftCellHit&) const;

      inline void* operator new(size_t);
      inline void  operator delete(void*);

      // methods from base class
      void Draw() override;
      void Print() override;

      // Setters
      void SetEventID(G4int id)           { fEventID = id; }
      void SetPDG(G4int pdg)              { fPDG = pdg; }
      void SetCharge(G4int q)             { fCharge = q; }
      void SetCellID(const CellID& id)    { fCellID = id; }
      void SetLocalPos(G4ThreeVector pos) { fLocalPos = pos; }
      void SetTimeDrift(G4double t)       { fTimeDrift = t; }
      void SetEnergyDeposit(G4double edep) { fEnergyDeposit = edep; }

      // Getters
      G4int GetEventID() const            { return fEventID; }
      G4int GetPDG() const                { return fPDG; }
      G4int GetCharge() const             { return fCharge; }
      CellID GetCellID() const            { return fCellID; }
      G4ThreeVector GetLocalPos() const   { return fLocalPos; }
      G4double GetTimeDrift() const       { return fTimeDrift; }
      G4double GetEnergyDeposit() const   { return fEnergyDeposit; }

  private:
      G4int fEventID;
      G4int fPDG;
      G4int fCharge;
      CellID fCellID;
      G4ThreeVector fLocalPos;
      G4double fTimeDrift;
      G4double fEnergyDeposit;  // ← Add this
};

// Define hits collection type
typedef G4THitsCollection<DriftCellHit> DriftCellHitsCollection;

// Memory allocation
extern G4ThreadLocal G4Allocator<DriftCellHit>* DriftCellHitAllocator;

inline void* DriftCellHit::operator new(size_t)
{
    if (!DriftCellHitAllocator) {
        DriftCellHitAllocator = new G4Allocator<DriftCellHit>;
    }
    return (void*)DriftCellHitAllocator->MallocSingle();
}

inline void DriftCellHit::operator delete(void* hit)
{
    DriftCellHitAllocator->FreeSingle((DriftCellHit*) hit);
}

}

#endif