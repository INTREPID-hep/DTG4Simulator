#ifndef RunAction_hh
#define RunAction_hh 1

#include "G4UserRunAction.hh"
#include "G4Types.hh"
#include "G4String.hh"
#include <vector>

class G4Run;
class G4GenericMessenger;

namespace DTSim
{

class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    ~RunAction() override;

    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;

    bool IsExtendedActive() const { return fExtendedOutput; }

  private:
    G4GenericMessenger* fMessenger;
    G4String fOutputFileName;
    G4bool fExtendedOutput;
    G4bool fNtupleAlreadyCreated;
    
  public:
    void DefineCommands();
    void CreateNtupleAndTree();
    // Vector storage for hit data (bound to ntuple columns)
    G4int fNSimHits_column;
    mutable std::vector<G4int> fHit_PDG;
    mutable std::vector<G4int> fHit_Charge;
    mutable std::vector<G4int> fHit_Wheel;
    mutable std::vector<G4int> fHit_Sector;
    mutable std::vector<G4int> fHit_Station;
    mutable std::vector<G4int> fHit_SuperLayer;
    mutable std::vector<G4int> fHit_Layer;
    mutable std::vector<G4int> fHit_Wire;
    mutable std::vector<G4double> fHit_XLocal;
    mutable std::vector<G4double> fHit_YLocal;
    mutable std::vector<G4double> fHit_ZLocal;
    mutable std::vector<G4double> fHit_Time;
    mutable std::vector<G4double> fHit_Edep;
    mutable std::vector<G4int> fHit_ProcessType;
    mutable std::vector<G4int> fHit_TrackID;
    mutable std::vector<G4int> fHit_ParentID;
    mutable std::vector<G4double> fHit_TrackLength;
    mutable std::vector<G4double> fHit_VertexKineticEnergy;
    mutable std::vector<G4double> fHit_VertexPosX;
    mutable std::vector<G4double> fHit_VertexPosY;
    mutable std::vector<G4double> fHit_VertexPosZ;
    
    // Vector storage for digi data (bound to ntuple columns)
    G4int fNDigis_column;
    mutable std::vector<G4int> fDigi_Wheel;
    mutable std::vector<G4int> fDigi_Sector;
    mutable std::vector<G4int> fDigi_Station;
    mutable std::vector<G4int> fDigi_SuperLayer;
    mutable std::vector<G4int> fDigi_Layer;
    mutable std::vector<G4int> fDigi_Wire;
    mutable std::vector<G4int> fDigi_TDC;
    mutable std::vector<G4int> fDigi_parentPDG;
    mutable std::vector<G4int> fDigi_TrackID;
    
    // Vector storage for generator-level data (bound to ntuple columns)
    G4int fNGen_column;
    mutable std::vector<G4int> fGen_PDG;
    mutable std::vector<G4int> fGen_Charge;
    mutable std::vector<G4double> fGen_Pt;
    mutable std::vector<G4double> fGen_Eta;
    mutable std::vector<G4double> fGen_Phi;
    mutable std::vector<G4double> fGen_RadEnergy;
    mutable std::vector<G4int> fGen_nSecondaries;
    
    // Vector storage for segment data (bound to ntuple columns)
    G4int fNSegments_column;
    mutable std::vector<G4int> fSeg_Wheel;
    mutable std::vector<G4int> fSeg_Sector;
    mutable std::vector<G4int> fSeg_Station;
    mutable std::vector<G4double> fSeg_LocalPosX;
    mutable std::vector<G4double> fSeg_LocalPosY;
    mutable std::vector<G4double> fSeg_LocalPosZ;
    mutable std::vector<G4double> fSeg_LocalDirX;
    mutable std::vector<G4double> fSeg_LocalDirY;
    mutable std::vector<G4double> fSeg_LocalDirZ;
    mutable std::vector<G4double> fSeg_GlobalPosX;
    mutable std::vector<G4double> fSeg_GlobalPosY;
    mutable std::vector<G4double> fSeg_GlobalPosZ;
    mutable std::vector<G4double> fSeg_GlobalDirX;
    mutable std::vector<G4double> fSeg_GlobalDirY;
    mutable std::vector<G4double> fSeg_GlobalDirZ;
};

}

#endif