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

  private:
    G4GenericMessenger* fMessenger;
    G4String fOutputFileName;
    
  public:
    // Vector storage for hit data (bound to ntuple columns)
    std::vector<G4int> fHit_PDG;
    std::vector<G4int> fHit_Charge;
    std::vector<G4int> fHit_Wheel;
    std::vector<G4int> fHit_Sector;
    std::vector<G4int> fHit_Station;
    std::vector<G4int> fHit_SuperLayer;
    std::vector<G4int> fHit_Layer;
    std::vector<G4int> fHit_Wire;
    std::vector<G4double> fHit_XLocal;
    std::vector<G4double> fHit_YLocal;
    std::vector<G4double> fHit_ZLocal;
    std::vector<G4double> fHit_Time;
    std::vector<G4double> fHit_Edep;
    std::vector<G4int> fHit_ProcessType;
    
    // Vector storage for digi data (bound to ntuple columns)
    std::vector<G4int> fDigi_Wheel;
    std::vector<G4int> fDigi_Sector;
    std::vector<G4int> fDigi_Station;
    std::vector<G4int> fDigi_SuperLayer;
    std::vector<G4int> fDigi_Layer;
    std::vector<G4int> fDigi_Wire;
    std::vector<G4int> fDigi_TDC;
    
    // Vector storage for generator-level data (bound to ntuple columns)
    std::vector<G4int> fGen_PDG;
    std::vector<G4int> fGen_Charge;
    std::vector<G4double> fGen_Pt;
    std::vector<G4double> fGen_Eta;
    std::vector<G4double> fGen_Phi;
    
    // Vector storage for segment data (bound to ntuple columns)
    std::vector<G4int> fSeg_Wheel;
    std::vector<G4int> fSeg_Sector;
    std::vector<G4int> fSeg_Station;
    std::vector<G4double> fSeg_LocalPosX;
    std::vector<G4double> fSeg_LocalPosY;
    std::vector<G4double> fSeg_LocalPosZ;
    std::vector<G4double> fSeg_LocalDirX;
    std::vector<G4double> fSeg_LocalDirY;
    std::vector<G4double> fSeg_LocalDirZ;
    std::vector<G4double> fSeg_GlobalPosX;
    std::vector<G4double> fSeg_GlobalPosY;
    std::vector<G4double> fSeg_GlobalPosZ;
    std::vector<G4double> fSeg_GlobalDirX;
    std::vector<G4double> fSeg_GlobalDirY;
    std::vector<G4double> fSeg_GlobalDirZ;
};

}

#endif