#ifndef RunAction_hh
#define RunAction_hh

#include "G4UserRunAction.hh"
#include <vector>

class G4Run;

namespace DTSim
{

class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    ~RunAction() = default;

    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;
    
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
};

}

#endif