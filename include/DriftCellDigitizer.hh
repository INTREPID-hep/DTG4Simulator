#ifndef DTSimDriftCellDigitizer_hh
#define DTSimDriftCellDigitizer_hh 1

#include "G4VDigitizerModule.hh"
#include "G4String.hh"
#include "G4Types.hh"

#include "DriftCellDigi.hh"

namespace DTSim
{

class DriftCellDigitizer : public G4VDigitizerModule
{
  public:
    DriftCellDigitizer(G4String name);
    ~DriftCellDigitizer() override;

    void Digitize() override;

  private:
    DriftCellDigiCollection* fDigiCollection = nullptr;
    G4int fHCID = -1;
};

}

#endif
