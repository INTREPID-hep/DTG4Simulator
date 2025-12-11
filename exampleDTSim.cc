//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file exampleDTSim.cc
/// \brief Main program of the basic DTSim example

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

#include "G4RunManagerFactory.hh"
#include "G4SteppingVerbose.hh"
#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "G4StepLimiterPhysics.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"


int main(int argc, char** argv)
{
  // Determine batch mode and macro file from command line
  G4String macroFileName;
  G4bool batchMode = false;

  if (argc > 1) {
    // First argument is macro file
    macroFileName = argv[1];
    batchMode = true;
  } else {
    // No arguments = interactive mode with default macro
    macroFileName = "macros/interactive/vis.mac";
  }

  // Setup UI for interactive mode
  G4UIExecutive* ui = nullptr;
  if (!batchMode) {
    ui = new G4UIExecutive(argc, argv);
  }


  // Use G4SteppingVerboseWithUnits
  G4int precision = 4;
  G4SteppingVerbose::UseBestUnit(precision);

  // Construct the default run manager
  //
  auto runManager = G4RunManagerFactory::CreateRunManager();
  
  // ------------- Mandatory user initialization classes ---------------
  // Physics list
  auto physicsList = new FTFP_BERT;
  physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  runManager->SetUserInitialization(physicsList);

  // Detector construction
  runManager->SetUserInitialization(new DTSim::DetectorConstruction);

  // User action initialization
  runManager->SetUserInitialization(new DTSim::ActionInitialization());

  // Visualization manager construction
  auto visManager = new G4VisExecutive;
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  auto UImanager = G4UImanager::GetUIpointer();

  // Execute macro file
  UImanager->ApplyCommand("/control/execute " + macroFileName);
  
  // Start interactive session if not in batch mode
  if (!batchMode && ui) {
    ui->SessionStart();
    delete ui;
  }

  // Cleanup
  delete visManager;
  delete runManager;
  
  return 0;
}
