///////////////////////////////////////////////////////////////////////
//// <CEPC>                                                        ////  
//// Wedge Geometry for Dual-reaout calorimter			   ////
//// Modification:                                                 ////
//// Option for Mokka   		fucd@ihep.ac.cn		   ////
////                                                               ////
//// Original Author: Mr.Jo Hyunsuk, Kyunpook National University  ////
//// E-Mail: hyunsuk.jo@cern.ch					   ////
//// 					                           ////
///////////////////////////////////////////////////////////////////////
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
// $Id: WGR16DetectorConstruction.hh 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file WGR16DetectorConstruction.hh
/// \brief Definition of the WGR16DetectorConstruction class

#ifndef WGR16DetectorConstruction_h
#define WGR16DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"
#include "G4FieldManager.hh"
#include <string.h>
#include <vector>
#include "dimensionB.hh"
#include "dimensionE.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"

#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Material.hh"
#include "G4VSensitiveDetector.hh"

using namespace std;
class WGR16MagneticField;

class G4VPhysicalVolume;
class G4Material;
class G4VSensitiveDetector;
class G4VisAttributes;
class G4GenericMessenger;
class G4VPhysicalVolume;


/// Detector construction

class WGR16DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    WGR16DetectorConstruction();
    virtual ~WGR16DetectorConstruction();
    
    virtual G4VPhysicalVolume* Construct();
#ifdef WGR16_STANDALONE
    virtual void ConstructSDandField();
#endif
  void Construct(G4LogicalVolume* worldLog, G4VSensitiveDetector* sd = 0);
    void ConstructMaterials();
	 void fiberBR(G4int i,G4double deltatheta_);
	 void fiberBL(G4int i,G4double deltatheta_);
	 void fiberER(G4int i,G4double deltatheta_);
	 void fiberEL(G4int i,G4double deltatheta_);
	 G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

private:
    void DefineCommands();
	
	 G4bool checkOverlaps;
    G4GenericMessenger* fMessenger;

#ifdef WGR16_STANDALONE    
    static G4ThreadLocal WGR16MagneticField* fMagneticField;
	 static G4ThreadLocal G4FieldManager* fFieldMgr;
#endif
	 G4VisAttributes* visAttrC;
         G4VisAttributes* visAttrS;
	 std::vector<G4VisAttributes*> fVisAttributes;
	 
	 G4double innerR;
	 G4double tower_height;

	 G4int NbOfBarrel;
	 G4int NbOfEndcap;
	 G4int NbOfZRot;
	 
	 G4double theta_unit;
	 G4double phi_unit;

	 G4double deltatheta;
	 G4double thetaofcenter;
	 G4double fulltheta;
	 G4double lastdeltatheta;
         //G4ThreeVector pt[8]={G4ThreeVector()};
         G4ThreeVector pt[8];

	 G4double PMTT;

	 dimensionB* dimB;
	 dimensionE* dimE;

	 char name[20];
	 G4Trap* tower;
	 G4Trap* pmtg;
	 G4Trap* pmtcath;
	 G4Tubs* fiber;
	 G4Tubs* fiber_S;
	 G4Tubs* fiber_C;
	 G4Tubs* fiberS;
	 G4Tubs* fiberC;
	 G4VSolid* intersect;
	 G4VSolid* intersect_;

  G4LogicalVolume* fiberCladCLog;
  G4LogicalVolume* fiberCladSLog;

	 G4LogicalVolume* towerLogicalBR[52];
	 G4LogicalVolume* towerLogicalBL[52];
//	 G4LogicalVolume* towerLogicalER[47];
//	 G4LogicalVolume* towerLogicalEL[47];
	 G4LogicalVolume* towerLogicalER[40];
	 G4LogicalVolume* towerLogicalEL[40];
	 
	 G4LogicalVolume* PMTGLogicalBR[52];
	 G4LogicalVolume* PMTGLogicalBL[52];
//	 G4LogicalVolume* PMTGLogicalER[47];
//	 G4LogicalVolume* PMTGLogicalEL[47];
	 G4LogicalVolume* PMTGLogicalER[40];
	 G4LogicalVolume* PMTGLogicalEL[40];

	 G4LogicalVolume* PMTCathLogicalBR[52];
	 G4LogicalVolume* PMTCathLogicalBL[52];
//	 G4LogicalVolume* PMTCathLogicalER[47];
//	 G4LogicalVolume* PMTCathLogicalEL[47];
	 G4LogicalVolume* PMTCathLogicalER[40];
	 G4LogicalVolume* PMTCathLogicalEL[40];

	 vector<G4LogicalVolume*> fiberLogical_BR[52];
	 vector<G4LogicalVolume*> fiberLogical_BR_[52];
	 vector<G4LogicalVolume*> fiberLogical_BL[52];
	 vector<G4LogicalVolume*> fiberLogical_BL_[52];
//	 vector<G4LogicalVolume*> fiberLogical_ER[47];
//	 vector<G4LogicalVolume*> fiberLogical_ER_[47];
//	 vector<G4LogicalVolume*> fiberLogical_EL[47];
//	 vector<G4LogicalVolume*> fiberLogical_EL_[47];
	 vector<G4LogicalVolume*> fiberLogical_ER[40];
	 vector<G4LogicalVolume*> fiberLogical_ER_[40];
	 vector<G4LogicalVolume*> fiberLogical_EL[40];
	 vector<G4LogicalVolume*> fiberLogical_EL_[40];
  
  map<int,G4LogicalVolume*> fiberCLog;
  map<int,G4LogicalVolume*> fiberSLog;
	 
	 G4VSensitiveDetector* PMTSDBR[52];
	 G4VSensitiveDetector* PMTSDBL[52];
//	 G4VSensitiveDetector* PMTSDER[47];
//	 G4VSensitiveDetector* PMTSDEL[47];
	 G4VSensitiveDetector* PMTSDER[40];
	 G4VSensitiveDetector* PMTSDEL[40];
	 
	 G4double clad_C_rMin; 
	 G4double clad_C_rMax; 
	 G4double clad_C_Dz  ; 
	 G4double clad_C_Sphi; 
	 G4double clad_C_Dphi; 

	 G4double core_C_rMin; 
	 G4double core_C_rMax; 
	 G4double core_C_Dz  ; 
	 G4double core_C_Sphi; 
	 G4double core_C_Dphi; 

	 G4double clad_S_rMin; 
	 G4double clad_S_rMax; 
	 G4double clad_S_Dz  ; 
	 G4double clad_S_Sphi; 
	 G4double clad_S_Dphi; 

	 G4double core_S_rMin; 
	 G4double core_S_rMax; 
	 G4double core_S_Dz  ; 
	 G4double core_S_Sphi; 
	 G4double core_S_Dphi; 

	 //---Materials for Cerenkov fiber---
	 G4Material *clad_C_Material;
	 G4Material *core_C_Material; 

	 //---Materials for Scintillation fiber---
	 G4Material *clad_S_Material;
	 G4Material *core_S_Material;

	 //--Material for PMT glass---
	 G4Material *Glass_Material;

	 //--- Material for PMT Photocathod ---
	 G4Material *PMTPC_Material; 


protected:
	 	G4LogicalVolume* fScoringVolume;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
