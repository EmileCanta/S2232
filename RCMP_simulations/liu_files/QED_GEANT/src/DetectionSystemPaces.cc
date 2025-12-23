#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemPaces.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units


DetectionSystemPaces::DetectionSystemPaces() :
	// Logical Volumes
	fAluminumHemisphereLog(0),
	fAluminumAnnulusTopLog(0),
	fAluminumAnnulusBotLog(0),
	fCanisterLog(0),
	fSiliconBlockLog(0),
	fSiliconDeadLayerLog(0),

	fTeflonAnnulusTopLog(0),
	fTeflonAnnulusBotLog(0),
	fDelrinHemisphereLog(0),
	fDelrinShellLog(0),
	fDelrinShell2Log(0)
	
{
	/*Measurement status*/
	//~   estimate
	//OK  confirmed measurement

	// Cut clearance
	fCutClearance = 0.01*mm;

	// Aluminum hemisphere
	fAluminumHemisphereInnerRadius =         (82.50 / 2)*mm; //OK
	fAluminumHemisphereOuterRadius =         (88.86 / 2)*mm; //OK
	fAluminumHemisphereBeamHoleRadius =     (34.50 / 2)*mm; //~
	fAluminumHemisphereBeamHoleRimHeight = (3.0)*mm; //~

	// Annulus that sits in front of silicon crystal
	fAluminumAnnulusTopInnerRadius =        (17.5 / 2)*mm; //OK
	fAluminumAnnulusTopOuterRadius =        (24. / 2)*mm; //OK
	fAluminumAnnulusTopThickness =           (1.5)*mm; //OK
	//fAluminumAnnulusTopThickness   =	0.001*mm;
	// Annulus that sits in back of silicon crystal
	fAluminumAnnulusBotInnerRadius =        (0./ 2)*mm; //~
	fAluminumAnnulusBotOuterRadius =        (24. / 2)*mm; //~
	fAluminumAnnulusBotThickness =           (1.5)*mm; //~

	// Canister holding crystal
	fCanisterInnerRadius =                    (15.86 / 2)*mm; //~
	fCanisterOuterRadius =                    (21.70 / 2)*mm; //~
	fCanisterThickness =                       (5.0)*mm; //~

	// Silicon crystal detector
	fSiliconBlockRadius =                     (15.86 / 2)*mm; //~
	fSiliconBlockThickness =                  (4.73-0.002*2)*mm; //~
	fSiliconDeadLayerThickness =             (0.0002)*mm; //~ Nominal Structure Stopping Power of Window - Ortec, 2u Si equivalent


	// Screws
	fScrewRadius =                             (1.00 / 2)*mm; //~
	fScrewPlacementRadius =                   (22.50 / 2)*mm; //~

	fConvert = 0.0254*m; //Needed to convert inches to meters


	// Teflon bases
	fTeflonAnnulusTopInnerRadius =          (21.00 / 2)*mm; //~
	fTeflonAnnulusTopOuterRadius =          (25.00 / 2)*mm; //~
	fTeflonAnnulusTopThickness =             (0.0)*mm; //OK
	fTeflonAnnulusBotInnerRadius =          (25.65 / 2)*mm; //~
	fTeflonAnnulusBotOuterRadius =          (30.50 / 2)*mm; //~
	fTeflonAnnulusBotThickness =             (0.00)*mm; //~

	// delrin hemisphere
	fDelrinHemisphereInnerRadius =           (84.4)*mm; //OK
	fDelrinHemisphereOuterRadius =           (89.4)*mm; //OK
	fDelrinHemisphereBeamHoleRadius =       (20.00 / 2)*mm; //OK
	
	
	fDelrinInnerRadius = 3.275*fConvert;
	fDelrinOuterRadius = 10.475*fConvert;
	fDelrin2InnerRadius = 8.90*cm;
	fDelrin2OuterRadius = 9.90*cm;
	fHevimetInnerRadius = 9.906*cm;
	fHevimetOuterRadius = 12.446*cm;
	fDelrinHoleRadius = 1.25*fConvert;

	// placement distances
	//fAluminumHemisphereDist =     (0.)*mm;
	//fDelrinHemisphereDist =       (0.)*mm;
	fAluminumAnnulusTopDist =    (fAluminumAnnulusTopThickness / 2);
	fCanisterDist =              (fAluminumAnnulusTopThickness) + (fCanisterThickness / 2) ;
	fSiliconBlockDist =          (fAluminumAnnulusTopThickness) + (fCanisterThickness-fSiliconBlockThickness)+ (fSiliconBlockThickness) / 2; 
	fSiliconDeadLayerFrontDist = (fAluminumAnnulusTopThickness) + (fSiliconDeadLayerThickness / 2) + (fCanisterThickness-fSiliconBlockThickness);
	fSiliconDeadLayerBackDist  = (fAluminumAnnulusTopThickness) + (fSiliconBlockThickness) + (fSiliconDeadLayerThickness/2)+ (fCanisterThickness-fSiliconBlockThickness);
	//fTeflonAnnulusTopDist =      (fAluminumAnnulusTopThickness) + (fSiliconBlockThickness) ;
	//fTeflonAnnulusBotDist =      (fTeflonAnnulusTopDist) + (fTeflonAnnulusTopThickness / 2) + (fTeflonAnnulusBotThickness / 2);
	fAluminumAnnulusBotDist =    ((fAluminumAnnulusTopThickness) + (fCanisterThickness)) + (fAluminumAnnulusBotThickness / 2);

	

	// CEDAR placement with respect to the center surface to center distance (41.5-((fAluminumAnnulusTopThickness) + (fCanisterThickness-fSiliconBlockThickness)) )
	fPacesPlacementDistance[0] = (41.5-((fAluminumAnnulusTopThickness) + (fCanisterThickness-fSiliconBlockThickness)) )*mm; 
	fPacesPlacementDistance[1] = (41.5-((fAluminumAnnulusTopThickness) + (fCanisterThickness-fSiliconBlockThickness)) )*mm;
	fPacesPlacementDistance[2] = (41.5-((fAluminumAnnulusTopThickness) + (fCanisterThickness-fSiliconBlockThickness)) )*mm;
	fPacesPlacementDistance[3] = (41.5-((fAluminumAnnulusTopThickness) + (fCanisterThickness-fSiliconBlockThickness)) )*mm;
	fPacesPlacementDistance[4] = (41.5-((fAluminumAnnulusTopThickness) + (fCanisterThickness-fSiliconBlockThickness)) )*mm;
	fPacesPlacementDistance[5] = (41.5-((fAluminumAnnulusTopThickness) + (fCanisterThickness-fSiliconBlockThickness)) )*mm;
	fPacesPlacementDistance[6] = (41.5-((fAluminumAnnulusTopThickness) + (fCanisterThickness-fSiliconBlockThickness)) )*mm;
	fPacesPlacementDistance[7] = (41.5-((fAluminumAnnulusTopThickness) + (fCanisterThickness-fSiliconBlockThickness)) )*mm;
	fPacesPlacementDistance[8] = (41.5-((fAluminumAnnulusTopThickness) + (fCanisterThickness-fSiliconBlockThickness)) )*mm;
	fPacesPlacementDistance[9] = (41.5-((fAluminumAnnulusTopThickness) + (fCanisterThickness-fSiliconBlockThickness)) )*mm;
	fPacesPlacementDistance[10] = (41.5-((fAluminumAnnulusTopThickness) + (fCanisterThickness-fSiliconBlockThickness)) )*mm;
	fPacesPlacementDistance[11] = (41.5-((fAluminumAnnulusTopThickness) + (fCanisterThickness-fSiliconBlockThickness)) )*mm;
	//fPacesPlacementDistance[12] = 44.5*mm;	
	//fPacesPlacementDistance[13] = 44.5*mm;
	//fPacesPlacementDistance[14] = 44.5*mm;
	//fPacesPlacementDistance[15] = 44.5*mm;
	
	

	// phi angle for top crystals
	//fPacesPlacementPhi[0] =   0.*deg; //OK
	fPacesPlacementPhi[0] =   67.5*deg; //OK
	fPacesPlacementPhi[1] =  157.5*deg; //OK
	fPacesPlacementPhi[2] = 247.5*deg; //OK
	fPacesPlacementPhi[3] = 337.5*deg; //OK
	
	fPacesPlacementPhi[4] = 22.5*deg; //OK
	fPacesPlacementPhi[5] = 112.5*deg; //OK
	fPacesPlacementPhi[6] = 202.5*deg; //OK
	fPacesPlacementPhi[7] = 292.5*deg; //OK
	
	fPacesPlacementPhi[8] = 67.5*deg; //OK
	fPacesPlacementPhi[9] = 157.5*deg; //OK
	fPacesPlacementPhi[10] = 247.5*deg; //OK
	fPacesPlacementPhi[11] = 337.5*deg; //OK
	
	
	
	
	
	
	//top crystals at theta angle of 108*deg
	//fPacesPlacementTheta[0] = 0.*deg; //~
	fPacesPlacementTheta[0] = 143.5*deg; //~
	fPacesPlacementTheta[1] = 143.5*deg; //~
	fPacesPlacementTheta[2] = 143.5*deg; //~
	fPacesPlacementTheta[3] = 143.5*deg; //~
	//shadowing LaBr3 at 125.3 deg
	fPacesPlacementTheta[4] = 125.3*deg; //~
	fPacesPlacementTheta[5] = 125.3*deg; //~
	fPacesPlacementTheta[6] = 125.3*deg; //~
	fPacesPlacementTheta[7] = 125.3*deg; //~
	
	fPacesPlacementTheta[8] = 112.2*deg; //~
	fPacesPlacementTheta[9] = 112.2*deg; //~
	fPacesPlacementTheta[10] = 112.2*deg; //~
	fPacesPlacementTheta[11] = 112.2*deg; //~
	




	
	
	//orientation of CEDAR detectors

	fPacesOrientationPhi[0] = 0.0*deg;
	fPacesOrientationPhi[1] = 0.0*deg;
	fPacesOrientationPhi[2] = 0.0*deg;
	fPacesOrientationPhi[3] = 0.0*deg;
	fPacesOrientationPhi[4] = 0.0*deg;
	fPacesOrientationPhi[5] = 0.0*deg;
	fPacesOrientationPhi[6] = 0.0*deg;
	fPacesOrientationPhi[7] = 0.0*deg;
	fPacesOrientationPhi[8] = 0.0*deg;
	fPacesOrientationPhi[9] = 0.0*deg;
	fPacesOrientationPhi[10] = 0.0*deg;
	fPacesOrientationPhi[11] = 0.0*deg;
	
	





	fPacesOrientationTheta[0] = 0.*deg; //OK
	fPacesOrientationTheta[1] = 0.*deg; //OK
	fPacesOrientationTheta[2] =-0.*deg; //OK
	fPacesOrientationTheta[3] = 0.*deg; //OK
	fPacesOrientationTheta[4] = 0.*deg; //OK
	fPacesOrientationTheta[5] = 0.*deg; //OK
	fPacesOrientationTheta[6] = 0.*deg; //OK
	fPacesOrientationTheta[7] =0.0*deg; //OK
	fPacesOrientationTheta[8] = 0.*deg; //OK
	fPacesOrientationTheta[9] = 0.*deg; //OK
	fPacesOrientationTheta[10] = 0.*deg; //OK
	fPacesOrientationTheta[11] =0.0*deg; //OK
	



	// assembly volume sizes
	fDetectorAssemblyRadius =    (34.50 / 2)*mm; //?
	fDetectorAssemblyThickness = (15.0)*mm; //?
}

DetectionSystemPaces::~DetectionSystemPaces() {
	// logical volumes
	//delete fAluminumHemisphereLog;
	delete fAluminumAnnulusTopLog;
	delete fAluminumAnnulusBotLog;
	delete fCanisterLog;
	delete fSiliconBlockLog;
	delete fSiliconDeadLayerLog;
	//delete fTeflonAnnulusTopLog;
	//delete fTeflonAnnulusBotLog;
	delete fDelrinHemisphereLog;
	
	delete fDelrinShellLog;
	delete fDelrinShell2Log;
}

G4int DetectionSystemPaces::Build() {

	// Build assembly volumes
	fAssembly = new G4AssemblyVolume();
	fAssemblyDetector = new G4AssemblyVolume();
	fAssemblySilicon = new G4AssemblyVolume();

	// Add silicon assembly
	AddSiliconBlock();
	AddSiliconDeadLayer();

	// Add detector assembly
	AddCanister();
	AddAluminumAnnulusTop();
	//AddTeflonAnnulusTop();
	//AddTeflonAnnulusBot();
	AddAluminumAnnulusBot();
	//AddScrews();
	//ConstructDelrinShell();

	// Add total assembly
	//AddAluminumHemisphere();

	// Assemble the ensemble
	
	CombineAssemblySilicon();
	CombineAssemblyDetector();
	//CombineAssembly();

	return 1;
}//end ::Build

G4int DetectionSystemPaces::PlaceDetector(G4LogicalVolume* expHallLog, G4int ndet) {
	//place detectors
	G4int dI;
	G4double* ptrPd = fPacesPlacementDistance;
	G4double* ptrPt = fPacesPlacementPhi;//t=phi, but p=theta???????
	G4double* ptrPp = fPacesPlacementTheta;
	G4double* ptrOt = fPacesOrientationPhi;
	G4double* ptrOp = fPacesOrientationTheta;
	/*G4double* ptrPp = fPacesPlacementPhi;//changed code, rotated by 90 degrees
	  G4double* ptrPt = fPacesPlacementTheta;
	  G4double* ptrOp = fPacesOrientationPhi;
	  G4double* ptrOt = fPacesOrientationTheta;*/
	if(ndet > 12 || ndet < 0) ndet = 12;
	G4double dDist, dPhi, dTheta, oriPhi, oriTheta;
	G4RotationMatrix* Ra;
	G4ThreeVector Ta, yprimeaxis;
	for(G4int i=0; i<ndet; i++) {
		dI = i;
		dDist = ptrPd[dI];
		dPhi = ptrPt[dI];
		dTheta = ptrPp[dI];
		Ra = new G4RotationMatrix;
		Ta.setX(dDist * cos(dPhi) * sin(dTheta));
		Ta.setY(dDist * sin(dPhi) * sin(dTheta));
		Ta.setZ(dDist *      1.0     * cos(dTheta));

		oriPhi = dPhi + ptrOt[dI] + M_PI/2; //plus 90 deg
		oriTheta = dTheta + ptrOp[dI];
		yprimeaxis = G4ThreeVector(cos(oriPhi), sin(oriPhi), 0);
		Ra->set(yprimeaxis, oriTheta);

		//G4cout<<"----------- dI = "<<dI<<G4endl;
		fAssemblySilicon->MakeImprint(expHallLog, Ta, Ra, dI+1);
		fAssemblyDetector->MakeImprint(expHallLog, Ta, Ra, dI*20);
	}

	//place hemisphere
	G4RotationMatrix* R0 = new G4RotationMatrix; G4ThreeVector T0;
	//fAssembly->MakeImprint(expHallLog, T0, R0, 5*40);

	return 1;
}//end ::PlaceDetector

G4int DetectionSystemPaces::CombineAssemblySilicon() {
	G4RotationMatrix* Ra = new G4RotationMatrix; G4ThreeVector Ta;
	Ta.setZ(fSiliconBlockDist);
	fAssemblySilicon->AddPlacedVolume(fSiliconBlockLog, Ta, Ra);

	return 1;
}//end ::CombineAssemblySilicon


G4int DetectionSystemPaces::CombineAssemblyDetector() {
	G4RotationMatrix* Ra = new G4RotationMatrix; G4ThreeVector Ta;
	Ta.setZ(fSiliconDeadLayerFrontDist);
	fAssemblyDetector->AddPlacedVolume(fSiliconDeadLayerLog, Ta, Ra);
	Ta.setZ(fSiliconDeadLayerBackDist);
	fAssemblyDetector->AddPlacedVolume(fSiliconDeadLayerLog, Ta, Ra);
	Ta.setZ(fCanisterDist);
	fAssemblyDetector->AddPlacedVolume(fCanisterLog, Ta, Ra);
	Ta.setZ(fAluminumAnnulusTopDist);
	fAssemblyDetector->AddPlacedVolume(fAluminumAnnulusTopLog, Ta, Ra);
	Ta.setZ(fAluminumAnnulusBotDist);
	fAssemblyDetector->AddPlacedVolume(fAluminumAnnulusBotLog, Ta, Ra);
	// Ta.setZ(fTeflonAnnulusTopDist);
	// fAssemblyDetector->AddPlacedVolume(fTeflonAnnulusTopLog, Ta, Ra);
	// Ta.setZ(fTeflonAnnulusBotDist);
	// fAssemblyDetector->AddPlacedVolume(fTeflonAnnulusBotLog, Ta, Ra);

	return 1;
}//end ::CombineAssemblyDetector

G4int DetectionSystemPaces::CombineAssembly() {
	//place aluminum hemisphere
	G4RotationMatrix* Ra = new G4RotationMatrix; G4ThreeVector Ta;
	//Ta.setZ(fAluminumHemisphereDist);
	//fAssembly->AddPlacedVolume(fAluminumHemisphereLog, Ta, Ra);
	//fAssembly->AddPlacedVolume(fDelrinShellLog, Ta, Ra);

	return 1;
}//end ::CombineAssembly

//Add silicon detector
G4int DetectionSystemPaces::AddSiliconBlock() {
	//material
	G4Material* material = G4Material::GetMaterial("Silicon");
	if(!material) {
		G4cout<<" ----> Material "<<"Silicon"<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	//vis attributes
	G4VisAttributes* siliconBlockVisAtt = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
	siliconBlockVisAtt->SetVisibility(true);

	//silicon block (without dead layer, smaller than actual)
	G4double deadLayerCut = fSiliconDeadLayerThickness + fCutClearance;
	G4double innerRadius = 0.*mm;
	G4double outerRadius = fSiliconBlockRadius;
	G4double halfLengthZ = fSiliconBlockThickness/2.0 - deadLayerCut;

	//primitive volume
	G4Tubs* siliconBlock = new G4Tubs("siliconBlock", innerRadius, outerRadius, halfLengthZ, 0*M_PI, 2*M_PI);

	//logical volume
	if(fSiliconBlockLog == nullptr) {
		fSiliconBlockLog = new G4LogicalVolume(siliconBlock, material, "pacesSiliconBlockLog", 0, 0, 0); // Renamed from "siliconBlockLog" to "pacesSiliconBlockLog"
		fSiliconBlockLog->SetVisAttributes(siliconBlockVisAtt);
	}

	return 1;
}//end ::end AddSiliconBlock

//Add detector dead layer
G4int DetectionSystemPaces::AddSiliconDeadLayer() {
	//material
	G4Material* material = G4Material::GetMaterial("Silicon"); // "Nominal Structure Stopping Power of Window-Ortec, 2u Si"
	if(!material) {
		G4cout<<" ----> Material "<<"Silicon"<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	//vis attributes
	G4VisAttributes* siliconDeadLayerVisAtt = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
	siliconDeadLayerVisAtt->SetVisibility(true);

	//dead layer
	G4double innerRadius = 0.*mm;
	G4double outerRadius = fSiliconBlockRadius;
	G4double halfLengthZ = fSiliconDeadLayerThickness / 2.0;

	//primitive volume
	G4Tubs* siliconDeadLayer = new G4Tubs("siliconDeadLayer", innerRadius, outerRadius, halfLengthZ, 0*M_PI, 2*M_PI);

	//logical volume
	if(fSiliconDeadLayerLog == nullptr) {
		fSiliconDeadLayerLog = new G4LogicalVolume(siliconDeadLayer, material, "siliconDeadLayerLog", 0, 0, 0);
		fSiliconDeadLayerLog->SetVisAttributes(siliconDeadLayerVisAtt);
	}

	return 1;
}//end ::end AddSiliconDeadLayer

//Add detector canister
G4int DetectionSystemPaces::AddCanister() {
	//material
	G4Material* material = G4Material::GetMaterial("Silicon");  // would be Silicon, Zhimin
	if(!material) {
		G4cout<<" ----> Material "<<"Silicon"<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	//vis attributes
	G4VisAttributes* canisterVisAtt = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	canisterVisAtt->SetVisibility(true);

	//main annulus
	G4double innerRadius = fCanisterInnerRadius + fCutClearance;
	G4double outerRadius = fCanisterOuterRadius;
	G4double halfLengthZ = fCanisterThickness/2.0;

	//primitive volume
	G4Tubs* canister = new G4Tubs("canister", innerRadius, outerRadius, halfLengthZ, 0*M_PI, 2*M_PI);

	//logical volume
	if(fCanisterLog == nullptr) {
		fCanisterLog = new G4LogicalVolume(canister, material, "canisterLog", 0, 0, 0);
		fCanisterLog->SetVisAttributes(canisterVisAtt);
	}

	return 1;
}//end ::end AddCanister

//Add annulus that is on top of silicon detector
G4int DetectionSystemPaces::AddAluminumAnnulusTop() {
	G4double z;
	G4double a = 14.00*g/mole;
	G4NistManager* man = G4NistManager::Instance();
	G4Element* Nit = man->FindOrBuildElement("N");
	G4Element* Oxy = man->FindOrBuildElement("O");
	G4Element* Alu = man->FindOrBuildElement("Al");
	
	
	G4int ncomp = 2;
	// total nucleon of Al(NO3)3 number is 294
	G4double density = 0.52068 *g/cm3;
	G4Material* Nitrate = new G4Material("Nitrate", density, ncomp);
	Nitrate->AddElement(Nit, 1);
	Nitrate->AddElement(Oxy, 3);
	G4double abundance, fracMass;
	
	density = 1.72 *g/cm3;
	G4Material* material = new G4Material("AlN", density, ncomp);
	G4int nAtoms;
	
	material->AddElement(Alu,fracMass=9.18*perCent);
	material->AddMaterial(Nitrate,fracMass=90.82*perCent);





	//material
	// G4Material* material = G4Material::GetMaterial("Aluminum");
	if(!material) {
		G4cout<<" ----> Material "<<"Aluminum Nitrate"<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	//vis attributes
	G4VisAttributes* aluminumAnnulusTopVisAtt = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
	aluminumAnnulusTopVisAtt->SetVisibility(true);

	//main annulus
	G4double innerRadius = fAluminumAnnulusTopInnerRadius;
	G4double outerRadius = fAluminumAnnulusTopOuterRadius;
	G4double halfLengthZ = fAluminumAnnulusTopThickness/2.0;
	//screw cuts
	G4double cutInnerRadius = 0.*mm;
	G4double cutOuterRadius = fScrewRadius;
	G4double cutHalfLengthZ = 20.*mm;

	//primitive volume
	G4Tubs* alAnToWithoutCuts = new G4Tubs("alAnToWithoutCuts", innerRadius, outerRadius, halfLengthZ, 0*M_PI, 2*M_PI);
	G4Tubs* cutScrewHole = new G4Tubs("cutScrewHole", cutInnerRadius, cutOuterRadius, cutHalfLengthZ, 0*M_PI, 2*M_PI);

	//cut out screw holes
	G4double cI;
	G4double cR = fScrewPlacementRadius;
	G4double cDt = 45.*deg;
	G4ThreeVector moveCut;
	//cI = 0; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	//G4SubtractionSolid* alAnTo0 = new G4SubtractionSolid("alAnTo0", alAnToWithoutCuts, cutScrewHole, 0, moveCut);
	//cI = 1; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	//G4SubtractionSolid* alAnTo1 = new G4SubtractionSolid("alAnTo1", alAnTo0, cutScrewHole, 0, moveCut);
	//cI = 2; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	//G4SubtractionSolid* alAnTo2 = new G4SubtractionSolid("alAnTo2", alAnTo1, cutScrewHole, 0, moveCut);
	//cI = 3; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
//	G4SubtractionSolid* alAnTo3 = new G4SubtractionSolid("alAnTo3", alAnTo2, cutScrewHole, 0, moveCut);
//	cI = 4; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
//	G4SubtractionSolid* alAnTo4 = new G4SubtractionSolid("alAnTo4", alAnTo3, cutScrewHole, 0, moveCut);
//	cI = 5; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
//	G4SubtractionSolid* alAnTo5 = new G4SubtractionSolid("alAnTo5", alAnTo4, cutScrewHole, 0, moveCut);
//	cI = 6; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
//	G4SubtractionSolid* alAnTo6 = new G4SubtractionSolid("alAnTo6", alAnTo5, cutScrewHole, 0, moveCut);
	//final cut
//	cI = 7; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	//G4SubtractionSolid* aluminumAnnulusTop = new G4SubtractionSolid("aluminumAnnulusTop", alAnTo6, cutScrewHole, 0, moveCut);

	//logical volume
	if(fAluminumAnnulusTopLog == nullptr) {
		fAluminumAnnulusTopLog = new G4LogicalVolume(alAnToWithoutCuts, material, "aluminumAnnulusTopLog", 0, 0, 0);
		fAluminumAnnulusTopLog->SetVisAttributes(aluminumAnnulusTopVisAtt);
	}

	return 1;
}//end ::end AddAluminumAnnulusTop

//Add annulus that is on bottom of silicon detector
G4int DetectionSystemPaces::AddAluminumAnnulusBot() {
	//material
	// G4Material* material = G4Material::GetMaterial("Aluminum");
	G4double z;
	G4double a = 14.00*g/mole;
	G4NistManager* man = G4NistManager::Instance();
	G4Element* Nit = man->FindOrBuildElement("N");
	G4Element* Oxy = man->FindOrBuildElement("O");
	G4Element* Alu = man->FindOrBuildElement("Al");

	
	G4int ncomp = 2;
	// total nucleon of Al(NO3)3 number is 294
	G4double density = 0.52068 *g/cm3;
	G4Material* Nitrate = new G4Material("Nitrate", density, ncomp);
	Nitrate->AddElement(Nit, 1);
	Nitrate->AddElement(Oxy, 3);
	G4double abundance, fracMass;
	
	density = 1.72 *g/cm3;
	G4Material* material = new G4Material("AlN", density, ncomp);
	G4int nAtoms;
	
	material->AddElement(Alu,fracMass=9.18*perCent);
	material->AddMaterial(Nitrate,fracMass=90.82*perCent);

	if(!material) {
		G4cout<<" ----> Material "<<"Aluminum"<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	//vis attributes
	G4VisAttributes* aluminumAnnulusBotVisAtt = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
	aluminumAnnulusBotVisAtt->SetVisibility(true);

	//main annulus
	G4double innerRadius = fAluminumAnnulusBotInnerRadius;
	G4double outerRadius = fAluminumAnnulusBotOuterRadius;
	G4double halfLengthZ = fAluminumAnnulusBotThickness/2.0;
	//screw cuts
	G4double cutInnerRadius = 0.*mm;
	G4double cutOuterRadius = fScrewRadius;
	G4double cutHalfLengthZ = 20.*mm;

	//primitive volume
	G4Tubs* alAnBoWithoutCuts = new G4Tubs("alAnBoWithoutCuts", innerRadius, outerRadius, halfLengthZ, 0*M_PI, 2*M_PI);
	G4Tubs* cutScrewHole = new G4Tubs("cutScrewHole", cutInnerRadius, cutOuterRadius, cutHalfLengthZ, 0*M_PI, 2*M_PI);

	//cut out screw holes
	G4double cI;
	G4double cR = fScrewPlacementRadius;
	G4double cDt = 45.*deg;
	G4ThreeVector moveCut;
//	cI = 0; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
//	G4SubtractionSolid* alAnBo0 = new G4SubtractionSolid("alAnBo0", alAnBoWithoutCuts, cutScrewHole, 0, moveCut);
//	cI = 1; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
//	G4SubtractionSolid* alAnBo1 = new G4SubtractionSolid("alAnBo1", alAnBo0, cutScrewHole, 0, moveCut);
//	cI = 2; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
//	G4SubtractionSolid* alAnBo2 = new G4SubtractionSolid("alAnBo2", alAnBo1, cutScrewHole, 0, moveCut);
//	cI = 3; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
//	G4SubtractionSolid* alAnBo3 = new G4SubtractionSolid("alAnBo3", alAnBo2, cutScrewHole, 0, moveCut);
//	cI = 4; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
//	G4SubtractionSolid* alAnBo4 = new G4SubtractionSolid("alAnBo4", alAnBo3, cutScrewHole, 0, moveCut);
//	cI = 5; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
//	G4SubtractionSolid* alAnBo5 = new G4SubtractionSolid("alAnBo5", alAnBo4, cutScrewHole, 0, moveCut);
//	cI = 6; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
//	G4SubtractionSolid* alAnBo6 = new G4SubtractionSolid("alAnBo6", alAnBo5, cutScrewHole, 0, moveCut);
	//final cut
//	cI = 7; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
//	G4SubtractionSolid* aluminumAnnulusBot = new G4SubtractionSolid("aluminumAnnulusBot", alAnBo6, cutScrewHole, 0, moveCut);

	//logical volume
	if(fAluminumAnnulusBotLog == nullptr) {
		fAluminumAnnulusBotLog = new G4LogicalVolume(alAnBoWithoutCuts, material, "aluminumAnnulusBotLog", 0, 0, 0);
		fAluminumAnnulusBotLog->SetVisAttributes(aluminumAnnulusBotVisAtt);
	}

	return 1;
}//end ::end AddAluminumAnnulusBot


G4int DetectionSystemPaces::ConstructDelrinShell() {
	G4Material* materialDelrin = G4Material::GetMaterial("Delrin");
	if(!materialDelrin) {
		G4cout<<" ----> Material "<<"Delrin"<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	G4VisAttributes* DelrinVisAtt = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
	DelrinVisAtt->SetVisibility(true);

	G4SubtractionSolid* shell = DelrinShell();

	fDelrinShellLog = new G4LogicalVolume(shell, materialDelrin, "DelrinShellLog", 0, 0, 0);
	fDelrinShellLog->SetVisAttributes(DelrinVisAtt);

	G4RotationMatrix* rotateNull = new G4RotationMatrix;
	G4ThreeVector moveNull(100.0,0.0,0.0);

	fAssembly->AddPlacedVolume(fDelrinShellLog, moveNull, rotateNull);

	return 1;
}




G4SubtractionSolid* DetectionSystemPaces::DelrinShell() {
	G4Sphere* sphere = new G4Sphere("sphere", fDelrinInnerRadius, fDelrinOuterRadius, 0, 2.0*M_PI, 0, M_PI);

	G4Tubs* chopTub = new G4Tubs("chopTub", 0, fDelrinHoleRadius, fDelrinOuterRadius +1.0, 0, 2.0*M_PI);

	G4RotationMatrix* rotateChopTub = new G4RotationMatrix;
	rotateChopTub->rotateX(-M_PI/2.0);

	G4ThreeVector moveNull(0.0,0.0,0.0);

	G4SubtractionSolid* delrinShell = new G4SubtractionSolid("DelrinShell", sphere, chopTub, rotateChopTub, moveNull);

	return delrinShell;
} //end ::DelrinShell




















//Add the aluminum hemisphere, housing for PACES
G4int DetectionSystemPaces::AddAluminumHemisphere() {
	//material
	// G4Material* material = G4Material::GetMaterial("Aluminum");
	G4double z;
	G4double a = 14.00*g/mole;
	G4NistManager* man = G4NistManager::Instance();
	G4Element* Nit = man->FindOrBuildElement("N");
	G4Element* Oxy = man->FindOrBuildElement("O");
	G4Element* Alu = man->FindOrBuildElement("Al");

	
	G4int ncomp = 2;
	// total nucleon of Al(NO3)3 number is 294
	G4double density = 0.52068 *g/cm3;
	G4Material* Nitrate = new G4Material("Nitrate", density, ncomp);
	Nitrate->AddElement(Nit, 1);
	Nitrate->AddElement(Oxy, 3);
	G4double abundance, fracMass;
	
	density = 1.72 *g/cm3;
	G4Material* material = new G4Material("AlN", density, ncomp);
	G4int nAtoms;
	
	material->AddElement(Alu,fracMass=9.18*perCent);
	material->AddMaterial(Nitrate,fracMass=90.82*perCent);

	if(!material) {
		G4cout<<" ----> Material "<<"Aluminum"<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	//vis attributes
	G4VisAttributes* aluminumHemisphereVisAtt = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
	aluminumHemisphereVisAtt->SetVisibility(true);

	//main hemisphere shell
	G4double innerRadius = fAluminumHemisphereInnerRadius;
	G4double outerRadius = fAluminumHemisphereOuterRadius;
	G4double thickness = std::abs(outerRadius - innerRadius);
	G4double startPhi = M_PI/2;
	G4double endPhi = M_PI - startPhi;
	//beam hole rim
	G4double rimInnerRadius = std::abs(innerRadius - fAluminumHemisphereBeamHoleRimHeight);
	G4double rimOuterRadius = innerRadius;
	G4double rimEndPhi = asin((fAluminumHemisphereBeamHoleRadius + thickness) / innerRadius);
	G4double rimStartPhi = M_PI - rimEndPhi;
	//beam hole cut
	G4double beamHoleInnerRadius = 0.*mm;
	G4double beamHoleOuterRadius = fAluminumHemisphereBeamHoleRadius;
	G4double beamHoleHalfLengthZ = outerRadius;
	//cuts where detectors go
	G4double detectorInnerRadius = 0.*mm;
	G4double detectorOuterRadius = fDetectorAssemblyRadius + fCutClearance;
	G4double detectorHalfLengthZ = fDetectorAssemblyThickness/2.0;

	//primitive volumes
	G4Sphere* alHeShell = new G4Sphere("alHeShell", innerRadius, outerRadius, 0*M_PI, 2*M_PI, startPhi, endPhi);
	G4Sphere* alHeRim = new G4Sphere("alHeRim", rimInnerRadius, rimOuterRadius, 0*M_PI, 2*M_PI, rimStartPhi, rimEndPhi);
	G4Tubs* cutBeamHole = new G4Tubs("cutBeamHole", beamHoleInnerRadius, beamHoleOuterRadius, beamHoleHalfLengthZ, 0*M_PI, 2*M_PI);
	G4Tubs* cutDetector = new G4Tubs("cutDetector", detectorInnerRadius, detectorOuterRadius, detectorHalfLengthZ, 0*M_PI, 2*M_PI);

	//add rim and beam hole
	G4RotationMatrix* R0 = new G4RotationMatrix; G4ThreeVector T0;
	G4UnionSolid* alHeShellRim = new G4UnionSolid("alHeShellRim", alHeShell, alHeRim, R0, T0);
	G4SubtractionSolid* alHeShellRimHole = new G4SubtractionSolid("alHeShellRimHole", alHeShellRim, cutBeamHole, R0, T0);
	//cut out detectors

	G4double* ptrPd = fPacesPlacementDistance;
	G4double* ptrPt = fPacesPlacementPhi;
	G4double* ptrPp = fPacesPlacementTheta;
	G4double* ptrOt = fPacesOrientationPhi;
	G4double* ptrOp = fPacesOrientationTheta;

	//  G4RotationMatrix* rotateCut[5];
	//  G4ThreeVector moveCut[5], yprimeaxis;
	//  G4double dDist, dPhi, dTheta, oriPhi, oriTheta;
	//
	//  for(int i=0; i<5; i++)
	//  {
	//    dDist = ptrPd[i] + detectorHalfLengthZ;
	//    dPhi = ptrPt[i];
	//    dTheta = ptrPp[i];
	//    rotateCut[i] = new G4RotationMatrix;
	//    moveCut[i].setX(dDist * cos(dPhi) * sin(dTheta));
	//    moveCut[i].setY(dDist * sin(dPhi) * sin(dTheta));
	//    moveCut[i].setZ(dDist *    1.0     * cos(dTheta));
	//    oriPhi = dPhi + ptrOt[i] - M_PI/2; //minus 90 deg
	//    oriTheta = dTheta + ptrOp[i];
	////    yprimeaxis = G4ThreeVector(cos(oriPhi), sin(oriPhi), 0);
	//	  yprimeaxis.set(cos(oriPhi), sin(oriPhi), 0);
	//    rotateCut[i]->set(yprimeaxis, oriTheta);
	//  }
	//  G4SubtractionSolid* alHe0 = new G4SubtractionSolid("alHe0", alHeShellRimHole, cutDetector, rotateCut[0], moveCut[0]);
	//  G4SubtractionSolid* alHe1 = new G4SubtractionSolid("alHe1", alHe0, cutDetector, rotateCut[1], moveCut[1]);
	//  G4SubtractionSolid* alHe2 = new G4SubtractionSolid("alHe2", alHe1, cutDetector, rotateCut[2], moveCut[2]);
	//  G4SubtractionSolid* alHe3 = new G4SubtractionSolid("alHe3", alHe2, cutDetector, rotateCut[3], moveCut[3]);
	//  G4SubtractionSolid* aluminumHemisphere = new G4SubtractionSolid("aluminumHemisphere", alHe3, cutDetector, rotateCut[4], moveCut[4]);
	//  for(int i=0; i<5; i++) delete rotateCut[i]; //safety

	G4ThreeVector moveCut, yprimeaxis;
	G4RotationMatrix* rotateCut = new G4RotationMatrix;
	G4double dDist, dPhi, dTheta, oriPhi, oriTheta;
	G4int dI;
	//one
	dI = 0;
	dDist = ptrPd[dI] + detectorHalfLengthZ;
	dPhi = ptrPt[dI];
	dTheta = ptrPp[dI];
	moveCut.setX(dDist * cos(dPhi) * sin(dTheta));
	moveCut.setY(dDist * sin(dPhi) * sin(dTheta));
	moveCut.setZ(dDist *      1.0     * cos(dTheta));
	oriPhi = dPhi + ptrOt[dI] - M_PI/2;
	oriTheta = dTheta + ptrOp[dI];
	yprimeaxis.set(cos(oriPhi), sin(oriPhi), 0);
	rotateCut->set(yprimeaxis, oriTheta);
	G4SubtractionSolid* alHe0 = new G4SubtractionSolid("alHe0", alHeShellRimHole, cutDetector, rotateCut, moveCut);
	//two
	dI = 1;
	dDist = ptrPd[dI] + detectorHalfLengthZ;
	dPhi = ptrPt[dI];
	dTheta = ptrPp[dI];
	moveCut.setX(dDist * cos(dPhi) * sin(dTheta));
	moveCut.setY(dDist * sin(dPhi) * sin(dTheta));
	moveCut.setZ(dDist *      1.0     * cos(dTheta));
	oriPhi = dPhi + ptrOt[dI] - M_PI/2;
	oriTheta = dTheta + ptrOp[dI];
	yprimeaxis.set(cos(oriPhi), sin(oriPhi), 0);
	rotateCut->set(yprimeaxis, oriTheta);
	G4SubtractionSolid* alHe1 = new G4SubtractionSolid("alHe1", alHe0, cutDetector, rotateCut, moveCut);
	//three
	dI = 2;
	dDist = ptrPd[dI] + detectorHalfLengthZ;
	dPhi = ptrPt[dI];
	dTheta = ptrPp[dI];
	moveCut.setX(dDist * cos(dPhi) * sin(dTheta));
	moveCut.setY(dDist * sin(dPhi) * sin(dTheta));
	moveCut.setZ(dDist *      1.0     * cos(dTheta));
	oriPhi = dPhi + ptrOt[dI] - M_PI/2;
	oriTheta = dTheta + ptrOp[dI];
	yprimeaxis.set(cos(oriPhi), sin(oriPhi), 0);
	rotateCut->set(yprimeaxis, oriTheta);
	G4SubtractionSolid* alHe2 = new G4SubtractionSolid("alHe2", alHe1, cutDetector, rotateCut, moveCut);
	//four
	dI = 3;
	dDist = ptrPd[dI] + detectorHalfLengthZ;
	dPhi = ptrPt[dI];
	dTheta = ptrPp[dI];
	moveCut.setX(dDist * cos(dPhi) * sin(dTheta));
	moveCut.setY(dDist * sin(dPhi) * sin(dTheta));
	moveCut.setZ(dDist *      1.0     * cos(dTheta));
	oriPhi = dPhi + ptrOt[dI] - M_PI/2;
	oriTheta = dTheta + ptrOp[dI];
	yprimeaxis.set(cos(oriPhi), sin(oriPhi), 0);
	rotateCut->set(yprimeaxis, oriTheta);
	G4SubtractionSolid* alHe3 = new G4SubtractionSolid("alHe3", alHe2, cutDetector, rotateCut, moveCut);
	//five
	dI = 4;
	dDist = ptrPd[dI] + detectorHalfLengthZ;
	dPhi = ptrPt[dI];
	dTheta = ptrPp[dI];
	moveCut.setX(dDist * cos(dPhi) * sin(dTheta));
	moveCut.setY(dDist * sin(dPhi) * sin(dTheta));
	moveCut.setZ(dDist *      1.0     * cos(dTheta));
	oriPhi = dPhi + ptrOt[dI] - M_PI/2;
	oriTheta = dTheta + ptrOp[dI];
	yprimeaxis.set(cos(oriPhi), sin(oriPhi), 0);
	rotateCut->set(yprimeaxis, oriTheta);
	G4SubtractionSolid* aluminumHemisphere = new G4SubtractionSolid("aluminumHemisphere", alHe3, cutDetector, rotateCut, moveCut);
	//  //end

	//logical volume
	if(fAluminumHemisphereLog == nullptr) {
		fAluminumHemisphereLog = new G4LogicalVolume(aluminumHemisphere, material, "aluminumHemisphereLog", 0, 0, 0);
		fAluminumHemisphereLog->SetVisAttributes(aluminumHemisphereVisAtt);
	}

	return 1;
}//end ::end AddAluminumHemisphere


//Add screws
G4int DetectionSystemPaces::AddScrews() {
	return 1;
}//end ::end AddScrews
