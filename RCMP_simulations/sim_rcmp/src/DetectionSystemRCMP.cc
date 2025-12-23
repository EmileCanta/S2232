#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"


#include "G4Material.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemRCMP.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units

#include "CADMesh.hh"

//RCMP Detection system constructor
DetectionSystemRCMP::DetectionSystemRCMP() :
	// LogicalVolumes
	fDSSSDpixelLog(0)
{ 	
	//BB7 DSSSD x,y length and thickness
	fXLength = 64.*mm;
	fYLength = 64.*mm;
	fDetectorThickness   = 1.*mm;
	fDetectorDistance1 = 41.93*mm;
	fDetectorDistance2 = 42.85*mm;
	fDetectorDistance3 = 48.71*mm;

	fXPosOffset = 0.;
	fYPosOffset = 0.;
	fZPosOffset = 0.;


	//number of channels per row and cols 32 and 32
	fPixelsXRow = (G4int)(32);
	fPixelsYRow = (G4int)(32);

	//total number of pixels 32*32
	fNumberOfPixels = fPixelsXRow*fPixelsYRow;
	
	
	
	//single pixel width
	fPixelWidth      = 2.*mm;
	
	
	
	



	
}


//RCMP detection system destructor

DetectionSystemRCMP::~DetectionSystemRCMP() {
	// LogicalVolumes
	delete fDSSSDpixelLog;
}



//Detector assembly
G4int DetectionSystemRCMP::Build() { 

	// Build assembly volume
	fAssembly = new G4AssemblyVolume(); 

    G4cout<<"BuildPixelVolume"<<G4endl;
	BuildPixelVolume();

	return 1;
}



//FUNCTION TO PLACE DETECTOR PIXEL INTO PHSYICAL VOLUME STATUS:::COMPLETE
G4int DetectionSystemRCMP::PlaceDetector(G4LogicalVolume* expHallLog) {

	//building detector assuming channel 1,1 starts at (-64,-64) This may change depends on detector rotations
	
	
	G4double posX;
	G4double posY;
	G4double posZ;

	G4int pixelNumber = 0;

	G4RotationMatrix* GriffinRotate = new G4RotationMatrix();
	G4RotationMatrix* test = new G4RotationMatrix();
	//GriffinRotate->rotateZ(-22.5*deg);   // orientation
	//GriffinRotate->rotateZ(-0.*deg);   // orientation
	GriffinRotate->rotateY(37.53*deg);   // orientation
	test->rotateY(37.4*deg);   // orientation

	// Using geant4 transform class to handle GRIFFIN rotation
	G4Transform3D rotateTransform(*GriffinRotate, G4ThreeVector()); // no translation
	G4Transform3D rotateTransformTest(*test, G4ThreeVector()); // no translation
	
	//G4ThreeVector FrameOffset(23.126, -41.55, -12.901);
	G4ThreeVector FrameOffset(0., 0., 0.);
	G4RotationMatrix* FrameRotate = new G4RotationMatrix();
	//FrameRotate->rotateY(225*deg);   // orientation

G4double startY = -fXLength / 2.0;   // rowX -> Y
G4double startZ = -fYLength / 2.0;   // rowY -> Z
G4double startX = -fYLength / 2.0;












	
	//32*32 pixel construction with row1 col(1-32) 0,1,2,3>>>>>>>31, row2 col(1-32)32>>>>>>>>63
	//First loop for Junction side in Ntuple with strip number increase along Y direction
	//Nested loop for Ohmic side in Ntuple with strip number increase along Z direction
	//construc detector in x-direction, x-distance is fixed, 
	for (G4int rowX = 0; rowX < fPixelsXRow; ++rowX) {
	
		
    		G4double posY = startX + fPixelWidth*rowX + fPixelWidth/2.0 + fXPosOffset;
		//posX = fDetectorDistance1;
    		for (G4int rowY = 0; rowY < fPixelsYRow; ++rowY) {
        		G4double posZ = startY + fPixelWidth*rowY + fPixelWidth/2.0 - 15*mm;
        		G4double posX = fDetectorDistance1 + fZPosOffset;
		//non-rotated coordinates using detector local positions
        	G4ThreeVector localPos(posX, posY, posZ);
		//construct detector in negative X at the same time
		G4ThreeVector localPos2(-posX,posY,posZ);
        // rotate the *position* about the global origin
        
	G4Point3D p(localPos);
	G4Point3D p2(localPos2);
    	G4Point3D pRot = rotateTransformTest * p;
	G4Point3D pRot2 = rotateTransformTest *p2; 
    	G4ThreeVector rotatedPos(pRot);
	G4ThreeVector rotatedPos2(pRot2);
        
        
        //Assiging this detector to RCMP1
        //fAssembly->MakeImprint(expHallLog, rotatedPos, test, pixelNumber);
        fAssembly->MakeImprint(expHallLog, rotatedPos, test, pixelNumber);

	++pixelNumber;
    }
}




	for (G4int rowX = 0; rowX < fPixelsXRow; ++rowX) {
	
		
    		G4double posY = startX + fPixelWidth*rowX + fPixelWidth/2.0 + fXPosOffset;
		//posX = fDetectorDistance1;
    		for (G4int rowY = 0; rowY < fPixelsYRow; ++rowY) {
        		G4double posZ = startY + fPixelWidth*rowY + fPixelWidth/2.0 + fYPosOffset;
        		G4double posX = fDetectorDistance1 + fZPosOffset;
		//non-rotated coordinates using detector local positions
        	G4ThreeVector localPos(posX, posY, posZ);
		//construct detector in negative X at the same time
		G4ThreeVector localPos2(-posX,posY,posZ);
        // rotate the *position* about the global origin
        
	G4Point3D p(localPos);
	G4Point3D p2(localPos2);
    	G4Point3D pRot = rotateTransform * p;
	G4Point3D pRot2 = rotateTransform *p2; 
    	G4ThreeVector rotatedPos(pRot);
	G4ThreeVector rotatedPos2(pRot2);
        
        
        //Assiging this detector to RCMP1
        //fAssembly->MakeImprint(expHallLog, rotatedPos, GriffinRotate, pixelNumber);

        //Assiging this detector to RCMP2
	//fAssembly->MakeImprint(expHallLog,rotatedPos2,GriffinRotate,pixelNumber+32*32);
	++pixelNumber;
    }
}












	//construction of inner detectors along y-axis 
pixelNumber = 32*32*2;
	GriffinRotate->rotateZ(74*deg);
	// Using geant4 transform class to handle GRIFFIN rotation
	//G4Transform3D rotateTransform2(*GriffinRotate, G4ThreeVector()); // no translation
	//First loop for Junction side in Ntuple with strip number increase along X direction
	//Nested loop for Ohmic side in Ntuple with strip number increase along Z direction
	//construc detector in y-direction, distance y is fixed
	for (G4int rowX = 0; rowX < fPixelsXRow; ++rowX) {
    		G4double posX = startX + fPixelWidth*rowX + fPixelWidth/2.0 + fXPosOffset;
		//posX = fDetectorDistance1;
    		for (G4int rowY = 0; rowY < fPixelsYRow; ++rowY) {
        		G4double posZ = startY + fPixelWidth*rowY + fPixelWidth/2.0 + fYPosOffset;
        		G4double posY = fDetectorDistance1 + fZPosOffset;
		//non-rotated coordinates using detector local positions
        	G4ThreeVector localPos3(posX, posY, posZ);
		//construct detector in negative X at the same time
		G4ThreeVector localPos4(posX,-posY,posZ);
        // rotate the *position* about the global origin
        
	G4Point3D p3(localPos3);
	G4Point3D p4(localPos4);
    	G4Point3D pRot3 = rotateTransform * p3;
	G4Point3D pRot4 = rotateTransform *p4; 
    	G4ThreeVector rotatedPos3(pRot3);
	G4ThreeVector rotatedPos4(pRot4);
        
        //fAssembly->MakeImprint(expHallLog, rotatedPos3, GriffinRotate, pixelNumber);
	//fAssembly->MakeImprint(expHallLog,rotatedPos4,GriffinRotate,pixelNumber+32*32);
	++pixelNumber;
    }
}



for (G4int rowX = 0; rowX < fPixelsXRow; ++rowX) {
    		G4double posX = startX + fPixelWidth*rowX + fPixelWidth/2.0 + fXPosOffset;
		//posX = fDetectorDistance1;
    		for (G4int rowY = 0; rowY < fPixelsYRow; ++rowY) {
        		G4double posZ = startY + fPixelWidth*rowY + fPixelWidth/2.0 + fYPosOffset;
        		G4double posY = fDetectorDistance1 + fZPosOffset;
		//non-rotated coordinates using detector local positions
        	G4ThreeVector localPos3(posX, posY, posZ);
		//construct detector in negative X at the same time
		G4ThreeVector localPos4(posX,-posY,posZ);
        // rotate the *position* about the global origin
        
	G4Point3D p3(localPos3);
	G4Point3D p4(localPos4);
    	G4Point3D pRot3 = rotateTransform * p3;
	G4Point3D pRot4 = rotateTransform *p4; 
    	G4ThreeVector rotatedPos3(pRot3);
	G4ThreeVector rotatedPos4(pRot4);
        
        //fAssembly->MakeImprint(expHallLog, rotatedPos3, GriffinRotate, pixelNumber);
	//fAssembly->MakeImprint(expHallLog,rotatedPos4,GriffinRotate,pixelNumber+32*32);
	++pixelNumber;
    }
}

 	//construction of inner detectors along y-axis 
    pixelNumber = 32*32*3;
    GriffinRotate->rotateX(90*deg);
// Using geant4 transform class to handle GRIFFIN rotation
//G4Transform3D rotateTransform2(*GriffinRotate, G4ThreeVector()); // no translation
//First loop for Junction side in Ntuple with strip number increase along X direction
//Nested loop for Ohmic side in Ntuple with strip number increase along Z direction
//construc detector in y-direction, distance y is fixed
for (G4int rowX = 0; rowX < fPixelsXRow; ++rowX) {
G4double posX = startX + fPixelWidth*rowX + fPixelWidth/2.0 + fXPosOffset;
//posX = fDetectorDistance1;
for (G4int rowY = 0; rowY < fPixelsYRow; ++rowY) {
G4double posY = startY + fPixelWidth*rowY + fPixelWidth/2.0 + fYPosOffset;
G4double posZ = fDetectorDistance1 + fZPosOffset;
//non-rotated coordinates using detector local positions
G4ThreeVector localPos5(posX, posY, posZ);
//construct detector in negative X at the same time
G4ThreeVector localPos6(posX,posY,-posZ);
// rotate the *position* about the global origin

G4Point3D p5(localPos5);
G4Point3D p6(localPos6);
G4Point3D pRot5 = rotateTransform * p5;
G4Point3D pRot6 = rotateTransform *p6; 
G4ThreeVector rotatedPos5(pRot5);
G4ThreeVector rotatedPos6(pRot6);

//fAssembly->MakeImprint(expHallLog, rotatedPos5, GriffinRotate, pixelNumber);
++pixelNumber;
}
}



for (G4int rowX = 0; rowX < fPixelsXRow; ++rowX) {
G4double posX = startX + fPixelWidth*rowX + fPixelWidth/2.0 + fXPosOffset;
//posX = fDetectorDistance1;
for (G4int rowY = 0; rowY < fPixelsYRow; ++rowY) {
G4double posY = startY + fPixelWidth*rowY + fPixelWidth/2.0 + fYPosOffset;
G4double posZ = fDetectorDistance1 + fZPosOffset;
//non-rotated coordinates using detector local positions
G4ThreeVector localPos5(posX, posY, posZ);
//construct detector in negative X at the same time
G4ThreeVector localPos6(posX,posY,-posZ);
// rotate the *position* about the global origin

G4Point3D p5(localPos5);
G4Point3D p6(localPos6);
G4Point3D pRot5 = rotateTransform * p5;
G4Point3D pRot6 = rotateTransform *p6; 
G4ThreeVector rotatedPos5(pRot5);
G4ThreeVector rotatedPos6(pRot6);

//fAssembly->MakeImprint(expHallLog, rotatedPos3, GriffinRotate, pixelNumber);
//fAssembly->MakeImprint(expHallLog,rotatedPos6,GriffinRotate,pixelNumber+32*32);
++pixelNumber;
}
}
    G4Material* FrameMaterial = G4Material::GetMaterial("G4_POLYETHYLENE");
	
    if(!FrameMaterial) 
    {
		G4cout << " ----> Material " << FrameMaterial << G4endl;
		return 0;
	}

	G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
	
    visAtt->SetVisibility(true);
    
    auto mesh = CADMesh::TessellatedMesh::FromSTL("../../frame_3d/frame.stl");

    auto solid = mesh->GetSolid();

    G4LogicalVolume* MeshLog = new G4LogicalVolume(solid, FrameMaterial, "MeshLog");

    //G4VPhysicalVolume* MeshPhys = new G4PVPlacement(0, FrameOffset, MeshLog, "MeshPhys", expHallLog, false, 0, true); 
    
    auto windows = CADMesh::TessellatedMesh::FromSTL("../../frame_3d/windows.stl");

    auto solidWind = windows->GetSolid();

    G4LogicalVolume* WindLog = new G4LogicalVolume(solidWind, FrameMaterial, "WindLog");

    G4VPhysicalVolume* WindPhys = new G4PVPlacement(0, FrameOffset, WindLog, "WindPhys", expHallLog, false, 0, true); 

	return 1;
}



//This function used for creating pixel volume:::::Material needs update

G4int DetectionSystemRCMP::BuildPixelVolume() {

	//need modification of silicon material
	G4Material* material = G4Material::GetMaterial("Silicon");
	if(!material) {
		G4cout<<" ----> Material "<<material<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
	visAtt->SetVisibility(true);


	G4Box* DSSSDpixel = BuildPixel();

	// Define rotation and movement objects, logic volume is defined at ceneter hence no rotation or movement
	G4ThreeVector direction 	= G4ThreeVector(0,0,0);
	G4double zPosition		= 0.0*mm;
	G4ThreeVector move 		= zPosition * direction;
	G4RotationMatrix* rotate  = new G4RotationMatrix;
	
		
	//Start by building one detector at distance of 47.6mm in z-direction
	//logical volume of detector pixel
	if(fDSSSDpixelLog == nullptr) {
		fDSSSDpixelLog = new G4LogicalVolume(DSSSDpixel, material, "DSSSDpixelLog", 0, 0, 0);
		fDSSSDpixelLog->SetVisAttributes(visAtt);

	}

	fAssembly->AddPlacedVolume(fDSSSDpixelLog, move, rotate);

	return 1;
}

///////////////////////////////////////////////////////////////////////
// BB7 Detector pixel construction 2mm*2mm*1mm :::::::FUNCTION COMPLETE
///////////////////////////////////////////////////////////////////////
G4Box* DetectionSystemRCMP::BuildPixel() {
	G4double halfLengthX = fPixelWidth/2.0;
	G4double halfLengthY = fPixelWidth/2.0;
	G4double halfLengthZ = fDetectorThickness/2.0;
	//Create BB7 DSSSD single pixel
	//logic volume along x axis
	G4Box* pixel = new G4Box("pixel", halfLengthZ, halfLengthY,halfLengthX );

	return pixel;
}//end ::BuildPixel
