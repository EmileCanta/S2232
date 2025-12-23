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

#include "DetectionSystemQED.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units

//QED Detection system constructor
DetectionSystemQED::DetectionSystemQED() :
	// LogicalVolumes
	fDSSSDpixelLog(0)
{ 	
	//BB7 DSSSD x,y length and thickness
	fXLength = 64.*mm;
	fYLength = 64.*mm;
	fDetectorThickness   = 1.*mm;
	fDetectorDistance = 47.6*mm;
	fDetectorOuterDistance =  57.2;

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


//QED detection system destructor

DetectionSystemQED::~DetectionSystemQED() {
	// LogicalVolumes
	delete fDSSSDpixelLog;
}



//Detector assembly
G4int DetectionSystemQED::Build() { 

	// Build assembly volume
	fAssembly = new G4AssemblyVolume(); 

	G4cout<<"BuildPixelVolume"<<G4endl;
	BuildPixelVolume();

	return 1;
}



//FUNCTION TO PLACE DETECTOR PIXEL INTO PHSYICAL VOLUME STATUS:::COMPLETE
G4int DetectionSystemQED::PlaceDetector(G4LogicalVolume* expHallLog) {

	//building detector assuming channel 1,1 starts at (-64,-64) This may change depends on detector rotations
	
	
	G4double posX;
	G4double posY;
	G4double posZ;

	G4int pixelNumber = 0;

	G4RotationMatrix* GriffinRotate = new G4RotationMatrix();
	GriffinRotate->rotateZ(-22.5*deg);   // orientation
	//GriffinRotate->rotateY(90*deg);
	// Using geant4 transform class to handle GRIFFIN rotation
	G4Transform3D rotateTransform(*GriffinRotate, G4ThreeVector()); // no translation
	
	

G4double startY = -fXLength / 2.0;   // rowX -> Y
G4double startZ = -fYLength / 2.0;   // rowY -> Z
G4double startX = -fYLength / 2.0;












	
	//32*32 pixel construction with row1 col(1-32) 0,1,2,3>>>>>>>31, row2 col(1-32)32>>>>>>>>63
	//First loop for Junction side in Ntuple with strip number increase along Y direction
	//Nested loop for Ohmic side in Ntuple with strip number increase along Z direction
	//construc detector in x-direction, x-distance is fixed, 
	for (G4int rowX = 0; rowX < fPixelsXRow; ++rowX) {
	
		
    		G4double posY = startX + fPixelWidth*rowX + fPixelWidth/2.0 + fXPosOffset;
		//posX = fDetectorDistance;
    		for (G4int rowY = 0; rowY < fPixelsYRow; ++rowY) {
        		G4double posZ = startY + fPixelWidth*rowY + fPixelWidth/2.0 + fYPosOffset;
        		G4double posX = fDetectorDistance + fZPosOffset;
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
        
        
        //Assiging this detector to QED1
        fAssembly->MakeImprint(expHallLog, rotatedPos, GriffinRotate, pixelNumber);

        //Assiging this detector to QED2
	//fAssembly->MakeImprint(expHallLog,rotatedPos2,GriffinRotate,pixelNumber+32*32);
	++pixelNumber;
    }
}




	for (G4int rowX = 0; rowX < fPixelsXRow; ++rowX) {
	
		
    		G4double posY = startX + fPixelWidth*rowX + fPixelWidth/2.0 + fXPosOffset;
		//posX = fDetectorDistance;
    		for (G4int rowY = 0; rowY < fPixelsYRow; ++rowY) {
        		G4double posZ = startY + fPixelWidth*rowY + fPixelWidth/2.0 + fYPosOffset;
        		G4double posX = fDetectorDistance + fZPosOffset;
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
        
        
        //Assiging this detector to QED1
        //fAssembly->MakeImprint(expHallLog, rotatedPos, GriffinRotate, pixelNumber);

        //Assiging this detector to QED2
	fAssembly->MakeImprint(expHallLog,rotatedPos2,GriffinRotate,pixelNumber+32*32);
	++pixelNumber;
    }
}












	//construction of inner detectors along y-axis 
pixelNumber = 32*32*2;
	GriffinRotate->rotateZ(90*deg);
	// Using geant4 transform class to handle GRIFFIN rotation
	//G4Transform3D rotateTransform2(*GriffinRotate, G4ThreeVector()); // no translation
	//First loop for Junction side in Ntuple with strip number increase along X direction
	//Nested loop for Ohmic side in Ntuple with strip number increase along Z direction
	//construc detector in y-direction, distance y is fixed
	for (G4int rowX = 0; rowX < fPixelsXRow; ++rowX) {
    		G4double posX = startX + fPixelWidth*rowX + fPixelWidth/2.0 + fXPosOffset;
		//posX = fDetectorDistance;
    		for (G4int rowY = 0; rowY < fPixelsYRow; ++rowY) {
        		G4double posZ = startY + fPixelWidth*rowY + fPixelWidth/2.0 + fYPosOffset;
        		G4double posY = fDetectorDistance + fZPosOffset;
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
        
        fAssembly->MakeImprint(expHallLog, rotatedPos3, GriffinRotate, pixelNumber);
	//fAssembly->MakeImprint(expHallLog,rotatedPos4,GriffinRotate,pixelNumber+32*32);
	++pixelNumber;
    }
}



for (G4int rowX = 0; rowX < fPixelsXRow; ++rowX) {
    		G4double posX = startX + fPixelWidth*rowX + fPixelWidth/2.0 + fXPosOffset;
		//posX = fDetectorDistance;
    		for (G4int rowY = 0; rowY < fPixelsYRow; ++rowY) {
        		G4double posZ = startY + fPixelWidth*rowY + fPixelWidth/2.0 + fYPosOffset;
        		G4double posY = fDetectorDistance + fZPosOffset;
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
	fAssembly->MakeImprint(expHallLog,rotatedPos4,GriffinRotate,pixelNumber+32*32);
	++pixelNumber;
    }
}






	//construction of outer layered detectors along y-axis



	pixelNumber = 32*32*4;
	//GriffinRotate->rotateY(90*deg);
	// Using geant4 transform class to handle GRIFFIN rotation
	//G4Transform3D rotateTransform2(*GriffinRotate, G4ThreeVector()); // no translation
	//32*32 pixel construction with row1 col(1-32) 0,1,2,3>>>>>>>31, row2 col(1-32)32>>>>>>>>63
	//hence all hit in row 1 is adding up copy number 0 to 31
	//all hit in col 1 is adding up 0,31,63,..............
	//construc detector in y-direction 
	for (G4int rowX = 0; rowX < fPixelsXRow; ++rowX) {
    		G4double posX = startX + fPixelWidth*rowX + fPixelWidth/2.0 + fXPosOffset;
		//posX = fDetectorDistance;
    		for (G4int rowY = 0; rowY < fPixelsYRow; ++rowY) {
        		G4double posZ = startY + fPixelWidth*rowY + fPixelWidth/2.0 + fYPosOffset;
        		G4double posY = fDetectorOuterDistance + fZPosOffset;
		//non-rotated coordinates using detector local positions
        	G4ThreeVector localPos5(posX, posY, posZ);
		//construct detector in negative X at the same time
		G4ThreeVector localPos6(posX,-posY,posZ);
        // rotate the *position* about the global origin
        
	G4Point3D p5(localPos5);
	G4Point3D p6(localPos6);
    	G4Point3D pRot5 = rotateTransform * p5;
	G4Point3D pRot6 = rotateTransform *p6; 
    	G4ThreeVector rotatedPos5(pRot5);
	G4ThreeVector rotatedPos6(pRot6);
       
        fAssembly->MakeImprint(expHallLog, rotatedPos5, GriffinRotate, pixelNumber);
	//fAssembly->MakeImprint(expHallLog,rotatedPos6,GriffinRotate,pixelNumber+32*32);
	 ++pixelNumber;
    }
}

for (G4int rowX = 0; rowX < fPixelsXRow; ++rowX) {
    		G4double posX = startX + fPixelWidth*rowX + fPixelWidth/2.0 + fXPosOffset;
		//posX = fDetectorDistance;
    		for (G4int rowY = 0; rowY < fPixelsYRow; ++rowY) {
        		G4double posZ = startY + fPixelWidth*rowY + fPixelWidth/2.0 + fYPosOffset;
        		G4double posY = fDetectorOuterDistance + fZPosOffset;
		//non-rotated coordinates using detector local positions
        	G4ThreeVector localPos5(posX, posY, posZ);
		//construct detector in negative X at the same time
		G4ThreeVector localPos6(posX,-posY,posZ);
        // rotate the *position* about the global origin
        
	G4Point3D p5(localPos5);
	G4Point3D p6(localPos6);
    	G4Point3D pRot5 = rotateTransform * p5;
	G4Point3D pRot6 = rotateTransform *p6; 
    	G4ThreeVector rotatedPos5(pRot5);
	G4ThreeVector rotatedPos6(pRot6);
       
        //fAssembly->MakeImprint(expHallLog, rotatedPos5, GriffinRotate, pixelNumber);
	fAssembly->MakeImprint(expHallLog,rotatedPos6,GriffinRotate,pixelNumber+32*32);
	 ++pixelNumber;
    }
}







	return 1;
}



//This function used for creating pixel volume:::::Material needs update

G4int DetectionSystemQED::BuildPixelVolume() {

	//need modification of silicon material
	G4Material* material = G4Material::GetMaterial("Silicon");
	if(!material) {
		G4cout<<" ----> Material "<<fPixelMaterial<<" not found, cannot build the detector shell! "<<G4endl;
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
G4Box* DetectionSystemQED::BuildPixel() {
	G4double halfLengthX = fPixelWidth/2.0;
	G4double halfLengthY = fPixelWidth/2.0;
	G4double halfLengthZ = fDetectorThickness/2.0;
	//Create BB7 DSSSD single pixel
	//logic volume along x axis
	G4Box* pixel = new G4Box("pixel", halfLengthZ, halfLengthY,halfLengthX );

	return pixel;
}//end ::BuildPixel
