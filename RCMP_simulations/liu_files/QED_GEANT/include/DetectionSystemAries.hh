#ifndef DETECTIONSYSTEMSARIES_HH
#define DETECTIONSYSTEMSARIES_HH

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4AssemblyVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"
#include "G4Para.hh"
#include "G4OpticalPhysics.hh"
#include "G4OpticalSurface.hh"



class DetectionSystemAries 	
{
public:
    DetectionSystemAries();
    ~DetectionSystemAries();

    void Build() ;
    void PlaceDetector(G4LogicalVolume* expHallLog);
    void ConstructSquareTile();
    void ConstructTriangleTile();
    void ConstructRectangleTile1();
    void ConstructRectangleTile2();
    void ConstructTestTile();
    void PhysicsDef();
   // void ConstructDome();
    void ConstructHold();
    G4LogicalVolume* scoringvol;
    G4LogicalVolume* scinttrack;

private:
// numbers
    G4double scintrectanglelength;
    G4double scintrectanglewidth;
    G4double scinttrianglelength;
    G4double scintsquarelength;
    G4double SeparateHemispheres;
    G4double gapbtwntiles;
    G4double innerradius;
    G4double outerradius;
    G4double smallSSL;
    G4double smallSST;
    G4double theta;
    G4double phi;
    G4double alpha;
    G4double beta;
    G4double gamma;
    G4double scintthickness;
    G4double silverthickness;
    G4double tiny;
    G4double zero;	
    G4double zeroangle;
    G4double halfrightangle; 
    G4double taperanglex;
    G4double taperangley;
    G4double thetaangle;
    G4double expo;
    G4double longslength;
    G4double longthickness;
    G4double longrlength;
    G4double trianglelong;
    G4double squarelong;
    G4double longerslength;
    G4double longerrectangle; 
    G4double phirect;
    G4double thetarect;
    G4double sl;
    G4double scmx;
    G4double scmz;    
    G4double sprotrude;
    G4double realthick;
    G4double t1;
    G4double t2;
    G4double t3;
    G4double tcmz;
    G4double tt1;
    G4double tt2;
    G4double tl;
    G4double ts;
    G4int n;
    G4double hrec;
    G4double drec;
    G4double rec1;
    G4double rec2;
    G4double tria;
    G4double trih;
    G4double trid;
    G4double rcma;
    G4double rcmt;
    G4double rcmb;
    G4double reclong;
    G4double rc1;
    G4double rc2;
    G4double rec3;
    G4double rec4;
    G4double rec5;
    G4double tri1;
    G4double tri2;
    G4double tri3;
    G4double tri4;
    G4double tri5;
    G4double tri6;
    G4double sipmConductThick;
    G4double sipmActiveThick;
    G4double sipmActiveDim;
    G4double sipmDim1;
    G4double sipmDim2;
    G4double sipmThick;



// logical volume
    G4LogicalVolume* fExpHallLog; 
//materials 
    G4Material* scintmat;
    G4Material* coatingmat;
    G4Material* uhv;
    G4Material* si;
    G4Material* delrin;
    G4Material* mylar;
    G4Material* copper;
    




//Vis attributes
    G4VisAttributes* coatingvis;
    G4VisAttributes* scintvis;
    G4VisAttributes* uhvis;
    G4VisAttributes* sivis;
    G4VisAttributes* shroudvis;
    G4VisAttributes* activevis;
    G4VisAttributes* conductvis;





// assembly volumes
    G4AssemblyVolume* fSquareTile;
    G4AssemblyVolume* ftriangle;
    G4AssemblyVolume* frectangle1;
    G4AssemblyVolume* frectangle2;  
    G4AssemblyVolume* ftest; 
    G4AssemblyVolume* fdome;
    G4AssemblyVolume* sSipm;

 
// rotation matrices
   G4RotationMatrix* rotateNull;
   G4RotationMatrix* rshroud;
   G4RotationMatrix* ractive;
   G4RotationMatrix* rconduct;
   G4RotationMatrix* rotate;
   G4RotationMatrix* rotatetest;
   G4RotationMatrix* rotright;
   G4RotationMatrix* rotleft;
   G4RotationMatrix* r1;
   G4RotationMatrix* r2;
   G4RotationMatrix* r3;
   G4RotationMatrix* r4;
   G4RotationMatrix* r5;
   G4RotationMatrix* r6;
   G4RotationMatrix* r7;
   G4RotationMatrix* r8;
   G4RotationMatrix* r9;
   G4RotationMatrix* r10;
   G4RotationMatrix* r11;
   G4RotationMatrix* r12;
   G4RotationMatrix* r13;
   G4RotationMatrix* r14;
   G4RotationMatrix* r15;
   G4RotationMatrix* r16;
   G4RotationMatrix* r17;
   G4RotationMatrix* r18;
   G4RotationMatrix* r19;
   G4RotationMatrix* r20;
   G4RotationMatrix* r21;
   G4RotationMatrix* r22;
   G4RotationMatrix* r23;
   G4RotationMatrix* r24;
   G4RotationMatrix* r25;
   G4RotationMatrix* r26;
   G4RotationMatrix* r27;
   G4RotationMatrix* r28;
   G4RotationMatrix* r29;
   G4RotationMatrix* r30;
   G4RotationMatrix* r31;
   G4RotationMatrix* r32;
   G4RotationMatrix* r33;
   G4RotationMatrix* r34;
   G4RotationMatrix* r35;
   G4RotationMatrix* r36;
   G4RotationMatrix* r37;
   G4RotationMatrix* r38;
   G4RotationMatrix* r39;
   G4RotationMatrix* r40;
   G4RotationMatrix* r41;
   G4RotationMatrix* r42;
   G4RotationMatrix* r43;
   G4RotationMatrix* r57;
   G4RotationMatrix* r44;
   G4RotationMatrix* r45;
   G4RotationMatrix* r46;
   G4RotationMatrix* r47;
   G4RotationMatrix* r48;
   G4RotationMatrix* r49;
   G4RotationMatrix* r50;
   G4RotationMatrix* r51;
   G4RotationMatrix* r52;
   G4RotationMatrix* r53;
   G4RotationMatrix* r54;
   G4RotationMatrix* r55;
   G4RotationMatrix* r56;
   G4RotationMatrix* r58;
   G4RotationMatrix* r59;
   G4RotationMatrix* r60;
   G4RotationMatrix* r61;
   G4RotationMatrix* r62;
   G4RotationMatrix* r63;
   G4RotationMatrix* r64;
   G4RotationMatrix* r65;
   G4RotationMatrix* r66;
   G4RotationMatrix* r67;
   G4RotationMatrix* r68;
   G4RotationMatrix* r69;
   G4RotationMatrix* r70;
   G4RotationMatrix* r71;
   G4RotationMatrix* r72;
   G4RotationMatrix* r73;
   G4RotationMatrix* r74;
   G4RotationMatrix* r75;
   G4RotationMatrix* r76;
   G4RotationMatrix* rotateTest;
   G4RotationMatrix* rotateWhole;
  
// translation vectors
   G4ThreeVector moveTest;
   G4ThreeVector mshroud;
   G4ThreeVector mactive;
   G4ThreeVector mconduct;
   G4ThreeVector move;
   G4ThreeVector moveNull;
   G4ThreeVector movesquare;
   G4ThreeVector movetriangle;
   G4ThreeVector moverectangle1;
   G4ThreeVector moverectangle2;   
   G4ThreeVector movetest; 
   G4ThreeVector subright;
   G4ThreeVector subleft;
   G4ThreeVector m1;
   G4ThreeVector m2;
   G4ThreeVector m3;
   G4ThreeVector m4;
   G4ThreeVector m5;
   G4ThreeVector m6;
   G4ThreeVector m7;
   G4ThreeVector m8;
   G4ThreeVector m9;
   G4ThreeVector m10;
   G4ThreeVector m11;
   G4ThreeVector m12;
   G4ThreeVector m13;
   G4ThreeVector m14;
   G4ThreeVector m15;
   G4ThreeVector m16;
   G4ThreeVector m17;
   G4ThreeVector m18;
   G4ThreeVector m19;
   G4ThreeVector m20;
   G4ThreeVector m21;
   G4ThreeVector m22;
   G4ThreeVector m23;
   G4ThreeVector m24;
   G4ThreeVector m25;
   G4ThreeVector m26;
   G4ThreeVector m27;
   G4ThreeVector m28;
   G4ThreeVector m29;
   G4ThreeVector m30;
   G4ThreeVector m31;
   G4ThreeVector m32;
   G4ThreeVector m33;
   G4ThreeVector m34;
   G4ThreeVector m35;
   G4ThreeVector m36;
   G4ThreeVector m37;
   G4ThreeVector m38;
   G4ThreeVector m39;
   G4ThreeVector m40;
   G4ThreeVector m41;
   G4ThreeVector m42;
   G4ThreeVector m43;
   G4ThreeVector m44;
   G4ThreeVector m45;
   G4ThreeVector m46;
   G4ThreeVector m47;
   G4ThreeVector m48;
   G4ThreeVector m49;
   G4ThreeVector m50;
   G4ThreeVector m51;
   G4ThreeVector m52;
   G4ThreeVector m53;
   G4ThreeVector m54;
   G4ThreeVector m55;
   G4ThreeVector m56;
      

   G4ThreeVector m57;
   G4ThreeVector m58;
   G4ThreeVector m59;
   G4ThreeVector m60;
   G4ThreeVector m61;
   G4ThreeVector m62;
   G4ThreeVector m63;
   G4ThreeVector m64;
   G4ThreeVector m65;
   G4ThreeVector m66;
   G4ThreeVector m67;
   G4ThreeVector m68;
   G4ThreeVector m69;
   G4ThreeVector m70;
   G4ThreeVector m71;
   G4ThreeVector m72;
   G4ThreeVector m73;
   G4ThreeVector m74;
   G4ThreeVector m75;
   G4ThreeVector m76;
   G4ThreeVector mr1shroud;
   G4ThreeVector mtest;
   G4ThreeVector mr2shroud;
   G4ThreeVector mtshroud;
   G4ThreeVector mtactive;
   G4ThreeVector mtconduct;










    //##################################################################
    //### methods to construct all of the components of the detector ###
    //##################################################################
    G4int Construct();
   

   
};
#endif
