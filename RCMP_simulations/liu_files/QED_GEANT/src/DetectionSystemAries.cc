#include "DetectionSystemAries.hh"
#include "G4Trap.hh" 
#include "G4SystemOfUnits.hh" 
#include "G4Material.hh" 
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4OpticalSurface.hh"
#include "G4PVPlacement.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh" 
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"
#include "G4Para.hh"
#include "cmath"
#include "iostream"
#include "G4Material.hh"
#include "G4Scintillation.hh"
#include "G4OpticalPhysics.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"

/////////////////////////////////////////////////////////
//////////////////Class Initialization///////////////////
/////////////////////////////////////////////////////////
DetectionSystemAries::DetectionSystemAries() {



    //Dimensions of the tiles
    SeparateHemispheres = 2.0*mm;
    gapbtwntiles        = 0.1*mm; 
    scintthickness      = 1.5*mm;
    silverthickness     = 300.001*nm;
    scinttrianglelength = 32*mm;
    trianglelong        = 33.18*mm;
    squarelong          = 16.62*mm;
    reclong             = 15.62*mm;
    scintsquarelength   = 16*mm;
    scintrectanglewidth = 15*mm;   
    innerradius         = 38.75*mm;
    outerradius         = 40.25*mm;

    //Values for making trapezoids
    longslength         = scintsquarelength/2 + silverthickness/2;
    longrlength         = scintrectanglewidth/2+silverthickness/2;
    longerrectangle     = 15.62/2*mm;
    sl			= scintsquarelength+silverthickness*2;
    longerslength       = squarelong/2 +silverthickness/2;
    longthickness       = scintthickness/2 + silverthickness/2;
    smallSSL            = scintsquarelength -1*nm;
    smallSST            = scintthickness -1*nm;
    zeroangle           = 0*deg;
    halfrightangle      = 5*deg;
    taperanglex         = atan(tan(22.5*deg)/sqrt(2))*rad;
    taperangley         = 45*deg;
    expo   = 2.5*cm;
    tiny   = 0.00001*nm;
    zero   = 0.0*m;
    thetarect=atan((scintsquarelength-squarelong)/(sqrt(2)*scintthickness))*rad ;
    phirect=45*deg;


    //Values for placing the tiles
    sprotrude           = squarelong-scintsquarelength;
    realthick           = scintthickness+silverthickness*2;
    scmx   =   (6*realthick*pow(sl,2)*(2*sl+sprotrude)+3*sprotrude*realthick*(pow(sprotrude,2)+4*sl*sprotrude+4*pow(sl,2)))/(8*realthick*(pow(sprotrude,2)+3*sl*sprotrude+3*pow(sl,2)))-sl/2  ;
   
    scmz   =   (3*(pow(sprotrude,1))+8*sl*sprotrude+6*pow(sl,2))*realthick/(pow(sprotrude,2)+3*sl*sprotrude+3*pow(sl,2))/4  ;
    t1     =   scmx;
    t2     =   scmz/sqrt(2)-scmx/sqrt(2);
    t3     =   scmz/sqrt(2)+scmx/sqrt(2);
    ts     =   scinttrianglelength+silverthickness*2;
    tl     =   trianglelong+silverthickness*2;
    tcmz   =   (realthick*(3*pow(tl,2)+2*ts*tl+pow(ts,2)))/(4*(pow(tl,2)+ts*tl+pow(ts,2)));
    tt1    =   sl*(sqrt(2)+2/3)+tcmz/sqrt(3);
    tt2    =   sl*(-1/3+sqrt(2))-tcmz/sqrt(3);
    rcmb   =   (6*(scintrectanglewidth+reclong)*(pow(scintsquarelength,2))+(squarelong-scintsquarelength)*((3*scintrectanglewidth+5*reclong)*scintsquarelength+squarelong*scintrectanglewidth+3*reclong*squarelong))/(4*((2*scintrectanglewidth+reclong)*scintsquarelength+squarelong+scintrectanglewidth+2*reclong*squarelong));
    rcmt   =   (scintthickness/2)*((scintrectanglewidth+reclong)*scintsquarelength+squarelong*scintrectanglewidth+3*reclong*squarelong)/((2*scintrectanglewidth+reclong)*scintsquarelength+squarelong*reclong+2*reclong*squarelong);
    rcma   =   (6*(scintsquarelength+squarelong)*(pow(scintrectanglewidth,2))+(reclong-scintrectanglewidth)*((3*scintsquarelength+5*squarelong)*scintrectanglewidth+reclong*scintsquarelength+3*squarelong*reclong))/(4*((2*scintsquarelength+squarelong)*scintrectanglewidth+reclong+scintsquarelength+2*squarelong*reclong));
    rc1    =   (1+sqrt(2))*scintsquarelength;
    rc2    =   (scintsquarelength-scintrectanglewidth);
    tria   =   scinttrianglelength/2+silverthickness/2;
    trih   =   tan(15*deg)*tria;
    trid   =   sin(45*deg)*scinttrianglelength+silverthickness;


    //fudge
    rec1   =   sin(45*deg)*sl/2+0.1091*mm;
    rec2   =   scintrectanglewidth/2+sin(45*deg)*sl/2+0.26375*mm;
    rec3   =   0.155*mm;
    rec4   =   0.159*mm;
    rec5   =   0.309*mm;
    hrec   =   sin(45*deg)*(longslength)+0.109*mm;
    drec   =   hrec+longslength+0.151175*mm;
    tri1   =   2.4*mm;
    tri2   =   2.35*mm;
    tri3   =   0.3*mm;
    tri5   =   5.9*mm;
    tri6   =   5.6*mm;

    //SiPM dimensions
    sipmConductThick = 0.01*mm;
    sipmActiveThick  = 0.20*mm;
    sipmActiveDim    = 1.00*mm;
    sipmDim1         = 1.5*mm;
    sipmDim2         = 1.8*mm;
    sipmThick        = 0.65*mm;

    //rotation matrices and translation vector 
    moveNull       = G4ThreeVector(0.0,0.0,0.0);
    rotateNull = new G4RotationMatrix;     
    rotateWhole= new G4RotationMatrix();
    rotateWhole->rotateZ(M_PI/8);  
    rotateWhole->rotateY(M_PI);
    rotateWhole->rotateX(0);
      
    //retrieving materials
    scintmat   = G4Material::GetMaterial("EJ228");
    coatingmat = G4Material::GetMaterial("Aluminum");
    uhv        = G4Material::GetMaterial("Vacuum");
    si         = G4Material::GetMaterial("Silicon");
    copper     = G4Material::GetMaterial("Copper");
    mylar      = G4Material::GetMaterial("Mylar");

    //vizualization
    coatingvis = new G4VisAttributes(G4Colour(0.9,0.9,0.9));
    coatingvis->SetVisibility(true);
    coatingvis->SetForceWireframe(true);    
    scintvis   = new G4VisAttributes(G4Colour(0,1.0,1.0));
    scintvis->SetVisibility(false);
    scintvis->SetForceWireframe(true);
    uhvis   = new G4VisAttributes(G4Colour(0,1.0,1.0));
    uhvis->SetVisibility(false);
    uhvis->SetForceWireframe(true);
    sivis    = new G4VisAttributes(G4Colour(1.0,0.64,0));
    sivis->SetVisibility(true);
//  shroudvis = new G4VisAttributes(G4Colour(1.0,0.64,0));
//  shroudvis->SetVisibility(true);
//  conductvis = new G4VisAttributes(G4Colour(0.196,0.804,0.196));
//  conductvis->SetVisibility(true);
//  activevis =new G4VisAttributes(G4Colour(1.0,0,0));
//  activevis->SetVisibility(true);
 

  

         
}





/////////////////////////////////////////////////////////////
//////////////////// Class destruction///////////////////////
/////////////////////////////////////////////////////////////
DetectionSystemAries::~DetectionSystemAries() {
	delete fSquareTile;
	delete ftriangle;
	delete frectangle1;
	delete frectangle2;
	delete ftest;	
        delete sSipm;
}





/////////////////////////////////////////////////////////////
/////////////////////Build class/////////////////////////////
/////////////////////////////////////////////////////////////
void DetectionSystemAries::Build(){
	fSquareTile = new G4AssemblyVolume();
	ftriangle   = new G4AssemblyVolume();
	frectangle1 = new G4AssemblyVolume();
        frectangle2 = new G4AssemblyVolume();
	ftest	    = new G4AssemblyVolume();
	sSipm	    = new G4AssemblyVolume();
	ConstructSquareTile();
	ConstructTriangleTile();
	ConstructRectangleTile1();
	ConstructRectangleTile2();
	
}

/////////////////////////////////////////////////////////////////////////////////
//////////////////Class for describing tile position in space////////////////////
///////////////////////////////////////////////////////////////////////////////// 
void DetectionSystemAries::PlaceDetector(G4LogicalVolume* expHallLog){


	

	//imprints of position and rotations of assemblies in the world volume	
	fSquareTile->MakeImprint(expHallLog, moveNull,rotateWhole,500);
	frectangle1->MakeImprint(expHallLog,moveNull,rotateWhole,200);
	ftriangle->MakeImprint(expHallLog,moveNull,rotateWhole,100);	
	frectangle2->MakeImprint(expHallLog,moveNull,rotateWhole,300);
	ftest->MakeImprint(expHallLog,moveNull,rotateNull);

}



///////////////////////////////////////////////////////////////////////
////////////////Class to construct the square tiles////////////////////
///////////////////////////////////////////////////////////////////////
void DetectionSystemAries::ConstructSquareTile(){


	//creating coating
	G4VSolid* squarecoatingsolid      = new G4Trap("squarecoatingsolid",longthickness,taperanglex,taperangley,longerslength,longerslength,longerslength,zeroangle,longslength,longslength,longslength,zeroangle);	
	G4LogicalVolume* squarecoatinglog = new G4LogicalVolume(squarecoatingsolid,coatingmat,"squarecoatinglog");

	//making and placing the scintillator tile
	G4VSolid* squarescintsolid        = new G4Trap("squarescintsolid",scintthickness/2,taperanglex,taperangley,squarelong/2,squarelong/2,squarelong/2,zeroangle,scintsquarelength/2,scintsquarelength/2,scintsquarelength/2,zeroangle);	
	G4LogicalVolume* squarescintlog   = new G4LogicalVolume(squarescintsolid,scintmat,"squarescintlog");

	G4VPhysicalVolume* squarescintphys = new G4PVPlacement(rotateNull,moveNull,squarescintlog,"squarescintphys",squarecoatinglog,false,1,false);      	

	//Sipm
	
	rshroud=new G4RotationMatrix();		/////////////////
	rshroud->rotateX(0);			/////Rotation////
	rshroud->rotateY(-M_PI/2);		/////////////////
	rshroud->rotateZ(0);

	mshroud=G4ThreeVector(-scintsquarelength/2+1*mm,-scintsquarelength/2+1*mm,-scintthickness/2);
    	mactive =G4ThreeVector(scintthickness/2+sipmThick/2,0,0);	////translations////
    	mconduct=G4ThreeVector(sipmThick/2,0,0);

	
	G4VSolid* sipmPlastic = new G4Box("sipmPlastic",sipmThick/2, sipmDim1/2, sipmDim2/2);
    	G4VSolid* sipmActive = new G4Box("sipmActive", sipmActiveThick/2, sipmActiveDim/2, sipmActiveDim/2);
    	G4VSolid* sipmConduct = new G4Box("sipmConduct", sipmConductThick/2, sipmDim1/2, sipmDim2/2);
    	G4VSolid* sipmShroud = new G4SubtractionSolid("sipmShroud", sipmPlastic, sipmActive, 0, moveNull);
	G4LogicalVolume* sSipmShroudLog = new G4LogicalVolume(sipmShroud, mylar, "sShroudLog", 0, 0, 0);
    	G4LogicalVolume* sSipmActiveLog = new G4LogicalVolume(sipmActive, si, "sActiveLog", 0, 0, 0);
    	G4LogicalVolume* sSipmConductLog = new G4LogicalVolume(sipmConduct, copper, "sConductLog", 0, 0, 0);
	//G4VPhysicalVolume* sshroudphys = new G4PVPlacement(rshroud,mshroud,sSipmShroudLog,"sshroudphys",squarecoatinglog,false,1,false);	
	//G4VPhysicalVolume* sactivephys = new G4PVPlacement(ractive,moveNull,sSipmActiveLog,"sactivephys",sSipmShroudLog,false,1,false);	
	//G4VPhysicalVolume* sconductphys = new G4PVPlacement(rconduct,mconduct,sSipmConductLog,"sconductpys",sSipmShroudLog,false,1,false);	


	
	//visualization	
	squarecoatinglog->SetVisAttributes(coatingvis);
	squarescintlog->SetVisAttributes(scintvis);
        sSipmShroudLog->SetVisAttributes(shroudvis);
	sSipmActiveLog->SetVisAttributes(activevis);
	sSipmConductLog->SetVisAttributes(conductvis);
	

	//placing the tiles
	//Downstream	
	r1=new G4RotationMatrix(0,0,90*deg);
	m1=G4ThreeVector(-sl/2-scmx,sl/2+scmx,-sl*(1+sqrt(2))-scmz);
	fSquareTile->AddPlacedVolume(squarecoatinglog,m1,r1);

	r2=new G4RotationMatrix(0,0,180*deg);
	m2=G4ThreeVector(sl/2+scmx,sl/2+scmx,-sl*(1+sqrt(2))-scmz);
	fSquareTile->AddPlacedVolume(squarecoatinglog,m2,r2);

	r3=new G4RotationMatrix(0,0,270*deg);
	m3=G4ThreeVector(sl/2+scmx,-sl/2-scmx,-sl*(1+sqrt(2))-scmz);	
	fSquareTile->AddPlacedVolume(squarecoatinglog,m3,r3);

	r4=new G4RotationMatrix(0,0,0);
	m4=G4ThreeVector(-sl/2-scmx,-sl/2-scmx,-sl*(1+sqrt(2))-scmz);
	fSquareTile->AddPlacedVolume(squarecoatinglog,m4,r4);

	r5=new G4RotationMatrix(0,315*deg,0);
	m5=G4ThreeVector(-sl/2-t1,sl/(4*cos(45*deg))+sl+t2,-sl*(1+sqrt(2))+sl/(4*cos(45*deg))-t3);
	fSquareTile->AddPlacedVolume(squarecoatinglog,m5,r5);

	r6=new G4RotationMatrix();
        r6->rotateZ(M_PI/2);
	r6->rotateX(M_PI/4);
	m6=G4ThreeVector(sl/2+scmx,sl/(4*cos(45*deg))+sl+scmz/sqrt(2)-scmx/sqrt(2),-sl*(1+sqrt(2))+sl/(4*cos(45*deg))-scmz/sqrt(2)-scmx/sqrt(2));
	fSquareTile->AddPlacedVolume(squarecoatinglog,m6,r6);

	r7=new G4RotationMatrix();
	r7->rotateZ(3*M_PI/2);
	r7->rotateY(7*M_PI/4);
	m7=G4ThreeVector(sl/(4*cos(45*deg))+sl+t2,sl/2+t1,-sl*(1+sqrt(2))+sl/(4*cos(45*deg))-t3);
	fSquareTile->AddPlacedVolume(squarecoatinglog,m7,r7);

	r8=new G4RotationMatrix();	
	r8->rotateY(7*M_PI/4);
	m8=G4ThreeVector(sl/(4*cos(45*deg))+sl+t2,-sl/2-t1,-sl*(1+sqrt(2))+sl/(4*cos(45*deg))-t3);
	fSquareTile->AddPlacedVolume(squarecoatinglog,m8,r8);

	r9=new G4RotationMatrix();	
	r9->rotateZ(M_PI);
	r9->rotateX(7*M_PI/4);
	m9=G4ThreeVector(sl/2+t1,-sl/(4*cos(45*deg))-sl-t2,-sl*(1+sqrt(2))+sl/(4*cos(45*deg))-t3);
	fSquareTile->AddPlacedVolume(squarecoatinglog,m9,r9);

	r10=new G4RotationMatrix();
	r10->rotateZ(3*M_PI/2);
	r10->rotateX(7*M_PI/4);
	m10=G4ThreeVector(-sl/2-t1,-sl/(4*cos(45*deg))-sl-t2,-sl*(1+sqrt(2))+sl/(4*cos(45*deg))-t3);
	fSquareTile->AddPlacedVolume(squarecoatinglog,m10,r10);

	r11=new G4RotationMatrix();
	r11->rotateZ(M_PI/2);
	r11->rotateY(M_PI/4);
	m11=G4ThreeVector(-sl/(4*cos(45*deg))-sl-t2,-sl/2-t1,-sl*(1+sqrt(2))+sl/(4*cos(45*deg))-t3);
	fSquareTile->AddPlacedVolume(squarecoatinglog,m11,r11);

	r12=new G4RotationMatrix();
	r12->rotateZ(M_PI);
	r12->rotateY(M_PI/4);	
	m12=G4ThreeVector(-sl/(4*cos(45*deg))-sl-t2,sl/2+t1,-sl*(1+sqrt(2))+sl/(4*cos(45*deg))-t3);
	fSquareTile->AddPlacedVolume(squarecoatinglog,m12,r12);

	r17=new G4RotationMatrix();
	r17->rotateZ(3*M_PI/2);
	r17->rotateX(M_PI/4);
	m17=G4ThreeVector(-sl/2-t1,3*(sl/(4*cos(45*deg)))+sl+t3,(-sl*(1+0.25*sqrt(2)))-t2);
	fSquareTile->AddPlacedVolume(squarecoatinglog,m17,r17);

	r18=new G4RotationMatrix();
        r18->rotateZ(M_PI);
	r18->rotateX(M_PI/4);
	m18=G4ThreeVector(sl/2+t1,3*(sl/(4*cos(45*deg)))+sl+t3,(-sl*(1+0.25*sqrt(2)))-t2);
	fSquareTile->AddPlacedVolume(squarecoatinglog,m18,r18);

	r19=new G4RotationMatrix();
	r19->rotateZ(M_PI);
	r19->rotateY(7*M_PI/4);
	m19=G4ThreeVector(3*(sl/(4*cos(45*deg)))+sl+t3,sl/2+t1,-sl*(1+0.25*sqrt(2))-t2);
	fSquareTile->AddPlacedVolume(squarecoatinglog,m19,r19);

	r20=new G4RotationMatrix();
	r20->rotateZ(M_PI/2);
	r20->rotateY(7*M_PI/4);
	m20=G4ThreeVector(3*(sl/(4*cos(45*deg)))+t3+sl,-sl/2-t1,-sl*(1+0.25*sqrt(2))-t2);
	fSquareTile->AddPlacedVolume(squarecoatinglog,m20,r20);

	r21=new G4RotationMatrix();
	r21->rotateZ(M_PI/2);
        r21->rotateX(7*M_PI/4);
	m21=G4ThreeVector(sl/2+t1,3*-sl/(4*cos(45*deg))-sl-t3,-sl*(1+0.25*sqrt(2))-t2);
	fSquareTile->AddPlacedVolume(squarecoatinglog,m21,r21);

	r22=new G4RotationMatrix();
        r22->rotateX(7*M_PI/4);
	m22=G4ThreeVector(-sl/2-t1,3*-sl/(4*cos(45*deg))-sl-t3,-sl*(1+0.25*sqrt(2))-t2);
	fSquareTile->AddPlacedVolume(squarecoatinglog,m22,r22);

	r23=new G4RotationMatrix();
	r23->rotateY(M_PI/4);
	m23=G4ThreeVector(-3*(sl/(4*cos(45*deg)))-sl-t3,-sl/2-t1,-sl*(1+0.25*sqrt(2))-t2);
	fSquareTile->AddPlacedVolume(squarecoatinglog,m23,r23);

	r24=new G4RotationMatrix();
	r24->rotateZ(3*M_PI/2);
	r24->rotateY(M_PI/4);
	m24=G4ThreeVector(-3*(sl/(4*cos(45*deg)))-sl-t3,+sl/2+t1,-sl*(1+0.25*sqrt(2))-t2);
	fSquareTile->AddPlacedVolume(squarecoatinglog,m24,r24);


	// //Upstream
	// r57=new G4RotationMatrix();
	// r57->rotateZ(0);
	// r57->rotateY(0);
	// r57->rotateX(3*M_PI/4);
	// m57=G4ThreeVector(-sl/2-t1,3*(sl/(4*cos(45*deg)))+sl+t3,sl*(1+0.25*sqrt(2))+t2);
	// fSquareTile->AddPlacedVolume(squarecoatinglog,m57,r57);

	// r58=new G4RotationMatrix();
    //     r58->rotateY(0);
	// r58->rotateZ(M_PI/2);
	// r58->rotateX(3*M_PI/4);
	// m58=G4ThreeVector(sl/2+t1,3*(sl/(4*cos(45*deg)))+sl+t3,(sl*(1+0.25*sqrt(2)))+t2);
	// fSquareTile->AddPlacedVolume(squarecoatinglog,m58,r58);

	// r59=new G4RotationMatrix();
	// r59->rotateZ(3*M_PI/2);
	// r59->rotateY(5*M_PI/4);
	// m59=G4ThreeVector(3*(sl/(4*cos(45*deg)))+sl+t3,sl/2+t1,sl*(1+0.25*sqrt(2))+t2);
	// fSquareTile->AddPlacedVolume(squarecoatinglog,m59,r59);

	// r60=new G4RotationMatrix();
	// r60->rotateZ(0);
	// r60->rotateY(5*M_PI/4);
	// m60=G4ThreeVector(3*(sl/(4*cos(45*deg)))+t3+sl,-sl/2-t1,sl*(1+0.25*sqrt(2))+t2);
	// fSquareTile->AddPlacedVolume(squarecoatinglog,m60,r60);

	// r61=new G4RotationMatrix();
	// r61->rotateY(M_PI);
	// r61->rotateZ(0);
    //     r61->rotateX(M_PI/4);
	// m61=G4ThreeVector(sl/2+t1,3*-sl/(4*cos(45*deg))-sl-t3,sl*(1+0.25*sqrt(2))+t2);
	// fSquareTile->AddPlacedVolume(squarecoatinglog,m61,r61);

	// r62=new G4RotationMatrix();
    //     r62->rotateY(M_PI);
	// r62->rotateZ(3*M_PI/2);
	// r62->rotateX(M_PI/4);
	// m62=G4ThreeVector(-sl/2-t1,3*-sl/(4*cos(45*deg))-sl-t3,sl*(1+0.25*sqrt(2))+t2);
	// fSquareTile->AddPlacedVolume(squarecoatinglog,m62,r62);

	// r63=new G4RotationMatrix();
	// r63->rotateZ(M_PI/2);
	// r63->rotateY(3*M_PI/4);
	// m63=G4ThreeVector(-3*(sl/(4*cos(45*deg)))-sl-t3,-sl/2-t1,sl*(1+0.25*sqrt(2))+t2);
	// fSquareTile->AddPlacedVolume(squarecoatinglog,m63,r63);

	// r64=new G4RotationMatrix();
	// r64->rotateZ(M_PI);
	// r64->rotateY(3*M_PI/4);
	// m64=G4ThreeVector(-3*(sl/(4*cos(45*deg)))-sl-t3,+sl/2+t1,sl*(1+0.25*sqrt(2))+t2);
	// fSquareTile->AddPlacedVolume(squarecoatinglog,m64,r64);

	// r69=new G4RotationMatrix();
	// r69->rotateY(0);
    //     r69->rotateZ(3*M_PI/2);
	// r69->rotateX(3*M_PI/4);
	// m69=G4ThreeVector(-sl/2-t1,sl/(4*cos(45*deg))+sl+t2,sl*(1+sqrt(2))-sl/(4*cos(45*deg))+t3);
	// fSquareTile->AddPlacedVolume(squarecoatinglog,m69,r69);

	// r70=new G4RotationMatrix();
    //     r70->rotateY(0);
	// r70->rotateZ(M_PI);
	// r70->rotateX(3*M_PI/4);
	// m70=G4ThreeVector(sl/2+scmx,sl/(4*cos(45*deg))+sl+scmz/sqrt(2)-scmx/sqrt(2),sl*(1+sqrt(2))-sl/(4*cos(45*deg))+scmz/sqrt(2)+scmx/sqrt(2));
	// fSquareTile->AddPlacedVolume(squarecoatinglog,m70,r70);


	// r71=new G4RotationMatrix();
	// r71->rotateZ(M_PI);
	// r71->rotateY(5*M_PI/4);
	// m71=G4ThreeVector(sl/(4*cos(45*deg))+sl+t2,sl/2+t1,+sl*(1+sqrt(2))-sl/(4*cos(45*deg))+t3);
	// fSquareTile->AddPlacedVolume(squarecoatinglog,m71,r71);

	// r72=new G4RotationMatrix();	
	// r72->rotateZ(M_PI/2);
	// r72->rotateY(5*M_PI/4);
	// m72=G4ThreeVector(sl/(4*cos(45*deg))+sl+t2,-sl/2-t1,+sl*(1+sqrt(2))-sl/(4*cos(45*deg))+t3);
	// fSquareTile->AddPlacedVolume(squarecoatinglog,m72,r72);


	// r73=new G4RotationMatrix();	
	// r73->rotateY(M_PI);
	// r73->rotateZ(M_PI/2);
	// r73->rotateX(M_PI/4);
	// m73=G4ThreeVector(sl/2+t1,-sl/(4*cos(45*deg))-sl-t2,sl*(1+sqrt(2))-sl/(4*cos(45*deg))+t3);
	// fSquareTile->AddPlacedVolume(squarecoatinglog,m73,r73);

	// r74=new G4RotationMatrix();
	// r74->rotateY(M_PI);
	// r74->rotateZ(M_PI);
	// r74->rotateX(M_PI/4);
	// m74=G4ThreeVector(-sl/2-t1,-sl/(4*cos(45*deg))-sl-t2,sl*(1+sqrt(2))-sl/(4*cos(45*deg))+t3);
	// fSquareTile->AddPlacedVolume(squarecoatinglog,m74,r74);


	// r75=new G4RotationMatrix();
	// r75->rotateZ(0);
	// r75->rotateY(3*M_PI/4);
	// m75=G4ThreeVector(-sl/(4*cos(45*deg))-sl-t2,-sl/2-t1,+sl*(1+sqrt(2))-sl/(4*cos(45*deg))+t3);
	// fSquareTile->AddPlacedVolume(squarecoatinglog,m75,r75);

	// r76=new G4RotationMatrix();
	// r76->rotateZ(3*M_PI/2);
	// r76->rotateY(3*M_PI/4);	
	// m76=G4ThreeVector(-sl/(4*cos(45*deg))-sl-t2,sl/2+t1,sl*(1+sqrt(2))-sl/(4*cos(45*deg))+t3);
	// fSquareTile->AddPlacedVolume(squarecoatinglog,m76,r76);



	


}





///////////////////////////////////////////////////////////////////////
///////////////Class to construct the triangle tiles///////////////////
///////////////////////////////////////////////////////////////////////
void DetectionSystemAries::ConstructTriangleTile(){
	

	//making coating using generic trapezoid definition
	G4VSolid* trianglecoatingsolid2     = new G4Trap("trianglecoatingsolid2",longthickness,zeroangle,zeroangle,(trianglelong/2+silverthickness/2)*(sqrt(3)/2),tiny,trianglelong/2+silverthickness/2,zeroangle,(scinttrianglelength/2+silverthickness/2)*(sqrt(3)/2),tiny,scinttrianglelength/2+silverthickness/2,zeroangle);	
	G4LogicalVolume* trianglecoatinglog = new G4LogicalVolume(trianglecoatingsolid2,coatingmat,"trianglecoatinglog");	

	//making and placing scint tile in the coating
	G4VSolid* trianglescint             = new G4Trap("trianglescint",scintthickness/2,zeroangle,zeroangle,trianglelong/2*(sqrt(3)/2),tiny,trianglelong/2,zeroangle,scinttrianglelength/2*(sqrt(3)/2),tiny,scinttrianglelength/2,zeroangle);
	G4LogicalVolume* scinttrianglelog   = new G4LogicalVolume(trianglescint,scintmat,"scinttrianglelog");
	G4VPhysicalVolume* trianglephys     = new G4PVPlacement(rotateNull,moveNull,scinttrianglelog,"trianglephys",trianglecoatinglog,false,1,false);

	//adding SiPm
	rshroud=new G4RotationMatrix();
	rshroud->rotateX(0);
	rshroud->rotateY(-M_PI/2);
	rshroud->rotateZ(0);

	mtshroud=G4ThreeVector(-4*mm,4*mm,-scintthickness/2);
    	mactive =G4ThreeVector(scintthickness/2+sipmThick/2,0,0);
    	mconduct=G4ThreeVector(sipmThick/2,0,0);




	G4VSolid* tsipmPlastic = new G4Box("tsipmPlastic",sipmThick/2, sipmDim1/2, sipmDim2/2);
    	G4VSolid* tsipmActive = new G4Box("tsipmActive", sipmActiveThick/2, sipmActiveDim/2, sipmActiveDim/2);
    	G4VSolid* tsipmConduct = new G4Box("tsipmConduct", sipmConductThick/2, sipmDim1/2, sipmDim2/2);
    	G4VSolid* tsipmShroud = new G4SubtractionSolid("tsipmShroud", tsipmPlastic, tsipmActive, 0, moveNull);
 
	G4LogicalVolume* tSipmShroudLog = new G4LogicalVolume(tsipmShroud, mylar, "tShroudLog", 0, 0, 0);
    	G4LogicalVolume* tSipmActiveLog = new G4LogicalVolume(tsipmActive, si, "tActiveLog", 0, 0, 0);
    	G4LogicalVolume* tSipmConductLog = new G4LogicalVolume(tsipmConduct, copper, "tConductLog", 0, 0, 0);
 
	//G4VPhysicalVolume* tshroudphys = new G4PVPlacement(rshroud,mtshroud,tSipmShroudLog,"tshroudphys",trianglecoatinglog,false,1,false);	

	//G4VPhysicalVolume* tactivephys = new G4PVPlacement(ractive,moveNull,tSipmActiveLog,"tactivephys",tSipmShroudLog,false,1,false);	

	//G4VPhysicalVolume* tconductphys = new G4PVPlacement(rconduct,mconduct,tSipmConductLog,"tconductpys",tSipmShroudLog,false,1,false);	



	//vizualization	
	trianglecoatinglog->SetVisAttributes(coatingvis);	
	scinttrianglelog->SetVisAttributes(scintvis);
	tSipmShroudLog->SetVisAttributes(shroudvis);
	tSipmActiveLog->SetVisAttributes(activevis);
	tSipmConductLog->SetVisAttributes(conductvis);
	

	//placing tiles
	r13=new G4RotationMatrix();
	r13->rotateZ(M_PI/3);
	r13->rotateX(M_PI/5.1);
	r13->rotateY(M_PI/4);
	m13=G4ThreeVector(-1.5*(sl/(4*cos(45*deg)))-sl+tri1,1.5*(sl/(4*cos(45*deg)))+sl-tri2,-sl*(1+0.5*sqrt(2))-tcmz+tri3);
	ftriangle->AddPlacedVolume(trianglecoatinglog,m13,r13);
		
	r14=new G4RotationMatrix();
	r14->rotateZ(M_PI/3);
	r14->rotateX(M_PI/5.1);
	r14->rotateY(-M_PI/4);
	m14=G4ThreeVector(-1*(-1.5*(sl/(4*cos(45*deg)))-sl)-tri1+tri5,1.5*(sl/(4*cos(45*deg)))+sl-tri2,-sl*(1+0.5*sqrt(2))-tcmz+tri3+tri6);
	ftriangle->AddPlacedVolume(trianglecoatinglog,m14,r14);




	r15=new G4RotationMatrix();
	r15->rotateZ(2*M_PI/3);
	r15->rotateX(-M_PI/5.1);
	r15->rotateY(-M_PI/4);
	m15=G4ThreeVector(-1*(-1.5*(sl/(4*cos(45*deg)))-sl)-tri1+tri5,-1*(1.5*(sl/(4*cos(45*deg)))+sl)+tri2,-sl*(1+0.5*sqrt(2))-tcmz+tri3+tri6);
	ftriangle->AddPlacedVolume(trianglecoatinglog,m15,r15);



	
	r16=new G4RotationMatrix();
	r16->rotateZ(2*M_PI/3);
	r16->rotateX(-M_PI/5.1);
	r16->rotateY(M_PI/4);
	m16=G4ThreeVector(-1.5*(sl/(4*cos(45*deg)))-sl+tri1,-1*(1.5*(sl/(4*cos(45*deg)))+sl)+tri2,-sl*(1+0.5*sqrt(2))-tcmz+tri3);
	ftriangle->AddPlacedVolume(trianglecoatinglog,m16,r16);




	// r65=new G4RotationMatrix();
	// r65->rotateZ(M_PI/3);
	// r65->rotateY(M_PI);
	// r65->rotateX(-M_PI/5.1);
	// r65->rotateY(-M_PI/4);
	// m65=G4ThreeVector(1*(-1.5*(sl/(4*cos(45*deg)))-sl)+tri1*mm-tri5,1.5*(sl/(4*cos(45*deg)))+sl-tri2*mm,-1*(-sl*(1+0.5*sqrt(2))-tcmz)-tri3*mm-tri6);
	// ftriangle->AddPlacedVolume(trianglecoatinglog,m65,r65);



	
	// r66=new G4RotationMatrix();
	// r66->rotateZ(M_PI/3);
	// r66->rotateY(M_PI);
	// r66->rotateX(-M_PI/5.1);
	// r66->rotateY(M_PI/4);
	// m66=G4ThreeVector(-1*(-1.5*(sl/(4*cos(45*deg)))-sl)-tri1*mm,1.5*(sl/(4*cos(45*deg)))+sl-tri2*mm,-1*(-sl*(1+0.5*sqrt(2))-tcmz)-tri3*mm);
	// ftriangle->AddPlacedVolume(trianglecoatinglog,m66,r66);




	
	// r67=new G4RotationMatrix();	
	// r67->rotateZ(2*M_PI/3);
	// r67->rotateY(M_PI);
	// r67->rotateX(M_PI/5.1);
	// r67->rotateY(M_PI/4);
	// m67=G4ThreeVector(-1*(-1.5*(sl/(4*cos(45*deg)))-sl)-tri1*mm,-1*(1.5*(sl/(4*cos(45*deg)))+sl-tri2*mm),-1*(-sl*(1+0.5*sqrt(2))-tcmz)-tri3*mm);
	// ftriangle->AddPlacedVolume(trianglecoatinglog,m67,r67);



	
	// r68=new G4RotationMatrix();
	// r68->rotateZ(2*M_PI/3);
	// r68->rotateY(M_PI);
	// r68->rotateX(M_PI/5.1);
	// r68->rotateY(-M_PI/4);
	// m68=G4ThreeVector(1*(-1.5*(sl/(4*cos(45*deg)))-sl)+tri1*mm-tri5,-1*(1.5*(sl/(4*cos(45*deg)))+sl-tri2*mm),-1*(-sl*(1+0.5*sqrt(2))-tcmz)-tri3*mm-tri6);
	// ftriangle->AddPlacedVolume(trianglecoatinglog,m68,r68);



}



///////////////////////////////////////////////////////////////////////
////////////////Class to construct rectangles type1////////////////////
///////////////////////////////////////////////////////////////////////
void DetectionSystemAries::ConstructRectangleTile1(){
	


	//creating coating
	G4VSolid* coatingrect1solid      = new G4Trap("coatingrect1solid",longthickness,thetarect,phirect-90*deg,longerrectangle,longerslength,longerslength,zeroangle,longrlength,longslength,longslength,zeroangle);	

	G4LogicalVolume* coatingrect1log = new G4LogicalVolume(coatingrect1solid,coatingmat,"coatingrect1log");

	//making the scintillator tile
	G4VSolid* rect1sol      = new G4Trap("rect1sol",scintthickness/2,thetarect,phirect-90*deg,longerrectangle-silverthickness/2,longerslength-silverthickness/2,longerslength-silverthickness/2,zeroangle,longrlength-silverthickness/2,longslength-silverthickness/2,longslength-silverthickness/2,zeroangle);	

	
//("rect1sol",scintthickness/2,thetarect,phirect,longerrectangle-silverthickness/2,longerslength-silverthickness/2,longerslength-silverthickness/2,zeroangle,longrlength-silverthickness/2,longslength-silverthickness/2,longslength-silverthickness/2,zeroangle);	




	G4LogicalVolume* rect1log = new G4LogicalVolume(rect1sol,coatingmat,"rect1log");

	G4VPhysicalVolume* rect1phys	= new G4PVPlacement(rotateNull,moveNull,rect1log,"rect1phys",coatingrect1log,false,1,false);


	//Adding the SiPMs to the tiles
	//adding SiPm
	rshroud=new G4RotationMatrix();
	rshroud->rotateX(0);
	rshroud->rotateY(-M_PI/2);
	rshroud->rotateZ(0);

	mr1shroud=G4ThreeVector(squarelong/2-1*mm,-reclong/2+1*mm,-scintthickness/2);
    	mactive =G4ThreeVector(scintthickness/2+sipmThick/2,0,0);
    	mconduct=G4ThreeVector(sipmThick/2,0,0);




	G4VSolid* r1sipmPlastic = new G4Box("r1sipmPlastic",sipmThick/2, sipmDim1/2, sipmDim2/2);
    	G4VSolid* r1sipmActive = new G4Box("r1sipmActive", sipmActiveThick/2, sipmActiveDim/2, sipmActiveDim/2);
    	G4VSolid* r1sipmConduct = new G4Box("r1sipmConduct", sipmConductThick/2, sipmDim1/2, sipmDim2/2);
    	G4VSolid* r1sipmShroud = new G4SubtractionSolid("r1sipmShroud", r1sipmPlastic, r1sipmActive, 0, moveNull);
 
	G4LogicalVolume* r1SipmShroudLog = new G4LogicalVolume(r1sipmShroud, mylar, "r1ShroudLog", 0, 0, 0);
    	G4LogicalVolume* r1SipmActiveLog = new G4LogicalVolume(r1sipmActive, si, "r1ActiveLog", 0, 0, 0);
    	G4LogicalVolume* r1SipmConductLog = new G4LogicalVolume(r1sipmConduct, copper, "r1ConductLog", 0, 0, 0);
 
//	G4VPhysicalVolume* r1shroudphys = new G4PVPlacement(rshroud,mr1shroud,r1SipmShroudLog,"r1shroudphys",coatingrect1log,false,1,false);	

//	G4VPhysicalVolume* r1activephys = new G4PVPlacement(ractive,moveNull,r1SipmActiveLog,"r1activephys",r1SipmShroudLog,false,1,false);	

//	G4VPhysicalVolume* r1conductphys = new G4PVPlacement(rconduct,mconduct,r1SipmConductLog,"r1conductpys",r1SipmShroudLog,false,1,false);	



	//Visualization
	coatingrect1log->SetVisAttributes(coatingvis);
	rect1log->SetVisAttributes(scintvis);
	r1SipmShroudLog->SetVisAttributes(shroudvis);
	r1SipmActiveLog->SetVisAttributes(activevis);
	r1SipmConductLog->SetVisAttributes(conductvis);
	



	//placing
	
	r26=new G4RotationMatrix();	m26=G4ThreeVector();
	r26->rotateZ(0);
	r26->rotateY(M_PI/4);
	r26->rotateX(M_PI/2);
	m26=G4ThreeVector(sl/2+t1-sl-drec-rec5,3*(sl/(4*cos(45*deg)))+sl+t3+(rec1)-hrec,(-sl*(1+0.25*sqrt(2)))-t2+(rec2));
	frectangle1->AddPlacedVolume(coatingrect1log,m26,r26);


	r28=new G4RotationMatrix();	
	r28->rotateZ(0);
	r28->rotateX(M_PI/2);
	m28=G4ThreeVector(sl/2+t1+rec3-rec4,3*(sl/(4*cos(45*deg)))+sl+t3+rec1,(-sl*(1+0.25*sqrt(2)))-t2+(rec2));
	frectangle1->AddPlacedVolume(coatingrect1log,m28,r28);


	r30=new G4RotationMatrix();
	r30->rotateZ(0);
	r30->rotateY(-M_PI/4);
	r30->rotateX(M_PI/2);
	m30=G4ThreeVector(3*(sl/(4*cos(45*deg)))+sl+t3+rec1-hrec,sl/2+t1+drec,-sl*(1+0.25*sqrt(2))-t2+rec2);

	frectangle1->AddPlacedVolume(coatingrect1log,m30,r30);



	r32=new G4RotationMatrix();
	r32->rotateZ(3*M_PI/2);
	r32->rotateY(3*M_PI/2);
	m32=G4ThreeVector(3*(sl/(4*cos(45*deg)))+t3+sl+rec1,-sl/2-t1,-sl*(1+0.25*sqrt(2))-t2+rec2);
	frectangle1->AddPlacedVolume(coatingrect1log,m32,r32);


	r34=new G4RotationMatrix();
	r34->rotateZ(M_PI);
	r34->rotateY(-M_PI/4);
        r34->rotateX(-M_PI/2);
	m34=G4ThreeVector(sl/2+t1+drec,3*-sl/(4*cos(45*deg))-sl-t3-rec1+hrec,-sl*(1+0.25*sqrt(2))-t2+rec2);
	frectangle1->AddPlacedVolume(coatingrect1log,m34,r34);


	r36=new G4RotationMatrix();
	r36->rotateZ(M_PI);
        r36->rotateX(3*M_PI/2);
	m36=G4ThreeVector(-(sl/2+t1),3*-sl/(4*cos(45*deg))-sl-t3-rec1,-sl*(1+0.25*sqrt(2))-t2+rec2);
	frectangle1->AddPlacedVolume(coatingrect1log,m36,r36);

	r38=new G4RotationMatrix();
	r38->rotateZ(M_PI);
	r38->rotateY(M_PI/4);
	r38->rotateX(-M_PI/2);
	m38=G4ThreeVector(-3*(sl/(4*cos(45*deg)))-sl-t3-rec1+hrec,-(sl/2+t1)-drec,-sl*(1+0.25*sqrt(2))-t2+rec2);
	frectangle1->AddPlacedVolume(coatingrect1log,m38,r38);


	r40=new G4RotationMatrix();
	r40->rotateZ(M_PI/2);
	r40->rotateY(M_PI/2);
	m40=G4ThreeVector(-3*(sl/(4*cos(45*deg)))-sl-t3-rec1,+sl/2+t1,-sl*(1+0.25*sqrt(2))-t2+rec2);
	frectangle1->AddPlacedVolume(coatingrect1log,m40,r40);


	// r41=new G4RotationMatrix();
	// r41->rotateZ(M_PI);
	// r41->rotateY(M_PI/4);
	// r41->rotateX(M_PI/2);
	// m41=G4ThreeVector(-3*(sl/(4*cos(45*deg)))-sl-t3-rec1+hrec,+sl/2+t1+drec,sl*(1+0.25*sqrt(2))+t2-rec2);
	// frectangle1->AddPlacedVolume(coatingrect1log,m41,r41);

	// r43=new G4RotationMatrix();
	// r43->rotateZ(M_PI);
	// r43->rotateX(M_PI/2);
	// //m43=G4ThreeVector(-8.5*mm,35*mm,6*mm);
	// m43=G4ThreeVector(sl/2+t1-sl-rec3-rec4,3*(sl/(4*cos(45*deg)))+sl+t3+(rec1),(sl*(1+0.25*sqrt(2)))+t2-rec2);
	// frectangle1->AddPlacedVolume(coatingrect1log,m43,r43);


	// r45=new G4RotationMatrix();
	// r45->rotateZ(M_PI);
	// r45->rotateY(-M_PI/4);
	// r45->rotateX(M_PI/2);
	// m45=G4ThreeVector(sl/2+t1+drec,3*(sl/(4*cos(45*deg)))+sl+t3+rec1-hrec,(sl*(1+0.25*sqrt(2)))+t2-(rec2));
	// frectangle1->AddPlacedVolume(coatingrect1log,m45,r45);



	// r47=new G4RotationMatrix();
	// r47->rotateZ(M_PI/2);
	// r47->rotateY(3*M_PI/2);
	// m47=G4ThreeVector(3*(sl/(4*cos(45*deg)))+sl+t3+rec1,sl/2+t1,sl*(1+0.25*sqrt(2))+t2-rec2);

	// frectangle1->AddPlacedVolume(coatingrect1log,m47,r47);






	// r49=new G4RotationMatrix();
	// r49->rotateZ(0);
	// r49->rotateY(-M_PI/4);
	// r49->rotateX(-M_PI/2);
	// m49=G4ThreeVector(3*(sl/(4*cos(45*deg)))+t3+sl+rec1-hrec,-sl/2-t1-drec,sl*(1+0.25*sqrt(2))+t2-rec2);
	// frectangle1->AddPlacedVolume(coatingrect1log,m49,r49);



	// r51=new G4RotationMatrix();
	// r51->rotateZ(0);
    //     r51->rotateX(3*M_PI/2);
	// m51=G4ThreeVector(sl/2+t1,3*-sl/(4*cos(45*deg))-sl-t3-rec1,sl*(1+0.25*sqrt(2))+t2-rec2);
	// frectangle1->AddPlacedVolume(coatingrect1log,m51,r51);



	// r53=new G4RotationMatrix();
	// r53->rotateZ(0);
	// r53->rotateY(M_PI/4);
    //     r53->rotateX(-M_PI/2);
	// m53=G4ThreeVector(-(sl/2+t1)-drec,3*-sl/(4*cos(45*deg))-sl-t3-rec1+hrec,sl*(1+0.25*sqrt(2))+t2-rec2);
	// frectangle1->AddPlacedVolume(coatingrect1log,m53,r53);



	// r55=new G4RotationMatrix();
	// r55->rotateZ(3*M_PI/2);
	// r55->rotateY(M_PI/2);
	// m55=G4ThreeVector(-3*(sl/(4*cos(45*deg)))-sl-t3-rec1,-(sl/2+t1),sl*(1+0.25*sqrt(2))+t2-rec2);
	// frectangle1->AddPlacedVolume(coatingrect1log,m55,r55);




}






///////////////////////////////////////////////////////////////////////
////////////////Class to construct rectangles type2////////////////////
///////////////////////////////////////////////////////////////////////
void DetectionSystemAries::ConstructRectangleTile2(){



	//creating coating
	G4VSolid* coatingrect2solid      = new G4Trap("coatingrect2solid",longthickness,thetarect,phirect,longerrectangle,longerslength,longerslength,zeroangle,longrlength,longslength,longslength,zeroangle);	

	G4LogicalVolume* coatingrect2log = new G4LogicalVolume(coatingrect2solid,coatingmat,"coatingrect2log");



	//making scintillator
	G4VSolid* rect2sol      = new G4Trap("rect2sol",scintthickness/2,thetarect,phirect,longerrectangle-silverthickness/2,longerslength-silverthickness/2,longerslength-silverthickness/2,zeroangle,longrlength-silverthickness/2,longslength-silverthickness/2,longslength-silverthickness/2,zeroangle);	


	G4LogicalVolume* rect2log = new G4LogicalVolume(rect2sol,coatingmat,"rect2log");

	G4VPhysicalVolume* rect2phys	= new G4PVPlacement(rotateNull,moveNull,rect2log,"rect2phys",coatingrect2log,false,1,false);




	
	//making SiPM
	rshroud=new G4RotationMatrix();
	rshroud->rotateX(0);
	rshroud->rotateY(-M_PI/2);
	rshroud->rotateZ(0);

	mr2shroud=G4ThreeVector(squarelong/2-1*mm,reclong/2-1*mm,-scintthickness/2);
    	mactive =G4ThreeVector(scintthickness/2+sipmThick/2,0,0);
    	mconduct=G4ThreeVector(sipmThick/2,0,0);




	G4VSolid* r2sipmPlastic = new G4Box("r2sipmPlastic",sipmThick/2, sipmDim1/2, sipmDim2/2);
    	G4VSolid* r2sipmActive = new G4Box("r2sipmActive", sipmActiveThick/2, sipmActiveDim/2, sipmActiveDim/2);
    	G4VSolid* r2sipmConduct = new G4Box("r2sipmConduct", sipmConductThick/2, sipmDim1/2, sipmDim2/2);
    	G4VSolid* r2sipmShroud = new G4SubtractionSolid("r2sipmShroud", r2sipmPlastic, r2sipmActive, 0, moveNull);
 
	G4LogicalVolume* r2SipmShroudLog = new G4LogicalVolume(r2sipmShroud, mylar, "r2ShroudLog", 0, 0, 0);
    	G4LogicalVolume* r2SipmActiveLog = new G4LogicalVolume(r2sipmActive, si, "r2ActiveLog", 0, 0, 0);
    	G4LogicalVolume* r2SipmConductLog = new G4LogicalVolume(r2sipmConduct, copper, "r2ConductLog", 0, 0, 0);
 
//	G4VPhysicalVolume* r2shroudphys = new G4PVPlacement(rshroud,mr2shroud,r2SipmShroudLog,"r2shroudphys",coatingrect2log,false,1,false);	

//	G4VPhysicalVolume* r2activephys = new G4PVPlacement(ractive,moveNull,r2SipmActiveLog,"r2activephys",r2SipmShroudLog,false,1,false);	

//	G4VPhysicalVolume* r2conductphys = new G4PVPlacement(rconduct,mconduct,r2SipmConductLog,"r2conductpys",r2SipmShroudLog,false,1,false);	





	//Visualization
	coatingrect2log->SetVisAttributes(coatingvis);
	rect2log->SetVisAttributes(scintvis);
	r2SipmShroudLog->SetVisAttributes(shroudvis);
	r2SipmActiveLog->SetVisAttributes(activevis);
	r2SipmConductLog->SetVisAttributes(conductvis);
	


	//placing detectors

	r25=new G4RotationMatrix();
	r25->rotateZ(M_PI);
	r25->rotateY(M_PI/4);
	r25->rotateX(M_PI/2);
	m25=G4ThreeVector(-3*(sl/(4*cos(45*deg)))-sl-t3-rec1+hrec,+sl/2+t1+drec,-sl*(1+0.25*sqrt(2))-t2+rec2);
	frectangle2->AddPlacedVolume(coatingrect2log,m25,r25);



	r27=new G4RotationMatrix();
	r27->rotateZ(M_PI);
	r27->rotateX(M_PI/2);
	//m27=G4ThreeVector(-rcmb,rc1+rcmt,-rc2-rcma);
	m27=G4ThreeVector(sl/2+t1-sl-rec3-rec4,3*(sl/(4*cos(45*deg)))+sl+t3+(rec1),(-sl*(1+0.25*sqrt(2)))-t2+(rec2));
	frectangle2->AddPlacedVolume(coatingrect2log,m27,r27);

	r29=new G4RotationMatrix();	
	r29->rotateZ(M_PI);
	r29->rotateY(-M_PI/4);
	r29->rotateX(M_PI/2);
	m29=G4ThreeVector(sl/2+t1+drec,3*(sl/(4*cos(45*deg)))+sl+t3+rec1-hrec,(-sl*(1+0.25*sqrt(2)))-t2+(rec2));
	frectangle2->AddPlacedVolume(coatingrect2log,m29,r29);


	r31=new G4RotationMatrix();
	r31->rotateZ(M_PI/2);
	r31->rotateY(3*M_PI/2);
	m31=G4ThreeVector(3*(sl/(4*cos(45*deg)))+sl+t3+rec1,sl/2+t1,-sl*(1+0.25*sqrt(2))-t2+rec2);

	frectangle2->AddPlacedVolume(coatingrect2log,m31,r31);



	r33=new G4RotationMatrix();
	r33->rotateZ(0);
	r33->rotateY(-M_PI/4);
	r33->rotateX(-M_PI/2);
	m33=G4ThreeVector(3*(sl/(4*cos(45*deg)))+t3+sl+rec1-hrec,-sl/2-t1-drec,-sl*(1+0.25*sqrt(2))-t2+rec2);
	frectangle2->AddPlacedVolume(coatingrect2log,m33,r33);


	r35=new G4RotationMatrix();	
	r35->rotateZ(0);
        r35->rotateX(3*M_PI/2);
	m35=G4ThreeVector(sl/2+t1,3*-sl/(4*cos(45*deg))-sl-t3-rec1,-sl*(1+0.25*sqrt(2))-t2+rec2);
	frectangle2->AddPlacedVolume(coatingrect2log,m35,r35);

	r37=new G4RotationMatrix();
	r37->rotateZ(0);
	r37->rotateY(M_PI/4);
        r37->rotateX(-M_PI/2);
	m37=G4ThreeVector(-(sl/2+t1)-drec,3*-sl/(4*cos(45*deg))-sl-t3-rec1+hrec,-sl*(1+0.25*sqrt(2))-t2+rec2);
	frectangle2->AddPlacedVolume(coatingrect2log,m37,r37);


	r39=new G4RotationMatrix();	
	r39->rotateZ(3*M_PI/2);
	r39->rotateY(M_PI/2);
	m39=G4ThreeVector(-3*(sl/(4*cos(45*deg)))-sl-t3-rec1,-(sl/2+t1),-sl*(1+0.25*sqrt(2))-t2+rec2);
	frectangle2->AddPlacedVolume(coatingrect2log,m39,r39);


	// r42=new G4RotationMatrix();
	// r42->rotateZ(0);
	// r42->rotateY(M_PI/4);
	// r42->rotateX(M_PI/2);
	// m42=G4ThreeVector(sl/2+t1-sl-drec-rec5,3*(sl/(4*cos(45*deg)))+sl+t3+(rec1)-hrec,(sl*(1+0.25*sqrt(2)))+t2-rec2);
	// frectangle2->AddPlacedVolume(coatingrect2log,m42,r42);

	// r44=new G4RotationMatrix();
	// r44->rotateZ(0);
	// r44->rotateX(M_PI/2);
	// m44=G4ThreeVector(sl/2+t1+rec3-rec4,3*(sl/(4*cos(45*deg)))+sl+t3+rec1,(sl*(1+0.25*sqrt(2)))+t2-(rec2));
	// frectangle2->AddPlacedVolume(coatingrect2log,m44,r44);

	// r46=new G4RotationMatrix();
	// r46->rotateZ(0);
	// r46->rotateY(-M_PI/4);
	// r46->rotateX(M_PI/2);
	// m46=G4ThreeVector(3*(sl/(4*cos(45*deg)))+sl+t3+rec1-hrec,sl/2+t1+drec,sl*(1+0.25*sqrt(2))+t2-rec2);

	// frectangle2->AddPlacedVolume(coatingrect2log,m46,r46);



	// r48=new G4RotationMatrix();
	// r48->rotateZ(3*M_PI/2);
	// r48->rotateY(3*M_PI/2);
	// m48=G4ThreeVector(3*(sl/(4*cos(45*deg)))+t3+sl+rec1,-sl/2-t1,sl*(1+0.25*sqrt(2))+t2-rec2);
	// frectangle2->AddPlacedVolume(coatingrect2log,m48,r48);


	// r50=new G4RotationMatrix();
	// r50->rotateZ(M_PI);
	// r50->rotateY(-M_PI/4);
    //     r50->rotateX(-M_PI/2);
	// m50=G4ThreeVector(sl/2+t1+drec,3*-sl/(4*cos(45*deg))-sl-t3-rec1+hrec,sl*(1+0.25*sqrt(2))+t2-rec2);
	// frectangle2->AddPlacedVolume(coatingrect2log,m50,r50);



	// r52=new G4RotationMatrix();
	// r52->rotateZ(M_PI);
    //     r52->rotateX(3*M_PI/2);
	// m52=G4ThreeVector(-(sl/2+t1),3*-sl/(4*cos(45*deg))-sl-t3-rec1,sl*(1+0.25*sqrt(2))+t2-rec2);
	// frectangle2->AddPlacedVolume(coatingrect2log,m52,r52);



	// r54=new G4RotationMatrix();
	// r54->rotateZ(M_PI);
	// r54->rotateY(M_PI/4);
	// r54->rotateX(-M_PI/2);
	// m54=G4ThreeVector(-3*(sl/(4*cos(45*deg)))-sl-t3-rec1+hrec,-(sl/2+t1)-drec,sl*(1+0.25*sqrt(2))+t2-rec2);
	// frectangle2->AddPlacedVolume(coatingrect2log,m54,r54);



	// r56=new G4RotationMatrix();
	// r56->rotateZ(M_PI/2);
	// r56->rotateY(M_PI/2);
	// m56=G4ThreeVector(-3*(sl/(4*cos(45*deg)))-sl-t3-rec1,+sl/2+t1,sl*(1+0.25*sqrt(2))+t2-rec2);
	// frectangle2->AddPlacedVolume(coatingrect2log,m56,r56);
}










///////////////////////////////////////////////////////////////////////
////////////////Class to construct the square tiles////////////////////
///////////////////////////////////////////////////////////////////////
//void DetectionSystemAries::ConstructDome(){
//}











//square 'full' of coating
   // G4VSolid* coatingsquarescintSolid = new G4Box("coatingsquarescientSolid", scintsquarelength + sivlerthickness, scintsquarelength + sivlerthickness, scintthickness + sivlerthickness);

   // G4LogicalVolume* coatingsquarescintLog  = new G4LogicalVolume( squarescintSolid, Ag,"coatingsquarescintLog");

//just the coating outer layer 'hollow'
   // G4VSolid* coatingtilesol =new G4SubtractionSolid subtract("coatingtilesol", &squarescientSolid , &coatingsquarescientSolid );  //still need to add rotationmatrix and positionvector 
    
    //G4LogicalVolume* justcoatingsquarescintLog  = new G4LogicalVolume( coatingtilesol, Ag,"justcoatingsquarescintLog");


//solid and logical volume for slightly smaller tile


//placing the tile in the coating (making it a daughter volume)
   // G4PhysicalVoume* coatedtilePhys = new G4Vplacement (rotateNull,moveNull,)   
    
//volume for the tile+its coating   
//still need to add rotationmatrix and positionvector 
   // G4VSolid* coatedtilesol = new G4UnionSolid("coatedtilesol",coatingtilesol,squarescintSolid);

//making it tapered by creating paralelipiped and substracting it to tile  
   // G4VSolid* batonsol = new G4Para("coatedtilesol",scintsquarelength,,,,);







//repeat process for triangle 



//Define SiPm using either G4SiPM online method or using max Winokan method (latter easier maybe less precise?)



//place SiPMs using subtractionvolume on the tile (subtract the layer of silver) all in corner of square 



//creating rotationmatrix and g4three vector to later place each tile with physical volume


//	G4RotationMatrix* rotation1=new G4RotationMatrix;
//	rotation1->rotateX
//	rotation1->rotateZ


 


  


  //create a for loop that places each tile using its number and rotation and postition number

    //G4VPhysicalVolume* simpletilePhys = new G4PVPlacement(0,G4ThreeVector(),squarescintLog,"simpletilePhys",worldLog,0,false); 

//	pos = G4ThreeVector(TransX(x,y,z,theta,phi),TransY(x,y,z,theta,phi),TransZ(x,z,theta));
//	x = scintsquarelength/2 +silverthickness; y = scintsquarelength/2 + silverthickness; z =scintthickness/2 + silverthickness;
///////////////////////////////////////////////////////////////////////
///////////////////Class to construct the tiles////////////////////////
///////////////////////////////////////////////////////////////////////
//void DetectionSystemAries::ConstructTile(){
//
//	//retrieving materials
//        G4Material* tilemat = G4Material::GetMaterial("BC422Q");
//	G4Material* coatingmat = G4Material::GetMaterial("Ag")
//    
//
//	//square scint uncoated tile definition(Volume and visualization)
//    	G4VSolid* squarescintSolid = new G4Box("squarescintSolid", scintsquarelength/2, scintsquarelength/2, scintthickness/2);
//    	G4LogicalVolume* squarescintLog = new G4LogicalVolume(squarescintSolid, tilemat,"squarescintLog");
//    	G4VisAttributes* tilevis = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
//        tilevis->SetVisibility(true);
//	squarescintLog->SetVisAttributes(tilevis);
//
//
//
//
//	//placing
//	ftiles->AddPlacedVolume(squarescintLog,moveNull,rotateNull);
//}
//

	
//	G4RotationMatrix* noRotate=new G4RotationMatrix;
//	noRotate->rotateX(M_PI/2.0);
//      noRotate->rotateX(0);
//      noRotate->rotateY(0);
//      noRotate->rotateZ(0);
//
//
//
//
//
////	// subtracting tile from block of coating
//	G4VSolid* justcoatingsol         = new G4SubtractionSolid("justcoatingsol",squarecoatingsolid,squarescintsolid,rotateNull,moveNull);
//	G4LogicalVolume* justcoatinglog  = new G4LogicalVolume(justcoatingsol,coatingmat,"justcoatinglog");
//	justcoatinglog->SetVisAttributes(coatingvis);
//
//G4VSolid* squarecoatingsolid      = new G4Trap("squarecoatingsolid",longthickness,taperanglex,taperangley,longslength,longslength,longslength,zeroangle,longslength-scintthickness/4,longslength-scintthickness/4,longslength-scintthickness/4,zeroangle);	
//	G4LogicalVolume* trianglecoatinglog2 =new G4LogicalVolume(trianglecoatingsolid2,coatingmat,"trianglecoatinglog2");	

	
