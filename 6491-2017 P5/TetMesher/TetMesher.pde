//  ******************* LITM: Layer-Interpolating Tet Mesh, 2017 ***********************
Boolean 
  animating=true, 
  PickedFocus=false, 
  center=true, 
  track=false, 
  showViewer=false, 
  showBalls=true, 
  showControl=true, 
  showCurve=true, 
  showPath=true, 
  showKeys=true, 
  showSkater=false, 
  scene1=false,
  solidBalls=false,
  showCorrectedKeys=true,
  showQuads=true,
  showVecs=true,
  showTube=true,
  showTriEdge=true,
  waterTight=false,
  waterAnimate=false,
  flipped = false;
float 
 h_floor=0, h_ceiling=600,  h=h_floor,
  t=0, 
  s=0,
  rb=40, rt=rb/2; // radius of the balls and tubes
  
int
  f=0, maxf=2*30, level=4, method=5;
String SDA = "angle";
float defectAngle=0;
pts P = new pts(); // polyloop in 3D
pts Q = new pts(); // second polyloop in 3D
pts R, S; 
    

void setup() {
  myFace = loadImage("data/pic.jpg");  // load image from file pic.jpg in folder data *** replace that file with your pic of your own face
  textureMode(NORMAL);          
  size(1600, 900, P3D); // P3D means that we will do 3D graphics
  P.declare(); Q.declare(); // P is a polyloop in 3D: declared in pts
  P.loadPts("data/pts");  
  Q.loadPts("data/pts2"); // loads saved models from file (comment out if they do not exist yet)
  noSmooth();
  frameRate(30);
  R=P; S=Q;
  }

int circumPts = 16;
int aniSpeed = 100;
float edgeStepLen;
float ballRadius = 0;

List<List<pt>> g_pointCloud;
List<Triangle> triangulation;
VoxelSpace voxelSpace;

void draw() {
  background(255);
  hint(ENABLE_DEPTH_TEST); 
  pushMatrix();   // to ensure that we can restore the standard view before writing on the canvas
  setView();  // see pick tab
  showFloor(h); // draws dance floor as yellow mat
  doPick(); // sets Of and axes for 3D GUI (see pick Tab)
  R.SETppToIDofVertexWithClosestScreenProjectionTo(Mouse()); // for picking (does not set P.pv)
  
  edgeStepLen = 2*PI*rt/float(circumPts);
  ballRadius = 0.7*max(edgeStepLen, TWO_PI*rt/float(circumPts));
  // Get Delaunay Triangulation from both planes
  Triangle[] P_DELAUNAY = getDelaunayTriangulation2D(P);
  Triangle[] Q_DELAUNAY = getDelaunayTriangulation2D(Q);
  
  // Find tetrahedrons
  Tetrahedron[] tets = getDelaunayTetra(Q_DELAUNAY, P_DELAUNAY, Q, P);
  //Tetrahedron[] tets = delaunayTetra(allPoints);  // Slower but proper delaunay
  
  // Beams
  List<Edge> g_beams = getBeams(tets);
  
  // Get sample points from Beams and beam ends
  if (waterTight) {
    waterAnimate = false;
    if (g_pointCloud == null) {
      g_pointCloud = samplePointCloud(g_beams, R, S);
      triangulation = ballRollPointCloud(g_pointCloud);
    }
  } else if (waterAnimate) {
    waterTight = false;
    if (g_pointCloud == null) {
      g_pointCloud = samplePointCloud(g_beams, R, S);
      ballRollAnimation(g_pointCloud);
    } else {
      for (int i = 0; i < aniSpeed; i++) {
        ballRollAnimationStep(); 
      }
    }
  } else {
    g_pointCloud = null;
    triangulation = null;
  }
  
  
  noFill();
  noStroke();  
  if(showTube) {
    fill(cyan);
    for (Edge e : g_beams) { e.draw(); }
  }
  
  if(showBalls) {
    fill(orange); P.drawBalls(rb);
    fill(green); Q.drawBalls(rb);  
    fill(red,100); R.showPicked(rb+5); 
  }
  
  //if (g_pointCloud != null) {
  //  for (List<pt> l : g_pointCloud) {
  //    for(pt p : l) {
  //      noStroke();
  //      fill(yellow);
  //      show(p, 1);
  //    }
  //  }
  //}
  
  if ((waterTight || waterAnimate) && triangulation != null) {
    for(Triangle t : triangulation) {
      fill(yellow);
      if (showTriEdge) {
        stroke(orange);
        strokeWeight(1);
      } else {
        noStroke();
      }
      t.drawShape();
    }
  }
    
  popMatrix(); // done with 3D drawing. Restore front view for writing text on canvas
  hint(DISABLE_DEPTH_TEST); // no z-buffer test to ensure that help text is visible
  // used for demos to show red circle when mouse/key is pressed and what key (disk may be hidden by the 3D model)
  if(mousePressed) {stroke(cyan); strokeWeight(3); noFill(); ellipse(mouseX,mouseY,20,20); strokeWeight(1);}
  if(keyPressed) {stroke(red); fill(white); ellipse(mouseX+14,mouseY+20,26,26); fill(red); text(key,mouseX-5+14,mouseY+4+20); strokeWeight(1); }
  if(scribeText) {fill(black); displayHeader();} // dispalys header on canvas, including my face
  if(scribeText && !filming) displayFooter(); // shows menu at bottom, only if not filming
  if(filming && (animating || change)) saveFrame("FRAMES/F"+nf(frameCounter++,4)+".tif");  // save next frame to make a movie
  change=false; // to avoid capturing frames when nothing happens (change is set uppn action)
  change=true;
  }