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
  waterTight=false,
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


List<List<pt>> g_pointCloud;
void draw() {
  println("loop");
  background(255);
  hint(ENABLE_DEPTH_TEST); 
  pushMatrix();   // to ensure that we can restore the standard view before writing on the canvas
  setView();  // see pick tab
  showFloor(h); // draws dance floor as yellow mat
  doPick(); // sets Of and axes for 3D GUI (see pick Tab)
  R.SETppToIDofVertexWithClosestScreenProjectionTo(Mouse()); // for picking (does not set P.pv)

  // Get Delaunay Triangulation from both planes
  Triangle[] P_DELAUNAY = getDelaunayTriangulation2D(P);
  Triangle[] Q_DELAUNAY = getDelaunayTriangulation2D(Q);
  
  // Find tetrahedrons
  Tetrahedron[] tets = getDelaunayTetra(Q_DELAUNAY, P_DELAUNAY, Q, P);
  //Tetrahedron[] tets = delaunayTetra(allPoints);  // Slower but proper delaunay
  
  // Beams
  List<Edge> g_beams = getBeams(tets);
  
  // Get sample points from Beams and beam ends
  List<Triangle> triangulation = null;
  if (waterTight) {
    g_pointCloud = samplePointCloud(g_beams, R, S);
    triangulation = ballRollPointCloud(g_pointCloud);
  } else {
    g_pointCloud = null;
  }
  println(R.G[0]);
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
  //  noStroke();
  //  fill(magenta);
  //  for (List<pt> l : g_pointCloud) {
  //    for (pt p : l) {
  //      show(p, 1);
  //    }
  //  }
  //}
  
  if (waterTight && triangulation != null) {
    for(Triangle t : triangulation) {
      fill(yellow);
      stroke(orange);
      strokeWeight(2);
      t.drawShape();
    }
  }
    
  popMatrix(); // done with 3D drawing. Restore front view for writing text on canvas
  hint(DISABLE_DEPTH_TEST); // no z-buffer test to ensure that help text is visible

  //*** TEAM: please fix these so that they provice the correct counts
  scribeHeader("Site count: "+3+" floor + "+7+" ceiling",1);
  scribeHeader("Beam count: "+3+" floor + "+7+" ceiling +"+6+" mixed",2);
  scribeHeader("Tet count: "+20,3);
 
  // used for demos to show red circle when mouse/key is pressed and what key (disk may be hidden by the 3D model)
  if(mousePressed) {stroke(cyan); strokeWeight(3); noFill(); ellipse(mouseX,mouseY,20,20); strokeWeight(1);}
  if(keyPressed) {stroke(red); fill(white); ellipse(mouseX+14,mouseY+20,26,26); fill(red); text(key,mouseX-5+14,mouseY+4+20); strokeWeight(1); }
  if(scribeText) {fill(black); displayHeader();} // dispalys header on canvas, including my face
  if(scribeText && !filming) displayFooter(); // shows menu at bottom, only if not filming
  if(filming && (animating || change)) saveFrame("FRAMES/F"+nf(frameCounter++,4)+".tif");  // save next frame to make a movie
  change=false; // to avoid capturing frames when nothing happens (change is set uppn action)
  change=true;
  }