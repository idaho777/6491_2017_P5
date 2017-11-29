import java.util.*;

// 2D Delaunay Triangulation ===========================================================================
Triangle[] getDelaunayTriangulation2D(pts points) {
  if (points.nv < 3) return new Triangle[0]; 
  float hite = points.G[0].z;
  List<Triangle> triangulation = new ArrayList<Triangle>();
  pt AA = P(-1000000, -1000000, hite);
  pt BB = P(0, 1000000, hite);
  pt CC = P(1000000, -1000000, hite);
  triangulation.add(TRI(AA, BB, CC));
  
  // insert each point
  for (int i = 0; i < points.nv; ++i) {
    pt curr = points.G[i];
    
    List<Triangle> insiders = getInsiders(triangulation, curr);         // Calculate which triangle circumcenters contain curr
    List<Edge> boundary = getBoundaryAroundCircumCenter(insiders);    // Get boundary points
    triangulation.removeAll(insiders);                                  // remove old triangles
    for (Edge e : boundary) { triangulation.add(TRI(e.A, e.B, curr)); } // add new triangles
  }
  
  // remove parent triangle
  int triInd = 0;
  while (triInd < triangulation.size()) {
    Triangle triToRem = triangulation.get(triInd);
    if (triToRem.contains(AA) || triToRem.contains(BB) || triToRem.contains(CC)) {
      triangulation.remove(triInd); 
    } else {
      triInd++;
    }
  }
  
  return triangulation.toArray(new Triangle[triangulation.size()]);
}


List<Triangle> getInsiders(List<Triangle> triangulation, pt C) {
  // Calculate which triangle circumcenters contain curr
  List<Triangle> insiders = new ArrayList<Triangle>();
  for (Triangle t : triangulation) {
    if (t.pointInsideCircum(C)) {
      insiders.add(t);
    }
  }
  
  return insiders;
}


List<Edge> getBoundaryAroundCircumCenter(List<Triangle> tri) {
  List<Edge> edges = new ArrayList<Edge>();
  
  for (int i = 0; i < tri.size(); ++i) {
    Triangle currTri = tri.get(i);
    List<Edge> currTriEdges = new ArrayList<Edge>();
    currTriEdges.add(EDGE(currTri.A, currTri.B));
    currTriEdges.add(EDGE(currTri.B, currTri.C));
    currTriEdges.add(EDGE(currTri.C, currTri.A));
    
    // Go through each triangle and attempt to remove all edges from the edge list
    for (int j = 0; j < tri.size(); ++j) {
      if (j == i) continue;
      
      int ind = -1;
      Triangle tempTri = tri.get(j);
      for (int ei = 0; ei < currTriEdges.size(); ei++) {
        if (tempTri.contains(currTriEdges.get(ei))) {
          ind = ei;
          break;
        }
      }
      
      if (ind >= 0) currTriEdges.remove(ind);
    }
    edges.addAll(currTriEdges);
  }
  
  return edges;
}


// Delauntay Triangulation 2nd method ================================================================================================================
// This triangulation is proper and doesn't utilize the floor/ceiling plane heuristic
Tetrahedron[] delaunayTetra(List<pt> points) {
  List<Tetrahedron> ret = new ArrayList<Tetrahedron>();
  
  vec V = V(10000, 0, -4000);
  
  pt AA = P(P(1000, 1000), V); V = R(V, TAU/3, V(1, 0, 0), V(0, 1, 0));
  pt BB = P(P(1000, 1000), V); V = R(V, TAU/3, V(1, 0, 0), V(0, 1, 0));
  pt CC = P(P(1000, 1000), V);
  pt DD = P(0, 0, 14000);
  ret.add(TETRA(AA, BB, CC, DD));
  
  for (pt p : points) {
    List<Tetrahedron> insiders = getInsidersTetra(ret, p);
    List<Triangle> boundary = getBoundaryAroundCircumCenter3D(insiders);
    ret.removeAll(insiders);
    for (Triangle t : boundary) { ret.add(TETRA(t.A, t.B, t.C, p)); }
  }
  
  //remove parent Tetrahedrons
  int tetInd = 0;
  while (tetInd < ret.size()) {
    Tetrahedron toRem = ret.get(tetInd);
    if (toRem.contains(AA) || toRem.contains(BB) || toRem.contains(CC) || toRem.contains(DD)) {
      ret.remove(tetInd); 
    } else {
      tetInd++; 
    }
  }
  
  return ret.toArray(new Tetrahedron[ret.size()]);
}


List<Tetrahedron> getInsidersTetra(List<Tetrahedron> tets, pt C) {
  // Calculate which triangle circumcenters contain curr
  List<Tetrahedron> insiders = new ArrayList<Tetrahedron>();
  for (Tetrahedron t : tets) {
    if (t.pointInsideCircum(C)) {
      insiders.add(t);
    }
  }
  return insiders;
}


List<Triangle> getBoundaryAroundCircumCenter3D(List<Tetrahedron> tets) {
  List<Triangle> boundary = new ArrayList<Triangle>();
  
  for (int i = 0; i < tets.size(); ++i) {
    Tetrahedron currTri = tets.get(i);
    List<Triangle> currTriFaces = new ArrayList<Triangle>();
    currTriFaces.add(TRI(currTri.A, currTri.B, currTri.C));
    currTriFaces.add(TRI(currTri.A, currTri.B, currTri.D));
    currTriFaces.add(TRI(currTri.A, currTri.D, currTri.C));
    currTriFaces.add(TRI(currTri.D, currTri.B, currTri.C));
    
    // Find which currTriFaces are shared with other tetrahedrons and remove them from boundary
    // In the end, boundary should only have non-sharing faces
    for (int j = 0; j < tets.size(); ++j) {
      if (i == j) continue;
      
      int ind = -1;
      Tetrahedron tempTet = tets.get(j);
      for (int ti = 0; ti < currTriFaces.size(); ti++) {
        if (tempTet.contains(currTriFaces.get(ti))) {
          ind = ti;
          break;
        }
      } 
      if (ind >= 0) { currTriFaces.remove(ind); }
    }
    boundary.addAll(currTriFaces);
  }
  
  return boundary;
}


// ===================================================================================================================
// Takes both sets of triangles and plane points
Tetrahedron[] getDelaunayTetra(Triangle[] topTri, Triangle[] botTri, pts topPlane, pts botPlane) {
  List<Tetrahedron> ret = new ArrayList<Tetrahedron>();
  
  // Calculate Tetrahedrons from bot up
  // Calculate Tetrahedrons from top down
  // Calculate tetrahedron from edge pairs
  List<Tetrahedron> botTriTopPoints = getDelaunayTetraPlanePoints(botTri, topPlane);
  List<Tetrahedron> topTriBotPoints = getDelaunayTetraPlanePoints(topTri, botPlane);
  List<Tetrahedron> pairEdge = getDelaunayTetraEdgeEdge(topTri, botTri, topPlane, botPlane);
  ret.addAll(botTriTopPoints);
  ret.addAll(topTriBotPoints);
  ret.addAll(pairEdge);
  
  Tetrahedron[] retArr = ret.toArray(new Tetrahedron[ret.size()]);
  return retArr;
}


List<Tetrahedron> getDelaunayTetraPlanePoints(Triangle[] tri, pts P) {
  List<Tetrahedron> tets = new ArrayList<Tetrahedron>();

  // For all triangles, create tetra with all opposite points
  // For each tetra, check if it excludes all other points.
  for (Triangle t : tri) {
    for (int i = 0; i < P.nv; ++i) {
      Tetrahedron currTet = TETRA(t.A, t.B, t.C, P.G[i]);
      boolean excludes = true;
      for (int j = 0; j < P.nv; ++j) {
        if (i == j) continue;
        if (currTet.pointInsideCircum(P.G[j])) {
          excludes = false;
          break;
        }
      }
      
      if (excludes) {
        tets.add(currTet); 
      }
    }
  }
  
  return tets;
}


List<Tetrahedron> getDelaunayTetraEdgeEdge(Triangle[] topTri, Triangle[] botTri, pts top, pts bot) {
  List<Tetrahedron> tets = new ArrayList<Tetrahedron>();
    
  List<Edge> topEdge = new ArrayList<Edge>();
  if (top.nv == 2) {
    topEdge.add(EDGE(top.G[0], top.G[1]));
  } else if (top.nv >= 3) {
    for (Triangle t : topTri) {
      Edge E;
      E = EDGE(t.A, t.B); if (!topEdge.contains(E)) topEdge.add(E);
      E = EDGE(t.B, t.C); if (!topEdge.contains(E)) topEdge.add(E);
      E = EDGE(t.C, t.A); if (!topEdge.contains(E)) topEdge.add(E);
    }
  }
  
  List<Edge> botEdge = new ArrayList<Edge>();
  if (bot.nv == 2) {
    botEdge.add(EDGE(bot.G[0], bot.G[1]));
  } else if (bot.nv >= 3) {
    for (Triangle t : botTri) {
      Edge E;
      E = EDGE(t.A, t.B); if (!botEdge.contains(E)) botEdge.add(E);
      E = EDGE(t.B, t.C); if (!botEdge.contains(E)) botEdge.add(E);
      E = EDGE(t.C, t.A); if (!botEdge.contains(E)) botEdge.add(E);
    }
  }
      
  for (Edge et : topEdge) {
    for (Edge eb : botEdge) {
      Tetrahedron currTet = TETRA(et.A, et.B, eb.A, eb.B);
      boolean excludes = true;
      for (int i = 0; i < top.nv; ++i) if (!currTet.contains(top.G[i]) && currTet.pointInsideCircum(top.G[i])) { excludes = false; break; }
      for (int i = 0; i < bot.nv; ++i) if (!currTet.contains(bot.G[i]) && currTet.pointInsideCircum(bot.G[i])) { excludes = false; break; }
      if (excludes) tets.add(currTet);
    }
  }
  
  return tets;
}


// ====================================================================================================================
// Get beam edges from tetrahedrons
// no unique edge beams
List<Edge> getBeams(Tetrahedron[] tets) {
  List<Edge> beams = new ArrayList<Edge>();
  for (Tetrahedron t : tets) {
    Edge E;
    //E = EDGE(t.A, t.B); if (!beams.contains(E)) beams.add(E);break;
    E = EDGE(t.A, t.B); if (!beams.contains(E)) beams.add(E);
    E = EDGE(t.B, t.C); if (!beams.contains(E)) beams.add(E);
    E = EDGE(t.C, t.A); if (!beams.contains(E)) beams.add(E);
    E = EDGE(t.A, t.D); if (!beams.contains(E)) beams.add(E);
    E = EDGE(t.B, t.D); if (!beams.contains(E)) beams.add(E);
    E = EDGE(t.C, t.D); if (!beams.contains(E)) beams.add(E);
  }
  
  return beams;
}

// Sample points on beams and end points
List<List<pt>> samplePointCloud(List<Edge> beams, pts R, pts S) {
  if (g_pointCloud == null) {
    g_pointCloud = new ArrayList<List<pt>>();
    g_pointCloud.addAll(sampleBeamCloud(beams));
    g_pointCloud.addAll(sampleBallCloud(R, S));
  }
  return g_pointCloud;
}


// Sample on beams
List<List<pt>> sampleBeamCloud(List<Edge> beams) {
  List<List<pt>> pointCloud = new ArrayList<List<pt>>();
  for (Edge b : beams) {
    pointCloud.add(b.samplePoints()); 
  }
  return pointCloud;
}


// Sample on Beam Ends
List<List<pt>> sampleBallCloud(pts A, pts B) {
  List<List<pt>> pointCloud = new ArrayList<List<pt>>();
  for (int i = 0; i < A.nv; ++i) { pointCloud.add(sampleBallPoints(A.G[i])); }
  for (int i = 0; i < B.nv; ++i) { pointCloud.add(sampleBallPoints(B.G[i])); }
  return pointCloud;
}


List<pt> sampleBallPoints(pt p) {
  List<pt> points = new ArrayList<pt>();
  
  vec I = V(1, 0, 0);
  vec J = V(0, 1, 0);
  vec K = V(0, 0, 1);
  vec rad = V(rb, K);
  
  
  float downSteps = PI*rb/edgeStepLen;
  float downAngle = PI/downSteps; 
  points.add(P(p, rad));
  points.add(P(p, M(rad)));
  boolean on = false;
  for (int i = 0; i < downSteps; ++i) { // Top Down
    rad = R(rad, downAngle, K, I);
    
    float r = sqrt(n2(rad) - sq(dot(rad, K)));
    float rStep = 2*PI*r/edgeStepLen;
    float rAngle = 2*PI/rStep;
    for (int j = 0; j < rStep; ++j) {
      points.add(P(p, R(rad, j*rAngle - ((on) ? rAngle/2 : 0), I, J)));
    }
    on =! on;
  }
  
  return points;
}


// Ball Rolling =========================================================================================================
List<Triangle> queue;
List<Triangle> ballRollPointCloud(List<List<pt>> pointCloud) {
  ballRollAnimation(pointCloud);
  
  // Main Loop.  Do BFS.  DFS contains holes
  while (!queue.isEmpty()) {
    ballRollAnimationStep();
  }
  println("Ball Roll Complete");
  return triangulation;
}


void ballRollAnimation(List<List<pt>> pointCloud) {
  println("BALL ROLL POINT CLOUD");
  println("This should take couple of seconds...");
  voxelSpace = setupVoxelSpace(pointCloud);
  triangulation = new ArrayList<Triangle>();
  
  int index = int(pointCloud.get(0).size()/2/circumPts);
  index *= circumPts;
  pt AA = pointCloud.get(0).get(index+4);
  pt BB = pointCloud.get(0).get(index+3);
  pt CC = pointCloud.get(0).get(index+4+circumPts);
  pt DD = pointCloud.get(0).get(index+3+circumPts);
  Triangle currTri = TRI(AA, BB, CC);
  
  if (d(AA, circumcenter2D(AA, BB, CC)) > d(AA, circumcenter2D(AA, BB, DD)) ){
    currTri = TRI(AA, BB, DD);
  }
  
  queue = new ArrayList<Triangle>();
  queue.add(currTri);
}


void ballRollAnimationStep() {
  if (queue.isEmpty()) { return; }
  Triangle currTri = queue.get(0); queue.remove(0);
  //Triangle currTri = queue.get(queue.size()-1); queue.remove(queue.size()-1);
    
  pt AA = currTri.A;
  pt BB = currTri.B;
  pt CC = currTri.C;
  
  if (triangulation.contains(currTri)) return;
  triangulation.add(currTri);

  Triangle leftTri = getBallRollTriangle(AA, CC, currTri, voxelSpace);
  if (leftTri  != null && !queue.contains(leftTri )) queue.add(leftTri);

  Triangle rightTri = getBallRollTriangle(CC, BB, currTri, voxelSpace);
  if (rightTri != null && !queue.contains(rightTri)) queue.add(rightTri);
}


Triangle getBallRollTriangle(pt L, pt R, Triangle tri, VoxelSpace voxelSpace) {
  List<pt> candidatePoints = new ArrayList<pt>();
  
  pt triCC2D = circumcenter2D(tri.A, tri.B, tri.C);
  pt triCC3D = P(triCC2D, sqrt(sq(ballRadius) - sq(d(triCC2D, tri.A))), tri.N());
  pt LRmid = P(L, 0.5, R);
  
  pt cc2DVoxIndex = P(voxelSpace.getX(triCC2D), voxelSpace.getY(triCC2D), voxelSpace.getZ(triCC2D));
  
  List<pt> voxelPoints = new ArrayList<pt>();
  voxelPoints.addAll(voxelSpace.getVoxel(triCC2D));
  for (int x = -1; x <=1; ++x) {
    for (int y = -1; y <=1; ++y) {
      for (int z = -1; z <=1; ++z) {
        pt b = P(triCC2D, V(x*ballRadius*3, y*ballRadius*3, z*ballRadius*3));
        pt bIndex = P(voxelSpace.getX(b), voxelSpace.getY(b), voxelSpace.getZ(b));
        if (!cc2DVoxIndex.equals(bIndex)) {
          voxelPoints.addAll(voxelSpace.getVoxel(b));
        }
      }
    }  
  }
    
  // Get candidate points
  // Find candidate points the ball can roll on
  for (pt p : voxelPoints) {
    if (tri.contains(p)) continue;
    pt cc2D = circumcenter2D(L, R, p);
    if (d(cc2D, LRmid) < d(triCC3D, LRmid)) {
      Triangle t = TRI(L, R, p);
      if (abs(dot(U(V(t.A, t.B)), U(V(t.B, t.C)))) == 1) continue;
      candidatePoints.add(p);
    }
  }
  
  // get points without other points inside 
  List<pt> betterPoints = new ArrayList<pt>();
  for (pt p : candidatePoints) {    
    // check if ball contains any other points
    Triangle t = TRI(L, R, p);
    pt cc2D = circumcenter2D(L, R, p);
    pt cc3D = P(cc2D, sqrt(sq(ballRadius) - sq(d(cc2D, L))), t.N());
    Tetrahedron tet = TETRA(L, R, p, P(cc3D, ballRadius, t.N()));
    
    boolean good = true;
    for (pt pp : candidatePoints) {
      if (tet.contains(pp)) continue;
      if (ballRadius > d(cc3D, pp)) {
        //println(d(cc3D, p), d(cc3D, pp));
        good = false;
        break;
      }
    }
    
    if (good) { 
      betterPoints.add(p);
    }
  }
  if (betterPoints.isEmpty()) return null;
  
  // For all candidate points, find which makes the ball roll the smallest angle
  vec refVec = U(V(LRmid, triCC3D));
  pt minPoint = null;;
  float minAngle = TAU;
  for (pt p : betterPoints) {
    Triangle t = TRI(L, R, p);
    pt cc2D = circumcenter2D(L, R, p);
    pt cc3D = P(cc2D, sqrt(sq(ballRadius) - sq(d(cc2D, p))), t.N());
    vec tempVec = U(V(LRmid, cc3D));
    
    float ang = angle(refVec, tempVec, V(L,R));
    if (dot(U(t.N()), U(tri.N())) == -1)  continue;// If this is flipped upside down
    if (dot(U(tempVec), U(refVec)) == -1) continue;
    
    if (ang < minAngle) {
      minAngle = ang;
      minPoint = p;
    }
  }
  if (minPoint == null) return null;
  if (minAngle >= PI) return null;

  return TRI(L, R, minPoint);
}