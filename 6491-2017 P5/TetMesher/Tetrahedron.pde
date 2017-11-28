      Tetrahedron TETRA(pt A, pt B, pt C, pt D) { return new Tetrahedron(A, B, C, D); }

class Tetrahedron {
  pt A, B, C, D;
  Tetrahedron() {}
  Tetrahedron(pt A, pt B, pt C, pt D) {
     this.A = A;
     this.B = B;
     this.C = C;
     this.D = D;
  }
  
  
  boolean pointInsideCircum(pt Z) {
    if (contains(Z)) return false;
    pt ccenter = circumcenter3D(A, B, C, D);
    return d(ccenter, Z) <= d(ccenter, A);
  }
  
  
  boolean contains(Triangle t) {
    return contains(t.A) && contains(t.B) && contains(t.C);
  }
  
  
  boolean contains(pt p) {
    return p.equals(A) || p.equals(B) || p.equals(C) || p.equals(D); 
  }
  
  void draw() {
    noStroke();
    beam(A, B, rt);
    beam(B, C, rt);
    beam(C, A, rt);
    beam(A, D, rt);
    beam(B, D, rt);
    beam(C, D, rt);
  }
  
  void drawSphere() {
    noStroke();
    pt cc = circumcenter3D(A, B, C, D);
    pushMatrix();
      translate(cc.x, cc.y, cc.z);
      fill(color(150, 100));
      sphere(d(A, cc));
    popMatrix(); 
  }  
  
  
  @Override
  public boolean equals(Object o) {
    if (!(o instanceof Tetrahedron)) return false;
    
    Tetrahedron obj = (Tetrahedron) o;
    return contains(obj.A) && contains(obj.B) && contains(obj.C) && contains(obj.D);
  }
}


pt circumcenter3D (pt A, pt B, pt C, pt D) {
  pt  CA = circumcenter2D(A, B, C);
  vec NA = cross(V(A, B), V(A, C)); 
  if (dot(NA, V(CA, D)) < 0) { NA = M(NA); }
  
  pt  CB = circumcenter2D(A, B, D);
  vec NB = cross(V(A, B), V(A, D)).normalize().mul(300);
  
  // rotation orthogonal axis
  vec NN = U(cross(NA, NB));
  vec I = U(NA);
  vec J = U(cross(I, NN));
  
  vec V1 = V(CB, CA);
  vec V2 = V(CB, P(CB, NB));
  vec V3 = R(NA, PI/2, I, J);
  
  float t = dot(V1, V3)/dot(V2, V3);
  return P(CB, t, P(CB, NB));
}