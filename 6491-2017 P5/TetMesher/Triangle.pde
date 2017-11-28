Triangle TRI(pt A, pt B, pt C) { return new Triangle(A, B, C); }

class Triangle {
  pt A, B, C;
  
  Triangle() {}
  Triangle(pt A, pt B, pt C) {
    this.A = A;
    this.B = B;
    this.C = C;
  }
  
  boolean pointInsideCircum(pt Z) {
    pt ccenter = circumcenter2D(A, B, C);
    float r = d(ccenter, A);
    return d(ccenter, Z) <= r;
  }
  
  boolean pointInsideTriangle(pt Z) {
    return false; 
  }
  
  boolean contains(Edge e) {
    return contains(e.A) && contains(e.B); 
  }
  
  boolean contains(pt Z) {
    return A.equals(Z) || B.equals(Z) || C.equals(Z); 
  }
  
  void drawShape() {
    beginShape();
      v(A);
      v(B);
      v(C);
    endShape(CLOSE);
  }
  
  void draw() {
    beam(A, B, rt);
    beam(B, C, rt);
    beam(C, A, rt);
    //fill(red); show(circumcenter2D2(A, B, C), rt);
    //fill(cyan); show(circumcenter2D(A, B, C), rt*2);
  }
  
  vec N() {
    return M(U(cross(V(A, B), V(B, C))));
  }
  
  @Override
  public boolean equals(Object o) {
    if (!(o instanceof Triangle)) return false;
    
    Triangle ot = (Triangle) o;
    return contains(ot.A) && contains(ot.B) && contains(ot.C);
  }
  
  @Override
  public String toString() {
    return "Triangle: " + A + ", " + B + ", " + C;
  }
}

//pt circumcenter2D (pt A, pt B, pt C) { // circumcenter to triangle (A,B,C)
//  vec AB = V(A, B);
//  float ab2 = dot(AB,AB);
//  vec AC = V(A, C); AC = R(AC);
//  float ac2 = dot(AC,AC);
//  float d = 2*dot(AB,AC);
//  AB = R(AB);
//  AB.rev(); AB.mul(ac2);
//  AC.mul(ab2);
//  AB.add(AC);
//  AB.div(d);
//  pt X = P(A);
//  X.add(AB);
//  return(X);
//};

pt circumcenter2D (pt A, pt B, pt C) {
  pt O = P(A, 0.5, B);
  
  vec planeNormal = U(cross(V(A, B), V(A, C)));
  vec I = U(V(A, B));
  vec J = U(cross(planeNormal, I));
  vec N = U(R(V(A, B), -PI/2, I, J));
  //vec N = U(R(V(A, B)));
  if (dot(N, V(A, C)) < 0) { N = M(N); }
  
  pt Q = P(A, 0.5, C);
  vec M = U(R(V(A, C), PI/2, I, J)).mul(300);
  //vec M = R(V(A, C)).normalize().mul(300);
  
  vec V1 = V(Q, O);
  vec V2 = V(Q, P(Q, M));
  vec V3 = R(N, PI/2, I, J);
  //vec V3 = R(N);
  
  float t = dot(V1, V3)/dot(V2, V3);
  return P(Q, t, P(Q, M));
}