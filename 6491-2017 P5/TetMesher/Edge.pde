Edge EDGE(pt A, pt B) { return new Edge(A, B); }


class Edge {
  pt A, B;
  
  Edge() {}
  Edge(pt A, pt B) {
    this.A = A;
    this.B = B;
  }
  
  void draw() {
    beam(A, B, rt);
  }
  
  public boolean contains(pt p) {
    return A.equals(p) || B.equals(p); 
  }
  
  
  List<pt> samplePoints() {
    float edgeStepLen = 2*PI*rt/circumPts;
    List<pt> list = new ArrayList<pt>();
    vec skewer = V(A, B);
    vec I = U(cross(skewer, R(skewer)));
    vec J = U(cross(skewer, I));
    
    pt curr = P(A);
    vec rad = V(curr, P(curr, rt, I));
    
    boolean on = false;
    while (d(curr, A) < d(A, B)) {
      for (int i = 0; i < circumPts; ++i) {
        list.add(P(curr, R(rad, i*TAU/circumPts + ((on) ? PI/circumPts : 0), I, J)));
      }
      on = !on;
      curr = P(curr, edgeStepLen, U(skewer));
    }
    
    if (d(curr, A) != d(A, B)) {
      for (int i = 0; i < circumPts; ++i) {
        list.add(P(B, R(rad, i*TAU/circumPts + ((on) ? PI/circumPts : 0), I, J)));
      }
    }
    
    return list;
  }
  
  
  @Override
  public boolean equals(Object o) {
    if (!(o instanceof Edge)) { return false; }
    Edge e = (Edge) o;    
    return contains(e.A) && contains(e.B);
  }
  
  @Override
  public String toString() {
    return "[" + A + ", " + B + "]";
  }
}