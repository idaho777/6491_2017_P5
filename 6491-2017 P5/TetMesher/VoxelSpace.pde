
public class VoxelSpace {
  List<pt>[][][] voxels;
  
  pt lowCorner, highCorner;
  float voxelSize;
  
  
  public VoxelSpace(pt lc, pt hc, float boxSize) {
    lowCorner = P(lc);
    highCorner = P(hc);
    voxelSize = boxSize;
    
    int xW = getX(highCorner) - getX(lowCorner) + 1;
    int yW = getY(highCorner) - getY(lowCorner) + 1;
    int zW = getZ(highCorner) - getZ(lowCorner) + 1;
    voxels = new ArrayList[xW][yW][zW];
    
    for (int x = 0; x < voxels.length; ++x) {
      for (int y = 0; y < voxels.length; ++y) {
        for (int z = 0; z < voxels.length; ++z) {
          voxels[x][y][z] = new ArrayList<pt>();
        }
      }
    }
  }
  
  void addPt(pt p) {
    if (p == null) return;
    getVoxel(p).add(p);
  }
  
  void removePt(pt p) {
    if (p == null) return;
    getVoxel(p).remove(p);
  }
  
  List<pt> getVoxel(pt p) {
    return voxels[getX(p)][getY(p)][getZ(p)];
  }
  
  int getX(pt p) {
    if (p.x <= lowCorner.x) return 0;
    if (p.x >= highCorner.x) return int(highCorner.x);
    return int(p.x/voxelSize);
  }
  
  int getY(pt p) {
    if (p.y <= lowCorner.y) return 0;
    if (p.y >= highCorner.y) return int(highCorner.y);
    return int(p.y/voxelSize);
  }
  
  int getZ(pt p) {
    if (p.z <= lowCorner.z) return 0;
    if (p.z >= highCorner.z) return int(highCorner.z);
    return int(p.z/voxelSize);
  }
}  



VoxelSpace setupVoxelSpace(List<List<pt>> pointCloud) {
  VoxelSpace vs = new VoxelSpace(P(0, 0, 0), P(2000, 2000, h_ceiling), 200);
  for (List<pt> l : pointCloud) {
    for (pt p : l) {
      vs.addPt(p);
    }
  }
  return vs;
}