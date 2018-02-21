#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;

#define PI 3.14159265



//classes

class Edge{
  vector<vector<float>> endPoints;
  float xMax;
  float yMax;
  float xMin;
  float yMin;
  vector<float> intersections;

  public:
    vector<float> calculateIntersection(vector<float>);
    void setEndpoints (vector<float>);
};

vector<float> Edge::calculateIntersection(vector<float> polyEdge){
    float dx, dy, c1, c2, m1, m2;
    vector<float> intersection;
    dx = polyEdge[2] - polyEdge[0];
    dy = polyEdge[3] - polyEdge[1];
 
    m1 = dy / dx;

    c1 = polyEdge[1] - m1 * polyEdge[0];

    dx = polyEdge[6] - polyEdge[4];
    dy = polyEdge[7] - polyEdge[5];
 
    m2 = dy / dx;
    c2 = polyEdge[7] - m2 * polyEdge[6]; 

    if( (m1 - m2) == 0)
        cout << "No Intersection between the lines\n";
    else
    {
        intersection.push_back((c2 - c1) / (m1 - m2));
        intersection.push_back(m1 * intersection[0] + c1);
    }

    return intersection;

}

void Edge::setEndpoints(vector<float> endpoint){
  endPoints.push_back(endpoint);
}


class Face {
    //index to use to find faces in 
    vector<float> faces;
    vector<float> coords;
    vector<float> viewport;
    vector<float> normals;
    vector<Edge> edges;

  public:
    void setCoords (vector<float>);
    void setFaces (vector<float>);
    void setViewport (vector<float>);
    void setEdges (Edge);
    int faceSize ();
    int coordsSize ();
    int viewSize ();
    float getFaceIndex(int);
    float getCoordIndex(int);
    float getViewIndex(int);
    Edge getEdges(int);
    void viewPushBack(float);
    void clearEdges ();
};

void Face::setCoords (vector<float> allCoords) {
  coords = allCoords;
}
void Face::setFaces (vector<float> allFaces) {
  faces = allFaces;
}
void Face::setViewport (vector<float> allView) {
  viewport = allView;
}
void Face::setEdges (Edge anEdge) {
  edges.push_back(anEdge);
}
int Face::faceSize(){
  int size;
  size = faces.size();

  return size;
}
int Face::coordsSize(){
  int size;
  size = coords.size();

  return size;
}
int Face::viewSize(){
  int size;
  size = viewport.size();

  return size;
}
float Face::getFaceIndex(int i){
  float value;
  value = faces[i];
  return value;
}
float Face::getCoordIndex(int i){
  float value;
  value = coords[i];
  return value;
}
float Face::getViewIndex(int i){
  float value;
  value = viewport[i];
  return value;
}
Edge Face::getEdges(int i){
  Edge value;
  value = edges[i];
  return value;
}
void Face::viewPushBack(float value){
  viewport.push_back(value);
}
void Face::clearEdges(){
  edges.clear();
}

//prototype
Face parsePS(string filename, int framew, int frameh, int loop);
void writeToBuffer(vector<int> pixelsX, vector<int> pixelsY, string filename, vector<vector<char>> xpm, int framew, int frameh, int loop);
vector<int> clipCords(int x1, int y1, int x2, int y2, int code1, int code2, int xl, int yl, int xh, int yh);
int getcode(int x,int y, int z);
vector<vector<char>> fillPoly(vector<vector<char>> xpm, vector<vector<float>> zBuffer, Face polygon);
float degToRad(int degrees);
vector<float> getSurfaceNormal(vector<vector<float>> p1, vector<vector<float>> p2, vector<vector<float>> p3);

//GLOBAL
static int LEFT=1, RIGHT=2, BOTTOM=4, TOP=8, NEAR=16, FAR=32;

uint64_t constexpr mix(char m, uint64_t s)
 {
  return ((s<<7) + ~(s>>3)) + ~m;
  }
constexpr unsigned str2int(const char* m)
{
    return (*m) ? mix(*m,str2int(m+1)) : 0;
}
const int MAX = 100;
 
 
// Function to multiply two matrices A[][] and B[][]
vector<vector<float>> multiplyMatrix(int row1, int col1, int row2, int col2, vector<vector<float>> A, vector<vector<float>> B)
{
    if (row2 != col1)
    {
      cout << "not possible" << endl;
    }
    vector<vector<float>> C;
    C.resize(row1, vector<float>(col2));

    for (int row = 0; row != row1; ++row) 
    {
      for (int col = 0; col != col2; ++col)
      {
        
        for (int inner = 0; inner != col1; ++inner)
          {
            C[row][col] += A[row][inner] * B[inner][col];
          }
         
      }
    }
  return C;
}

int main(int argc, char* argv[]){
  vector<int> drawX, drawY;
  int px, py, a, b, c, rotation, transX, transY, transZ, xl = 0, yl = 0, 
  xh = 500, yh = 500, lowXView, lowYView, highXView, highYView;
  float x1, y1, z1, x2, y2, x3, y3, preX1, preY1, preZ1, preX2, preY2, preZ2,dy, dx, m, rx, ry, dex, dey, yInter, scale, prpX = 0, prpY = 0, prpZ = 1, vrpX = 0, vrpY = 0, vrpZ = 0, vpnX = 0, 
  vpnY = 0, vpnZ = -1, vuvX = 0, vuvY = 1, vuvZ = 0, vrcUMin = -0.7, vrcVMin = -0.7, vrcUMax = .7, vrcVMax = .7;
  string filename1, filename2, filename3;
  vector<string> allFiles;
  Face coordinates;
  bool projectionType = 0;
  float front, back;
  front = .6;
  back = -.6;

  transX = 0;
  transY = 0;
  transZ = 0;
  scale = 1;
  rotation = 0;
  int loop = 1;
  int loopTotal = 1;

  argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present

  //switch statements for args
  for (int i = 0; i < argc; i+=2) {
    switch (str2int(argv[i])) {
      case str2int("-f"):
        //filename
        cout << argv[i+1] << endl;
        filename1 = argv[i+1];
        allFiles.push_back(filename1);
        break;
      case str2int("-j"):
        //The next argument is an integer lower bound in the x dimension of the viewport (0)
        cout << argv[i+1] << endl;
        lowXView = atof(argv[i+1]);
        break;
      case str2int("-k"):
        //  The next argument is an integer lower bound in the y dimension of the viewport (0)
        cout << argv[i+1] << endl;
        lowYView = atoi(argv[i+1]);
        break;
      case str2int("-o"):
        // The next argument is an integer upper bound in the x dimension of the viewport window (499)
        cout << argv[i+1] << endl;
        highXView = atoi(argv[i+1]);
        break;
      case str2int("-p"):
        //  The next argument is an integer upper bound in the y dimension of the viewport window (499)
        cout << argv[i+1] << endl;
        highYView = atoi(argv[i+1]);
        break;
      case str2int("-x"):
        // The next argument is a floating point value for the x of the Projection Reference Point (PRP) in VRC coordinates. (0.0)
        cout << argv[i+1] << endl;
        prpX = atoi(argv[i+1]);
        break;
      case str2int("-y"):
        // The next argument is a floating point value for the y of the Projection Reference Point (PRP) in VRC coordinates. (0.0)
        cout << argv[i+1] << endl;
        prpY = atoi(argv[i+1]);
        break;
      case str2int("-z"):
        //The next argument is a floating point value for the z of the Projection Reference Point (PRP) in VRC coordinates. (1.0)
        cout << argv[i+1] << endl;
        prpZ = atoi(argv[i+1]);
        break;
      case str2int("-X"):
        //The next argument is a floating point value for the x of the View Reference Point (VRP) in VRC coordinates. (0.0)
        cout << argv[i+1] << endl;
        vrpX = atoi(argv[i+1]);
      case str2int("-Y"):
        // The next argument is a floating point value for the y of the View Reference Point (VRP) in VRC coordinates. (0.0)
        cout << argv[i+1] << endl;
        vrpY = atoi(argv[i+1]);
        break;
      case str2int("-Z"):
        // The next argument is a floating point value for the z of the View Reference Point (VRP) in VRC coordinates. (0.0)
        cout << argv[i+1] << endl;
        vrpZ = atoi(argv[i+1]);
        break;
      case str2int("-q"):
        // The next argument is a floating point value for the x of the View Plane Normal (VPN) in VRC coordinates. (0.0)
        cout << argv[i+1] << endl;
        vpnX = atoi(argv[i+1]);
        break;
      case str2int("-r"):
        // The next argument is a floating point value for the y of the View Plane Normal (VPN) in VRC coordinates. (0.0)
        cout << argv[i+1] << endl;
        vpnY = atoi(argv[i+1]);
        break;
      case str2int("-w"):
        //The next argument is a floating point value for the z of the View Plane Normal (VPN) in VRC coordinates. (-1.0)
        cout << argv[i+1] << endl;
        vpnZ = atoi(argv[i+1]);
        break;
      case str2int("-Q"):
        //The next argument is a floating point value for the x of the View Up Vector (VUP) in VRC coordinates. (0.0)
        cout << argv[i+1] << endl;
        vuvX = atoi(argv[i+1]);
        break;
      case str2int("-R"):
        // The next argument is a floating point value for the y of the View Up Vector (VUP) in VRC coordinates. (1.0)
        cout << argv[i+1] << endl;
        vuvY = atoi(argv[i+1]);
        break;
      case str2int("-W"):
        //The next argument is a floating point value for the z of the View Up Vector (VUP) in VRC coordinates. (0.0)
        cout << argv[i+1] << endl;
        vuvZ = atoi(argv[i+1]);
        break;
      case str2int("-u"):
        // The next argument is a floating point value for the u min of the VRC window in VRC coordinates. (-0.7)
        cout << argv[i+1] << endl;
        vrcUMin = atoi(argv[i+1]);
        break;
      case str2int("-v"):
        // The next argument is a floating point value for the v min of the VRC window in VRC coordinates. (-0.7)
        cout << argv[i+1] << endl;
        vrcVMin = atoi(argv[i+1]);
        break;
      case str2int("-U"):
        // The next argument is a floating point value for the u max of the VRC window in VRC coordinates. (0.7)
        cout << argv[i+1] << endl;
        vrcUMax = atoi(argv[i+1]);
        break;
      case str2int("-V"):
        //The next argument is a floating point value for the v max of the VRC window in VRC coordinates. (0.7)
        cout << argv[i+1] << endl;
        vrcVMax = atoi(argv[i+1]);
        break;
      case str2int("-P"):
        // Is a flag to indicate to use parallel projection. If not present, use perspective projection
        projectionType = 1;
        break;
      case str2int("-F"):
        // Near clipping plane
        front = atoi(argv[i+1]);
        break;
      case str2int("-B"):
        // far clipping plane
        back = atoi(argv[i+1]);
        break;
      case str2int("-g"):
        // another model input that is green
        cout << argv[i+1] << endl;
        loopTotal++;
        filename2 = argv[i+1];
        allFiles.push_back(filename2);
        break;
      case str2int("-i"):
        // another model input that is blue
        cout << argv[i+1] << endl;
        loopTotal++;
        filename3 = argv[i+1];
        allFiles.push_back(filename3);
        break;
  }
}

  //set matricies
  vector<vector<float>> prp, vrp, vpn, negVRP, vuv, unitV, viewCoordU, viewCoordV;
  prp.resize(4,vector<float>(1));
  vrp.resize(4,vector<float>(1));
  vpn.resize(4,vector<float>(1));
  negVRP.resize(4,vector<float>(1));
  vuv.resize(4,vector<float>(1));
  unitV.resize(4,vector<float>(1));
  viewCoordU.resize(4,vector<float>(1));
  viewCoordV.resize(4,vector<float>(1));
  for(int row = 0; row < 4; row++){
    for(int col = 0; col < 1; col++){
      prp[row][col] = 0;
      vrp[row][col] = 0;
      vpn[row][col] = 0;
      negVRP[row][col] = 0;
      vuv[row][col] = 0;
      unitV[row][col] = 0;
      viewCoordU[row][col] = 0;
      viewCoordV[row][col] = 0;
    }
  }
  prp[0][0] = prpX;
  prp[1][0] = prpY;
  prp[2][0] = prpZ;
  prp[3][0] = 1;

  vrp[0][0] = vrpX;
  vrp[1][0] = vrpY;
  vrp[2][0] = vrpZ;
  vrp[3][0] = 1;

  negVRP[0][0] = -1*vrpX;
  negVRP[1][0] = -1*vrpY;
  negVRP[2][0] = -1*vrpZ;
  negVRP[3][0] = 1;


  vpn[0][0] = vpnX;
  vpn[1][0] = vpnY;
  vpn[2][0] = vpnZ;
  vpn[3][0] = 1;

  vuv[0][0] = vuvX;
  vuv[1][0] = vuvY;
  vuv[2][0] = vuvZ;
  vuv[3][0] = 1;

  unitV[0][0] = 1;
  unitV[1][0] = 1;
  unitV[2][0] = 1;
  unitV[3][0] = 1;

  //get u,v,n for camera coordinate system
  //n = VPN
  //u = VUP x VPN
  //v = n x u

  float uA = vuv[1][0] * vpn[2][0] - vuv[2][0] * vpn[1][0];
  float uB = vuv[2][0] * vpn[0][0] - vuv[0][0] * vpn[2][0];
  float uC = vuv[0][0] * vpn[1][0] - vuv[1][0] * vpn[0][0];
  viewCoordU[0][0] = uA;
  viewCoordU[1][0] = uB;
  viewCoordU[2][0] = uC;
  viewCoordU[3][0] = 1;

  float vA = vpn[1][0] * viewCoordU[2][0] - vpn[2][0] * viewCoordU[1][0] ;
  float vB = vpn[2][0] * viewCoordU[0][0] - vpn[0][0] * viewCoordU[2][0] ;
  float vC = vpn[0][0] * viewCoordU[1][0] - vpn[1][0] * viewCoordU[0][0] ;
  viewCoordV[0][0] = vA;
  viewCoordV[1][0] = vB;
  viewCoordV[2][0] = vC;
  viewCoordV[3][0] = 1;

  //get window center
  vector<float> windowCenter;
  windowCenter.push_back((vrcUMax+vrcUMin)/2);
  windowCenter.push_back((vrcVMax+vrcVMin)/2);

  //initialize scale matrix
  vector<vector<float> > scaleM3;
  scaleM3.resize(4, vector<float>(4));
  scaleM3.at(0).at(0) = scale;
  scaleM3.at(0).at(1) = 0;
  scaleM3.at(0).at(2) = 0;
  scaleM3.at(0).at(3) = 0;

  scaleM3.at(1).at(0) = 0; 
  scaleM3.at(1).at(1) = scale;
  scaleM3.at(1).at(2) = 0;
  scaleM3.at(1).at(3) = 0;

  scaleM3.at(2).at(0) = 0;
  scaleM3.at(2).at(1) = 0;
  scaleM3.at(2).at(2) = scale;
  scaleM3.at(2).at(3) = 0;

  scaleM3.at(3).at(0) = 0;
  scaleM3.at(3).at(1) = 0;
  scaleM3.at(3).at(2) = 0;
  scaleM3.at(3).at(3) = 1;

  //___________________________________________________

  //3D TRANSFORM MATRIX

  //___________________________________________________

  vector<vector<float> > translationM3;
  translationM3.resize(4, vector<float>(4));
  translationM3.at(0).at(0) = 1;
  translationM3.at(0).at(1) = 0;
  translationM3.at(0).at(2) = 0;
  translationM3.at(0).at(3) = negVRP[0][0];

  translationM3.at(1).at(0) = 0; 
  translationM3.at(1).at(1) = 1;
  translationM3.at(1).at(2) = 0;
  translationM3.at(1).at(3) = negVRP[1][0];

  translationM3.at(2).at(0) = 0;
  translationM3.at(2).at(1) = 0;
  translationM3.at(2).at(2) = 1;
  translationM3.at(2).at(3) = negVRP[2][0];

  translationM3.at(3).at(0) = 0;
  translationM3.at(3).at(1) = 0;
  translationM3.at(3).at(2) = 0;
  translationM3.at(3).at(3) = 1;


  //initialize rotation matrix for each axis
  vector<vector<float> > rotationXM3, rotationYM3, rotationZM3;
  rotationXM3.resize(4, vector<float>(4));
  rotationYM3.resize(4, vector<float>(4));
  rotationZM3.resize(4, vector<float>(4));

  rotationXM3.at(0).at(0) = 1;
  rotationXM3.at(0).at(1) = 0;
  rotationXM3.at(0).at(2) = 0;
  rotationXM3.at(0).at(3) = 0;

  rotationXM3.at(1).at(0) = 0; 
  rotationXM3.at(1).at(1) = cos(degToRad(rotation));
  rotationXM3.at(1).at(2) = -1*sin(degToRad(rotation));
  rotationXM3.at(1).at(3) = 0;

  rotationXM3.at(2).at(0) = 0;
  rotationXM3.at(2).at(1) = sin(degToRad(rotation));
  rotationXM3.at(2).at(2) = cos(degToRad(rotation));
  rotationXM3.at(2).at(3) = 0;

  rotationXM3.at(3).at(0) = 0;
  rotationXM3.at(3).at(1) = 0;
  rotationXM3.at(3).at(2) = 0;
  rotationXM3.at(3).at(3) = 1;

  rotationYM3.at(0).at(0) = cos(degToRad(rotation));
  rotationYM3.at(0).at(1) = 0;
  rotationYM3.at(0).at(2) = sin(degToRad(rotation));
  rotationYM3.at(0).at(3) = 0;

  rotationYM3.at(1).at(0) = 0; 
  rotationYM3.at(1).at(1) = 1;
  rotationYM3.at(1).at(2) = 0;
  rotationYM3.at(1).at(3) = 0;

  rotationYM3.at(2).at(0) = -1*sin(degToRad(rotation));
  rotationYM3.at(2).at(1) = 0;
  rotationYM3.at(2).at(2) = cos(degToRad(rotation));
  rotationYM3.at(2).at(3) = 0;

  rotationYM3.at(3).at(0) = 0;
  rotationYM3.at(3).at(1) = 0;
  rotationYM3.at(3).at(2) = 0;
  rotationYM3.at(3).at(3) = 1;

  rotationZM3.at(0).at(0) = cos(degToRad(rotation));
  rotationZM3.at(0).at(1) = -1*sin(degToRad(rotation));
  rotationZM3.at(0).at(2) = 0;
  rotationZM3.at(0).at(3) = 0;

  rotationZM3.at(1).at(0) = sin(degToRad(rotation)); 
  rotationZM3.at(1).at(1) = cos(degToRad(rotation));
  rotationZM3.at(1).at(2) = 0;
  rotationZM3.at(1).at(3) = 0;

  rotationZM3.at(2).at(0) = 0;
  rotationZM3.at(2).at(1) = 0;
  rotationZM3.at(2).at(2) = 1;
  rotationZM3.at(2).at(3) = 0;

  rotationZM3.at(3).at(0) = 0;
  rotationZM3.at(3).at(1) = 0;
  rotationZM3.at(3).at(2) = 0;
  rotationZM3.at(3).at(3) = 1;
  
  cout << "high x world " << xh << endl;
  cout << "low x world: " << xl << endl;
  cout << "high y world: " << yh << endl;
  cout << "low y world: " << yl << endl;

  //create composite matrix for world view to viewport
  //transform to origin * scale
  //then that * transform back
  xh = 500;
  yh = 500;
  xl = 0;
  yl = 0;


  int frameh = (yh - yl);
  int framew = (xh - xl);
  vector<vector<char>> xpm;
  vector<vector<float>> zBuffer;
  xpm.resize(framew,vector<char>(frameh));
  zBuffer.resize(framew,vector<float>(frameh));

  //initialize xpm with all X's
  // cout << "the size of the vector of vectors " << xpm.size() << " " << xpm[0].size() << endl;
  for(int row = 0; row < framew; row++){
    for(int col = 0; col < frameh; col++){
      xpm[row][col] = 'X';
      zBuffer[row][col] = back;
    }
  }

  //get coordinates from ps file
  vector<vector<float>> orthoProj, perspProj;
  //orthographic projection matrix
  orthoProj.resize(4,vector<float>(4));
  for(int row = 0; row < 4; row++){
    for(int col = 0; col < 4; col++){
      orthoProj[row][col] = 0;
    }
  }
  //first row
  orthoProj[0][0] = 2/(vrcUMax - vrcUMin);
  orthoProj[0][2] = ((vrcUMax + vrcUMin ) - 2* prp[0][0]) / ((vrcUMax - vrcUMin )* prp[2][0]);
  orthoProj[0][3] = -1*(vrcUMax + vrcUMin ) / (vrcUMax - vrcUMin );
  //second row
  orthoProj[1][1] = 2/(vrcVMax - vrcVMin);
  orthoProj[1][2] = ((vrcVMax + vrcVMin ) - 2* prp[1][0]) / ((vrcVMax - vrcVMin )* prp[2][0]);
  orthoProj[1][3] = -1*(vrcVMax + vrcVMin ) / (vrcVMax - vrcVMin );;

  //third row
  orthoProj[2][2] = 1/(front - back);
  orthoProj[2][3] = -1* (front/(front - back));

  //last row
  orthoProj[3][3] = 1;
  cout << "orthographic projection matrix: " << endl;
  for(int i = 0; i <   4; ++i)
  {
    for(int j = 0; j < 4; ++j)
    {
      cout << orthoProj[i][j] << " ";
      if(j == 4 - 1)
        cout << endl;
    }
  }

  //perspective projection matrix
  perspProj.resize(4,vector<float>(4));
  for(int row = 0; row < 4; row++){
    for(int col = 0; col < 4; col++){
      perspProj[row][col] = 0;
    }
  }

  perspProj[0][0] = (2 * prp[2][0])/((vrcUMax - vrcUMin)* (prp[2][0] - back)) ;
  perspProj[0][2] = ((vrcUMax + vrcUMin ) - 2* prp[0][0]) / ((vrcUMax - vrcUMin )* (prp[2][2] - back));
  perspProj[0][3] = -1*((vrcUMax + vrcUMin ) - 2* prp[0][0]) / ((vrcUMax - vrcUMin )* (prp[2][0] - back));
  //second row
  perspProj[1][1] = (2 * prp[2][0])/((vrcVMax - vrcVMin)* (prp[2][0] - back));
  perspProj[1][2] = ((vrcVMax + vrcVMin ) - 2* prp[0][0]) / ((vrcVMax - vrcVMin )* (prp[2][0] - back));
  perspProj[1][3] = -1*((vrcVMax + vrcVMin ) - 2* prp[0][0]) / ((vrcVMax - vrcVMin )* (prp[2][0] - back));

  //third row
  perspProj[2][2] = 1/(prp[2][0] - back);
  perspProj[2][3] = -1* (prp[2][0]/(prp[2][0] - back));

  //last row
  perspProj[3][3] = 1;
  cout << "perspective projection matrix: " << endl;
  for(int i = 0; i <   4; ++i)
  {
    for(int j = 0; j < 4; ++j)
    {
      cout << perspProj[i][j] << " ";
      if(j == 4 - 1)
        cout << endl;
    }
  }


  //lastly create view matrix
  vector<vector<float>> finalRot;
  vector<vector<float>> viewM, finalTrans;
  viewM.resize(4, vector<float>(4));
  finalRot.resize(4, vector<float>(4));

  //first translate
  cout << "translationM3: " << endl;
  for(int i = 0; i <   4; ++i)
  {
    for(int j = 0; j < 4; ++j)
    {
      cout << translationM3[i][j] << " ";
      if(j == 4 - 1)
        cout << endl;
    }
  }
  //then create rotation matrix for final view
  finalRot.at(0).at(0) = viewCoordU[0][0];
  finalRot.at(0).at(1) = viewCoordU[1][0];
  finalRot.at(0).at(2) = viewCoordU[2][0];
  finalRot.at(0).at(3) = 0;

  finalRot.at(1).at(0) = viewCoordV[0][0]; 
  finalRot.at(1).at(1) = viewCoordV[1][0];
  finalRot.at(1).at(2) = viewCoordV[2][0];
  finalRot.at(1).at(3) = 0;

  finalRot.at(2).at(0) = vpnX;
  finalRot.at(2).at(1) = vpnY;
  finalRot.at(2).at(2) = vpnZ;
  finalRot.at(2).at(3) = 0;

  finalRot.at(3).at(0) = 0;
  finalRot.at(3).at(1) = 0;
  finalRot.at(3).at(2) = 0;
  finalRot.at(3).at(3) = 1;
  cout << "rotate camera to view reference coordinates: " << endl;
  for(int i = 0; i <   4; ++i)
  {
    for(int j = 0; j < 4; ++j)
    {
      cout << finalRot[i][j] << " ";
      if(j == 4 - 1)
        cout << endl;
    }
  }

  //Create projection times view matrix this will be used to multiply each set of world coordinates to get connonical coords
  viewM = multiplyMatrix(4, 4, 4, 4, finalRot, translationM3);
  cout << "View matrix: " << endl;
  for(int i = 0; i <   4; ++i)
  {
    for(int j = 0; j < 4; ++j)
    {
      cout << viewM[i][j] << " ";
      if(j == 4 - 1)
        cout << endl;
    }
  }

  vector<vector<float>> viewTimesProj;
  viewTimesProj.resize(4, vector<float>(4));

  //multiply view matrix and projection matrix. if -P flag use parallel projection else use perspective
  if(projectionType){
    viewTimesProj = multiplyMatrix(4, 4, 4, 4, orthoProj, viewM);
    cout << "orthographic projection" << endl;
  }
  
  else{
    viewTimesProj = multiplyMatrix(4, 4, 4, 4, perspProj, viewM);
    cout << "perspective projection" << endl;
  }

  cout << "viewTimesProj: " << endl;
  for(int i = 0; i <   4; ++i)
  {
    for(int j = 0; j < 4; ++j)
    {
      cout << viewTimesProj[i][j] << " ";
      if(j == 4 - 1)
        cout << endl;
  
    }
  }
  while(loop <= loopTotal){
    coordinates = parsePS(allFiles.at(loop-1), framew, frameh, loop);
  // cout << "did the coordinates actually get set to Face? " << coordinates.getCoordIndex(1) << endl;


  // cout << "we've reached past parsing" << endl;

  //start calculating matricies
  

  //for each line of new coordinates...convert 3d coordinates to 2d viewport coords
  for(int i = 0; i != coordinates.coordsSize(); i += 3){

    
    vector<vector<float>> canonical1(4, vector<float>(1)), canonical2(4, vector<float>(1)), pre1(4, vector<float>(1)), pre2(4, vector<float>(1));

    //calculate coordinates from world coordinates and view matrix x projection matrix
    preX1 = coordinates.getCoordIndex(i);
    preY1 = coordinates.getCoordIndex(i+1);
    preZ1 = coordinates.getCoordIndex(i+2);
    pre1[0][0] = preX1;
    pre1[1][0] = preY1;
    pre1[2][0] = preZ1;
    pre1[3][0] = 1;

    // preX2 = coordinates.getCoordIndex(i+3);
    // preY2 = coordinates.getCoordIndex(i+4);
    // preZ2 = coordinates.getCoordIndex(i+5);
    // pre2[0][0] = preX2;
    // pre2[1][0] = preY2;
    // pre2[2][0] = preZ2;
    // pre2[3][0] = 1;

    // cout << "preX1 coords: " << preX1 << endl;
    // cout << "preY1 coords: " << preY1 << endl;
    // cout << "preX2 coords: " << preX2 << endl;
    // cout << "preY2 coords: " << preY2 << endl;

    canonical1 = multiplyMatrix(4, 4, 4, 1, viewTimesProj, pre1); 
    // cout << "canonical1: " << endl;
    // for(int i = 0; i <   4; ++i)
    // {
    //   for(int j = 0; j < 1; ++j)
    //   {
    //     cout << canonical1[i][j] << " ";
    //     if(j == 1 - 1)
    //       cout << endl;
    //   }
    // }

    // canonical2 = multiplyMatrix(4, 4, 4, 1, viewTimesProj, pre2);
    // cout << "canonical2: " << endl;
    // for(int i = 0; i <   4; ++i)
    // {
    // for(int j = 0; j < 1; ++j)
    // {
    //   cout << canonical2[i][j] << " ";
    //   if(j == 1 - 1)
    //     cout << endl;
    // }
    // }

    if(!projectionType){
      x1 = canonical1[0][0]/(-1*canonical1[2][0]);
      y1 = canonical1[1][0]/(-1*canonical1[2][0]);
      z1 = canonical1[2][0];
      // x2 = canonical2[0][0]/(-1*canonical2[2][0]);
      // y2 = canonical2[1][0]/(-1*canonical2[2][0]);  
      
    }
    else{
      x1 = canonical1[0][0];
      y1 = canonical1[1][0];
      z1 = canonical1[2][0];
      // x2 = canonical2[0][0];
      // y2 = canonical2[1][0];    
    }

    // cout << "x1 canonical coords: " << x1 << endl;
    // cout << "y1 canonical coords: " << y1 << endl;
    // cout << "x2 canonical coords: " << x2 << endl;
    // cout << "y2 canonical coords: " << y2 << endl;

    // vector<vector<float>>canonical2D1, canonical2D2;
    // canonical2D1.resize(3, vector<float>(1));
    // //canonical2D2.resize(3, vector<float>(1));

    // canonical2D1[0][0] = x1;
    // canonical2D1[1][0] = y1;
    // canonical2D1[2][0] = 1;

    // canonical2D2[0][0] = x2;
    // canonical2D2[1][0] = y2;
    // canonical2D2[2][0] = 1;




    //change to viewport coordinate system
    //translate T(-(-1),-(-1))
    //then scale S((umax - umin)/(1-(-1)) . (vmax - vmin)/(1-(-1)))
    //then translate again T(umin, vmin)
    vector<vector<float>> step1, worldView, viewport1, viewport2;
    step1.resize(4, vector<float>(4));
    worldView.resize(4, vector<float>(4));
    viewport1.resize(4, vector<float>(1));
    //viewport2.resize(3, vector<float>(1));



    vector<vector<float> > translationViewport1;
    translationViewport1.resize(4, vector<float>(4));
    translationViewport1.at(0).at(0) = 1;
    translationViewport1.at(0).at(1) = 0;
    translationViewport1.at(0).at(2) = 0;
    translationViewport1.at(0).at(3) = 1;

    translationViewport1.at(1).at(0) = 0; 
    translationViewport1.at(1).at(1) = 1;
    translationViewport1.at(1).at(2) = 0;
    translationViewport1.at(1).at(3) = 1;

    translationViewport1.at(2).at(0) = 0;
    translationViewport1.at(2).at(1) = 0;
    translationViewport1.at(2).at(2) = 1;
    translationViewport1.at(2).at(3) = 1;

    translationViewport1.at(3).at(0) = 0;
    translationViewport1.at(3).at(1) = 0;
    translationViewport1.at(3).at(2) = 0;
    translationViewport1.at(3).at(3) = 1;

    vector<vector<float> > translationViewport2;
    translationViewport2.resize(4, vector<float>(4));
    translationViewport2.at(0).at(0) = 1;
    translationViewport2.at(0).at(1) = 0;
    translationViewport2.at(0).at(2) = 0;
    translationViewport2.at(0).at(3) = 0;


    translationViewport2.at(1).at(0) = 0; 
    translationViewport2.at(1).at(1) = 1;
    translationViewport2.at(1).at(2) = 0;
    translationViewport2.at(1).at(3) = 0;


    translationViewport2.at(2).at(0) = 0;
    translationViewport2.at(2).at(1) = 0;
    translationViewport2.at(2).at(2) = 1;
    translationViewport2.at(2).at(3) = 0;


    translationViewport2.at(3).at(0) = 0;
    translationViewport2.at(3).at(1) = 0;
    translationViewport2.at(3).at(2) = 0;
    translationViewport2.at(3).at(3) = 1;




    vector<vector<float> > scaleViewport;
    scaleViewport.resize(4, vector<float>(4));
    scaleViewport.at(0).at(0) = (xh - xl)/2;
    scaleViewport.at(0).at(1) = 0;
    scaleViewport.at(0).at(2) = 0;
    scaleViewport.at(0).at(3) = 0;

    scaleViewport.at(1).at(0) = 0; 
    scaleViewport.at(1).at(1) = (yh - yl)/2;
    scaleViewport.at(1).at(2) = 0;
    scaleViewport.at(1).at(3) = 0;

    scaleViewport.at(2).at(0) = 0; 
    scaleViewport.at(2).at(1) = (yh - yl)/2;
    scaleViewport.at(2).at(2) = 0;
    scaleViewport.at(2).at(3) = 0;

    scaleViewport.at(3).at(0) = 0;
    scaleViewport.at(3).at(1) = 0;
    scaleViewport.at(3).at(2) = 0;
    scaleViewport.at(3).at(3) = 0;

    step1 = multiplyMatrix(4,4,4,4, translationViewport2, scaleViewport);
    worldView = multiplyMatrix(4,4,4,4, step1, translationViewport1);

    viewport1 = multiplyMatrix(4,4,4,1, worldView, canonical1);
    //viewport2 = multiplyMatrix(3,3,3,1, worldView, canonical2D2);

    x1 = viewport1[0][0];
    y1 = viewport1[1][0];
    z1 = viewport1[2][0];
    // x2 = viewport2[0][0];
    // y2 = viewport2[1][0];


    // cout << "x1 viewport coords: " << x1 << endl;
    // cout << "y1 viewport coords: " << y1 << endl;
    // cout << "x2 viewport coords: " << x2 << endl;
    // cout << "y2 viewport coords: " << y2 << endl;

    coordinates.viewPushBack(x1); 
    coordinates.viewPushBack(y1);
    coordinates.viewPushBack(z1);
    // coordinates.viewPushBack(x2);
    // coordinates.viewPushBack(y2);
  }

  vector<float> faceCoords(8);
  // cout << "2d coordinate vector is like : "<< coordinates.faceSize() << endl;
  //go through faces and draw the connecting vertecies
  for(int i = 0; i != coordinates.faceSize(); i += 3){
    vector<int> clip(6);
    vector<float> translated(3), normals(3);
    vector<float> rotated(3);
    vector<float> scaled(3);
    vector<vector<float>> bf1, bf2, bf3;

    bf1.resize(3, vector<float>(1));
    bf2.resize(3, vector<float>(1));
    bf3.resize(3, vector<float>(1));


      // // cout << "x1 clip: " << clip.at(0) << endl;
      // // cout << "y1 clip: " << clip.at(1) << endl;
      // // cout << "x2 clip: " << clip.at(2) << endl;
      // // cout << "y2 clip: " << clip.at(3) << endl;
      // //swap coordinates if x1 > x2

      // //compute helper variables
      // x1 = clip.at(0);
      // y1 = clip.at(1);
      // x2 = clip.at(2);
      // y2 = clip.at(3);
      
      int vert1, vert2, vert3;
      vert1 = coordinates.getFaceIndex(i);
      vert2 = coordinates.getFaceIndex(i+1);
      vert3 = coordinates.getFaceIndex(i+2);

      //for each triangle decide if back face culling is necessary.  getCoordsIndex for 3d coordinates and then call getSurfaceNormal
      //If surface normal z value is < 0 discard face. else if the normal in world coordinates times the VPN is greater than 0 also discard face
      //1 -> 0,1,2
      //2 -> 3,4,5
      //3 -> 6,7,8
      //4 -> 9,10,11
      //3D coordinates of Face
      bf1[0][0] = coordinates.getCoordIndex(vert1*3 - 3);
      bf1[1][0] = coordinates.getCoordIndex(vert1*3 - 2);
      bf1[2][0] = coordinates.getCoordIndex(vert1*3 - 1);

      bf2[0][0] = coordinates.getCoordIndex(vert2*3 - 3);
      bf2[1][0] = coordinates.getCoordIndex(vert2*3 - 2);
      bf2[2][0] = coordinates.getCoordIndex(vert2*3 - 1);

      bf3[0][0] = coordinates.getCoordIndex(vert3*3 - 3);
      bf3[1][0] = coordinates.getCoordIndex(vert3*3 - 2);
      bf3[2][0] = coordinates.getCoordIndex(vert3*3 - 1);
      normals = getSurfaceNormal(bf1, bf2, bf3);


      //1 -> 0,1
      //2 -> 2,3
      //3 -> 4,5
      //4 -> 6,7
      //5 -> 8,9
      //6 -> 10,11
      //2D Coordinates of Face
      if(normals[2] > 0){
        faceCoords.push_back(coordinates.getViewIndex(vert1*3 - 3));
        faceCoords.push_back(coordinates.getViewIndex(vert1*3 - 2));
        cout << "the z value of vert1: " <<  coordinates.getViewIndex(vert1*3 - 1) << endl;
        if(zBuffer[coordinates.getViewIndex(vert1*3 - 3)][coordinates.getViewIndex(vert1*3 - 2)] < coordinates.getViewIndex(vert1*3 - 1)){
          zBuffer[coordinates.getViewIndex(vert1*3 - 3)][coordinates.getViewIndex(vert1*3 - 2)] = coordinates.getViewIndex(vert1*3 - 1);
        }

        faceCoords.push_back(coordinates.getViewIndex(vert2*3 - 3));
        faceCoords.push_back(coordinates.getViewIndex(vert2*3 - 2));
        cout << "the z value of vert2: " <<  coordinates.getViewIndex(vert2*3 - 1) << endl;
        if(zBuffer[coordinates.getViewIndex(vert2*3 - 3)][coordinates.getViewIndex(vert2*3 - 2)] < coordinates.getViewIndex(vert2*3 - 1)){
          zBuffer[coordinates.getViewIndex(vert2*3 - 3)][coordinates.getViewIndex(vert2*3 - 2)] = coordinates.getViewIndex(vert2*3 - 1);
        }

        faceCoords.push_back(coordinates.getViewIndex(vert3*3 - 3));
        faceCoords.push_back(coordinates.getViewIndex(vert3*3 - 2));
        cout << "the z value of vert3: " <<  coordinates.getViewIndex(vert3*3 - 1) << endl;
        if(zBuffer[coordinates.getViewIndex(vert3*3 - 3)][coordinates.getViewIndex(vert3*3 - 2)] < coordinates.getViewIndex(vert3*3 - 1)){
          zBuffer[coordinates.getViewIndex(vert3*3 - 3)][coordinates.getViewIndex(vert3*3 - 2)] = coordinates.getViewIndex(vert3*3 - 1);
        }
        faceCoords.push_back(coordinates.getViewIndex(vert1*3 - 3));
        faceCoords.push_back(coordinates.getViewIndex(vert1*3 - 2));
        // cout << "faceCoords is like : "<< faceCoords.size() << endl;

        //add z value to z buffer if it is less negative than the current z value (meaning it is closer to the camera)

      }
      
      Edge polygonEdges;

      for(int j = 0; j < faceCoords.size(); j+=2){

        if(j + 2 < faceCoords.size()){
          x1 = faceCoords.at(j);
          y1 = 500 - faceCoords.at(j+1);
          x2 = faceCoords.at(j+2);
          y2 = 500 - faceCoords.at(j+3);

          vector<float> endpoints;
          endpoints.push_back(x1);
          endpoints.push_back(y1);
          endpoints.push_back(x2);
          endpoints.push_back(y2);

          polygonEdges.setEndpoints(endpoints);

          // cout << "x1 plotting face: " << x1 << endl;
          // cout << "y1 plotting face: " << y1 << endl;
          // cout << "x2 plotting face: " << x2 << endl;
          // cout << "y2 plotting face: " << y2 << endl;

          dy = (y2-y1);
          dx = (x2-x1);

          if(x1 > x2){
            float tempX = x1;
            float tempY = y1;
            x1 = x2;
            x2 = tempX;
            y1 = y2;
            y2 = tempY;
            cout << "swapped" << endl;
          }



          // cout << "clip size: " << clip.size() << endl;
          // clip.clear();
          // cout << "x1 clear: " << clip.at(0) << endl;
          // cout << "y1 clear: " << clip.at(1) << endl;
          // cout << "x2 clear: " << clip.at(2) << endl;
          // cout << "y2 clear: " << clip.at(3) << endl;

          // cout << "x1: " << clipped[i] << endl;
          // cout << "y1: " << clipped[i+1] << endl;
          // cout << "x2: " << clipped[i+2] << endl;
          // cout << "y2: " << clipped[i+3] << endl;

          // cout << "dx: " << dx << endl;
          // cout << "dy: " << dy << endl;
          cout << "line from (" << (int)x1 << ","<< (int)y1 << ") to ("  << (int)x2 << "," << (int)y2 << ")." << endl;



          if(dx != 0){
            m = (float)(dy/dx);
            cout << "m: " << m << endl << endl;
            yInter = (y1 - x1*m);
            // cout << "yInter: " << yInter << endl; 
          }

          else{
            if(y1 > y2){
              float tempX = x1;
              float tempY = y1;
              x1 = x2;
              x2 = tempX;
              y1 = y2;
              y2 = tempY;
              // cout << "swapped" << endl;
            }
            m = NAN;
            cout << "m: " << m << endl << endl; 
             for(int j = 0; y1 < y2; j+=4){

              drawX.push_back(x1);
              drawY.push_back(y1);
              if(loop == 1){
                xpm.at(x1).at(y1) = '-';
              }
              if(loop == 2){
                xpm.at(x1).at(y1) = '+';
              }
              if(loop == 3){
                xpm.at(x1).at(y1) = '*';
              }
              

              y1++;
              }
          // cout << endl << "line created dx = 0" << endl << endl;

          }
          if(m == 0){  
            //if 0 <= m <= 1
            if(x1 > x2){
              float tempX = x1;
              float tempY = y1;
              x1 = x2;
              x2 = tempX;
              y1 = y2;
              y2 = tempY;
              // cout << "swapped" << endl;
            }
            while(x1 < x2){
                drawX.push_back(x1 + 1);
                drawY.push_back(y1);
                
                if(loop == 1){
                xpm.at(x1+ 1).at(y1) = '-';
              }
              if(loop == 2){
                xpm.at(x1+ 1).at(y1) = '+';
              }
              if(loop == 3){
                xpm.at(x1+ 1).at(y1) = '*';
              }

                x1++;
              }
            // cout << endl << "line created m == 0" << endl << endl;
          }

          if(m > 0 && m <= 1){  
            //if 0 <= m <= 1
            if(x1 > x2){
              float tempX = x1;
              float tempY = y1;
              x1 = x2;
              x2 = tempX;
              y1 = y2;
              y2 = tempY;
              // cout << "swapped" << endl;
              // dy = (y2-y1);
              // dx = (x2-x1);
              // cout << "new dx: " << dx << endl;
              // cout << "new dy: " << dy << endl;
            }
            rx = x2;
            ry = y2;
            
            while(x1 < x2){
              rx = (y1-yInter)/m;
              ry = y1;

              dex = rx - x1;
              dey = ry - y1;
              int decisionC = ((2*dey) - dex + (2*(dex*ry - dey*rx)));
              int decision = (2*dey*x1) - (2*dex*y1) + decisionC;
              // cout << "decision " << decision << endl;
              // cout << "x1 " << x1 << endl;
              // cout << "y1 " << y1 << endl;
              
              if(decision >= 0){
                drawX.push_back(x1 + 1);
                drawY.push_back(y1 + 1);
                
              if(loop == 1){
                xpm.at(x1+1).at(y1+1) = '-';
              }
              if(loop == 2){
                xpm.at(x1+1).at(y1+1) = '+';
              }
              if(loop == 3){
                xpm.at(x1+1).at(y1+1) = '*';
              }
                // draw.push_back(1);
                y1++;
                x1++;
              }
              else{
                drawX.push_back(x1 + 1);
                drawY.push_back(y1);
                
                if(loop == 1){
                  xpm.at(x1+1).at(y1) = '-';
                }
                if(loop == 2){
                  xpm.at(x1+1).at(y1) = '+';
                }
                if(loop == 3){
                  xpm.at(x1+1).at(y1) = '*';
                }
                  // draw.push_back(1);
                  x1++;
              }

            }
          }
          
            // cout << endl << "line created 0 < m < 1" << endl << endl;
          }

          if( m > 1 ){
            //m > 1
            if(y1 > y2){
              float tempX = x1;
              float tempY = y1;
              x1 = x2;
              x2 = tempX;
              y1 = y2;
              y2 = tempY;
              // cout << "swapped" << endl;
              dy = (y2-y1);
              dx = (x2-x1);
              // cout << "new dx: " << dx << endl;
              // cout << "new dy: " << dy << endl;
            }
            rx = x1;
            ry = y1;
            

            while(y1 < y2){
              rx = (y1-yInter)/m;
              ry = y1;

              dex = rx - x1;
              dey = ry - y1;
              int decisionC = dey - 2*dex + 2*(dex*ry - dey*rx);
              int decision = (2*dey*x1) - (2*dex*y1) + decisionC;

              if(decision <= 0){
                drawX.push_back(x1 + 1);
                drawY.push_back(y1 + 1);
                
                if(loop == 1){
                  xpm.at(x1+1).at(y1+1) = '-';
                }
                if(loop == 2){
                  xpm.at(x1+1).at(y1+1) = '+';
                }
                if(loop == 3){
                  xpm.at(x1+1).at(y1+1) = '*';
                }
                  // draw.push_back(1);
                  y1++;
                  x1++;
              }
              else{
                drawX.push_back(x1);
                drawY.push_back(y1 + 1);
                
                if(loop == 1){
                xpm.at(x1).at(y1+1) = '-';
              }
              if(loop == 2){
                xpm.at(x1).at(y1+1) = '+';
              }
              if(loop == 3){
                xpm.at(x1).at(y1+1) = '*';
              }
                // draw.push_back(1);
                y1++;
              }

            }
              // cout << endl << "line created m > 1" << endl << endl;
          }

          if(m >= -1 && m < 0){
            //-1 < m < 0
            if(x1 > x2){
              float tempX = x1;
              float tempY = y1;
              x1 = x2;
              x2 = tempX;
              y1 = y2;
              y2 = tempY;
              // cout << "swapped" << endl;
            }
            rx = x2;
            ry = y2;
            
            while(x1 < x2){
              rx = (y1-yInter)/m;
              ry = y1;

              dex = rx - x1;
              dey = ry - y1;
              int decisionC = (2*dey + dex +2*(dex*ry - dey*rx));
              int decision = 2*dey*x1 - 2*dex*y1 + decisionC;
              // cout << "decision: " << decision << endl << endl;

              if(decision <= 0){
                drawX.push_back(x1 + 1);
                drawY.push_back(y1 - 1);
                
                if(loop == 1){
                xpm.at(x1+1).at(y1-1) = '-';
              }
              if(loop == 2){
                xpm.at(x1+1).at(y1-1) = '+';
              }
              if(loop == 3){
                xpm.at(x1+1).at(y1-1) = '*';
              }
                // draw.push_back(1);
                y1--;
                x1++;
              }
              else{
                drawX.push_back(x1 + 1);
                drawY.push_back(y1);
                
                if(loop == 1){
                xpm.at(x1+1).at(y1) = '-';
              }
              if(loop == 2){
                xpm.at(x1+1).at(y1) = '+';
              }
              if(loop == 3){
                xpm.at(x1+1).at(y1) = '*';
              }
                // draw.push_back(1);
                x1++;
              }
            }
            // cout << endl << "line created -1 < m < 0" << endl << endl;
          }

          if(m < -1){
            // m < -1
            if(y1 < y2){
              float tempX = x1;
              float tempY = y1;
              x1 = x2;
              x2 = tempX;
              y1 = y2;
              y2 = tempY;
              // cout << "swapped" << endl;
              // dy = (y1-y2);
              // dx = (x1-x2);
              // cout << "new dx: " << dx << endl;
              // cout << "new dy: " << dy << endl;
            }
            rx = x1;
            ry = y1;


            

            while(y1 > y2){

              rx = (y1-yInter)/m;
              ry = y1;

              dex = rx - x1;
              dey = ry - y1;

              // cout << "rx: " << rx << endl;
              // cout << "ry: " << ry<< endl;
              // cout << "dex: " << dex << endl;
              // cout << "dey: " << dey << endl;
              // cout << "x1: " << x1 << endl;
              // cout << "y1: " << y1 << endl;


              int decisionC = 2*dey - dex + 2*(dex*ry - dey*rx);
              // cout << "decision Constant: " << decisionC << endl;
              // int decision = decY + (2 * decX);
              int decision = (2*dey*x1) - (2*dex*y1) + decisionC;

              // cout << "decision: " << decision << endl;

              if(decision <=  0){
                // drawX.push_back(x1++);
                // drawY.push_back(y1--);
                
                if(loop == 1){
                  xpm.at(x1++).at(y1--) = '-';
                }
                if(loop == 2){
                  xpm.at(x1++).at(y1--) = '+';
                }
                if(loop == 3){
                  xpm.at(x1++).at(y1--) = '*';
                }
                // draw.push_back(1);
              }
             else{
                // drawX.push_back(x1);
                // drawY.push_back(y1--);
                
                if(loop == 1){
                  xpm.at(x1).at(y1--) = '-';
                }
                if(loop == 2){
                  xpm.at(x1).at(y1--) = '+';
                }
                if(loop == 3){
                  xpm.at(x1).at(y1--) = '*';
                }
                  // draw.push_back(1);
              }

          }

        // cout << endl << "line created m < -1" << endl << endl;
          

      }
      coordinates.setEdges(polygonEdges);
      // xpm = fillPoly(xpm, zBuffer, coordinates);
      coordinates.clearEdges();
    }
    
    faceCoords.clear();


      

    //write to xpm file
    
  }
      if(loop <= loopTotal){
      cout << "loop total " << loopTotal << endl;
      cout << "loop " << loop << endl;
      loop++;
    }
  }
  

    // cout << "we're setting xpm array from drawX.size(): " << xpm.size() << endl;
    // cout << "we're setting xpm array from drawY.size(): " << xpm[0].size() << endl;
   //set xpm array to have - for pixel locations

  // for(int i = 0; i < drawX.size(); i++){
  //   if(drawX[i] < xpm.size() && drawY[i] < xpm[0].size()){
  //       // cout << "how the f did this happen?" << endl;
  //       // cout << "drawX[i] drawY[i]: " << drawX[i] << " " << drawY[i] << endl;
  //       // vector<float> p1(3), p2(3), p3(3);
  //       // p1.push_back(coordinates.getViewIndex(vert1*3 - 3));
  //       // p1.push_back(coordinates.getViewIndex(vert1*3 - 2));
  //       // p1.push_back(coordinates.getViewIndex(vert1*3 - 1));

  //       // p2.push_back(coordinates.getViewIndex(vert2*3 - 3));
  //       // p2.push_back(coordinates.getViewIndex(vert2*3 - 2));
  //       // p2.push_back(coordinates.getViewIndex(vert2*3 - 1));

  //       // p3.push_back(coordinates.getViewIndex(vert3*3 - 3));
  //       // p3.push_back(coordinates.getViewIndex(vert3*3 - 2));
  //       // p3.push_back(coordinates.getViewIndex(vert3*3 - 1));

  //       // fillPoly(p1, p2, p3, xpm);
  //       xpm.at(drawY[i]).at(drawX[i]) = '-';
  //   }
    
  // }
  // cout << "displaying z buffer" << endl;

  // for(int i = 0; i < zBuffer.size(); i++){
  //       cout << "some values in zbuffer: " << zBuffer[i][i] << " " << endl;
  // }

  //Fill polygons
  // xpm = fillPoly(xpm, zBuffer, coordinates);

  writeToBuffer(drawX, drawY, filename1, xpm, framew, frameh, loopTotal);
  return 0;
}

//PARSE PS FILE
//edited to take in new ps format
Face parsePS(string filename, int framew, int frameh, int loop){

  // cout << "filename: " << filename << endl;
  ifstream infile;
  infile.open(filename.c_str());
  string line;
  vector<float> faces;
  vector<float> coords;
  Face allFaces;
  vector<Face> faceInfo;

  while(getline(infile, line))
  {
    istringstream iss(line);
    float n, m, u;
    string instruct, next;

    bool nums = static_cast<bool>(iss >> instruct >>  n >> m >> u );
    // cout << nums << endl;
    

      // if(word){
      //   cout << instruct<< endl;
      //   word = false;
      // }
      if(nums){
        // cout << n << " " << m << " " << u << " " << instruct << endl;
        if(instruct == "v"){
          coords.push_back(n);
          coords.push_back(m);
          coords.push_back(u);
        }
        if(instruct == "f"){
          faces.push_back(n);
          faces.push_back(m);
          faces.push_back(u);
        }
        
        nums = false;
      }

  }
  allFaces.setFaces(faces);
  allFaces.setCoords(coords);
  // cout << "did the coordinates actually get set to Face? " << allFaces.getCoordIndex(1) << endl;
  return allFaces;
}

vector<float> getSurfaceNormal(vector<vector<float>> p1, vector<vector<float>> p2, vector<vector<float>> p3){
  // Set Vector U to (Triangle.p2 minus Triangle.p1)
  // Set Vector V to (Triangle.p3 minus Triangle.p1)

  // Set Normal.x to (multiply U.y by V.z) minus (multiply U.z by V.y)
  // Set Normal.y to (multiply U.z by V.x) minus (multiply U.x by V.z)
  // Set Normal.z to (multiply U.x by V.y) minus (multiply U.y by V.x)

  // Returning normal

  vector<float> u(3), v(3), normal(3);

  u[0] = p2[0][0] - p1[0][0]; 
  u[1] = p2[1][0] - p1[1][0]; 
  u[2] = p2[2][0] - p1[2][0]; 

  v[0] = p3[0][0] - p1[0][0]; 
  v[1] = p3[1][0] - p1[1][0]; 
  v[2] = p3[2][0] - p1[2][0]; 

  normal[0] = (u[1]*v[2]) - (u[2]*v[1]); 
  normal[1] = (u[2]*v[0]) - (u[0]*v[2]); 
  normal[2] = (u[0]*v[1]) - (u[1]*v[0]); 


  return normal;
}

//LINE DRAWING CALCULATIONS
//clip coordinates to viewing plane
// vector<int> clipCords(int x1, int y1, int x2, int y2, int code1, int code2, int code3, int xl, int yl, int xh, int yh){
//   int accept;
//   vector<int> clipped(6);
//   accept = 0;   //decides if line is to be drawn
    
//     while(1){
//       float m =(float)(y2-y1)/(x2-x1);
//       //Both points inside. Accept line
//       cout << "code1 && code2 == 0: " << code1 << " " << code2 << endl;
//       if(code1==0 && code2==0 && code3==0){ 
//         // cout << "points are inside " << code1 << " AND " << code2 << endl;
//         // cout << "x1: " << x1 << endl;
//         // cout << "y1: " << y1 << endl;
//         // cout << "x2: " << x2 << endl;
//         // cout << "y2: " << y2 << endl;
//         clipped[0] = x1;
//         clipped[1] = y1;
//         clipped[2] = x2;
//         clipped[3] = y2;
//         clipped[4] = z1;
//         clipped[5] = z2;
//         // cout << "x1: " << clipped[0] << endl;
//         // cout << "y1: " << clipped[1] << endl;
//         // cout << "x2: " << clipped[2] << endl;
//         // cout << "y2: " << clipped[3] << endl;
//         // cout << "clip index 2 and 3: " << clipped[2] << " " << clipped[3] << endl;
//         accept = 1;
//         break;
//       }
//       //AND of both codes != 0.Line is outside. Reject line
//       else if((code1 & code2)!=0){
//         cout << " points are outside " << code1 << " AND " << code2 << endl;
//         clipped[0] = 0;
//         clipped[1] = 0;
//         clipped[2] = 0;
//         clipped[3] = 0;
//         break;
//       }else{
//         int x,y;
//         int temp;
//         //Decide if point1 is inside, if not, calculate intersection
//         if(code1==0){
//           temp = code2;
//           if(temp & TOP){
//             cout << "x2 y2 clips top edge" << endl;
//             // cout << "x1: " << x1 << endl;
//             // cout << "y1: " << y1 << endl;
//             // cout << "x2: " << x2 << endl;
//             // cout << "y2: " << y2 << endl;
//             x = x2 + (yh-y2)/m;
//             cout << "m:           " << m << endl;
//             y = yh;
//             cout << "clipped x and y: " << x << " " << y << endl;
//           }
//           else if(temp & BOTTOM){   //Line clips bottom edge
//             cout << "line clips bottom edge" << endl;
//             x = x2 + (yl-y2)/m;
//             y = yl;
//           }else if(temp & LEFT){  //Line clips left edge
//             cout << "line clips left edge" << endl;
//             x = xl;
//             y = y2 + m*(xl-x2);
//           }else if(temp & RIGHT){   //Line clips right edge
//             cout << "line clips right edge" << endl;
//             x = xh;
//             y = y2 + m*(xh-x2);
//           }
//           //Check which point we had selected earlier as temp, and replace its co-ordinates
//           if(temp == code1){ 
//             // cout << "temp was x2 y2" << endl;
//             x2 = x;
//             y2 = y;
//             code1 = getcode(x,y, xl, yl, xh, yh);
//           }else{
//             // cout << "temp was x2 y2" << endl;
//             x2 = x;
//             y2 = y;
//             // cout << "clip index 2 and 3: " << clipped[2] << " " << clipped[3] << endl;
//             code2 = getcode(x,y, xl, yl, xh, yh);
//           }
//         }
//         else{
//           temp = code1;
//           if(temp & TOP){
//             cout << "x1 y1 line clips top edge" << endl;
//             // cout << "x1: " << x1 << endl;
//             // cout << "y1: " << y1 << endl;
//             // cout << "x2: " << x2 << endl;
//             // cout << "y2: " << y2 << endl;
//             x = x1 + (yh-y1)/m;
//             cout << "m:           " << m << endl;
//             y = yh;
//             cout << "clipped x and y: " << x << " " << y << endl;
//           }
//           else if(temp & BOTTOM){   //Line clips bottom edge
//             cout << "line clips bottom edge" << endl;
//             x = x1 + (yl-y1)/m;
//             y = yl;
//           }else if(temp & LEFT){  //Line clips left edge
//             cout << "line clips left edge" << endl;
//             x = xl;
//             y = y1 + m*(xl-x1);
//           }else if(temp & RIGHT){   //Line clips right edge
//             cout << "line clips right edge" << endl;
//             x = xh;
//             y = y1 + m*(xh-x1);
//           }
//           //Check which point we had selected earlier as temp, and replace its co-ordinates
//           if(temp == code1){ 
//             // cout << "temp was x1 y1" << endl;
//             x1 = x;
//             y1 = y;
//             code1 = getcode(x,y, xl, yl, xh, yh);
//           }else{
//             // cout << "temp was x2 y2" << endl;
//             x2 = x;
//             y2 = y;
//             // cout << "clip index 2 and 3: " << clipped[2] << " " << clipped[3] << endl;
//             code2 = getcode(x,y, xl, yl, xh, yh);
//           }
//         }
          
//       }
//         //Line clips top edge
//     }
//   // for(int z = 0; z < clipped.size(); z++){
//     // if(clipped[z] > 500){
//     //   cout << "this isn't clipped " << clipped[z] << endl;
//     //       cout << "x1: " << clipped[0] << endl;
//     //       cout << "y1: " << clipped[1] << endl;
//     //       cout << "x2: " << clipped[2] << endl;
//     //       cout << "y2: " << clipped[3] << endl;
//     // }
//   // }
//     return clipped;
// }

// int getcode(int x,int y, int z){
//   int code = 0;

//     // cout << "x: " << x << endl;
//     // cout << "y: " << y << endl;
//   //Perform Bitwise OR to get outcode
//   if(y > 1) {
//     code |=TOP;
//     // cout << "code: " << code << endl;
//   }
//   if(y < -1) {
//     code |=BOTTOM;
//     // cout << "code: " << code << endl;
//   }
//   if(x < -1) {
//     code |=LEFT;
//     // cout << "code: " << code << endl;
//   }
//   if(x > 1) {
//     code |=RIGHT;
//     // cout << "code: " << code << endl;
//   }
//   if(z < -1) {
//     code |=NEAR;
//     // cout << "code: " << code << endl;
//   }
//   if(z > 0) {
//     code |=FAR;
//     // cout << "code: " << code << endl;
//   }
  
//   return code;
// }

//WRITE TO FRAME BUFFER (XPM FILE)
void writeToBuffer(vector<int> pixelsX, vector<int> pixelsY, string filename, vector<vector<char>> xpm, int framew, int frameh, int loop){
  ofstream ofs;
  string prefix;
  prefix = filename.erase(filename.size() - 3, filename.size());
  cout << prefix << endl;
  string output = prefix + "xpm";
  ofs.open(output.c_str(), std::ofstream::out | std::ofstream::app);
  int xIter, yIter;
  xIter = 0;
  yIter = 0;

  ofs << "/* XPM /*" << endl << "static char *sco100[] = {" << endl << endl << "/* width height num_colors chars_per_pixel/*" << endl << 
  "\"" << framew  << " " << frameh << " 4 1\"," << endl << endl << "/* colors /*"  << endl << "\"- c #ff0000\"," << endl << "\"X c #ffffff\"," << endl
  << "\"+ c #00ff00\"," << endl << "\"* c #0000ff\"," << endl << endl << "/* pixels /*" << endl <<
  "\"";

  //print xpm values row by row
  for(int i = 0; i < framew; i++){
    for(int j = 0; j < frameh; j++){
      // cout << "i j: " << i << " " << j << endl;
      ofs << xpm.at(j).at(i);
    }
    ofs << "\"," << endl << "\"";
  }



  ofs.close();
}

//need to store the x,y values of the whole line segment that I'm passing in
vector<vector<char>> fillPoly(vector<vector<char>> xpm, vector<vector<float>> zBuffer, Face polygon){
  int za, zb, zp, xa, xb;
  //start scan fill by checking each pixel in clipping plane for number of intersections
  //if isFilled (odd number of extrema crossed) then add points to xpm array
  bool isFilled = false;
  int parity;
  for(int n = 0; n < 500; n++){
    parity = 0;
    for(int m = 0; m < 500; m++){
      //using current line check against each edge in polygon for intersection
      //if there is an intersection store that as an extrema
      Edge edge1, edge2, edge3;
      edge1 = polygon.getEdges(0);
      edge2 = polygon.getEdges(1);
      edge3 = polygon.getEdges(2);

      

      
      //test pixel (n,m) for intersection
      if(xpm[n][m] == '-'){
        parity++;
        xa = m;
      }
      if(xpm[n][m] == '-'){
        parity++;
      }
      isFilled = parity % 2;

      if(isFilled){
        // calculate extrema points
      // za = p1[2] - (p1[2] - p2[2]) * (p1[1] - n)/(p1[1] - p2[1]);
      // zb = p1[2] - (p1[2] - p3[2]) * (p1[1] - n)/(p1[1] - p3[1]);
      // zp = zb - (zb - za) * (xb - m)/(xb - xa);
      // cout << "zp: " << zp << endl;
        xpm[n][m] = '-';
      
      }
    }
  }
  return xpm;
}


float degToRad(int degrees){
  float rads;
  rads = degrees * PI/180;
  // cout << "radians: " << rads << endl;
  return rads;
}

