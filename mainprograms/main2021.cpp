

#include <iostream>
#include <math.h>
#include "IntRule.h"
#include "IntRule1d.h"
#include "IntRuleQuad.h"
#include "IntRuleTetrahedron.h"
#include "IntRuleTriangle.h"
#include "Topology1d.h"
#include "TopologyTriangle.h"
#include "TopologyQuad.h"
#include "TopologyTetrahedron.h"
#include "DataTypes.h"
#include "Analysis.h"
#include "VTKGeoMesh.h"
#include "ReadGmsh.h"

using std::cout;
using std::endl;
using std::cin;

int main (){

    ReadGmsh readingGmsh;
    GeoMesh mesh;
    std::string filename = "/home/jdasilva/test-gmsh.msh";
    std::string outputFile = "/home/jdasilva/mesh.vtk";
    VTKGeoMesh vtk;

    cout << "Reading Gmsh file: " << filename << " and converting to VTK." << endl;
    readingGmsh.Read(mesh, filename);
    vtk.PrintGMeshVTK(&mesh, outputFile);

    return 0;
}
