
///\cond
#include <iostream>
#include <math.h>
///\endcond
#include "Topology1d.h"
#include "TopologyTriangle.h"
#include "TopologyQuad.h"
#include "TopologyTetrahedron.h"
#include "Analysis.h"
#include "VTKGeoMesh.h"
#include "GeoElement.h"
#include "ReadGmsh.h"
#include "CompMesh.h"
#include "CompElement.h"
#include "Poisson.h"
#include "L2Projection.h"
#include "Assemble.h"

int main (){

    GeoMesh gmesh;
    ReadGmsh readingGmsh;
    std::string filename = "/home/jdasilva/triangle.msh";
    std::string outputFile = "/home/jdasilva/triangle-mesh.vtk";

    std::cout << "Reading Gmsh file: " << filename << 
        " and converting to VTK." << std::endl;
    readingGmsh.Read(gmesh, filename);
    VTKGeoMesh plotmesh;
    plotmesh.PrintGMeshVTK(&gmesh, outputFile);
    
    // Coputational mesh:
    CompMesh cmesh(&gmesh);
    MatrixDouble perm(3,3);
    perm.setZero();
    perm(0,0) = 1.0;
    perm(1,1) = 1.0;
    perm(2,2) = 1.0;
    Poisson *mat1 = new Poisson(3, perm);
    MatrixDouble proj(1,1), val1(1,1), val2(1,1);
    proj.setZero();
    val1.setZero();
    val2.setZero();
    L2Projection *bc_linha = new L2Projection(0, 2, proj, val1, val2);
    L2Projection *bc_point = new L2Projection(0, 1, proj, val1, val2);
    std::vector<MathStatement *> mathvec = {0, bc_point, bc_linha, mat1};
    cmesh.SetMathVec(mathvec);
    cmesh.AutoBuild();
    cmesh.Resequence();

    Assemble assemble(&cmesh);
    auto neq = assemble.NEquations();
    MatrixDouble globmat(neq, neq), rhs(neq,1);
    assemble.Compute(globmat, rhs);
    
    return 0;
}
