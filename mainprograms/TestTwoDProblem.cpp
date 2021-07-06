//
//  TestOneDProblem.cpp MODIFICADO DO ORIGINAL
//  FemSC
//
//  Created by Eduardo Ferri on 08/17/15.
//
//
//TestOneDProblem cpp
/*
 Os testes foram preparados com um proposito educacional,
 recomenda-se que o aluno entenda a funcionalidade de cada
 teste e posteriormente use com seu cÛdigo caso a caso
 */
//      Obs: O xmax e xmin estao tomados como 4 e 0, respectivamente,
//      caso estes valores sejam alterados, editar o teste TestNodes.
//
//
#include <iostream>
#include <math.h>

#include "GeoMesh.h"
#include "ReadGmsh.h"
#include "CompMesh.h"
#include "Poisson.h"
#include "L2Projection.h"
#include "Analysis.h"
#include "PostProcessTemplate.h"

using std::cout;
using std::endl;
using std::cin;

int main ()
{
    
    GeoMesh gmesh;
    ReadGmsh readingGmsh;
    std::string filename = "/home/jdasilva/Desktop/triangle.msh";
    //std::string outputFile = "quads-mesh.vtk";
    
    cout << "Reading Gmsh file: " << filename <<  endl;
    readingGmsh.Read(gmesh, filename);

    // Computational mesh:
    CompMesh cmesh(&gmesh);
    MatrixDouble perm(3,3);
    perm.setZero();
    perm(0,0) = 1.0;
    perm(1,1) = 1.0;
    perm(2,2) = 1.0;
    Poisson *mat1 = new Poisson(1, perm);
    mat1->SetDimension(2);

    //Force vector:
    auto force = [](const VecDouble &x, VecDouble &res){
        res[0] = 2.*(1.-x[0])*x[0] + 2.*(1.-x[1])*x[1]; //Laplacian
    };
    mat1->SetForceFunction(force);

    MatrixDouble proj(1,1), val1(1,1), val2(1,1);
    proj.setZero();
    val1.setZero();
    val2.setZero();
    L2Projection *bc_linha = new L2Projection(0, 2, proj, val1, val2);
    L2Projection *bc_point = new L2Projection(0, 3, proj, val1, val2);
    std::vector<MathStatement *> mathvec = {0, bc_point, bc_linha, mat1};
    cmesh.SetMathVec(mathvec);
    cmesh.SetDefaultOrder(2);
    cmesh.AutoBuild();
    cmesh.Resequence();
    
    Analysis locAnalysis(&cmesh);
    locAnalysis.RunSimulation();
    PostProcessTemplate<Poisson> postprocess;
    auto exact = [](const VecDouble &x, VecDouble &val, MatrixDouble &deriv){
        val[0] = (1.-x[0])*x[0]*(1.-x[1])*x[1];
        deriv(0,0) = (1.-2.*x[0])*(1.-x[1])*x[1];
        deriv(1,0) = (1.-2.*x[1])*(1.-x[0])*x[0];
    };
    postprocess.AppendVariable("Sol");
    postprocess.AppendVariable("DSol");
    postprocess.AppendVariable("Flux");
    postprocess.AppendVariable("Force");
    postprocess.AppendVariable("SolExact");
    postprocess.AppendVariable("DSolExact");
    postprocess.SetExact(exact);
    mat1->SetExactSolution(exact);
    locAnalysis.PostProcessSolution("triangle.vtk",postprocess);

    VecDouble errvec;
    errvec = locAnalysis.PostProcessError(std::cout, postprocess);

    return 0;
}