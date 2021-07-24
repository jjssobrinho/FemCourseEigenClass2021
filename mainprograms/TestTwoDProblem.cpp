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

int main (int argc, char *argv[])
{
    
    GeoMesh gmesh;
    ReadGmsh readingGmsh;
    std::string dirName = "../../tests/";
    std::string filepath = dirName + argv[1];
    std::string outputFile = filepath.substr(0, filepath.find_last_of("."));
    outputFile = outputFile +".vtk";

    std::cout << outputFile << std::endl;
    
    cout << "Reading Gmsh file: " << filepath <<  endl;
    readingGmsh.Read(gmesh, filepath);

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
        //res[0] = 2.*(1.-x[0])*x[0] + 2.*(1.-x[1])*x[1]; //Laplacian
        res[0] = (-200.*(5.*(-1. + 2*x[0])*(-1 + x[1])*x[1]*cos(5*x[0])*cos(4*x[1]) + 
                sin(5*x[0])*(-((-1 + x[1])*(-x[1] - 2*x[0]*(1. + 10*x[1]) + std::pow(x[0],2)*(2. + 20*x[1]))*
                cos(4*x[1])) + 4.*(-1. + x[0])*x[0]*(1. - 3.*x[1] + std::pow(x[1],2))*sin(4.*x[1]))))*exp(-x[1]);
    };

    mat1->SetForceFunction(force);

    MatrixDouble proj(1,1), val1(1,1), val2(1,1);
    proj.setZero();
    val1.setZero();
    val2.setZero();
    L2Projection *bc_linha = new L2Projection(0, 2, proj, val1, val2);
    L2Projection *bc_point = new L2Projection(0, 3, proj, val1, val2);
    std::vector<MathStatement *> mathvec = {0, mat1, bc_linha, bc_point};
    cmesh.SetMathVec(mathvec);
    cmesh.SetDefaultOrder(1);
    cmesh.AutoBuild();
    cmesh.Resequence();
    
    Analysis locAnalysis(&cmesh);
    locAnalysis.RunSimulation();
    PostProcessTemplate<Poisson> postprocess;
    auto exact = [](const VecDouble &x, VecDouble &val, MatrixDouble &deriv){
        val[0] = 100.*(sin(5.*x[0])*exp(-x[1])*cos(4.*x[1]))*(x[0]*(1.-x[0])*x[1]*(1.-x[1]));
        deriv(0,0) = (100.*(-1. + x[1])*x[1]*cos(4.*x[1])*(5.*(-1. + x[0])*x[0]*cos(5.*x[0]) + 
                    (-1. + 2*x[0])*sin(5.*x[0])))*exp(-x[1]);
        deriv(1,0) = (-100.*(-1. + x[0])*x[0]*sin(5.*x[0])*((1. - 3.*x[1] + std::pow(x[1],2))*cos(4.*x[1]) + 
                    4.*(-1. + x[1])*x[1]*sin(4.*x[1])))*exp(-x[1]);
    };
    postprocess.AppendVariable("Sol");
    postprocess.AppendVariable("DSol");
    postprocess.AppendVariable("Flux");
    postprocess.AppendVariable("Force");
    postprocess.AppendVariable("SolExact");
    postprocess.AppendVariable("DSolExact");
    postprocess.SetExact(exact);
    mat1->SetExactSolution(exact);
    locAnalysis.PostProcessSolution(outputFile ,postprocess);

    VecDouble errvec;
    errvec = locAnalysis.PostProcessError(std::cout, postprocess);

    return 0;
}