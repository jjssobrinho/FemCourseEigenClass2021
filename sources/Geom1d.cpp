/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "Geom1d.h"

Geom1d::Geom1d() {
}

Geom1d::~Geom1d() {
}

Geom1d::Geom1d(const Geom1d &copy) {
    fNodeIndices = copy.fNodeIndices;
}

Geom1d& Geom1d::operator=(const Geom1d& copy) {
    fNodeIndices = copy.fNodeIndices;
    return *this;
}

void Geom1d::Shape(const VecDouble &xi, VecDouble &phi, MatrixDouble &dphi) {
    if(xi.size() <= 0 || xi.size() > Dimension) DebugStop();

    phi.resize(2);
    dphi.resize(1,2);

    phi[0]  = (1.-xi[0])*0.5;
    phi[1]  = (1.+xi[0])*0.5;
    dphi(0,0) = -0.5; 
    dphi(0,1) = 0.5; 
}

void Geom1d::X(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x) {
    if (NodeCo.rows() <= 0 || NodeCo.rows() > 3) DebugStop();
    if (NodeCo.cols() <= 0 || NodeCo.cols() > NumNodes()) DebugStop();
    // if (xi.size() != Dimension) DebugStop();
    int dim = NodeCo.rows();
    if (x.size() <= 0 || x.size() > dim) DebugStop();
    
    int nrow = NodeCo.rows();
    
    VecDouble phi;
    MatrixDouble dphi;
    Shape(xi, phi, dphi);
    x.resize(dim);
    x.setZero();
    int nnodes = NumNodes();

    for (int i = 0; i < nrow; i++){
        for (int j = 0; j < nnodes; j++) {
            x[i] += NodeCo(i,j)*phi[j];
        }
    }
}

void Geom1d::GradX(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x, 
                   MatrixDouble &gradx) {
    if (NodeCo.rows() <= 0 || NodeCo.rows() > 3) DebugStop();
    if (NodeCo.cols() <= 0 || NodeCo.cols() > NumNodes()) DebugStop();
    int dim = NodeCo.rows();
    if (x.size() <= 0 || x.size() > dim) DebugStop();
    
    VecDouble phi;
    MatrixDouble dphi;
    Shape(xi, phi, dphi);
    int nnodes = NumNodes();
    int masterdim = Dimension;
    x.resize(dim);
    x.setZero();
    gradx.resize(dim,masterdim);
    gradx.setZero();

    for (int k = 0; k < nnodes; k++){
        for (int i = 0; i < dim; i++){
            x[i] += phi[k]*NodeCo(i, k);
            for (int j = 0; j < masterdim; j++)gradx(i,j) += NodeCo(i, k)*dphi(j, k);
        }
    }
}

void Geom1d::SetNodes(const VecInt &nodes) {
    if(nodes.rows() != 2) DebugStop();
    fNodeIndices = nodes;
}

void Geom1d::GetNodes(VecInt &nodes) const {
    nodes = fNodeIndices;
}

int Geom1d::NodeIndex(int node) const {
    if(node<0 || node > 2) DebugStop();
    return fNodeIndices[node];
}

int Geom1d::NumNodes(){
    return nCorners;    
}

GeoElementSide Geom1d::Neighbour(int side) const{
    if(side <0 || side>2) DebugStop();
    return fNeighbours[side];
}

void Geom1d::SetNeighbour(int side, const GeoElementSide &neighbour) {
    if(side < 0 || side > 2) DebugStop();
    fNeighbours[side]=neighbour;
}
