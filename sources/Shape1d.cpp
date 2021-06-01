//
//  Shape1d.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include <cmath>
#include <math.h>
#include "tpanic.h"
#include "DataTypes.h"
#include "Shape1d.h"

void Shape1d::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, MatrixDouble &dphi){
    
    if (orders[0] < 0 || orders[1] < 0 || orders[2] < 0) {
        std::cout << "Shape1d::Shape: Invalid dimension for arguments: order\n";
        DebugStop();
    }
    if (orders[0] > 1 || orders[1] > 1) {
        std::cout << "Shape1d::Shape: Invalid dimension for arguments: order\n";
        DebugStop();
    }
    if (orders[2] > 2) {
        std::cout << "Shape1d::Shape: Please implement it for order > 2\n";
        DebugStop();
    }
    
    auto nshape = NShapeFunctions(orders);
    //std::cout << "nshape = " << nshape << std::endl;
    // outputs were: 2, 2, 3, 3 for each test called
    phi.resize(nshape);
    dphi.resize(1,nshape);
        
    //std::cout << "Please implement me\n";
    double csi = xi[0];
    phi[0] = 0.5*(1.-csi);
    phi[1] = 0.5*(1.+csi);
    dphi(0,0) = -0.5;
    dphi(0,1) = 0.5;

    if (orders[2] == 2){
        phi[2] = 1. - csi*csi;
        dphi(0,2) = -2.*csi;
    }

}

/// returns the number of shape functions associated with a side
int Shape1d::NShapeFunctions(int side, int order){

    if(order < 1 || order >2) DebugStop();
    switch (side)
    {
    case 0:
        return 1;
        break;
    case 1:
        return 1;
        break;
    case 2:
        return order-1;
        break;
    
    default:
        std::cout << "Shape1d::NShapeFunctions : Wrong side " << side << "\n";
        DebugStop();
        return -1;
        break;
    }
    return -1;
}

/// returns the total number of shape functions
int Shape1d::NShapeFunctions(VecInt &orders) {
    
    int nsf_tot = 0;
    for (int is=0; is<3; is++) {
        nsf_tot += NShapeFunctions(is, orders[is]);
    }
    
    return nsf_tot;
}
