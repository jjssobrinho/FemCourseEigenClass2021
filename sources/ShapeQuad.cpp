//
//  ShapeQuad.cpp
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "Shape1d.h"
#include "ShapeQuad.h"

// computes the shape functions in function of the coordinate in parameter
// space and orders of the shape functions (size of orders is number of 
// sides of the element topology)
void ShapeQuad::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, 
                      MatrixDouble &dphi){
    
    for (int i = 0; i < orders.size(); i++) {
        if (orders[i] < 0) {
            std::cout << "ShapeQuad::Shape: Invalid dimension for arguments: order\n";
            DebugStop();
        }
    }
    if (orders[0] > 1 || orders[1] > 1 || orders[2] > 1 || orders[3] > 1) {
        std::cout << "ShapeQuad::Shape: Invalid dimension for arguments: order\n";
        DebugStop();
    }

    auto nf = NShapeFunctions(orders);

    if (orders[nf-1] > 2) {
        std::cout << "ShapeQuad::Shape, only implemented until order = 2" << std::endl;
        DebugStop();
    }

    phi.resize(nf);
    dphi.resize(2,nf);

    double csi = xi[0];
    double eta = xi[1];

    phi[0] = 0.25*(1.-csi)*(1.-eta);
    phi[1] = 0.25*(1.+csi)*(1.-eta);
    phi[2] = 0.25*(1.+csi)*(1.+eta);
    phi[3] = 0.25*(1.-csi)*(1.+eta);

    dphi(0,0) = -0.25*(1.-eta);
    dphi(1,0) = -0.25*(1.-csi);
    dphi(0,1) = 0.25*(1.-eta);
    dphi(1,1) = -0.25*(1.+csi);
    dphi(0,2) = 0.25*(1.+eta);
    dphi(1,2) = 0.25*(1.+csi);
    dphi(0,3) = -0.25*(1.+eta);
    dphi(1,3) = 0.25*(1.-csi);

    int count = 4;
    if ( nf > nCorners) {
        if (orders[4] == 2) {
            phi[4] = 4.*phi[0]*(phi[1] + phi[2]);
            dphi(0,4) = 4.*(phi[1] + phi[2])*dphi(0,0) + 4.*phi[0]*
                        (dphi(0,1) + dphi(0,2));
            dphi(1,4) = 4.*(phi[1] + phi[2])*dphi(1,0) + 4.*phi[0]*
                        (dphi(1,1) + dphi(1,2));
            count++;
        } else if(orders[4] != 1) DebugStop();

        if (orders[5] == 2) {
            phi[5] = 4.*phi[1]*(phi[2] + phi[3]);
            dphi(0,5) = 4.*(phi[1] + phi[2])*dphi(0,1) + 4.*phi[1]*
                        (dphi(0,2) + dphi(0,3));
            dphi(1,5) = 4.*(phi[1] + phi[2])*dphi(1,1) + 4.*phi[1]*
                        (dphi(1,2) + dphi(1,3));
            count++;
        } else if(orders[5] != 1) DebugStop();

        if (orders[6] == 2) {
            phi[6] = 4.*phi[2]*(phi[3] + phi[0]);
            dphi(0,6) = 4.*(phi[3] + phi[0])*dphi(0,2) + 4.*phi[2]*
                        (dphi(0,3) + dphi(0,0));
            dphi(1,6) = 4.*(phi[3] + phi[0])*dphi(1,2) + 4.*phi[2]*
                        (dphi(1,3) + dphi(1,0));
            count++;
        } else if(orders[6] != 1) DebugStop();

        if (orders[7] == 2) {
            phi[7] = 4.*phi[3]*(phi[0] + phi[1]);
            dphi(0,7) = 4.*(phi[0] + phi[1])*dphi(0,3) + 4.*phi[3]*
                        (dphi(0,0) + dphi(0,1));
            dphi(1,7) = 4.*(phi[0] + phi[1])*dphi(1,3) + 4.*phi[3]*
                        (dphi(1,0) + dphi(1,1));
            count++;
        } else if(orders[7] != 1) DebugStop();

        if (orders[8] >= 2) {
            phi[8] = 16.*phi[0]*phi[2];
            dphi(0,8) = 16.*dphi(0,0)*phi[2] + 16.*dphi(0,2)*phi[0];
            dphi(1,8) = 16.*dphi(1,0)*phi[2] + 16.*dphi(1,2)*phi[0];
            count++;
        } else if(orders[8] != 1) DebugStop();
    }
    
    if (count != nf) DebugStop();
}

/// returns the number of shape functions associated with a side
int ShapeQuad::NShapeFunctions(int side, int order){
    if(order < 1 || order >2) DebugStop();
    if(side<4)
        return 1;//0 a 4
    else if(side<8)
        return (order-1);//6 a 14
    else if(side==8)
        return ((order-1)*(order-1));
    
    std::cout << "ShapeQuad::NShapeFunctions, bad parameter side " << side << std::endl;
    DebugStop();
    
    return 0;
}

/// returns the total number of shape functions
int ShapeQuad::NShapeFunctions(VecInt &orders){
    
    int res=4;
    for(int in=4; in<orders.size(); in++) {
        res += NShapeFunctions(in, orders[in]);
    }
    
    return res;
}
