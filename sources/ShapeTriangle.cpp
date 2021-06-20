//
//  ShapeTriangle.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "ShapeTriangle.h"
#include "Shape1d.h"

/// computes the shape functions in function of the coordinate in parameter space and orders of the shape functions (size of orders is number of sides of the element topology)
void ShapeTriangle::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, MatrixDouble &dphi){
    
    for (int i = 0; i < orders.size(); i++)
    {
        if (orders[i] < 0) {
            std::cout << "ShapeTriangle::Shape: Invalid dimension for arguments: order\n";
            DebugStop();
        }
    }
    if (orders[0] > 1 || orders[1] > 1 || orders[2] > 1) {
        std::cout << "ShapeTriangle::Shape: Invalid dimension for arguments: order\n";
        DebugStop();
    }

    int nshape = NShapeFunctions(orders);
    phi.resize(nshape);
    dphi.resize(2,nshape);

    if (orders[nshape-1] > 2) {
        std::cout << "ShapeTriangle::Shape, only implemented until order = 2" << std::endl;
        DebugStop();
    }
    if (xi.size() < 2) DebugStop();
    
    // Linear order
    phi[0] =  1.-xi[0]-xi[1];
    phi[1] =  xi[0];
    phi[2] =  xi[1];

    dphi(0,0) = -1.;
    dphi(1,0) = -1.;
    dphi(0,1) =  1.;
    dphi(1,1) =  0.;
    dphi(0,2) =  0.;
    dphi(1,2) =  1.;

    int count = 3;
    if (nshape > nCorners){
        int is;
        for (is = 3; is < 6; is++){ // sides loop
            if (orders[is] == 2){
                int is1 = SideNodeLocIndex(is,0);
                int is2 = SideNodeLocIndex(is,1);
                phi[is] = 4.*phi[is1]*phi[is2];
                dphi(0,is) = 4.*(dphi(0, is1) * phi[is2] + phi[is1] * dphi(0,is2));
                dphi(1,is) = 4.*(dphi(1, is1) * phi[is2] + phi[is1] * dphi(1,is2));
                count++;
            }
        }
    }
    //if (orders[6] == 3) {
    //    phi[6] = 27.*phi[0]*phi[1]*phi[2];
    //    dphi(0,6) = 27.*dphi(0,0)*phi[1]*phi[2] + 27.*dphi(0,1)*phi[0]*phi[2] +
    //                27.*dphi(0,2)*phi[0]*phi[1];
    //    dphi(1,6) = 27.*dphi(1,0)*phi[1]*phi[2] + 27.*dphi(1,1)*phi[0]*phi[2] +
    //                27.*dphi(1,2)*phi[0]*phi[1];
    //    count++;
    //}
    if (count != nshape) {
        std::cout << "Error: count /= nshape " << std::endl;
        DebugStop();
    }
   for(int is = 3 ; is< nSides; is++) if(orders[is] != 1 && orders[is] != 2)
   DebugStop();
}

/// returns the number of shape functions associated with a side
int ShapeTriangle::NShapeFunctions(int side, int order){
    switch(side) {
        case 0:
        case 1:
        case 2:
            return 1;
        case 3:
        case 4:
        case 5:
            return order-1;
        case 6:
            return 0;
    }
    
    DebugStop();
    std::cout << "ShapeTriangle::NShapeFunctions, bad parameter side " << std::endl;
    return 0;
}

/// returns the total number of shape functions
int ShapeTriangle::NShapeFunctions(VecInt &orders){
    
    int res=3;
    for(int in=3; in<orders.size(); in++) {
        res += NShapeFunctions(in, orders[in]);
    }
    
    return res;
    
}
