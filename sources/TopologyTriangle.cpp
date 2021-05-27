//
//  TopologyTriangle.cpp
//  FemSC
//
//  Created by Karolinne Oliveira Coelho on 4/28/18.
//
//

#include <stdio.h>
#include "TopologyTriangle.h"
#include "tpanic.h"

// Number of sides associated with triangle elements elements
const int TopologyTriangle::nSides;

// Number of corner nodes associated with triangle elements
const int TopologyTriangle::nCorners;

// Dimension of triangle elements
const int TopologyTriangle::Dimension;


int TopologyTriangle::NSideNodes(int side)
{
    if (side>6) {
        std::cout << "TopologyTriangle::NSideNodes: Bad parameter side" << std::endl;
        DebugStop();
        return EXIT_FAILURE;
    }
    
    int nsidenodes[7] = {1,1,1,2,2,2,3};
    return nsidenodes[side];
}

// local node index of a node associated with a side
int TopologyTriangle::SideNodeLocIndex(int side, int node)
{
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    std::cout << "side= " << side << std::endl;
    std::cout << "node= " << node << std::endl;
    if(side<3 && node==0){
        std::cout << "returns= " << side << "\n" << std::endl;
        return side;
    }
    if(side>=3 && side<6 && node<2) {
        std::cout << "returns= " << (side+node)%3 << "\n" << std::endl;
        return (side+node)%3;
    }
    if(side==6 && node <3) {
        std::cout << "returns= " << node << "\n" <<std::endl;
        return node;
    }
    
    std::cout << "TopologyTriangle::SideNodeIndex inconsistent side or node" << std::endl;
    DebugStop();
    return EXIT_FAILURE;
}

// return the enumerated element type
MElementType TopologyTriangle::Type(){
    return ETriangle;
}
