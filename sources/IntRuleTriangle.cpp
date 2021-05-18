/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <iostream> 
#include "IntRuleTriangle.h"
#include "tpanic.h"

#define Perm3(a)	a, a, a
#define Dup3(w)		(0.5L*w)
#define Perm21(a)	a, a, 1.L - 2.0L*a, a, 1.L- 2.L*a, a, 1.L- 2.L*a, a, a
#define Dup21(w)	Dup3(w), Dup3(w), Dup3(w)

IntRuleTriangle::IntRuleTriangle(){

}

IntRuleTriangle::IntRuleTriangle(int order) {
    SetOrder(order);
}

void IntRuleTriangle::SetOrder(int order) {
    fOrder = order;
    if (order < 0 || order > MaxOrder()) DebugStop();

    switch (order)
    {
    case 0:
    case 1:
        fPoints.resize(1,Dimension());
        fWeights.resize(1);
        fPoints(0,0) = 1./3.;
        fPoints(0,1) = 1./3.;
        fWeights(0) = 1./2.;
        break;
    case 2:
        fPoints.resize(3,Dimension());
        fWeights.resize(3);
        fPoints(0,0) = 1./6.;
        fPoints(0,1) = 1./6.;
        fPoints(1,0) = 2./3.;
        fPoints(1,1) = 1./6.;
        fPoints(2,0) = 1./6.;
        fPoints(2,1) = 2./3.;
        fWeights(0) = 1./6.;
        fWeights(1) = 1./6.;
        fWeights(2) = 1./6.;
        break;

    case 3:
        fPoints.resize(4,Dimension());
        fWeights.resize(4);
        fPoints(0,0) = 1./3.;
        fPoints(0,1) = 1./3.;
        fPoints(1,0) = 1./5.;
        fPoints(1,1) = 1./5.;
        fPoints(2,0) = 3./5.;
        fPoints(2,1) = 1./5.;
        fPoints(3,0) = 1./5.;
        fPoints(3,1) = 3./5.;
        fWeights(0) = -27./96.;
        fWeights(1) = 25./96.;
        fWeights(2) = 25./96.;
        fWeights(3) = 25./96.;
        break;

    case 4:
    case 5:
        fPoints.resize(7,Dimension());
        fWeights.resize(7);
        fPoints(0,0) = (6.-sqrt(15.))/21.;
        fPoints(0,1) = (6.-sqrt(15.))/21.;
        fPoints(1,0) = (6.-sqrt(15.))/21.;
        fPoints(1,1) = 1.-2.*(6.-sqrt(15.))/21.;
        fPoints(2,0) = 1.-2.*(6.-sqrt(15.))/21.;
        fPoints(2,1) = (6.-sqrt(15.))/21.;
        fWeights(0)  = 0.5*((155. - sqrt(15.)) / 1200.);
        fWeights(1)  = 0.5*((155. - sqrt(15.)) / 1200.);
        fWeights(2)  = 0.5*((155. - sqrt(15.)) / 1200.);
        
        fPoints(3,0) = (6.+sqrt(15.))/21.;
        fPoints(3,1) = (6.+sqrt(15.))/21.;
        fPoints(4,0) = (6.+sqrt(15.))/21.;
        fPoints(4,1) = 1.-2.*(6.+sqrt(15.))/21.;
        fPoints(5,0) = 1.-2.*(6.+sqrt(15.))/21.;
        fPoints(5,1) = (6.+sqrt(15.))/21.;
        fWeights(3)  = 0.5*((155. + sqrt(15.)) / 1200.);
        fWeights(4)  = 0.5*((155. + sqrt(15.)) / 1200.);
        fWeights(5)  = 0.5*((155. + sqrt(15.)) / 1200.);

        fPoints(6,0) = 1./3.;
        fPoints(6,1) = 1./3.;
        fWeights(6) = 0.5*(9./40.);
        break;
    
    default:
        DebugStop();
        break;
    }
}
