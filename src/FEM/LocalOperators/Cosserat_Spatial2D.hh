//
//  CosseratSpatial.hh
//
//
//  Created by Tim Dodwell on 10/11/2015.
//
//

#ifndef CosseratSpatial_h
#define CosseratSpatial_h

#include <dune/pdelab/localoperator/idefault.hh>

template<typename PARAM, int dispdofel>
class CosseratSpatial2D :
    public Dune::PDELab::Cosserat2D<PARAM,dispdofel>, // Inherit Cosserat Linear Operator from Stationary Problem
    public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double> // default methods
{
    PARAM& param;

public:
    // Constructor for Cosserat Spartial - Wrap onto stationary Cosserat operator
    CosseratSpatial2D(PARAM& param_, unsigned int intorder_=1)
    :Dune::PDELab::Cosserat2D<PARAM,dispdofel>(param_,intorder_),param(param_) {

    }

    //
    void preStep(double time, double dt, int stages){
        // CAN ADD LINE TO CHANGE BOUNDARY CONDITIONS

        param.setTime(time + dt);

        Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>::preStep(time,dt,stages);

    }

};




#endif /* CosseratSpatial_h */
