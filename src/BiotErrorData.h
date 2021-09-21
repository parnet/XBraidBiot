//
// Created by parnet on 2021-09-21.
//

#ifndef UG_PLUGIN_XBRAIDBIOT_BIOTERRORDATA_H
#define UG_PLUGIN_XBRAIDBIOT_BIOTERRORDATA_H

#include <math.h>
#include "lib_disc/function_spaces/integrate.h"
#include "lib_disc/function_spaces/grid_function.h"



template<typename TDomain, typename TAlgebra>
class BiotErrorData {
public:
    typedef ug::GridFunction<TDomain, TAlgebra> TGridFunction;
    typedef SmartPtr<TGridFunction> SPGridFunction;


    double l2_norm_p = 0.0;
    double l2_norm_ux = 0.0;
    double l2_norm_uy = 0.0;

    double h1_norm_ux = 0.0;
    double h1_norm_uy = 0.0;


    int porder = 2;
    int uorder = 4;

    void compute(SPGridFunction u) {

        this->l2_norm_p = ug::L2Norm(*u.get(), "p", porder);;
        this->l2_norm_ux = ug::L2Norm(*u.get(), "ux", uorder);
        this->l2_norm_uy = ug::L2Norm(*u.get(), "uy", uorder);

        this->h1_norm_ux = ug::H1SemiNorm(*u.get(), "ux", uorder);
        this->h1_norm_uy = ug::H1SemiNorm(*u.get(), "uy", uorder);

    }
};



#endif //UG_PLUGIN_XBRAIDBIOT_BIOTERRORDATA_H
