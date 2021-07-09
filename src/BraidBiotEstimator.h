//
// Created by parnet on 09.06.21.
//

#ifndef UG_PLUGIN_XBRAIDBIOT_BRAIDBIOTESTIMATOR_H
#define UG_PLUGIN_XBRAIDBIOT_BRAIDBIOTESTIMATOR_H

#include <math.h>
#include "lib_disc/function_spaces/integrate.h"
#include "../../XBraidForUG4/src/interface/spatial_norm.h"


template<typename TDomain, typename TAlgebra>
class BiotBraidSpatialNorm : public BraidSpatialNorm<TDomain,TAlgebra> {
public:
    typedef ug::GridFunction<TDomain, TAlgebra> TGridFunction;
    typedef SmartPtr<TGridFunction> SPGridFunction;

    int m_uorder = 2;
    int m_porder = 4;
    double u_factor = 1;
    double p_factor = 1;
    BiotBraidSpatialNorm() : BraidSpatialNorm<TDomain,TAlgebra>(){}

    ~BiotBraidSpatialNorm(){}

    void set_order(int uorder, int porder){
        this->m_uorder = uorder;
        this->m_porder = porder;
    }

    void set_parameter(double alpha, double lambda,double mu){
        u_factor = lambda + 2*mu;
        p_factor = alpha;
    }
    /*void set_problem(SPBiotProblem problem){

    }*/

    double norm(SPGridFunction u) override {
        double unorm_x =  ug::H1SemiNorm(*u.get(), "ux", this->m_uorder);
        double unorm_y =  ug::H1SemiNorm(*u.get(), "uy", this->m_uorder);
        double unorm_p =  ug::L2Norm(*u.get(), "p", this->m_porder);

        double total_norm = sqrt(u_factor*(unorm_x*unorm_x  + unorm_y*unorm_y) + p_factor*unorm_p*unorm_p);
        return total_norm;
    }
};


#endif //UG_PLUGIN_XBRAIDBIOT_BRAIDBIOTESTIMATOR_H
