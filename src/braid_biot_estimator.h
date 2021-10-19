//
// Created by parnet on 09.06.21.
//

#ifndef UG_PLUGIN_XBRAIDBIOT_BRAID_BIOT_ESTIMATOR_H
#define UG_PLUGIN_XBRAIDBIOT_BRAID_BIOT_ESTIMATOR_H

#include <math.h>
#include "lib_disc/function_spaces/integrate.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "../../XBraidForUG4/src/interface/spatial_norm.h"

namespace ug {

    namespace XBraidBiot {
        template<typename TDomain, typename TAlgebra>
        class BiotBraidSpatialNorm : public ug::XBraidForUG4::BraidSpatialNorm<TDomain, TAlgebra> {
        public:
            typedef ug::GridFunction<TDomain, TAlgebra> TGridFunction;
            typedef SmartPtr<TGridFunction> SPGridFunction;

            int m_uorder = 4;
            int m_porder = 2;
            double u_factor = 1;
            double p_factor = 1;

            BiotBraidSpatialNorm() : ug::XBraidForUG4::BraidSpatialNorm<TDomain, TAlgebra>() {}

            ~BiotBraidSpatialNorm() {}

            void set_order(int uorder, int porder) {
                this->m_uorder = uorder;
                this->m_porder = porder;
            }

            void set_parameter(double alpha, double lambda, double mu) {
                u_factor = lambda + 2 * mu;
                p_factor = alpha;
            }

            /*void set_problem(SPBiotProblem problem){

            }*/

            double norm(SPGridFunction u) override {
                double norm_x = ug::H1SemiNorm(*u.get(), "ux", this->m_uorder);
                double norm_y = ug::H1SemiNorm(*u.get(), "uy", this->m_uorder);
                double norm_p = ug::L2Norm(*u.get(), "p", this->m_porder);


                double pnorm = p_factor * norm_p * norm_p; //p_factor*unorm_p*unorm_p;
                double unorm = u_factor * (norm_x * norm_x + norm_y * norm_y);

                double total_norm = sqrt(unorm + pnorm);

                return total_norm;
            }
        };
    }}

#endif //UG_PLUGIN_XBRAIDBIOT_BRAID_BIOT_ESTIMATOR_H
