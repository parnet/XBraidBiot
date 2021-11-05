//
// Created by parnet on 09.06.21.
//

#ifndef UG_PLUGIN_XBRAIDBIOT_BRAID_BIOT_ESTIMATOR_H
#define UG_PLUGIN_XBRAIDBIOT_BRAID_BIOT_ESTIMATOR_H

#include <math.h>
#include "lib_disc/function_spaces/integrate.h"
#include "lib_disc/function_spaces/grid_function.h"

#include "../../XBraidForUG4/src/interface/spatial_norm.h"
#include "../../XBraidForUG4/src/util/paralog.h"

#include "biot_error_data.h"


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


        template<typename TDomain, typename TAlgebra>
        class BiotBraidDisplacementNorm : public ug::XBraidForUG4::BraidSpatialNorm<TDomain, TAlgebra> {
        public:
            typedef ug::GridFunction<TDomain, TAlgebra> TGridFunction;
            typedef SmartPtr<TGridFunction> SPGridFunction;

            typedef ug::XBraidForUG4::Paralog TParalog;
            typedef SmartPtr<TParalog> SPParalog;


            BiotBraidDisplacementNorm() : ug::XBraidForUG4::BraidSpatialNorm<TDomain, TAlgebra>() {}

            ~BiotBraidDisplacementNorm() {}

            int count = 0;

            SPParalog m_log;

            void set_log(SPParalog log) {
                this->m_log = log;
            }

            double norm(SPGridFunction u) override {
                BiotErrorData<TDomain,TAlgebra> errdata = BiotErrorData<TDomain,TAlgebra>();
                errdata.compute(u);

                m_log->o << ">R> rnorm idx=" << count << std::endl;
                m_log->o << std::setw(10) << ">R>  l2(p)"
                         << std::setw(20) << errdata.l2_norm_p
                         << std::endl;

                m_log->o << std::setw(10) << ">R> l2(ux)"
                         << std::setw(20) << errdata.l2_norm_ux
                         << std::endl;

                m_log->o << std::setw(10) << ">R> l2(uy)"
                         << std::setw(20) << errdata.l2_norm_uy
                         << std::endl;

                m_log->o << std::setw(10) << ">R> h1(ux)"
                         << std::setw(20) << errdata.h1_norm_ux
                         << std::endl;

                m_log->o << std::setw(10) << ">R> h1(uy)"
                         << std::setw(20) << errdata.h1_norm_uy
                         << std::endl;
                count ++ ;
                return errdata.l2_norm_ux;
            }
        };

    }}

#endif //UG_PLUGIN_XBRAIDBIOT_BRAID_BIOT_ESTIMATOR_H
