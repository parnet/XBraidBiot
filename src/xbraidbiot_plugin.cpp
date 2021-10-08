//
// Created by parnet on 09.06.21.
//

//
// Created by parnet on 09.06.21.
//

//  zzzzz
// Created by parnet on 12.05.19.
//


#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"

#include "braid_biot_estimator.h"
#include "braid_biot_control.h"

#include "braid_biot_precomputed.h"


using namespace std;
using namespace ug::bridge;

namespace ug {

    namespace XBraidBiot {

        struct Functionality {


            template<typename TDomain, typename TAlgebra>
            static void DomainAlgebra(Registry &reg, string grp) {

                string suffix = GetDomainAlgebraSuffix<TDomain, TAlgebra>();
                string tag = GetDomainAlgebraTag<TDomain, TAlgebra>();

                // Braid Time Integrator
                {
                    typedef BiotBraidSpatialNorm<TDomain, TAlgebra> TBiotSpatialNorm;
                    typedef ug::XBraidForUG4::BraidSpatialNorm<TDomain, TAlgebra> TSpatialNorm;
                    string name_gf = string("BiotBraidSpatialNorm").append(suffix);
                    reg.add_class_<TBiotSpatialNorm, TSpatialNorm>(name_gf, grp)
                            .add_constructor()
                            .add_method("norm", &TBiotSpatialNorm::norm, "None", "verbose",
                                        "set the level of verbose (true / false)")
                            .add_method("set_order", &TBiotSpatialNorm::set_order, "None", "verbose",
                                        "set the level of verbose (true / false)")
                            .add_method("set_parameter", &TBiotSpatialNorm::set_parameter, "None", "verbose",
                                        "set the level of verbose (true / false)")

                                    //.add_method("set_adapt_convergence", &TBraidIntegrator::set_adapt_conv, "None", "initial time","set t0 as initial time")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name_gf, "BiotBraidSpatialNorm", tag);

                }

                // Braid Time Integrator
                {
                    typedef BraidBiotCheckPrecomputed<TDomain, TAlgebra> TBraidBiotCheckPrecomputed;
                    typedef ug::XBraidForUG4::Scriptor<TDomain, TAlgebra> TScriptor;
                    string name_gf = string("BraidBiotCheckPrecomputed").append(suffix);
                    reg.add_class_<TBraidBiotCheckPrecomputed, TScriptor>(name_gf, grp)
                            .add_constructor()
                            .add_method("compare_norms", &TBraidBiotCheckPrecomputed::compare_norms, "None", "verbose",
                                        "set the level of verbose (true / false)")
                            .add_method("set_num_ref", &TBraidBiotCheckPrecomputed::set_num_ref, "None", "verbose",
                                        "set the level of verbose (true / false)")
                            .add_method("set_max_index", &TBraidBiotCheckPrecomputed::set_max_index, "None", "verbose",
                                        "set the level of verbose (true / false)")
                            .add_method("set_base_path", &TBraidBiotCheckPrecomputed::set_base_path, "None", "verbose",
                                        "set the level of verbose (true / false)")
                            .add_method("set_c_factor", &TBraidBiotCheckPrecomputed::set_c_factor, "None", "verbose",
                                        "set the level of verbose (true / false)")
                            .add_method("set_solution_name", &TBraidBiotCheckPrecomputed::set_solution_name, "None",
                                        "verbose", "set the level of verbose (true / false)")
                            .add_method("set_diff_name", &TBraidBiotCheckPrecomputed::set_diff_name, "None", "verbose",
                                        "set the level of verbose (true / false)")
                            .add_method("lua_write", &TBraidBiotCheckPrecomputed::lua_write, "None", "verbose",
                                        "set the level of verbose (true / false)")
                            .add_method("set_log", &TBraidBiotCheckPrecomputed::set_log, "None", "verbose",
                                        "set the level of verbose (true / false)")
                            .add_method("set_vtk_write_mode", &TBraidBiotCheckPrecomputed::set_vtk_write_mode, "None",
                                        "verbose", "set the level of verbose (true / false)")
                            .add_method("set_io_write_mode", &TBraidBiotCheckPrecomputed::set_io_write_mode, "None",
                                        "verbose", "set the level of verbose (true / false)")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name_gf, "BraidBiotCheckPrecomputed", tag);
                }


            }


            template<typename TDomain>
            static void Domain(Registry &reg, string grp) {
                /*
                string suffix = GetDomainSuffix<TDomain>();
                string tag = GetDomainTag<TDomain>();
                 */

            }

            template<int dim>
            static void Dimension(Registry &reg, string grp) {
                /*
                string suffix = GetDimensionSuffix<dim>();
                string tag = GetDimensionTag<dim>();
                 */
            }

            template<typename TAlgebra>
            static void Algebra(Registry &reg, string grp) {
                string suffix = GetAlgebraSuffix<TAlgebra>();
                string tag = GetAlgebraTag<TAlgebra>();
            }

            static void Common(Registry &reg, string grp) {}

        };
        // end namespace XBraidIntegrator


        struct FunctionalityFor2D {
            template<typename TDomain, typename TAlgebra>
            static void DomainAlgebra(Registry &reg, string grp) {
                //	useful defines
                const string suffix = GetDomainAlgebraSuffix<TDomain, TAlgebra>();
                const string tag = GetDomainAlgebraTag<TDomain, TAlgebra>();

                // Braid Time Integrator
                {
                    typedef BraidBiotCheck<TDomain, TAlgebra> TBiotSpatialNorm;
                    typedef ug::XBraidForUG4::Scriptor<TDomain, TAlgebra> TSpatialNorm;
                    string name_gf = string("BraidBiotCheck").append(suffix);
                    reg.add_class_<TBiotSpatialNorm, TSpatialNorm>(name_gf, grp)
                            .add_constructor()
                            .add_method("set_problem", &TBiotSpatialNorm::set_problem, "None", "verbose",
                                        "set the level of verbose (true / false)")
                            .add_method("set_napprox", &TBiotSpatialNorm::set_napprox, "None", "verbose",
                                        "set the level of verbose (true / false)")

                                    //.add_method("set_adapt_convergence", &TBraidIntegrator::set_adapt_conv, "None", "initial time","set t0 as initial time")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name_gf, "BraidBiotCheck", tag);
                }


            }
        };

    }



    extern "C" void
    InitUGPlugin_XBraidBiot(Registry *reg, string grp) {
        using namespace XBraidBiot;
        grp.append("XBraidBiot");
        // Space Time Communicator


        try {
            RegisterCommon<Functionality>(*reg, grp);
            RegisterDimensionDependent<Functionality>(*reg, grp);
            RegisterDomainDependent<Functionality>(*reg, grp);
            RegisterAlgebraDependent<Functionality>(*reg, grp);
            RegisterDomainAlgebraDependent<Functionality>(*reg, grp);
            RegisterDomain2dAlgebraDependent<FunctionalityFor2D>(*reg,grp);
        }
        UG_REGISTRY_CATCH_THROW(grp);

    }
}