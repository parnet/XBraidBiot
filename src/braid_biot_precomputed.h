//
// Created by parnet on 2021-09-21.
//

#ifndef UG_PLUGIN_XBRAIDBIOT_BRAIDBIIOTPRECOMPUTED_H
#define UG_PLUGIN_XBRAIDBIOT_BRAIDBIIOTPRECOMPUTED_H
//
// Created by parnet on 14.06.21.
//

#include "../../XBraidForUG4/src/interface/scriptor.h"
#include "../../XBraidForUG4/src/util/scriptor.h"
#include "../../XBraidForUG4/src/util/paralog.h"

#include "../../XBraidUtil/src/io_gridfunction.h"
#include "common/math/math_vector_matrix/math_vector_functions.h"
#include "lib_algebra/vector_interface/vec_functions.h"
#include "braid_biot_estimator.h"
#include "biot_error_data.h"

#include "../../XBraidPoroelasticity/src/biot_tools.h"
#include "../../XBraidPoroelasticity/src/barry_mercer.h"

using namespace ug::XBraidForUG4;

namespace ug {

    namespace XBraidBiot {
        template<typename TDomain, typename TAlgebra>
        class BraidBiotCheckPrecomputed : public Scriptor<TDomain, TAlgebra> {
        public:

            typedef ug::GridFunction<TDomain, TAlgebra> TGridFunction;
            typedef SmartPtr<TGridFunction> SPGridFunction;

            typedef ug::VTKOutput<TDomain::dim> TVTKOutput;
            typedef SmartPtr<TVTKOutput> SPVTKOutput;

            typedef VTKScriptor<TDomain, TAlgebra> TVTKScriptor;
            typedef SmartPtr<TVTKScriptor> SPVTKScriptor;

            typedef IOScriptor<TDomain, TAlgebra> TIOScriptor;
            typedef SmartPtr<TIOScriptor> SPIOScriptor;

            typedef Paralog TParalog;
            typedef SmartPtr<TParalog> SPParalog;

            std::vector<int> index_level;

            int max_index = 512;
            int max_index_precomputed = 512;

            bool write_solution = false;
            bool write_error = false;

            bool io_write_solution = false;
            bool io_write_error = false;

            BiotErrorData<TDomain, TAlgebra> err_u;
            BiotErrorData<TDomain, TAlgebra> err_sol;
            BiotErrorData<TDomain, TAlgebra> err_udiffsol;

            SPVTKScriptor m_out_solution;
            SPVTKScriptor m_out_diff;
            SPIOScriptor m_ioout_solution;
            SPIOScriptor m_ioout_diff;
            SPParalog m_log;

            std::string base_path = "/home/mparnet/analyticsolution";


            typedef std::tuple<int, int, int> TKey;
            std::map<TKey, int> map;

            int num_ref = 3;

            BraidBiotCheckPrecomputed() : Scriptor<TDomain, TAlgebra>() {
                index_level = std::vector<int>();
                err_u = BiotErrorData<TDomain, TAlgebra>();
                err_sol = BiotErrorData<TDomain, TAlgebra>();
                err_udiffsol = BiotErrorData<TDomain, TAlgebra>();
            };

            ~BraidBiotCheckPrecomputed() = default;

            void set_log(SPParalog log) {
                this->m_log = log;
            }

            void set_base_path(std::string path) {
                this->base_path = path;
            }

            void set_solution_name(SPVTKOutput vtk, const char *fname) {
                this->m_out_solution = make_sp(new TVTKScriptor(vtk, fname));
                this->m_ioout_solution = make_sp(new TIOScriptor(fname));

            }

            void set_diff_name(SPVTKOutput vtk, const char *fname) {
                this->m_out_diff = make_sp(new TVTKScriptor(vtk, fname));
                this->m_ioout_diff = make_sp(new TIOScriptor(fname));
            }

            void compare_norms(int index, double time, int iteration, int level, int c, bool done) {
                m_log->o << ">> norms idx=" << index
                         << " t=" << time
                         << " iter=" << iteration
                         << " level=" << level
                         << " c=" << c
                         << " done=" << done
                         << std::endl;

                m_log->o << std::setw(10) << ">> norm"
                         << std::setw(20) << "solution"
                         << std::setw(20) << "error"
                         << std::setw(20) << "relative"
                         << std::endl;

                m_log->o << std::setw(10) << ">> l2(p)"
                         << std::setw(20) << err_u.l2_norm_p
                         << std::setw(20) << err_udiffsol.l2_norm_p
                         << std::setw(20) << (err_udiffsol.l2_norm_p / err_sol.l2_norm_p)
                         << std::endl;

                m_log->o << std::setw(10) << ">> l2(ux)"
                         << std::setw(20) << err_u.l2_norm_ux
                         << std::setw(20) << err_udiffsol.l2_norm_ux
                         << std::setw(20) << (err_udiffsol.l2_norm_ux / err_sol.l2_norm_ux)
                         << std::endl;

                m_log->o << std::setw(10) << ">> l2(uy)"
                         << std::setw(20) << err_u.l2_norm_uy
                         << std::setw(20) << err_udiffsol.l2_norm_uy
                         << std::setw(20) << (err_udiffsol.l2_norm_uy / err_sol.l2_norm_uy)
                         << std::endl;


                m_log->o << std::setw(10) << ">> h1(ux)"
                         << std::setw(20) << err_u.h1_norm_ux
                         << std::setw(20) << err_udiffsol.h1_norm_ux
                         << std::setw(20) << (err_udiffsol.h1_norm_ux / err_sol.h1_norm_ux)
                         << std::endl;

                m_log->o << std::setw(10) << ">> h1(uy)"
                         << std::setw(20) << err_u.h1_norm_uy
                         << std::setw(20) << err_udiffsol.h1_norm_uy
                         << std::setw(20) << (err_udiffsol.h1_norm_uy / err_sol.h1_norm_uy)
                         << std::endl << std::endl;
            }

            void set_num_ref(int ref) {
                this->num_ref = ref;
            }

            void set_max_index(int precomputed, int problem) {
                this->max_index_precomputed = precomputed;
                this->max_index = problem;
                index_level.resize(1);
                index_level[0] = this->max_index;
            }

            void set_c_factor(int level, int factor) {
                index_level.resize(level + 2);
                index_level[level + 1] = index_level[level] / factor;
                this->m_log->o << "level: " << index_level[level] << "\t" << factor << "\t" << index_level[level + 1]
                               << std::endl;
            }

            void set_vtk_write_mode(bool solution, bool error) {
                this->write_solution = solution;
                this->write_error = error;
            }

            void set_io_write_mode(bool solution, bool error) {
                this->io_write_solution = solution;
                this->io_write_error = error;
            }

            bool lua_write(SPGridFunction u, int index, double time) {
                return this->write(u, index, time);
            }

            bool write(SPGridFunction u, int index, double time) override {
                int zidx = (index * this->max_index_precomputed) / index_level[0];
                int rem = (index * this->max_index_precomputed) % index_level[0];
                if (rem == 0 && index != 0) {
                    SPGridFunction sol = u->clone_without_values();
                    SPGridFunction udiffsol = u->clone();

                    // write vtk output
                    if (this->write_solution) {
                        m_out_solution->write(u, index, time);
                    }

                    if (this->io_write_solution) {
                        m_ioout_solution->write(u, index, time);
                    }


                    // load gridfunction file (ref solution)
                    ug::XBraidUtil::IOGridFunction <TDomain, TAlgebra> io = ug::XBraidUtil::IOGridFunction<TDomain, TAlgebra>();
                    std::stringstream ss_ref;
                    ss_ref << this->base_path << "/num_ref_" << this->num_ref << "/BarryMercer2D_" << zidx
                           << ".gridfunction";
                    io.read(sol, ss_ref.str().c_str());

                    // substract
                    ::VecAdd(1.0, *udiffsol.get(), -1.0, *sol.get());

                    // write vtk error
                    if (this->write_error) {
                        m_out_diff->write(udiffsol, index, time);
                    }

                    if (this->io_write_error) {
                        m_ioout_diff->write(udiffsol, index, time);
                    }

                    // compute norms
                    err_u.compute(u->clone());
                    err_sol.compute(sol->clone());
                    err_udiffsol.compute(udiffsol->clone());

                    // write norms
                    compare_norms(index, time, 0, 0, 0, true);
                }
                return false; // no error
            };

            bool write(SPGridFunction u, int index, double time, int iteration, int level) override {

                int count = 0;
                auto tuple = std::make_tuple(index, iteration, level);
                auto it = map.find(tuple);
                if (it != map.end()) {
                    count = it->second;
                    count += 1;
                    map[tuple] = count;
                } else {
                    count = 0;
                    map.emplace(tuple, 0);
                }


                int zidx = (index * this->max_index_precomputed) / index_level[level];
                int rem = (index * this->max_index_precomputed) % index_level[level];

                if (rem == 0 && index != 0) {
                    SPGridFunction sol = u->clone_without_values();
                    SPGridFunction udiffsol = u->clone();

                    // write vtk output
                    if (this->write_solution) {
                        m_out_solution->write(u, index, time, iteration, level);
                    }
                    if (this->io_write_solution) {
                        m_ioout_solution->write(u, index, time, iteration, level);
                    }

                    // load gridfunction file (ref solution)
                    ug::XBraidUtil::IOGridFunction <TDomain, TAlgebra> io = ug::XBraidUtil::IOGridFunction<TDomain, TAlgebra>();
                    std::stringstream ss_ref;
                    ss_ref << this->base_path << "/num_ref_" << this->num_ref << "/BarryMercer2D_" << zidx
                           << ".gridfunction";
                    std::cout << index << "\t";
                    std::cout << ss_ref.str().c_str() << std::endl;
                    io.read(sol, ss_ref.str().c_str());
                    // substract
                    ::VecAdd(1.0, *udiffsol.get(), -1.0, *sol.get());

                    // write vtk error
                    if (this->write_solution) {
                        m_out_diff->write(udiffsol, index, time, iteration, level);
                    }
                    if (this->io_write_solution) {
                        m_ioout_diff->write(udiffsol, index, time, iteration, level);
                    }

                    // compute norms
                    err_u.compute(u->clone());
                    err_sol.compute(sol->clone());
                    err_udiffsol.compute(udiffsol->clone());

                    // write norms
                    compare_norms(index, time, iteration, level, count, false);
                }

                return false; // no error
            };
        };
    }}
#endif //UG_PLUGIN_XBRAIDBIOT_BRAIDBIOTPRECOMPUTED_H
