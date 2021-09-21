//
// Created by parnet on 2021-09-21.
//

#ifndef UG_PLUGIN_XBRAIDBIOT_BRAIDBIIOTPRECOMPUTED_H
#define UG_PLUGIN_XBRAIDBIOT_BRAIDBIIOTPRECOMPUTED_H
//
// Created by parnet on 14.06.21.
//

#include "../../XBraidForUG4/src/interface/scriptor.h"
#include "../../XBraidUtil/src/IOGridFunction.h"

#include "BraidBiotEstimator.h"
#include "BiotErrorData.h"

template<typename TDomain, typename TAlgebra>
class BraidBiotCheckPrecomputed : public Scriptor<TDomain, TAlgebra> {
public:

    typedef ug::GridFunction <TDomain, TAlgebra> TGridFunction;
    typedef SmartPtr <TGridFunction> SPGridFunction;

    typedef ug::VTKOutput <TDomain::dim> TVTKOutput;
    typedef SmartPtr <TVTKOutput> SPVTKOutput;

    typedef SmartPtr <TProblem> SPProblem;

    typedef Paralog TParalog;
    typedef SmartPtr<TParalog> SPParalog;

    SPProblem m_problem;
    std::vector<int> index_level;

    int max_index = 512;
    int max_index_precomputed = 512;

    BiotErrorData err_u;
    BiotErrorData err_sol;
    BiotErrorData err_udiffsol;

    const char *m_sol_filename;
    SPVTKOutput m_out_solution;
    const char *m_diff_filename;
    SPVTKOutput m_out_diff;
    SPParalog m_log;

    std::string path_ref_3 = "/home/mparnet/analyticsolution/num_ref_3/BarryMercer2D_";
    std::string path_ref_4 = "/home/mparnet/analyticsolution/num_ref_4/BarryMercer2D_";

    int num_ref = 3;

    BraidBiotCheck() : Scriptor<TDomain, TAlgebra>() {
        index_level = std::vector<int>()
        err_u = BiotErrorData();
        err_sol = BiotErrorData();
        err_udiffsol = BiotErrorData();
    };

    ~BraidBiotCheck() = default;


    void set_vtk_solution(SPVTKOutput vtk, const char * fname){
        this->m_sol_filename = fname;
        this->m_out_solution = vtk;
    }

    void set_vtk_diff(SPVTKOutput vtk, const char * fname){
        this->m_out_diff = fname;
        this->m_diff_filename = vtk;
    }

    void compare_norms(int index, double time, int iteration, int level) {
        m_log->o << "norms idx=" <<index
                 <<" t=" << time
                 << " iter=" <<iteration
                 <<" level=" <<level
                 << std::endl

        m_log->o << std::setw(10) << "norm"
                  << std::setw(20) << "solution"
                  << std::setw(20) << "error"
                  << std::setw(20) << "relative"
                  << std::endl;

        m_log->o << std::setw(10) << "l2(p)"
                  << std::setw(20) << err_u.l2_norm_p
                  << std::setw(20) << err_udiffsol.l2_norm_p
                  << std::setw(20) << (err_udiffsol.l2_norm_p / err_sol.l2_norm_p)
                  << std::endl;

        m_log->o << std::setw(10) << "l2(ux)"
                  << std::setw(20) << err_u.l2_norm_ux
                  << std::setw(20) << err_udiffsol.l2_norm_ux
                  << std::setw(20) << (err_udiffsol.l2_norm_ux / err_sol.l2_norm_ux)
                  << std::endl;

        m_log->o << std::setw(10) << "l2(uy)"
                  << std::setw(20) << err_u.l2_norm_uy
                  << std::setw(20) << err_udiffsol.l2_norm_uy
                  << std::setw(20) << (err_udiffsol.l2_norm_uy / err_sol.l2_norm_uy)
                  << std::endl;


        m_log->o << std::setw(10) << "h1(ux)"
                  << std::setw(20) << err_u.h1_norm_ux
                  << std::setw(20) << err_udiffsol.h1_norm_ux
                  << std::setw(20) << (err_udiffsol.h1_norm_ux / err_sol.h1_norm_ux)
                  << std::endl;

        m_log->o << std::setw(10) << "h1(uy)"
                  << std::setw(20) << err_u.h1_norm_uy
                  << std::setw(20) << err_udiffsol.h1_norm_uy
                  << std::setw(20) << (err_udiffsol.h1_norm_uy / err_sol.h1_norm_uy)
                  << std::endl;
    }


    void set_problem(SPProblem problem) {
        this->m_problem = problem;
    }

    void set_num_ref(int ref){
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
        this->m_log->o << "level: " << index_level[level] << "\t" << factor << "\t" << index_level[level+1];
    }


    bool write(SPGridFunction u, int index, double time) override {
        int z = (index * this->max_index_precomputed) / index_level[0];
        int rem = (index * this->max_index_precomputed) % index_level[0];
        if (rem == 0 && index != 0) {
            SPGridFunction sol = u->clone_without_values();
            SPGridFunction udiffsol = u->clone();

            // write vtk output
            m_out_solution->print(this->m_sol_filename, *u, index, time);

            // load gridfunction file (ref solution)
            IOGridFunction<TDomain,TAlgebra> io = IOGridFunction<TDomain,TAlgebra>();
            std::stringstream ss_ref;
            if(this->num_ref == 3){
                ss_ref << this->path_ref_3;
            } else {
                ss_ref << this->path_ref_4;
            }
            ss_ref << zidx << ".gridfunction";
            io.read(udiffsol, ss_ref.str().c_str())

            // substract
            VecAdd(1, *udiffsol.get(), -1, *sol.get());

            // write vtk error
            m_out_diff->print(this->m_diff_filename, *udiffsol, index, time);

            // compute norms
            err_u.compute(u->clone());
            err_sol.compute(sol->clone());
            err_udiffsol.compute(udiffsol->clone());

            // write norms
            compare_norms(index, time, 0, 0)
        }
        return false; // no error
    };

    bool write(SPGridFunction u, int index, double time, int iteration, int level) override {


        int zidx = (index * this->max_index_precomputed) / index_level[0];
        int rem = (index * this->max_index_precomputed) % index_level[0];

        if (rem == 0 && index != 0) {
            SPGridFunction sol = u->clone_without_values();
            SPGridFunction udiffsol = u->clone();

            // write vtk output
            std::stringstream ss_solution;
            ss_solution << m_sol_filename << "_k" << iteration << "_l" << level << "_c" << count;
            m_out_solution->print(ss_solution.str().c_str(), *u, index, time);

            // load gridfunction file (ref solution)
            IOGridFunction<TDomain,TAlgebra> io = IOGridFunction<TDomain,TAlgebra>();
            std::stringstream ss_ref;
            if(this->num_ref == 3){
                ss_ref << this->path_ref_3;
            } else {
                ss_ref << this->path_ref_4;
            }
            ss_ref << zidx << ".gridfunction";
            io.read(udiffsol, ss_ref.str().c_str())
            // substract
            VecAdd(1, *udiffsol.get(), -1, *sol.get());

            // write vtk error
            std::stringstream ss_diff;
            ss_diff << m_diff_filename << "_k" << iteration << "_l" << level << "_c" << count;
            m_out_diff->print(ss_diff.str().c_str(), *udiffsol, index, time);

            // compute norms
            err_u.compute(u->clone());
            err_sol.compute(sol->clone());
            err_udiffsol.compute(udiffsol->clone());

            // write norms
            compare_norms(index, time, 0, 0)
        }

        return false; // no error
    };
};

#endif //UG_PLUGIN_XBRAIDBIOT_BRAIDBIOTPRECOMPUTED_H
