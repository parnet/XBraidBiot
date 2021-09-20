//
// Created by parnet on 14.06.21.
//

#ifndef UG_PLUGIN_XBRAIDBIOT_BRAIDBIOTCONTROL_H
#define UG_PLUGIN_XBRAIDBIOT_BRAIDBIOTCONTROL_H

#include "../../XBraidForUG4/src/interface/scriptor.h"

#include "../../Poroelasticity/src/biot_tools.h"
#include "../../Poroelasticity/src/barry_mercer.h"

template<typename TDomain, typename TAlgebra>
class BraidBiotCheck : public Scriptor<TDomain, TAlgebra>{
public:

    typedef ug::GridFunction <TDomain, TAlgebra> TGridFunction;
    typedef SmartPtr <TGridFunction> SPGridFunction;

    typedef ug::Poroelasticity::BarryMercerProblem<TDomain,TAlgebra> TProblem;
    typedef SmartPtr<TProblem> SPProblem;

    SPProblem m_problem;
    int napprox = 16;

    BraidBiotCheck() : Scriptor<TDomain, TAlgebra>(){

    };

    ~BraidBiotCheck() = default;

    void set_problem(SPProblem problem){
        this->m_problem = problem;
    }
    void set_napprox(int napprox){
        this->napprox = napprox;
    }


    virtual bool write(SPGridFunction u, int index, double time) {
        m_problem->m_errData.napprox = this->napprox;
        m_problem->m_errData.iteration = -1;
        m_problem->m_errData.level = -1;
        m_problem->post_processing(u,index, time);
        return false;
    };

    virtual bool write(SPGridFunction u, int index, double time, int iteration, int level) {
        m_problem->m_errData.napprox = this->napprox;
        m_problem->m_errData.iteration = iteration;
        m_problem->m_errData.level = level;

        m_problem->post_processing(u,index, time);

        m_problem->m_errData.iteration = -1;
        m_problem->m_errData.level = -1;
        return false;
    } ;

};


#endif //UG_PLUGIN_XBRAIDBIOT_BRAIDBIOTCONTROL_H
