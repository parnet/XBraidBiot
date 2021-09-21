//
// Created by parnet on 2021-09-21.
//

#ifndef UG_PLUGIN_XBRAIDBIOT_BIOTERRORDATA_H
#define UG_PLUGIN_XBRAIDBIOT_BIOTERRORDATA_H

class BiotErrorData {
public:
    double l2_norm_p = 0.0;
    double l2_norm_ux = 0.0;
    double l2_norm_uy = 0.0;

    double h1_norm_ux = 0.0;
    double h1_norm_uy = 0.0;


    int porder = 2;
    int uorder = 4;

    void compute(SPGridFunction u) {
        auto *uref = u.get();

        this->l2_norm_p = L2Norm(uref, "p", porder);;
        this->l2_norm_ux = L2Norm(uref, "ux", uorder);
        this->l2_norm_uy = L2Norm(uref, "uy", uorder);

        this->h1_norm_ux = H1SemiNorm(uref, "ux", uorder);
        this->h1_norm_uy = H1SemiNorm(uref, "uy", uorder);

    }
};



#endif //UG_PLUGIN_XBRAIDBIOT_BIOTERRORDATA_H
