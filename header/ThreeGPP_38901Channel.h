#pragma once
#include <iostream>
#include <vector>
#include <complex>
#include <algorithm>

#include "TDLChannel.h"

using namespace std;

class ThreeGPP_38901Channel : public TDL_Model
{
public:
    ThreeGPP_38901Channel();
    ~ThreeGPP_38901Channel();

    vector<double> ThreeGppPowerScaledB(vector<double> &pModel, double &kDesire);

    double ThreeGppcomputeDopplerMax(const double &msV, const double &fc, const double &alpha);

    virtual vector<complex<double>> getPsdFilter(const double &tap_fm, const double &taps_fs,
                                                 const int &N_fading,
                                                 const PsdSpecturm &psdType) const;

    virtual void generate();

    void ThreeGpp_38901_TDL_A(const double &ds_desired, const double &msV,
                              const double &fc, const double &Os,
                              const double &N_fading, const double &Nsamples,
                              const double &alpha);

private:
    double _fm;
};