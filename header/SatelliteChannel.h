#pragma once
#include <iostream>
#include <vector>
#include <complex>
#include <algorithm>
#include "TDLChannel.h"

using namespace std;

class SatelliteChannel : public TDL_Model
{
public:
    SatelliteChannel();
    ~SatelliteChannel();

    // 设置参数方法
    void setElevation(const double &elevation);
    void setSatV(const double &satV);
    void setEarthR(const double &earthR);
    void setSatH(const double &satH);

    vector<double> satPowerScaledB(vector<double> &pModel, double &kDesire);

    double computeDopplerMax(const double &msv, const double &alpha,
                             const double &fc);

    virtual vector<complex<double>> getPsdFilter(const double &tap_fm, const double &taps_fs,
                                                 const int &N_fading,
                                                 const PsdSpecturm &psdType) const;

    virtual void generate();

    void NTN_38811_TDL_A(const double &ds_desired, const double &alpha_model,
                         const double &msV, const double &fc,
                         const double &Os, const double &N_fading,
                         const double &Nsamples);
    void NTN_38811_TDL_B(const double &ds_desired, const double &alpha_model,
                         const double &msV, const double &fc,
                         const double &Os, const double &N_fading,
                         const double &Nsamples);

    void NTN_38811_TDL_C(const double &ds_desired, const double &alpha_model,
                         const double &msV, const double &fc,
                         const double &Os, const double &N_fading,
                         const double &Nsamples, const double &kDesired);

    void NTN_38811_TDL_D(const double &ds_desired, const double &alpha_model,
                         const double &msV, const double &fc,
                         const double &Os, const double &N_fading,
                         const double &Nsamples, const double &kDesired);

private:
    double _elevation; // 定义卫星地面仰角
    double _satV;      // 定义卫星移动速度
    double _earthR;    // 定义地球半径
    double _satH;      // 定义卫星高度
    double _fdShift;   // 卫星运动引起的多普勒
    double _satFm;
};