#include "D:\my_programm\TDL\header\SatelliteChannel.h"
#include "D:\my_programm\TDL\header\FFT.h"
#include <random>
#include <time.h>

SatelliteChannel::SatelliteChannel()
{
    // 默认构造参数
    _satH = 35786 * pow(10, 3); // 同步卫星
    _earthR = 6371 * pow(10, 3);
    _satV = 3007;
    _elevation = 0;
}

SatelliteChannel::~SatelliteChannel()
{
}

void SatelliteChannel::setElevation(const double &elevation)
{
    this->_elevation = elevation;
}

void SatelliteChannel::setSatV(const double &satV)
{
    this->_satV = satV;
}

void SatelliteChannel::setEarthR(const double &earthR)
{
    this->_earthR = earthR;
}

void SatelliteChannel::setSatH(const double &satH)
{
    this->_satH = satH;
}

double SatelliteChannel::computeDopplerMax(const double &msv, const double &alpha,
                                           const double &fc)
{
    this->_fdShift = this->_satV / 3 / pow(10, 8) * (this->_earthR / (this->_earthR + this->_satH)) * cos(this->_elevation * 3.14159 / 180) * fc;
    this->_satFm = msv / 3 / pow(10, 8) * (fc + this->_fdShift) * cos(alpha * 3.14159 / 180);
    return this->_satFm;
}

vector<complex<double>> SatelliteChannel::getPsdFilter(const double &tap_fm, const double &taps_fs,
                                                       const int &N_fading,
                                                       const PsdSpecturm &psdType) const
{
    return this->TDL_Model::getPsdFilter(tap_fm, taps_fs, N_fading, psdType);
}

void SatelliteChannel::generate()
{
    this->TDL_Model::generate();
}

// 功率缩放
vector<double> SatelliteChannel::satPowerScaledB(vector<double> &pModel, double &kDesire)
{
    double K_model = 0;
    double sum = 0;
    vector<double> P_scaled;
    P_scaled.resize(pModel.size());
    for (int i = 0; i < pModel.size(); i++)
    {
        sum = pow(10, pModel[i] * 0.1) + sum;
    }
    K_model = pModel[0] - 10 * log10(sum);
    for (int j = 0; j < pModel.size(); j++)
    {
        P_scaled[j] = pModel[j] - kDesire + K_model;
    }
}

void SatelliteChannel::NTN_38811_TDL_A(const double &ds_desired, const double &alpha_model,
                                       const double &msV, const double &fc,
                                       const double &Os, const double &N_fading,
                                       const double &Nsamples)
{
    this->setMsV(msV);
    this->setCenterFrequencyHz(fc);
    this->setOs(Os);
    this->setTdlCir(N_fading, Nsamples, 0);
    this->computeDopplerMax(this->_msV, alpha_model, this->_centerFrequencyHz);
    this->setMaxDopplerFrequencyHz(this->_satFm);
    this->addTaps(0 * ds_desired, 0, this->_satFm, this->FLAT);
    this->addTaps(1.0811 * ds_desired, -4.675, this->_satFm, this->FLAT);
    this->addTaps(2.8416 * ds_desired, -6.482, this->_satFm, this->FLAT);
    this->generate();
}

void SatelliteChannel::NTN_38811_TDL_B(const double &ds_desired, const double &alpha_model,
                                       const double &msV, const double &fc,
                                       const double &Os, const double &N_fading,
                                       const double &Nsamples)
{
    this->setMsV(msV);
    this->setCenterFrequencyHz(fc);
    this->setOs(Os);
    this->setTdlCir(N_fading, Nsamples, 0);
    this->computeDopplerMax(this->_msV, alpha_model, this->_centerFrequencyHz);
    this->setMaxDopplerFrequencyHz(this->_satFm);
    this->addTaps(0 * ds_desired, 0, this->_satFm, this->FLAT);
    this->addTaps(0.7249 * ds_desired, -1.973, this->_satFm, this->FLAT);
    this->addTaps(0.7410 * ds_desired, -4.332, this->_satFm, this->FLAT);
    this->addTaps(5.7392 * ds_desired, -11.914, this->_satFm, this->FLAT);
    this->generate();
}

void SatelliteChannel::NTN_38811_TDL_C(const double &ds_desired, const double &alpha_model,
                                       const double &msV, const double &fc,
                                       const double &Os, const double &N_fading,
                                       const double &Nsamples, const double &kDesired)
{
    this->setMsV(msV);
    this->setCenterFrequencyHz(fc);
    this->setOs(Os);
    this->setTdlCir(N_fading, Nsamples, 0);
    this->computeDopplerMax(this->_msV, alpha_model, this->_centerFrequencyHz);
    this->setMaxDopplerFrequencyHz(this->_satFm);
    double k_linear = pow(10, kDesired * 0.1);
    if (kDesired != 10.224)
        this->addTaps_Los(0 * ds_desired, 20 * log10(sqrt(1 / (k_linear + 1))), this->_satFm, this->FLAT, N_fading, alpha_model, kDesired, 20 * log10(sqrt(k_linear / (k_linear + 1))));
    else
        this->addTaps_Los(0 * ds_desired, 20 * log10(sqrt(1 / (k_linear + 1))), this->_satFm, this->FLAT, N_fading, alpha_model, 10.224, 20 * log10(sqrt(k_linear / (k_linear + 1))));
    this->addTaps(14.8124 * ds_desired, -23.373, this->_satFm, this->FLAT);
    this->generate();
}

void SatelliteChannel::NTN_38811_TDL_D(const double &ds_desired, const double &alpha_model,
                                       const double &msV, const double &fc,
                                       const double &Os, const double &N_fading,
                                       const double &Nsamples, const double &kDesired)
{
    this->setMsV(msV);
    this->setCenterFrequencyHz(fc);
    this->setOs(Os);
    this->setTdlCir(N_fading, Nsamples, 0);
    this->computeDopplerMax(this->_msV, alpha_model, this->_centerFrequencyHz);
    this->setMaxDopplerFrequencyHz(this->_satFm);
    double k_linear = pow(10, kDesired * 0.1);
    if (kDesired != 11.707)
        this->addTaps_Los(0 * ds_desired, 20 * log10(sqrt(1 / (k_linear + 1))), this->_satFm, this->FLAT, N_fading, alpha_model, kDesired, 20 * log10(sqrt(k_linear / (k_linear + 1))));
    else
        this->addTaps_Los(0 * ds_desired, 20 * log10(sqrt(1 / (k_linear + 1))), this->_satFm, this->FLAT, N_fading, alpha_model, 11.707, 20 * log10(sqrt(k_linear / (k_linear + 1))));
    this->addTaps(0.5596 * ds_desired, -9.887, this->_satFm, this->FLAT);
    this->addTaps(7.3340 * ds_desired, -16.771, this->_satFm, this->FLAT);
    this->generate();
}
