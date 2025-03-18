#include "D:\my_programm\TDL\header\ThreeGPP_38901Channel.h"
#include "D:\my_programm\TDL\header\FFT.h"
#include <random>
#include <time.h>

ThreeGPP_38901Channel::ThreeGPP_38901Channel()
{
}

ThreeGPP_38901Channel::~ThreeGPP_38901Channel()
{
}

vector<double> ThreeGPP_38901Channel::ThreeGppPowerScaledB(vector<double> &pModel, double &kDesire)
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

vector<complex<double>> ThreeGPP_38901Channel::getPsdFilter(const double &tap_fm, const double &taps_fs,
                                                            const int &N_fading,
                                                            const PsdSpecturm &psdType) const
{
    return this->TDL_Model::getPsdFilter(tap_fm, taps_fs, N_fading, psdType);
}

void ThreeGPP_38901Channel::generate()
{
    this->TDL_Model::generate();
}

double ThreeGPP_38901Channel::ThreeGppcomputeDopplerMax(const double &msV, const double &fc, const double &alpha)
{
    this->_fm = (msV / 3 / pow(10, 8) * (fc)*cos(alpha * 3.14159 / 180));

    return this->_fm;
}

void ThreeGPP_38901Channel::ThreeGpp_38901_TDL_A(const double &ds_desired, const double &msV,
                                                 const double &fc, const double &Os,
                                                 const double &N_fading, const double &Nsamples,
                                                 const double &alpha)
{
    this->setMsV(msV);
    this->setCenterFrequencyHz(fc);
    this->setOs(Os);
    this->setTdlCir(N_fading, Nsamples, 0);
    this->ThreeGppcomputeDopplerMax(msV, fc, alpha);
    this->setMaxDopplerFrequencyHz(this->_fm);
    this->addTaps(0 * ds_desired, -13.4, this->_fm, this->FLAT);
    this->addTaps(0.3819 * ds_desired, 0, this->_fm, this->FLAT);
    this->addTaps(0.4025 * ds_desired, -2.2, this->_fm, this->FLAT);
    this->addTaps(0.5868 * ds_desired, -4, this->_fm, this->FLAT);
    this->addTaps(0.4620 * ds_desired, -6, this->_fm, this->FLAT);
    this->addTaps(0.5375 * ds_desired, -8.2, this->_fm, this->FLAT);
    this->addTaps(0.6708 * ds_desired, -9.9, this->_fm, this->FLAT);
    this->addTaps(0.5750 * ds_desired, -10.5, this->_fm, this->FLAT);
    this->addTaps(0.7618 * ds_desired, -7.5, this->_fm, this->FLAT);
    this->addTaps(1.5375 * ds_desired, -15.9, this->_fm, this->FLAT);
    this->addTaps(1.8978 * ds_desired, -6.6, this->_fm, this->FLAT);
    this->addTaps(2.2242 * ds_desired, -16.7, this->_fm, this->FLAT);
    this->addTaps(2.1718 * ds_desired, -12.4, this->_fm, this->FLAT);
    this->addTaps(2.4942 * ds_desired, -15.2, this->_fm, this->FLAT);
    this->addTaps(2.5119 * ds_desired, -10.8, this->_fm, this->FLAT);
    this->addTaps(3.0582 * ds_desired, -11.3, this->_fm, this->FLAT);
    this->addTaps(4.0810 * ds_desired, -12.7, this->_fm, this->FLAT);
    this->addTaps(4.4579 * ds_desired, -16.2, this->_fm, this->FLAT);
    this->addTaps(4.5695 * ds_desired, -18.3, this->_fm, this->FLAT);
    this->addTaps(4.7966 * ds_desired, -18.9, this->_fm, this->FLAT);
    this->addTaps(5.0066 * ds_desired, -16.6, this->_fm, this->FLAT);
    this->addTaps(5.3043 * ds_desired, -19.9, this->_fm, this->FLAT);
    this->addTaps(9.6586 * ds_desired, -29.7, this->_fm, this->FLAT);
    this->generate();
}