#include "D:\my_programm\TDL\header\TDLChannel.h"
#include "D:\my_programm\TDL\header\FFT.h"
#include <random>
#include <time.h>

TDL_Model::TDL_Model()
{
    // 默认构造参数
    _OS = 2;
    _centerFrequencyHz = 5e9;
    _maxDopplerFrequencyHz = 0;

    _hCir._nTaps = 0;
    _hCir._nFading = 0;
    _hCir._nSamples = 0;
    _hCir._intervalTimeOfSamples = 0;
    _hCir._nDelays.resize(_hCir._nTaps);
    _hCir._hCoeff.reserve(_hCir._nTaps * _hCir._nSamples);
}

TDL_Model::~TDL_Model()
{
}

void TDL_Model::setOs(const double &Nos)
{
    this->_OS = Nos;
}

void TDL_Model::setCenterFrequencyHz(const double &Fre)
{
    this->_centerFrequencyHz = Fre;
}

void TDL_Model::setMaxDopplerFrequencyHz(const double &fm)
{
    this->_maxDopplerFrequencyHz = fm;
}

void TDL_Model::setMsV(const double &v)
{
    this->_msV = v;
}

const double TDL_Model::getCenterFrequencyHz()
{
    return this->_centerFrequencyHz;
}

const double TDL_Model::getOS()
{
    return this->_OS;
}

const double TDL_Model::getMaxDopplerFrequencyHz()
{
    return this->_maxDopplerFrequencyHz;
}

const TDL_Model::TdlCir TDL_Model::getTdlCir()
{
    return this->_hCir;
}

void TDL_Model::setTdlCir(const int N_fading, const double N_samples,
                          const int nTaps)
{
    this->_hCir._nFading = N_fading;
    this->_hCir._nSamples = N_samples;
    this->_hCir._intervalTimeOfSamples = 3 * pow(10, 8) / _centerFrequencyHz / _OS / _msV;
    this->_hCir._fs = 1 / this->_hCir._intervalTimeOfSamples;
    this->_hCir._nTaps = nTaps;
    this->_hCir._timeOfSamples = this->_hCir._nSamples * this->_hCir._intervalTimeOfSamples;
}

vector<complex<double>> TDL_Model::getPsdFilter(const double &tap_fm, const double &taps_fs,
                                                const int &N_fading,
                                                const PsdSpecturm &psdType) const
{
    vector<complex<double>> psd_filter;
    vector<complex<double>> temp;
    vector<double> freq;

    double ts = 1 / taps_fs;
    double fs2 = ceil(taps_fs / 2);

    // 频率轴设置
    double i = 0;
    int j = 0;
    if (fs2 * 2 < N_fading)
    {
        freq.resize(N_fading);
        for (i = -N_fading / 2 + 0.5, j = 0; i < N_fading / 2; i++, j++)
        {
            freq[j] = i;
        }
    }
    else
    {
        freq.resize(fs2 * 2);
        for (i = -fs2 + 0.5, j = 0; i < fs2; i++, j++)
            freq[j] = i;
    }

    // 获得频域psd_filter系数
    psd_filter.resize(freq.size());
    for (int k = 0; k < freq.size(); k++)
    {
        switch (psdType)
        {
        case 0:
            // JAKES谱
            if (fabs(freq[k]) < tap_fm)
            {
                psd_filter[k].real(1 / sqrt(1 - pow(freq[k] / tap_fm, 2)));
                psd_filter[k].imag(0);
            }
            else
            {
                psd_filter[k].real(0);
                psd_filter[k].imag(0);
            }
            break;
        case 1:
            // RICE谱
            if (fabs(freq[k]) < tap_fm)
            {
                if (round(freq[k]) != round(0.7 * tap_fm))
                    psd_filter[k].real(0.41 / sqrt(1 - pow(freq[k] / tap_fm, 2)) / 2 / 3.14159);
                else
                    psd_filter[k].real(0.41 / sqrt(1 - pow(freq[k] / tap_fm, 2)) / 2 / 3.14159 + 0.91);
                psd_filter[k].imag(0);
            }
            else
            {
                psd_filter[k].real(0);
                psd_filter[k].imag(0);
            }
            break;
        default:
            cout << "psdtype error!" << endl;
        }
    }

    // for (int h = 0; h < psd_filter.size(); h++)
    // {
    //     cout << psd_filter[h].real();
    //     if (psd_filter[h].imag() >= 0)
    //     {
    //         cout << "+" << psd_filter[h].imag() << "i";
    //     }
    //     else
    //     {
    //         cout << psd_filter[h].imag() << "i";
    //     }
    //     cout << " ";
    // }
    // cout << endl;

    // 进行shift
    temp.resize(freq.size());
    for (int n = 0; n < freq.size(); n++)
    {
        if (n > freq.size() / 2 - 1)
        {
            temp[n].real(psd_filter[n - freq.size() / 2].real());
            temp[n].imag(0);
        }
        else
        {
            temp[n].real(psd_filter[freq.size() / 2 + n].real());
            temp[n].imag(0);
        }
    }

    // for (int h = 0; h < temp.size(); h++)
    // {
    //     cout << temp[h].real();
    //     if (temp[h].imag() >= 0)
    //     {
    //         cout << "+" << temp[h].imag() << "i";
    //     }
    //     else
    //     {
    //         cout << temp[h].imag() << "i";
    //     }
    //     cout << " ";
    // }
    // cout << endl;

    // 进行IFFT运算
    FFT ifft(N_fading);
    vector<complex<double>> XK;
    XK = ifft(temp, ifft._P, ifft._N, ifft._N);

    //     for (int h = 0; h < XK.size(); h++)
    // {
    //     cout << XK[h].real();
    //     if (XK[h].imag() >= 0)
    //     {
    //         cout << "+" << XK[h].imag() << "i";
    //     }
    //     else
    //     {
    //         cout << XK[h].imag() << "i";
    //     }
    //     cout << " ";
    // }
    // cout << endl;

    // 需要进行一步shift运算
    ifft.fftshift(XK);
    // 进行归一化
    int m = 0;
    double sum = 0;
    for (int l = 0; l < XK.size(); l++)
    {
        sum = sum + pow(XK[l].real(), 2);
    }
    while (m < XK.size())
    {
        XK[m].imag(0);
        double t;
        t = XK[m].real();
        XK[m].real(t / sqrt(sum));
        m++;
    }

    // for (int h = 0; h < XK.size(); h++)
    // {
    //     cout << XK[h].real();
    //     if (XK[h].imag() >= 0)
    //     {
    //         cout << "+" << XK[h].imag() << "i";
    //     }
    //     else
    //     {
    //         cout << XK[h].imag() << "i";
    //     }
    //     cout << " ";
    // }
    // cout << endl;

    return XK;
}

void TDL_Model::addTaps(const double &delayNs, const double &powerdB,
                        const double &tapFd, const PsdSpecturm &psdType)
{
    // 抽头数据
    TapCharacters tap;
    tap._tapDelays = delayNs;
    tap._tapPowerdB = pow(10.0, powerdB * 0.1); // 转化为线性值
    tap._tapDopplerFrequencyHz = tapFd;
    tap._psdType = psdType;
    tap.los_flag = 0;

    // 将抽头添加到抽头容器中
    this->_taps.push_back(tap);
    this->_hCir._nTaps++;
    this->_hCir._nDelays.push_back(tap._tapDelays);
    this->_hCir._hCoeff.reserve(this->_hCir._nTaps * this->_hCir._nSamples);
}

void TDL_Model::addTaps_Los(const double &delayNs, const double &powerdB,
                            const double &tapFd, const PsdSpecturm &psdType,
                            const double &Nsamples, const double &alpha,
                            const double &K_dB, const double &losPdB)
{
    // 抽头数据
    TapCharacters tap;
    tap._tapDelays = delayNs;
    tap._tapPowerdB = pow(10.0, powerdB * 0.1); // 转化为线性值
    tap.K_dB = pow(10.0, K_dB * 0.1);
    tap._tapDopplerFrequencyHz = tapFd;
    tap._psdType = psdType;
    tap.los_flag = 1;
    tap.los_p = pow(10.0, losPdB * 0.1);
    tap.alpha = alpha * 3.14169 / 180;

    // 将抽头添加到抽头容器中
    this->_taps.push_back(tap);
    this->_hCir._nTaps++;
    this->_hCir._nDelays.push_back(tap._tapDelays);
    this->_hCir._hCoeff.reserve(this->_hCir._nTaps * this->_hCir._nSamples);
}

void TDL_Model::generate()
{
    for (int h = 0; h < this->_hCir._nTaps; h++)
    {

        // 生成IQ两路数据
        vector<complex<double>> G;
        vector<double> G1;
        vector<double> G2;

        G.resize(this->_hCir._nSamples, 0.0);
        G1.resize(this->_hCir._nSamples + this->_hCir._nFading - 1);
        G2.resize(this->_hCir._nSamples + this->_hCir._nFading - 1);

        // 多普勒系数--先测试一个的时候
        vector<complex<double>> psd;
        psd.resize(this->_hCir._nFading);
        psd = this->getPsdFilter(_taps[h]._tapDopplerFrequencyHz, _hCir._fs,
                                 _hCir._nFading, _taps[h]._psdType);
        // 定义高斯分布随机生成器
        const double mean = 0.0;                          // 均值
        const double stddev = 1;                          // 标准差
        std::default_random_engine generator(time(NULL)); // 增加随机性
        // std::default_random_engine generator;
        std::normal_distribution<double> dist(mean, stddev);

        // 将高斯分布添加到IQ两路数据
        for (auto &g : G)
        {
            g.real(sqrt(_taps[h]._tapPowerdB) * sqrt(0.5) * dist(generator));
            g.imag(sqrt(_taps[h]._tapPowerdB) * sqrt(0.5) * dist(generator));
        }

        int n = psd.size(), m = G.size();
        int k = n + m - 1;
        int i, j;
        for (i = 0; i < k; i++)
        {
            for (j = max(0, i + 1 - n); j <= min(i, m - 1); j++)
            {
                G1[i] += G[j].real() * psd[i - j].real();
                G2[i] += G[j].imag() * psd[i - j].real();
            }
        }

        G.resize(k - n / 2 - n + 1 + n / 2);
        for (int l = n / 2 - 1, r = 0; l < k - n / 2; l++, r++)
        {
            G[r].real(G1[l]);
            G[r].imag(G2[l]);
        }

        // G.resize(k -2*n+2);
        // for (int l = n-1, r = 0; l < k-n-1; l++, r++)
        // {
        //     G[r].real(G1[l]);
        //     G[r].imag(G2[l]);
        // }

        // 如果是直射径，在滤波后的系数中直接叠加直射信号
        if (_taps[h].los_flag == 1)
        {

            for (int kk = 0; kk < G.size(); kk++)
            {
                G[kk].real(G[kk].real() + sqrt(_taps[0].los_p) * cos(2 * 3.14159 * _taps[0]._tapDopplerFrequencyHz * cos(_taps[h].alpha) * kk * this->_hCir._intervalTimeOfSamples));
                G[kk].imag(G[kk].imag() + sqrt(_taps[0].los_p) * sin(2 * 3.14159 * _taps[0]._tapDopplerFrequencyHz * cos(_taps[h].alpha) * kk * this->_hCir._intervalTimeOfSamples));
            }
        }
        this->_hCir._hCoeff.push_back(G);
    }
}

bool TDL_Model::run()
{
    return true;
}
