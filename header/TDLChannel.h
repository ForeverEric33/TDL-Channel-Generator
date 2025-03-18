#pragma once
#include <iostream>
using namespace std;
#include <vector>
#include <complex>
#include <algorithm>
#include <time.h>

using namespace std;

class TDL_Model
{
    // 所需的TDL信道数据结构
public:
    // 多普勒谱
    enum PsdSpecturm
    {
        FLAT = 0,
        RICE = 1
    };
    // 单个抽头数据
    struct TapCharacters
    {
        double _tapDelays;
        double _tapPowerdB;
        double _tapDopplerFrequencyHz;
        double alpha;
        double K_dB;
        int los_flag;
        double los_p;
        PsdSpecturm _psdType;
    };
    // 信道冲击响应数据结构
    struct TdlCir
    {
        int _nTaps;
        int _nFading;
        double _nSamples;
        double _fs;
        double _intervalTimeOfSamples;
        double _timeOfSamples;
        // 信道时延存储矩阵
        vector<double> _nDelays;
        // 信道系数矩阵 _nTaps*_nSamples
        vector<vector<complex<double>>> _hCoeff;
    };

public:
    TDL_Model();
    virtual ~TDL_Model();
    // 设置参数方法
    void setTdlCir(const int N_fading, const double N_samples,
                   const int nTaps);                 // 设置信道CIR数据参数
    void setCenterFrequencyHz(const double &Fre);    // 设置中心频率
    void setOs(const double &Nos);                   // 设置过采样因子
    void setMaxDopplerFrequencyHz(const double &fd); // 设置最大多普勒频移
    void setMsV(const double &v);                    // 设置终端移动速度

    // 获取信道参数
    const TdlCir getTdlCir();
    const double getCenterFrequencyHz();
    const double getOS();
    const double getMaxDopplerFrequencyHz();

    // 获取时域多普勒滤波系数
    virtual vector<complex<double>> getPsdFilter(const double &tap_fm, const double &taps_fs,
                                                 const int &N_fading,
                                                 const PsdSpecturm &psdType) const;

    // 添加抽头
    void addTaps(const double &delayNs, const double &powerdB,
                 const double &tapFd, const PsdSpecturm &psdType);

    // 添加直射径------这里采用的是class直接加直射信号，其实也可以采用RICE谱
    void addTaps_Los(const double &delayNs, const double &powerdB,
                     const double &tapFd, const PsdSpecturm &psdType,
                     const double &Nsamples, const double &alpha,
                     const double &K_dB, const double &losPdB);

    // 生成信道CIR的核心算法
    virtual void generate();

    // 运行仿真
    bool run();

protected:
    TdlCir _hCir;                  // 信道响应输出
    vector<TapCharacters> _taps;   // 抽头数据
    double _centerFrequencyHz;     // 中心频率
    double _OS;                    // 过采样因子
    double _maxDopplerFrequencyHz; // 最大多普勒频移
    double _msV;                   // 终端移动速度
};