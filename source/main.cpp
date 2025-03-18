#include "D:\my_programm\TDL\header\FFT.h"
#include "D:\my_programm\TDL\header\TDLChannel.h"
#include "D:\my_programm\TDL\header\SatelliteChannel.h"
#include "D:\my_programm\TDL\header\ThreeGPP_38901Channel.h"
#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include <fstream>
using namespace std;
// 存入txt
void save_txt(vector<vector<complex<double>>> hCoeff, TDL_Model &tdl)
{
    int nTaps = hCoeff.size();
    int nSamples = hCoeff[0].size();
    ofstream ofs;
    ofs.open("D:/my_programm/TDL/text.txt", ios::out);
    for (int i = 0; i < nTaps; i++)
    {
        // ofs << "Tap" << i << ": ";
        // ofs << endl;
        for (int j = 0; j < nSamples; j++)
        {
            ofs << hCoeff[i][j].real();
            if (hCoeff[i][j].imag() > 0)
            {
                ofs << "+" << hCoeff[i][j].imag() << "i";
            }
            else
            {
                ofs << hCoeff[i][j].imag() << "i";
            }
            ofs << " ";
            ofs << tdl.getTdlCir()._intervalTimeOfSamples * j;
            ofs << endl;
            // ofs << " ";
        }
        // ofs << endl;
    }
    // ofs << "All tap data has been read!";
    ofs.close();
}

void psd_filter_test()
{
    double taps_fs = 240;
    double tap_fm = 30;
    vector<complex<double>> psd_filter;
    vector<complex<double>> temp;
    vector<double> freq;

    double ts = 1 / taps_fs;
    double fs2 = ceil(taps_fs / 2);
    // srand((unsigned int)time(NULL)); // 生成随机数种子
    //  频率轴设置
    double i = 0;
    int j = 0;
    freq.resize(fs2 * 2);
    for (i = -fs2 + 0.5, j = 0; i < fs2; i++, j++)
    {
        freq[j] = i;
    }
    // 获得频域psd_filter系数
    psd_filter.resize(freq.size());
    for (int k = 0; k < freq.size(); k++)
    {
        if (fabs(freq[k]) < tap_fm - 2)
        {
            psd_filter[k].real(1 / sqrt(1 - pow(freq[k] / tap_fm, 2)));
            psd_filter[k].imag(0);
        }
        else
        {
            psd_filter[k].real(0);
            psd_filter[k].imag(0);
        }
        // cout << psd_filter[k] << " ";
    }
    cout << endl;
    // 这里是对频域进行shift
    temp.resize(freq.size());
    for (int n = 0; n < freq.size(); n++)
    {
        temp[n].real(psd_filter[freq.size() / 2 + n].real());
        temp[n].imag(0);
        if (n > freq.size() / 2)
        {
            temp[n].real(psd_filter[n - freq.size() / 2].real());
            temp[n].imag(0);
        }
        // cout << temp[n].real() << " ";
    }
    cout << endl;
    // 进行ifft转换成时域序列
    FFT ifft(256);
    vector<complex<double>> XK;
    XK = ifft(temp, ifft._P, ifft._N, ifft._N);
    // 进行fftshift
    ifft.fftshift(XK);
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
        cout << XK[m].real() << " ";
        m++;
    }
    cout << endl;
}
void fft_test()
{
    vector<double> xn = {1, 2, 3, 4, 5, 6, 7, 8};
    int length = xn.size();
    FFT fft(256);
    vector<complex<double>> XK;
    XK = fft(xn, fft._P, fft._N);
    int m = 0;
    while (m < fft._N)
    {
        cout << XK[m].real() << "+" << XK[m].imag() << "i"
             << " ";
        m++;
    }
    cout << endl;
    cout << "hello" << endl;

    vector<complex<double>> X;
    // 这里输入的是FFT变换后的X序列，进行的是ifft运算
    X = fft(XK, fft._P, fft._N, length);
    m = 0;
    while (m < X.size())
    {
        cout << X[m] << " ";
        m++;
    }
    cout << endl;
}

void gauss_test()
{
    // 随机数种子
    srand((unsigned int)time(NULL)); // 生成随机数种子

    std::vector<double> data(5, 1.0); // 模拟的多普勒系数
    std::vector<double> G1;
    std::vector<double> G2;
    vector<complex<double>> G;
    G.resize(20, 0.0);
    G1.resize(G.size() + data.size() - 1);
    G2.resize(G.size() + data.size() - 1);

    // Define random generator with Gaussian distribution
    const double mean = 0.0; // 均值
    const double stddev = 1; // 标准差
    std::default_random_engine generator;
    std::normal_distribution<double> dist(mean, stddev);

    // 将高斯噪声添加到IQ两路数据
    for (auto &g : G)
    {
        g.real(dist(generator));
        g.imag(dist(generator));
        cout << g.real() << " ";
    }
    cout << endl;

    // 进行卷积运算
    int n = data.size(), m = G.size();
    int k = m + n - 1;
    int i, j;
    for (i = 0; i < k; i++)
    {
        for (j = max(0, i + 1 - n); j <= min(i, m - 1); j++)
        {
            G1[i] += G[j].real() * data[i - j];
            G2[i] += G[j].imag() * data[i - j];
        }
        cout << G1[i] << " ";
    }
    cout << endl;

    // 进行截取
    G.resize(k - 2 * n);
    for (int l = n - 1, r = 0; l < k - n - 1; l++, r++)
    {
        G[r].real(G1[l]);
        G[r].imag(G2[l]);
        cout << G[r] << " ";
    }

    // Add Gaussian noise
    for (auto &x : data)
    {
        x = x + dist(generator);
        // cout << x << " ";
    }
    cout << endl;
}
void sat_test()
{
    // 基本配置
    SatelliteChannel sat;
    double ds_desired = 10;
    double alpha_mode = 60;
    double msV = 10;
    double fc = 4 * pow(10, 9);
    double Os = 7.5;
    double N_fading = 1024;
    double Nsamples = 5000;
    sat.setSatV(7000);
    sat.NTN_38811_TDL_D(ds_desired, alpha_mode, msV, fc, Os, N_fading, Nsamples, 11.707);
    // 存入文件
    save_txt(sat.getTdlCir()._hCoeff, sat);
    cout << sat.getTdlCir()._hCoeff.at(0).size() << endl;
}
void ThreeGpp_test()
{
    // 基本配置
    ThreeGPP_38901Channel ThreeGpp;
    double ds_desired = 10;
    double alpha_mode = 0;
    double msV = 10;
    double fc = 4 * pow(10, 9);
    double Os = 7.5;
    double N_fading = 1024;
    double Nsamples = 5000;
    ThreeGpp.ThreeGpp_38901_TDL_A(ds_desired, msV, fc, Os, N_fading, Nsamples, alpha_mode);
    // 存入文件
    save_txt(ThreeGpp.getTdlCir()._hCoeff, ThreeGpp);
    cout << ThreeGpp.getTdlCir()._hCoeff.at(0).size() << endl;
}
void tdl_test()
{
    // 基本配置
    TDL_Model tdl;
    tdl.setMsV(10);
    tdl.setCenterFrequencyHz(4 * pow(10, 9));
    tdl.setOs(7.5);
    tdl.setTdlCir(1024, 5000, 0);
    double fm;
    fm = 20.0 / 300.0 * 900.0;
    tdl.setMaxDopplerFrequencyHz(fm);

    // 添加抽头
    tdl.addTaps(0, 0, fm, tdl.FLAT);

    tdl.addTaps(0.7249, -1.973, fm, tdl.FLAT);

    tdl.addTaps(0.7410, -4.332, fm, tdl.FLAT);

    tdl.addTaps(5.7392, -11.914, fm, tdl.FLAT);

    tdl.generate();
    // 存入文件
    save_txt(tdl.getTdlCir()._hCoeff, tdl);
    cout << tdl.getTdlCir()._hCoeff.at(0).size() << endl;
}
int main()
{

    ThreeGpp_test();
    //  sat_test();
    //  psd_filter_test();
    //  tdl_test();
    return 0;
}

// 需要进行一步shift运算
// for (int n = 0; n < XK.size(); n++)
// {

//     temp[n].real(XK[XK.size() / 2 + n].real());
//     temp[n].imag(XK[XK.size() / 2 + n].imag());
//     if (n > XK.size() / 2)
//     {
//         temp[n].real(XK[n - XK.size() / 2].real());
//         temp[n].imag(XK[n - XK.size() / 2].imag());
//     }
//     //  cout << temp[n] << " ";
// }

// 归一化操作
// int m = 0;
// double sum = 0;
// for (int l = 0; l < temp.size(); l++)
// {
//     sum = sum + pow(temp[l].real(), 2);
// }
// while (m < temp.size())
// {
//     temp[m].imag(0);
//     double t;
//     t = temp[m].real();
//     temp[m].real(t / sqrt(sum));
//     cout << temp[m].real() << " ";
//     m++;
// }
// cout << endl;
