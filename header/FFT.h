#pragma once
#include <iostream>
#include <vector>
#include <complex>
#include <algorithm>

using namespace std;

class FFT
{
public:
    FFT(int N);
    vector<complex<double>> operator()(vector<double> &xn, const int p, const int N);                       // FFT运算
    vector<complex<double>> operator()(vector<complex<double>> &Xk, const int p, const int N, const int L); // IFFT运算

    void fftshift(vector<complex<double>> &XK);
    unsigned int reverse_bit(int num, const int p);
    complex<double> Wn(const int m, const int N);
    void First_level_butterfly_operation(int i, int j, complex<double> Wn);
    int calfftNum(int L);

public:
    int _P; // 蝶形运算级数
    int _N; // FFT运算点数

private:
    vector<complex<double>> _fftSerise;
    vector<complex<double>> _xn;
};
