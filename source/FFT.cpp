#include "D:\my_programm\TDL\header\FFT.h"

FFT::FFT(int N)
{
    _P = calfftNum(N);
    _N = N;
    _fftSerise.resize(_N);
}

int FFT::calfftNum(int L)
{
    int p = 0;
    while (L > pow(2, p))
        p++;
    return p;
}
unsigned int FFT::reverse_bit(int num, const int p)
{
    int i;
    int bit;
    unsigned new_num = 0;
    for (i = 0; i < p; i++)
    {
        bit = num & 1;
        new_num <<= 1;
        new_num = new_num | bit;
        num >>= 1;
    }
    return new_num;
}
complex<double> FFT::Wn(const int m, const int N)
{
    complex<double> wn(cos(2 * acos(-1) * m / N), sin(-2 * acos(-1) * m / N));
    return wn;
}
void FFT::First_level_butterfly_operation(int i, int j, complex<double> Wn)
{
    complex<double> temp = _fftSerise[i];
    _fftSerise[i] = temp + Wn * _fftSerise[j];
    _fftSerise[j] = temp - Wn * _fftSerise[j];
}
vector<complex<double>> FFT::operator()(vector<double> &xn, const int p, const int N)
{
    int length = xn.size();
    if (length < N)
    {
        xn.resize(N);
        for (int n = xn.size(); n < N; n++)
        {
            xn[n] = 0;
        }
    }
    for (int j = 0; j < N; j++)
    {
        _fftSerise[j].real(xn[reverse_bit(j, p)]);
        _fftSerise[j].imag(0);
    }

    for (int k = 1; k <= p; k++)
    {
        // 每一层需要进行的蝶形运算
        for (int l = 0; l < N / pow(2, k); l++)
        {
            // 每一层的蝶形运算
            for (int n = 0; n < pow(2, k - 1); n++)
            {
                First_level_butterfly_operation(l * pow(2, k) + n, l * pow(2, k) + n + pow(2, k - 1), Wn(n, pow(2, k)));
            }
        }
    }
    return _fftSerise;
}

// 进行IFFT运算
vector<complex<double>> FFT::operator()(vector<complex<double>> &Xk, const int p, const int N, const int L)
{
    // 对输入序列Xk进行变换,取共轭，幅值*1/N
    int i = 0;
    while (i < Xk.size())
    {
        complex<double> temp = 0;
        temp.real(conj(Xk[i]).real() / N);
        temp.imag(conj(Xk[i]).imag() / N);
        Xk[i] = temp;
        i++;
    }
    // 进行FFT同样的运算
    int length = Xk.size();
    if (length < N)
    {
        Xk.resize(N);
        for (int n = Xk.size(); n < N; n++)
        {
            Xk[n] = 0;
        }
    }
    for (int j = 0; j < N; j++)
    {
        _fftSerise[j].real(Xk[reverse_bit(j, p)].real());
        _fftSerise[j].imag(Xk[reverse_bit(j, p)].imag());
    }

    for (int k = 1; k <= p; k++)
    {
        // 每一层需要进行的蝶形运算
        for (int l = 0; l < N / pow(2, k); l++)
        {
            // 每一层的蝶形运算
            for (int n = 0; n < pow(2, k - 1); n++)
            {
                First_level_butterfly_operation(l * pow(2, k) + n, l * pow(2, k) + n + pow(2, k - 1), Wn(n, pow(2, k)));
            }
        }
    }

    this->_xn.resize(L);
    for (int m = 0; m < L; m++)
    {
        _xn[m].real(_fftSerise[m].real());
        _xn[m].imag(_fftSerise[m].imag());
    }
    return _xn;
}

void FFT::fftshift(vector<complex<double>> &XK)
{
    vector<complex<double>> temp;
    complex<double> t;
    t.real(XK[0].real());
    t.imag(XK[0].imag());
    temp.resize(XK.size());

    for (int n = 0; n < XK.size(); n++)
    {
        if (n > XK.size() / 2 - 1)
        {
            temp[n].real(XK[n - XK.size() / 2].real());
            temp[n].imag(XK[n - XK.size() / 2].imag());
        }
        else
        {
            temp[n].real(XK[XK.size() / 2 + n].real());
            temp[n].imag(XK[XK.size() / 2 + n].imag());
        }
    }

    for (int i = 0; i < XK.size(); i++)
    {
        XK[i].real(temp[i].real());
        XK[i].imag(temp[i].imag());
    }
    XK[XK.size() / 2].real(t.real());
    XK[XK.size() / 2].imag(t.imag());
}