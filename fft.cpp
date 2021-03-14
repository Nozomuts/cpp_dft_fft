#include <complex>
#include <iostream> // for cout

#define N 64            // 分割数
#define Fs (double)8000 // サンプリング周波数
#define A (double)1     // 振幅
#define F0 (double)440  // 周波数
#define phi (double)0   // 初期位相

using namespace std; // cout, endl, string, to_string

void fft(complex<double> *x, int _N, int q) {
    if (_N > 1) {
        for (int k = 0; k < _N / 2; k++) {
            const complex<double> t = exp(complex<double>(0, -2 * M_PI * k / _N));
            const complex<double> a = x[q + k];
            const complex<double> b = x[q + k + _N / 2];
            x[q + k] = a + b;
            x[q + k + _N / 2] = (a - b) * t;
        }
        // 再帰
        fft(x, _N / 2, q);
        fft(x, _N / 2, _N / 2 + q);
    }
}

// ビット反転並べ替え
void bit_reverse(complex<double> *x, int _N) {
    for (int i = 0, j = 1; j < _N - 1; j++) {
        for (int k = _N >> 1; k > (i ^= k); k >>= 1)
            ;
        if (i < j)
            std::swap(x[i], x[j]); // x[i]とx[j]を交換する
    }
}

int main() {
    complex<double> x_out[N]; // 元データ

    // 元データ作成
    for (int i = 0; i < N; i++) {
        // x(t) = A * sin(2 * pi * F0 * t + phi) ( 0 <= t < 0.008 )
        x_out[i] = complex<double>(A * sin(2 * M_PI * F0 * i / Fs + phi), 0); // t = i / Fs
    }

    fft(x_out, N, 0);
    bit_reverse(x_out, N);
    for (int i = 0; i < N; i++) {
        cout << x_out[i] << endl;
    }
    return 0;
}
