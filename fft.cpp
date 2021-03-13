#include <complex>
#include <math.h>
#include <stdio.h>

#define N 80            // 分割数
#define Fs (double)8000 // サンプリング周波数
#define A (double)1     // 振幅
#define F0 (double)440  // 周波数

using namespace std; // cout, endl, string, to_string

void fft(std::complex<double> *x, int _N) {
    // Nが1以下で早期リターン
    if (_N <= 1) {
        return;
    }

    // 偶奇に分割
    std::complex<double> even[_N / 2]; // 偶数
    std::complex<double> odd[_N / 2];  // 奇数
    for (int i = 0; i < _N / 2; i++) {
        even[i] = x[i * 2];
        odd[i] = x[i * 2 + 1];
    }

    // 再帰
    fft(even, _N / 2);
    fft(odd, _N / 2);

    // DFT計算
    for (int k = 0; k < _N / 2; k++) {
        std::complex<double> t = exp(std::complex<double>(0, -2 * M_PI * k / N)) * odd[k];
        x[k] = even[k] + t;
        x[N / 2 + k] = even[k] - t;
    }
}

int main() {
    std::complex<double> x_out[N]; // 元データ

    // 元データ作成
    for (int i = 0; i < N; i++) {
        x_out[i] = std::complex<double>(A * sin(2 * M_PI * F0 * i / Fs), 0.0); // t = i / Fs
    }

    fft(x_out, N);
    return 0;
}
