#include <complex>
#include <iostream> // for cout
#include <math.h>   // for sin(), cos()

#define N 80            // 分割数
#define Fs (double)8000 // サンプリング周波数
#define A (double)1     // 振幅
#define F0 (double)440  // 周波数

using namespace std; // cout, endl, string, to_string

void dft(std::complex<double> x[], int _N) {
    double SRC_re[N]; // 元データの実部
    double SRC_im[N]; // 元データの虚部
    double DFT_re[N]; // DFTの実部
    double DFT_im[N]; // DFTの虚部

    // DFT計算
    for (int k = 0; k < N; k++) {
        DFT_re[k] = 0.0;
        DFT_im[k] = 0.0;
        for (int n = 0; n < N; n++) {
            DFT_re[k] += SRC_re[n] * (cos((2 * M_PI / N) * k * n)) + SRC_im[n] * (sin((2 * M_PI / N) * k * n));
            DFT_im[k] += SRC_re[n] * (-sin((2 * M_PI / N) * k * n)) + SRC_im[n] * (cos((2 * M_PI / N) * k * n));
        }
        cout << k / Fs << " " << (Fs / N) * k << " " << SRC_re[k] << " " << SRC_im[k] << " " << DFT_re[k] << " "
             << DFT_im[k] << endl;
    }
}

// 離散フーリエ変換 f(t) = A * sin(2 * pi * F0 * t) ( 0 <= t < 0.01 )
int main() {
    double DFT_re[N]; // DFTの実部
    double DFT_im[N]; // DFTの虚部

    // ヘッダ出力 ( 時間, 角周波数, 元データ(実部), 元データ(虚部), DFT(実部), DFT(虚部) )
    printf("t, f, x_re, x_im, X_re, X_im\n");

    std::complex<double> x_out[N]; // 元データ

    // 元データ作成
    for (int i = 0; i < N; i++) {
        x_out[i] = std::complex<double>(A * sin(2 * M_PI * F0 * i / Fs), 0.0); // t = i / Fs
    }
    dft(x_out, N);

    return 0;
}