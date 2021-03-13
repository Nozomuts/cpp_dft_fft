#include <iostream> // for cout
#include <math.h>   // for sin(), cos()
#include <stdio.h>  // for printf()

#define N 80            // 分割数
#define Fs (double)8000 // サンプリング周波数
#define A (double)1     // 振幅
#define F0 (double)440  // 周波数

using namespace std; // cout, endl, string, to_string

// 実部と虚部を合わせて、フォーマット
string format(double re, double im) {

    if (re == 0 && im == 0)
        return "0";
    else if (re == 0)
        return to_string(im) + "i";
    else if (im == 0)
        return to_string(re);
    else
        return to_string(re) + (im > 0 ? "+" + to_string(im) : to_string(im)) + "i";
}

// 離散フーリエ変換 f(t) = A * sin(2 * pi * F0 * t) ( 0 <= t < 0.01 )
int main() {
    double SRC_re[N]; // 元データの実部
    double SRC_im[N]; // 元データの虚部
    double DFT_re[N]; // DFTの実部
    double DFT_im[N]; // DFTの虚部

    // ヘッダ出力 ( 時間, 角周波数, 元データ(実部)+元データ(虚部), DFT(実部)+DFT(虚部) )
    printf("t,f,x(t),X(ω)\n");

    // 元データ作成
    for (int i = 0; i < N; i++) {
        SRC_re[i] = A * sin(2 * M_PI * F0 * i / Fs); // t = i / Fs
        SRC_im[i] = 0.0;
    }

    // 計算・結果出力
    for (int k = 0; k < N; k++) {
        DFT_re[k] = 0.0;
        DFT_im[k] = 0.0;
        for (int n = 0; n < N; n++) {
            DFT_re[k] += SRC_re[n] * (cos((2 * M_PI / N) * k * n)) + SRC_im[n] * (sin((2 * M_PI / N) * k * n));
            DFT_im[k] += SRC_re[n] * (-sin((2 * M_PI / N) * k * n)) + SRC_im[n] * (cos((2 * M_PI / N) * k * n));
        }
        printf("%lf,%lf,%s,%s\n", k / Fs, (Fs / N) * k, format(SRC_re[k], SRC_im[k]).c_str(),
               format(DFT_re[k], DFT_im[k]).c_str());
    }
    return 0;
}