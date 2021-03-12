/*********************************************
 * 離散フーリエ変換                             *
 *   f(t) = A * sin(2 * pi * F0 * t)         *
 *          ( 0 <= t < 0.01 )                *
 *********************************************/
#include <iostream> // for cout
#include <math.h>   // for sin(), cos()
#include <stdio.h>  // for printf()

#define N 800             // 分割数
#define A 1               // 振幅
#define F0 440            // 周波数
#define CSV_DFT "DFT.csv" // 出力ファイル (DFT)

int main() {
    double SRC_re[N]; //元データの実部
    double SRC_im[N]; //元データの虚部
    double DFT_re[N]; // DFTの実部
    double DFT_im[N]; // DFTの虚部

    try {
        // 元データ作成
        for (int i = 0; i < N; i++) {
            SRC_re[i] = A * sin((2 * M_PI * F0 / N) * i);
            SRC_im[i] = 0.0;

            // ファイル構造体のポインタを作成
            FILE *pf;

            // 出力ファイルOPEN(書き込み専用)
            pf = fopen(CSV_DFT, "w");

            // ヘッダ出力 ( k, 角周波数, 元データ(実部), 元データ(虚部), DFT(実部), DFT(虚部) )
            fprintf(pf, "k,f,x_re,x_im,X_re,X_im\n");

            // 計算・結果出力
            for (int k = 0; k < N; k++) {
                DFT_re[k] = 0.0;
                DFT_im[k] = 0.0;
                for (int n = 0; n < N; n++) {
                    DFT_re[k] += SRC_re[n] * (cos((2 * M_PI / N) * k * n)) + SRC_im[n] * (sin((2 * M_PI / N) * k * n));
                    DFT_im[k] += SRC_re[n] * (-sin((2 * M_PI / N) * k * n)) + SRC_im[n] * (cos((2 * M_PI / N) * k * n));
                }
                fprintf(pf, "%d,%lf,%lf,%lf,%lf,%lf\n", k, (2 * M_PI / N) * k, SRC_re[k], SRC_im[k], DFT_re[k],
                        DFT_im[k]);
            }

            // 出力ファイルCLOSE
            fclose(pf);
        }

    } catch (...) {
        std::cout << "ERROR" << std::endl;
        return -1;
    }

    return 0;
}