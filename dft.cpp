/*********************************************
 * 離散フーリエ変換                          *
 *   f(t) = 2 * sin(4 * t) + 3 * cos(2 * t)  *
 *          ( 0 <= t < 2 * pi )              *
 *********************************************/
#include <iostream>  // for cout
#include <math.h>    // for sin(), cos()
#include <stdio.h>   // for printf()

#define N 100                // 分割数
#define CSV_DFT  "DFT.csv"   // 出力ファイル (DFT)
#define CSV_IDFT "IDFT.csv"  // 出力ファイル (IDFT)

using namespace std;

/*
 * 計算クラス
 */
class Calc
{
    double SRC_re[N];   //元データの実部
    double SRC_im[N];   //元データの虚部
    double DFT_re[N];   //DFTの実部
    double DFT_im[N];   //DFTの虚部
    double IDFT_re[N];  //IDFTの実部
    double IDFT_im[N];  //IDFTの虚部

    public:
        void makeSourceData();  // 元データ作成
        void executeDFT();      // 離散フーリエ変換
        void executeIDFT();     // 逆離散フーリエ変換

    private:
        double calcTerm(int n, double x);  //各項計算
};

/*
 * 元データ作成
 */
void Calc::makeSourceData()
{
    int i;

    for (i = 0; i < N; i++) {
        SRC_re[i] = 2 * sin(4 * (2 * M_PI / N) * i)
                  + 3 * cos(2 * (2 * M_PI / N) * i);
        SRC_im[i] = 0.0;
    }
}

/*
 * 離散フーリエ変換
 */
void Calc::executeDFT()
{
    int k, n;  // LOOPインデックス
    FILE *pf;  // ファイルポインタ

    // 出力ファイルOPEN
    pf = fopen(CSV_DFT, "w");

    // ヘッダ出力 ( k, 角周波数, 元データ(実部), 元データ(虚部), DFT(実部), DFT(虚部) )
    fprintf(pf, "k,f,x_re,x_im,X_re,X_im\n");

    // 計算・結果出力
    for (k = 0; k < N; k++) {
        DFT_re[k] = 0.0;
        DFT_im[k] = 0.0;
        for (n = 0; n < N; n++) {
            DFT_re[k] += SRC_re[n] * ( cos((2 * M_PI / N) * k * n))
                       + SRC_im[n] * ( sin((2 * M_PI / N) * k * n));
            DFT_im[k] += SRC_re[n] * (-sin((2 * M_PI / N) * k * n))
                       + SRC_im[n] * ( cos((2 * M_PI / N) * k * n));
        }
        fprintf(pf, "%d,%lf,%lf,%lf,%lf,%lf\n",
            k, (2 * M_PI / N) * k, SRC_re[k], SRC_im[k], DFT_re[k], DFT_im[k]);
    }

    // 出力ファイルCLOSE
    fclose(pf);
}

/*
 * メイン処理
 */
int main()
{
    try
    {
        // 計算クラスインスタンス化
        Calc objCalc;

        // 元データ作成
        objCalc.makeSourceData();

        // 離散フーリエ変換
        objCalc.executeDFT();
    }
    catch (...) {
        cout << "ERROR" << endl;
        return -1;
    }

    // 正常終了
    return 0;
}