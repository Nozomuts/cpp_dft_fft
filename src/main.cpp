#include <chrono>
#include <fstream>
#include <iostream> // for cout
#include "main.h"

#define N 1024  // 分割数
#define Fs 8000 // サンプリング周波数
#define A 1     // 振幅
#define F0 440  // 周波数
#define phi 0   // 初期位相
#define DIVISOR 1000

using namespace std;    // cout, endl, swap, ios, complex
using namespace chrono; // system_clock, duration_cast, microseconds, ofstream

int main() {
    short x_r[N], x_i[N];
    double y_r[N], y_i[N]; // x_r,x_iは元データ兼fftのデータ
    double sum = 0;

    srand(10);
    // 元データ作成
    for (int i = 0; i < N; i++) {
        // x(t) = A * sin(2 * pi * F0 * t + phi) ( 0 <= t < 0.008 )
        double rand1 = (double)rand() / RAND_MAX;
        double rand2 = (double)rand() / RAND_MAX;
        x_r[i] = rand1 * DIVISOR; // t = i / Fs
        x_i[i] = rand2 * DIVISOR;
        y_r[i] = rand1;
        y_i[i] = rand2;
    }

    create_table_short();
    create_table();

    fft(y_r, y_i);
    fft_short(x_r, x_i);
    bit_reverse(y_r, y_i);
    bit_reverse_short(x_r, x_i);

    // 結果をファイルに出力
    ofstream fft_ofs("int.csv");
    ofstream fft2_ofs("int2.csv");
    ofstream fft3_ofs("int3.csv");
    for (int i = 0; i < N; i++) {
        fft_ofs << y_r[i] - x_r[i] / (double)DIVISOR << "," << y_i[i] - x_i[i] / (double)DIVISOR << endl;
        fft2_ofs << y_r[i] << "," << y_i[i] << endl;
        fft3_ofs << x_r[i] / (double)DIVISOR << "," << x_i[i] / (double)DIVISOR << endl;
    }
    fft_ofs.close();
    fft2_ofs.close();
    fft3_ofs.close();

    for (int i = 0; i < N; i++) {
        // √(実部^2+虚部^2)^2 (差分の2乗)
        sum += ((y_r[i] - x_r[i] / (double)DIVISOR) * (y_r[i] - x_r[i] / (double)DIVISOR) +
                (y_i[i] - x_i[i] / (double)DIVISOR) * (y_i[i] - x_i[i] / (double)DIVISOR));
    }

    // dftとfftの結果の比較(差分の2乗平均の平方根)
    cout << sqrt(sum / N) << endl;

    return 0;
}
