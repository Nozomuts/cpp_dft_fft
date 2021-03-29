#include <chrono>
#include <complex>
#include <iostream> // for cout

#define N 64            // 分割数
#define Fs (double)8000 // サンプリング周波数
#define A (double)1     // 振幅
#define F0 (double)440  // 周波数
#define phi (double)0   // 初期位相

using namespace std;    // cout, endl
using namespace chrono; // system_clock, duration_cast, microseconds

void dft(complex<double> x[], int _N) {
    complex<double> dft_out[N];

    // 出力用の配列初期化
    for (int i = 0; i < N; i++) {
        dft_out[i] = complex<double>(0, 0);
    }

    // DFT計算
    for (int k = 0; k < N; k++) {
        for (int n = 0; n < N; n++) {
            dft_out[k] += x[n] * exp(complex<double>(0, -2 * M_PI * k * n / N)); // (実部,虚部)
        }
        cout << dft_out[k] << endl;
    }
}

// 離散フーリエ変換
int main() {
    complex<double> x_out[N]; // 元データ
    system_clock::time_point start, end;

    // 元データ作成
    for (int i = 0; i < N; i++) {
        // x(t) = A * sin(2 * pi * F0 * t + phi) ( 0 <= t < 0.008 )
        x_out[i] = complex<double>(A * sin(2 * M_PI * F0 * i / Fs + phi), 0); // t = i / Fs
    }
    start = system_clock::now(); // 計測スタート時刻を保存
    dft(x_out, N);
    end = system_clock::now(); // 計測終了時刻を保存
    double time = static_cast<double>(duration_cast<microseconds>(end - start).count() / 1000.0); // 要した時間を計算
    // 要した時間をミリ秒（1/1000秒）に変換して表示
    cout << time << " ms \n";

    return 0;
}