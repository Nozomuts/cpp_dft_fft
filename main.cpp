#include <chrono>
#include <complex>
#include <fstream>
#include <iostream> // for cout

#define N 64            // 分割数
#define Fs (double)8000 // サンプリング周波数
#define A (double)1     // 振幅
#define F0 (double)440  // 周波数
#define phi (double)0   // 初期位相

using namespace std;    // cout, endl, swap
using namespace chrono; // system_clock, duration_cast, microseconds, ofstream

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
            swap(x[i], x[j]); // x[i]とx[j]を交換する
    }
}

int main() {
    complex<double> x_out[N]; // 元データ
    system_clock::time_point start, end;

    // 元データ作成
    for (int i = 0; i < N; i++) {
        // x(t) = A * sin(2 * pi * F0 * t + phi) ( 0 <= t < 0.008 )
        x_out[i] = complex<double>(A * sin(2 * M_PI * F0 * i / Fs + phi), 0); // t = i / Fs
    }

    cout << "- - - DFT - - -" << endl;
    start = system_clock::now(); // 計測スタート時刻を保存
    dft(x_out, N);
    end = system_clock::now(); // 計測終了時刻を保存
    // 要した時間を計算
    double dft_time = static_cast<double>(duration_cast<microseconds>(end - start).count() / 1000.0);

    cout << "- - - FFT - - -" << endl;
    start = system_clock::now();
    fft(x_out, N, 0);
    bit_reverse(x_out, N);
    end = system_clock::now();
    double fft_time = static_cast<double>(duration_cast<microseconds>(end - start).count() / 1000.0);

    for (int i = 0; i < N; i++) {
        cout << x_out[i] << endl;
    }

    // 要した時間をミリ秒（1/1000秒）に変換して表示
    cout << "DFT: " << dft_time << " ms" << endl;
    cout << "FFT: " << fft_time << " ms" << endl;

    ofstream dft_ofs("dft.txt");
    dft_ofs << dft_time;
    dft_ofs.close();
    ofstream fft_ofs("fft.txt");
    fft_ofs << fft_time;
    fft_ofs.close();

    return 0;
}
