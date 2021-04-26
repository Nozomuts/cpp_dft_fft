#include <chrono>
#include <fstream>
#include <iostream> // for cout

#define N 64            // 分割数
#define Fs (double)8000 // サンプリング周波数
#define A (double)1     // 振幅
#define F0 (double)440  // 周波数
#define phi (double)0   // 初期位相

using namespace std;    // cout, endl, swap, ios, complex
using namespace chrono; // system_clock, duration_cast, microseconds, ofstream

void dft(double x_r[N], double x_i[N], double *dft_r, double *dft_i) {
    for (int k = 0; k < N; k++) {
        for (int n = 0; n < N; n++) {
            dft_r[k] += x_r[n] * cos(2 * M_PI / N * k * n) + x_i[n] * sin(2 * M_PI / N * k * n);
            dft_i[k] += x_r[n] * (-sin(2 * M_PI / N * k * n)) + x_i[n] * cos(2 * M_PI / N * k * n);
        }
    }
}

void fft(double *x_r, double *x_i) {
    int m = N;
    while (m > 1) {
        for (int i = 0; i < N / m; i++) {
            for (int j = 0; j < m / 2; j++) {
                double a_r = x_r[i * m + j];
                double a_i = x_i[i * m + j];
                double b_r = x_r[i * m + j + m / 2];
                double b_i = x_i[i * m + j + m / 2];
                x_r[i * m + j] = a_r + b_r;
                x_i[i * m + j] = a_i + b_i;
                x_r[i * m + j + m / 2] = (a_r - b_r) * cos(2 * M_PI / m * j) + (a_i - b_i) * sin(2 * M_PI / m * j);
                x_i[i * m + j + m / 2] = (a_r - b_r) * (-sin(2 * M_PI / m * j)) + (a_i - b_i) * cos(2 * M_PI / m * j);
            }
        }
        m /= 2;
    }
}

// ビット反転並べ替え
void bit_reverse(double *x_r, double *x_i) {
    for (int i = 0, j = 1; j < N; j++) {
        for (int k = N >> 1; k > (i ^= k); k >>= 1)
            ;
        if (i < j) {
            swap(x_r[i], x_r[j]); // 入れ替え
            swap(x_i[i], x_i[j]);
        }
    }
}

int main() {
    double x_r[N], x_i[N], dft_r[N], dft_i[N]; // x_r,x_iは元データ兼fftのデータ
    double sum = 0;
    system_clock::time_point start, end;

    // 元データ作成
    for (int i = 0; i < N; i++) {
        // x(t) = A * sin(2 * pi * F0 * t + phi) ( 0 <= t < 0.008 )
        x_r[i] = A * sin(2 * M_PI * F0 * i / Fs + phi); // t = i / Fs
        x_i[i] = 0;
    }

    dft(x_r, x_i, dft_r, dft_i);

    fft(x_r, x_i);
    bit_reverse(x_r, x_i);

    // 結果をファイルに出力
    ofstream dft_ofs("dft.csv");
    for (int i = 0; i < N; i++) {
        dft_ofs << dft_r[i] << "," << dft_i[i] << endl;
    }
    dft_ofs.close();

    ofstream fft_ofs("fft.csv");
    for (int i = 0; i < N; i++) {
        fft_ofs << x_r[i] << "," << x_i[i] << endl;
    }
    fft_ofs.close();

    for (int i = 0; i < N; i++) {
        // √(実部^2+虚部^2)^2 (差分の2乗)
        sum += ((x_r[i] - dft_r[i]) * (x_r[i] - dft_r[i]) + (x_i[i] - dft_i[i]) * (x_i[i] - dft_i[i]));
    }

    // dftとfftの結果の比較(差分の2乗平均の平方根)
    cout << sqrt(sum / N) << endl;

    return 0;
}
