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

double sin_table[] = {
    0,       0.0980171, 0.19509,  0.290285, 0.382683, 0.471397, 0.55557,  0.634393, 0.707107,
    0.77301, 0.83147,   0.881921, 0.92388,  0.95694,  0.980785, 0.995185, 1,
};

double use_table_sin(int i) {
    int n = i % N;
    if (n <= N / 4) {
        return sin_table[n];
    } else if (n <= N / 2) {
        return sin_table[N / 2 - n];
    } else if (n <= 3 * N / 4) {
        return -sin_table[n - N / 2];
    } else if (n <= N) {
        return -sin_table[N - n];
    } else {
        printf("Error!");
        return -1;
    }
}

double use_table_cos(int i) {
    i += N / 4;
    return use_table_sin(i);
}

void dft(double x_r[N], double x_i[N], double *dft_r, double *dft_i) {
    for (int k = 0; k < N; k++) {
        for (int n = 0; n < N; n++) {
            dft_r[k] += x_r[n] * use_table_cos(k * n) + x_i[n] * use_table_sin(k * n);
            dft_i[k] += x_r[n] * (-use_table_sin(k * n)) + x_i[n] * use_table_cos(k * n);
        }
    }
}

int main() {
    double x_r[N], x_i[N], dft_r[N], dft_i[N]; // x_r,x_iは元データ兼fftのデータ

    // 元データ作成
    for (int i = 0; i < N; i++) {
        // x(t) = A * sin(2 * pi * F0 * t + phi) ( 0 <= t < 0.008 )
        x_r[i] = A * sin(2 * M_PI * F0 * i / Fs + phi); // t = i / Fs
        x_i[i] = 0;
    }

    dft(x_r, x_i, dft_r, dft_i);

    // 結果をファイルに出力
    ofstream dft_ofs("dft.csv");
    for (int i = 0; i < N; i++) {
        dft_ofs << dft_r[i] << "," << dft_i[i] << endl;
    }
    dft_ofs.close();

    return 0;
}
