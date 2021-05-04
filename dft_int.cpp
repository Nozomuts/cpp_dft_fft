#include <chrono>
#include <fstream>
#include <iostream> // for cout

#define N 64    // 分割数
#define Fs 8000 // サンプリング周波数
#define A 1     // 振幅
#define F0 440  // 周波数
#define phi 0   // 初期位相

using namespace std;    // cout, endl, swap, ios, complex
using namespace chrono; // system_clock, duration_cast, microseconds, ofstream

int sin_table[] = {
    0,      9802,  19509, 29029, 38268, 47140, 55557, 63439,  70711,
    77301, 83147, 88192, 92388, 95694, 98079, 99519, 100000,
};

int use_table_sin(int i) {
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

int use_table_cos(int i) {
    i += N / 4;
    return use_table_sin(i);
}

void dft(int x_r[N], int x_i[N], int *dft_r, int *dft_i) {
    for (int k = 0; k < N; k++) {
        for (int n = 0; n < N; n++) {
            dft_r[k] += x_r[n] * use_table_cos(k * n) + x_i[n] * use_table_sin(k * n);
            dft_i[k] += x_r[n] * (-use_table_sin(k * n)) + x_i[n] * use_table_cos(k * n);
        }
    }
}

int main() {
    int x_r[N], x_i[N];
    int dft_r[N], dft_i[N]; // x_r,x_iは元データ兼fftのデータ

    // 元データ作成
    for (int i = 0; i < N; i++) {
        // x(t) = A * sin(2 * pi * F0 * t + phi) ( 0 <= t < 0.008 )
        x_r[i] = A * sin(2 * M_PI * F0 * i / Fs + phi) * 1000; // t = i / Fs
        x_i[i] = 0;
    }

    dft(x_r, x_i, dft_r, dft_i);

    // 結果をファイルに出力
    ofstream dft_ofs("dft_int.csv");
    for (int i = 0; i < N; i++) {
        dft_ofs << dft_r[i] << "," << dft_i[i] << endl;
    }
    dft_ofs.close();

    return 0;
}
