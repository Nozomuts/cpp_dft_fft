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

int sin_table_int[N / 4 + 1];
double sin_table[N / 4 + 1];

// テーブルを用いたsin関数
double use_table_sin(int n) {
    n %= N;
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

// テーブルを用いたcos関数
double use_table_cos(int n) {
    n += N / 4;
    return use_table_sin(n);
}

int use_table_sin_int(int i) {
    int n = i % N;
    if (n <= N / 4) {
        return sin_table_int[n];
    } else if (n <= N / 2) {
        return sin_table_int[N / 2 - n];
    } else if (n <= 3 * N / 4) {
        return -sin_table_int[n - N / 2];
    } else if (n <= N) {
        return -sin_table_int[N - n];
    } else {
        printf("Error!");
        return -1;
    }
}

int use_table_cos_int(int i) {
    i += N / 4;
    return use_table_sin_int(i);
}

// テーブル作成
void create_table_int() {
    sin_table_int[0] = 0;

    for (int i = 1; i <= N / 4; i++) {
        sin_table_int[i] = sin(2 * M_PI / N * i) * 10000;
    }
}

// テーブル作成
void create_table() {
    sin_table[0] = 0;

    for (int i = 1; i <= N / 4; i++) {
        sin_table[i] = sin(2 * M_PI / N * i);
    }
}

// 高速フーリエ変換
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
                x_r[i * m + j + m / 2] =
                    (a_r - b_r) * use_table_cos(N / m * j) + (a_i - b_i) * use_table_sin(N / m * j);
                x_i[i * m + j + m / 2] =
                    (a_r - b_r) * (-use_table_sin(N / m * j)) + (a_i - b_i) * use_table_cos(N / m * j);
            }
        }
        m /= 2;
    }
}

// 高速フーリエ変換
void fft_int(int *x_r, int *x_i) {
    int m = N;
    while (m > 1) {
        for (int i = 0; i < N / m; i++) {
            for (int j = 0; j < m / 2; j++) {
                int a_r = x_r[i * m + j];
                int a_i = x_i[i * m + j];
                int b_r = x_r[i * m + j + m / 2];
                int b_i = x_i[i * m + j + m / 2];
                x_r[i * m + j] = a_r + b_r;
                x_i[i * m + j] = a_i + b_i;
                x_r[i * m + j + m / 2] = (a_r - b_r) * use_table_cos_int(N / m * j) / 10000 +
                                         (a_i - b_i) * use_table_sin_int(N / m * j) / 10000;
                x_i[i * m + j + m / 2] = (a_r - b_r) * (-use_table_sin_int(N / m * j)) / 10000 +
                                         (a_i - b_i) * use_table_cos_int(N / m * j) / 10000;
            }
        }
        m /= 2;
    }
}

// ビット反転並べ替え
void bit_reverse_int(int *x_r, int *x_i) {
    for (int i = 0, j = 1; j < N; j++) {
        for (int k = N >> 1; k > (i ^= k); k >>= 1)
            ;
        if (i < j) {
            swap(x_r[i], x_r[j]); // 入れ替え
            swap(x_i[i], x_i[j]);
        }
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
    int x_r[N], x_i[N];
    double y_r[N], y_i[N]; // x_r,x_iは元データ兼fftのデータ

    // 元データ作成
    for (int i = 0; i < N; i++) {
        // x(t) = A * sin(2 * pi * F0 * t + phi) ( 0 <= t < 0.008 )
        x_r[i] = A * sin(2 * M_PI * F0 * i / Fs + phi) * 10000; // t = i / Fs
        x_i[i] = 0;
        y_r[i] = A * sin(2 * M_PI * F0 * i / Fs + phi); // t = i / Fs
        y_i[i] = 0;
    }

    create_table_int();
    create_table();

    fft(y_r, y_i);
    fft_int(x_r, x_i);
    bit_reverse(y_r, y_i);
    bit_reverse_int(x_r, x_i);

    // 結果をファイルに出力
    ofstream fft_ofs("int.csv");
    for (int i = 0; i < N; i++) {
        fft_ofs << y_r[i] - x_r[i] / (double)10000 << "," << y_i[i] - x_i[i] / (double)10000 << endl;
    }
    fft_ofs.close();

    return 0;
}
